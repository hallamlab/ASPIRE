import importlib.util
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform


SCRIPT = Path(__file__).parents[1] / "processes/community_predictor_comparison/community_predictor_comparison.py"


def load_module():
    spec = importlib.util.spec_from_file_location("community_predictor_comparison", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_nested_distance_model_detects_added_signal():
    module = load_module()
    rng = np.random.default_rng(4)
    n = 36
    frame = pd.DataFrame({
        "Season": np.tile(["winter", "summer", "spring"], 12),
        "signal": np.linspace(-2, 2, n),
    })
    response = np.column_stack([
        frame["signal"].to_numpy() + rng.normal(0, 0.1, n),
        -frame["signal"].to_numpy() + rng.normal(0, 0.1, n),
    ])
    distance = squareform(pdist(response, metric="euclidean"))
    stats = module.nested_distance_stats(
        distance, frame, ["Season"], ["signal"], {"signal"}
    )
    assert stats["incremental_r2"] > 0.8
    assert stats["adjusted_r2"] > stats["baseline_r2"]


def test_cross_validation_uses_held_out_blocks_and_compares_models():
    module = load_module()
    rng = np.random.default_rng(7)
    blocks = np.repeat(["A", "B", "C", "D"], 8)
    predictor = np.tile(np.linspace(-1, 1, 8), 4)
    frame = pd.DataFrame({"block": blocks, "predictor": predictor})
    response = np.column_stack([
        predictor + rng.normal(0, 0.05, len(predictor)),
        predictor * 0.5 + rng.normal(0, 0.05, len(predictor)),
    ])
    folds, summary = module.multivariate_cv(
        response,
        frame,
        {
            "intercept": ([], set()),
            "predictor": (["predictor"], {"predictor"}),
        },
        "block",
    )
    assert folds["fold"].nunique() == 4
    scores = summary.set_index("model")["cv_r2_mean"]
    assert scores["predictor"] > scores["intercept"]
    difference = module.paired_difference(folds, "predictor", "intercept", 42)
    assert difference.loc[0, "mean_difference_a_minus_b"] > 0


def test_cruise_groups_are_compared_on_the_matched_cruise_cohort(tmp_path, monkeypatch):
    module = load_module()
    rows = []
    counts = {}
    for cruise in range(1, 13):
        year = 2020 + (cruise - 1) // 3
        group = "cruise_group_1" if cruise % 2 else "cruise_group_2"
        for bottle in range(2):
            sample = f"s{cruise}_{bottle}"
            rows.append({
                "sampleID": sample, "Cruise": cruise, "Year": year,
                "Season": ["winter", "spring", "summer"][cruise % 3],
                "date": f"{year}-{(cruise % 12) + 1:02d}-01", "Depth": 10 + bottle * 50,
                "pea_J_m3": 100 + cruise * 4, "depth_centroid_distance": 20 + cruise,
                "cruise_group": group, "cruise_group_max_prob": 0.95,
                "cruise_group_uncertain": False,
                "renewal_phase": (
                    "renewal" if cruise % 4 == 0
                    else "post-renewal" if cruise % 4 == 1
                    else "baseline"
                ),
                "o2_compartment": "oxic" if bottle == 0 else "dysoxic",
                "gmm_component": bottle, "o2_subcompartment_final": f"g{bottle}",
                "o2_subcompartment_final_assignment_source": "observed",
            })
            counts[sample] = [80 if group.endswith("1") else 10, 10 if group.endswith("1") else 80, 10 + cruise]
    metadata = tmp_path / "metadata.tsv"
    count_path = tmp_path / "counts.tsv"
    outdir = tmp_path / "results"
    pd.DataFrame(rows).to_csv(metadata, sep="\t", index=False)
    pd.DataFrame(counts, index=["asv1", "asv2", "asv3"]).to_csv(count_path, sep="\t")
    monkeypatch.setattr(sys, "argv", [
        "community_predictor_comparison.py", "--metadata", str(metadata),
        "--asv-counts", str(count_path), "--outdir", str(outdir),
        "--permutations", "19", "--formats", "png",
    ])
    module.main()
    comparison = pd.read_csv(outdir / "tables/cruise_group_model_comparison.tsv", sep="\t")
    assert set(comparison["model"]) == {
        "Season baseline", "PEA", "Depth centroid distance", "PEA + centroid", "Cruise group"
    }
    assert (outdir / "tables/cruise_group_paired_performance.tsv").exists()
    assert (outdir / "plots/cruise_group_ordination.png").exists()
    renewal = pd.read_csv(outdir / "tables/renewal_group_model_comparison.tsv", sep="\t")
    assert set(renewal["model"]) == {
        "Season baseline", "PEA", "Depth centroid distance", "PEA + centroid",
        "Renewal phase",
    }
    assert (outdir / "tables/renewal_group_paired_performance.tsv").exists()
    assert (outdir / "plots/renewal_group_ordination.png").exists()
    stratified = pd.read_csv(
        outdir / "tables/cruise_group_compartment_model_comparison.tsv", sep="\t"
    )
    assert set(stratified["cruise_group"]) == {"cruise_group_1", "cruise_group_2"}
    assert set(stratified["model"]) == {"Depth + season baseline", "Hybrid compartment"}
    assert (outdir / "tables/cruise_group_compartment_paired_performance.tsv").exists()
    paired_compartments = pd.read_csv(
        outdir / "tables/sample_compartment_paired_performance.tsv", sep="\t"
    )
    assert len(paired_compartments) == 3
    assert "q_value" in paired_compartments
    refinement = pd.read_csv(
        outdir / "tables/sample_hybrid_parent_refinement.tsv", sep="\t"
    )
    assert set(refinement["parent_classification"]) == {"Legacy O2", "GMM"}
    assert {"q_value", "cv_q_value", "partial_r2"} <= set(refinement)
    pairwise = pd.read_csv(outdir / "tables/cruise_group_compartment_pairwise.tsv", sep="\t")
    assert {"q_value", "permdisp_q_value", "compartment_a", "compartment_b"} <= set(pairwise)
    assert (outdir / "plots/cruise_group_compartment_ordination.png").exists()
    renewal_stratified = pd.read_csv(
        outdir / "tables/renewal_group_compartment_model_comparison.tsv", sep="\t"
    )
    assert set(renewal_stratified["renewal_phase"]) == {
        "baseline", "renewal", "post-renewal"
    }
    assert set(renewal_stratified["model"]) == {
        "Depth + season baseline", "Hybrid compartment"
    }
    assert (outdir / "tables/renewal_group_compartment_pairwise.tsv").exists()
    assert (outdir / "plots/renewal_group_compartment_ordination.png").exists()
