from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "processes" / "genome_cooccurrence" / "prepare_genome_cooccurrence.py"
SPEC = importlib.util.spec_from_file_location("prepare_genome_cooccurrence", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)

INFERENCE_SCRIPT = ROOT / "processes" / "genome_cooccurrence" / "infer_genome_proportionality.py"
INFERENCE_SPEC = importlib.util.spec_from_file_location("infer_genome_proportionality", INFERENCE_SCRIPT)
INFERENCE = importlib.util.module_from_spec(INFERENCE_SPEC)
assert INFERENCE_SPEC and INFERENCE_SPEC.loader
INFERENCE_SPEC.loader.exec_module(INFERENCE)


def _qc_row(genome: str, species: str, tier: str, length: int, completeness: float = 90.0):
    return {
        "Genome_Id": genome,
        "Species": species,
        "Domain": "Bacteria",
        "Phylum": "Pseudomonadota",
        "Class": "Gammaproteobacteria",
        "Order": "Order",
        "Family": "Family",
        "Genus": "Genus",
        "mimag_tier": tier,
        "sum_len": length,
        "Completeness": completeness,
        "Contamination": 1.0,
        "N50": 10000,
        "fasta_path": f"/{genome}.fna",
    }


def test_species_representative_hierarchy_prefers_mimag_then_length():
    qc = pd.DataFrame([
        _qc_row("medium_long", "Species A", "medium", 4_000_000, 99.0),
        _qc_row("high_short", "Species A", "high", 2_000_000, 80.0),
        _qc_row("b_short", "Species B", "medium", 2_000_000, 99.0),
        _qc_row("b_long", "Species B", "medium", 3_000_000, 70.0),
    ])
    selected = MODULE.select_representatives(qc, expected_n=2)
    winners = set(selected.loc[selected["selected_species_representative"], "Genome_Id"])
    assert winners == {"high_short", "b_long"}


def test_modality_matrix_closes_each_usable_sample(tmp_path: Path):
    qc = pd.DataFrame([
        _qc_row("g1", "Species A", "high", 2_000_000),
        _qc_row("g2", "Species B", "medium", 2_000_000),
    ])
    selected = MODULE.select_representatives(qc, expected_n=2)
    abundance = pd.DataFrame([
        {"sample_id": sample, "genome_id": genome, "fpkm": value}
        for sample, values in {
            "s1": (1.0, 3.0), "s2": (2.0, 2.0), "s3": (4.0, 0.0),
            "s4": (1.0, 1.0), "s5": (3.0, 2.0),
        }.items()
        for genome, value in zip(("g1", "g2"), values)
    ])
    matrix, audit = MODULE.prepare_modality(
        abundance,
        selected,
        sample_col="sample_id",
        genome_col="genome_id",
        value_col="fpkm",
        closure_total=1_000_000.0,
        input_scale="closed_abundance",
        input_feature_level="genome",
        label="metagenome",
        outdir=tmp_path,
    )
    assert matrix.shape == (2, 5)
    assert np.allclose(matrix.sum(axis=0).to_numpy(), 1_000_000.0)
    assert np.allclose(audit["prepared_matrix_sum"], 1_000_000.0)


def test_raw_recruitment_counts_are_not_renormalized(tmp_path: Path):
    qc = pd.DataFrame([
        _qc_row("g1", "Species A", "high", 2_000_000),
        _qc_row("g2", "Species B", "medium", 2_000_000),
    ])
    selected = MODULE.select_representatives(qc, expected_n=2)
    abundance = pd.DataFrame([
        {"sample_id": sample, "genome_id": genome, "read_count": value}
        for sample, values in {
            "s1": (10, 30), "s2": (20, 20), "s3": (40, 0),
            "s4": (10, 10), "s5": (30, 20),
        }.items()
        for genome, value in zip(("g1", "g2"), values)
    ])
    matrix, audit = MODULE.prepare_modality(
        abundance,
        selected,
        sample_col="sample_id",
        genome_col="genome_id",
        value_col="read_count",
        closure_total=1_000_000.0,
        input_scale="raw_counts",
        input_feature_level="genome",
        label="metagenome",
        outdir=tmp_path,
    )
    assert np.issubdtype(matrix.to_numpy().dtype, np.floating)
    assert matrix["s1"].sum() == 40
    assert audit.loc[audit["sample_id"].eq("s1"), "prepared_matrix_sum"].iloc[0] == 40
    assert audit["closure_total"].isna().all()


def test_species_sample_read_floor_precedes_transformation(tmp_path: Path):
    qc = pd.DataFrame([
        _qc_row("g1", "Species A", "high", 2_000_000),
        _qc_row("g2", "Species B", "medium", 2_000_000),
    ])
    selected = MODULE.select_representatives(qc, expected_n=2)
    abundance = pd.DataFrame([
        {"sample_id": sample, "genome_id": genome, "read_count": value}
        for sample, values in {
            "s1": (99, 100), "s2": (101, 40), "s3": (150, 150),
            "s4": (200, 0), "s5": (0, 250),
        }.items()
        for genome, value in zip(("g1", "g2"), values)
    ])
    matrix, audit = MODULE.prepare_modality(
        abundance, selected,
        sample_col="sample_id", genome_col="genome_id", value_col="read_count",
        closure_total=1_000_000.0, input_scale="raw_counts",
        input_feature_level="genome", label="metagenome", outdir=tmp_path,
        minimum_read_count=100,
    )
    assert matrix.loc["g1", "s1"] == 0
    assert matrix.loc["g2", "s1"] == 100
    assert matrix.loc["g1", "s2"] == 101
    assert matrix.loc["g2", "s2"] == 0
    floor = pd.read_csv(tmp_path / "metagenome_genome_sample_read_floor_audit.tsv", sep="\t")
    assert floor["set_to_zero_by_read_floor"].sum() == 2
    assert audit["minimum_read_count"].eq(100).all()


def test_genome_counts_are_averaged_within_species_including_zeros(tmp_path: Path):
    qc = pd.DataFrame([
        _qc_row("g1", "Species A", "high", 2_000_000),
        _qc_row("g1b", "Species A", "medium", 1_500_000),
        _qc_row("g2", "Species B", "medium", 2_000_000),
    ])
    selected = MODULE.select_representatives(qc, expected_n=2)
    abundance = pd.DataFrame([
        {"sample_id": sample, "genome_id": genome, "read_count": value}
        for sample, values in {
            "s1": (200, 400, 300), "s2": (200, 0, 300),
            "s3": (300, 500, 300), "s4": (400, 600, 300),
            "s5": (500, 700, 300),
        }.items()
        for genome, value in zip(("g1", "g1b", "g2"), values)
        if value != 0
    ])
    matrix, _ = MODULE.prepare_modality(
        abundance, selected,
        sample_col="sample_id", genome_col="genome_id", value_col="read_count",
        closure_total=1_000_000.0, input_scale="raw_counts",
        input_feature_level="genome", label="metagenome", outdir=tmp_path,
        minimum_read_count=100,
    )
    assert matrix.loc["g1", "s1"] == 300
    assert matrix.loc["g1", "s2"] == 100
    assert matrix.loc["g2", "s2"] == 300


def test_species_recruitment_counts_map_to_representative_node_ids(tmp_path: Path):
    qc = pd.DataFrame([
        _qc_row("g1", "Species A", "high", 2_000_000),
        _qc_row("g2", "Species B", "medium", 2_000_000),
    ])
    selected = MODULE.select_representatives(qc, expected_n=2)
    abundance = pd.DataFrame([
        {"sample_id": sample, "taxon": taxon, "read_count": value}
        for sample, values in {
            "s1": (100, 300), "s2": (200, 200), "s3": (400, 10),
            "s4": (100, 100), "s5": (300, 200),
        }.items()
        for taxon, value in zip(
            ("Pseudomonadota | Family A | Species A", "Bacteroidota | Family B | Species B"),
            values,
        )
    ])
    matrix, audit = MODULE.prepare_modality(
        abundance,
        selected,
        sample_col="sample_id",
        genome_col="taxon",
        value_col="read_count",
        closure_total=1_000_000.0,
        input_scale="raw_counts",
        input_feature_level="species",
        label="metagenome",
        outdir=tmp_path,
    )
    assert list(matrix.index) == ["g1", "g2"]
    assert matrix.loc["g1", "s1"] == 100
    assert matrix.loc["g2", "s1"] == 300
    assert audit["input_feature_level"].eq("species").all()


def test_proportionality_and_multiplicative_replacement_are_auditable():
    counts = np.asarray([[10, 20, 40, 80], [5, 10, 20, 40], [0, 3, 0, 6]], dtype=float)
    composition, audit = INFERENCE.multiplicative_replacement(counts, 0.65)
    assert np.allclose(composition.sum(axis=0), 1.0)
    assert (composition > 0).all()
    rho, vlr = INFERENCE.rho_p_matrix(INFERENCE.clr(composition))
    assert np.allclose(np.diag(rho), 1.0)
    assert np.allclose(np.diag(vlr), 0.0)
    assert rho[0, 1] > 0.99
    assert audit["closed_sum"].between(0.999999, 1.000001).all()


def test_network_outputs_preserve_disconnected_species(tmp_path: Path):
    selected = pd.DataFrame([{
        "feature_1": "g1", "feature_2": "g2", "rho_p": 0.8,
        "direction": "positive", "bootstrap_recovery": 1.0,
        "sign_consistency": 1.0, "q_value": 0.01,
    }])
    INFERENCE.network_outputs(["g1", "g2", "g3"], selected, tmp_path, "test", 42)
    graph = __import__("networkx").read_graphml(tmp_path / "test_network_all.graphml")
    nodes = pd.read_csv(tmp_path / "test_node_features.csv")
    clusters = pd.read_csv(tmp_path / "test_proportionality_clusters.tsv", sep="\t")
    assert set(graph.nodes()) == {"g1", "g2", "g3"}
    assert graph.degree("g3") == 0
    assert nodes.loc[nodes["Taxon"].eq("g3"), "EigenCentral"].iloc[0] == 0.0
    assert set(clusters["proportionality_cluster"].str.startswith("GN-PC")) == {True}


def test_workflow_uses_proportionality_and_retains_genome_isolates():
    workflow = (ROOT / "workflow" / "modules" / "networks_reporting.nf").read_text()
    graph_script = (ROOT / "processes" / "graph_network" / "graph_network.py").read_text()
    config = (ROOT / "examples" / "si.local.yml").read_text()
    assert "path(metagenome_abundance, stageAs: 'metagenome_abundance_long.tsv')" in workflow
    assert "path(metatranscriptome_abundance, stageAs: 'metatranscriptome_abundance_long.tsv')" in workflow
    assert "path(sample_metadata, stageAs: 'genome_cooccurrence_metadata.tsv')" in workflow
    assert '--input-scale "${genomeCooccurrenceInputScale}"' in workflow
    assert "metagenome_value_col: read_count" in config
    assert "metatranscriptome_value_col: read_count" in config
    assert "input_scale: raw_counts" in config
    assert "input_feature_level: genome" in config
    assert "genome_abundance_long.tsv" in config
    assert "minimum_read_count: 0" in config
    assert '--input-feature-level "${genomeCooccurrenceInputFeatureLevel}"' in workflow
    assert '--minimum-read-count ${genomeCooccurrenceMinimumReadCount}' in workflow
    assert "min_abs_rho: 0.30" in config
    assert "min_bootstrap_recovery: 0.80" in config
    assert "min_sign_consistency: 0.90" in config
    assert "permutations: 1000" in config
    assert '--min-abs-rho ${genomeCooccurrenceMinAbsRho}' in workflow
    assert '--min-bootstrap-recovery ${genomeCooccurrenceMinBootstrapRecovery}' in workflow
    assert '--min-sign-consistency ${genomeCooccurrenceMinSignConsistency}' in workflow
    assert 'infer_genome_proportionality.py' in (ROOT / "asv_pipeline.nf").read_text()
    assert "Channel.value(file(genomeCooccurrenceMetadata))" in (ROOT / "asv_pipeline.nf").read_text()
    assert 'Rscript "${spieceasiScriptPath}"' not in workflow.split("process GENOME_COOCCURRENCE", 1)[1].split("process GRAPH_NETWORK", 1)[0]
    assert "--module-subnetworks" in workflow
    assert "uniform module zoom without node displacement" in graph_script
    assert r'--pair-statistics "\${outdir}/\${prefix}_all_pair_statistics.tsv"' in workflow
    assert "proportionality_heatmap_all_pairs" in (ROOT / "processes" / "genome_cooccurrence" / "plot_genome_cooccurrence.py").read_text()
    assert "network_propagated_microbial_compartment" in workflow
    assert "network_propagated_renewal_associated_modules" in workflow
    assert "network_propagated_renewal_post_renewal_vs_stagnation" in workflow
    assert "proportionality_clusters.tsv" in workflow
    assert "module-renewal-profiles" in workflow
