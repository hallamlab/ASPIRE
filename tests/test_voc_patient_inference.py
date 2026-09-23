"""Small deterministic checks of patient-level VOC inference (no patient data)."""
import importlib.util
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata

SCRIPT = Path(__file__).parents[1] / "processes/voc_correlation/plot_voc_corr.py"
spec = importlib.util.spec_from_file_location("voc_patient_test", SCRIPT)
voc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(voc)


def fixture():
    idx = [f"s{i}" for i in range(8)]
    counts = pd.DataFrame({"a": [1, 3, 2, 4, 5, 6, 8, 9],
                           "b": [9, 7, 8, 6, 5, 4, 2, 1], "zero": [0]*8}, index=idx)
    chem = pd.DataFrame({"v": [1, 3, 2, 4, 5, 6, 8, 9]}, index=idx)
    meta = pd.DataFrame({"normalized_sample_id": idx,
                         "patient_id": ["p1", "p1", "p2", "p3", "p4", "p5", "p6", "p7"],
                         "case_status": ["Control"]*8})
    return counts, chem, meta


def run(counts, chem, meta, candidates=("a", "zero"), **kwargs):
    return voc.patient_correlation_results(counts, chem, meta, candidates,
                                          n_permutations=99, **kwargs)


def test_patient_averaging_normalization_and_reproducibility():
    c, v, m = fixture()
    r, rel, clr, summary = run(c, v, m)
    assert len(rel) == 7 and summary["patients"] == 7
    assert np.isclose(rel.loc["p1", "a"], .2)
    assert np.allclose(rel.sum(axis=1), 1)
    assert np.allclose(clr.sum(axis=1), 0)
    pd.testing.assert_frame_equal(r, run(c, v, m)[0])
    a = r.query("asv == 'a' and normalization == 'relative_abundance'").iloc[0]
    assert np.isclose(a.rho, spearmanr([.2,.2,.4,.5,.6,.8,.9], [2,2,4,5,6,8,9]).statistic)
    assert a.n_patients == 7
    assert r.query("asv == 'zero'").status.eq("too_few_detected_patients").all()


def test_primary_is_invariant_to_sample_depth_and_uses_full_denominator():
    c, v, m = fixture()
    r, rel, _, _ = run(c, v, m, candidates=["a"])
    scaled = c.mul(np.arange(1, 9), axis=0)
    r2, rel2, _, _ = run(scaled, v, m, candidates=["a"])
    pd.testing.assert_frame_equal(rel, rel2)
    pd.testing.assert_frame_equal(r.query("normalization == 'relative_abundance'"),
                                  r2.query("normalization == 'relative_abundance'"))
    assert (rel.a < 1).all()  # b remains in denominator despite not being tested


def test_duplicate_all_replicates_does_not_increase_patient_n():
    c, v, m = fixture()
    c2, v2, m2 = c.copy(), v.copy(), m.copy()
    c2.index = c2.index + "rep"
    v2.index = v2.index + "rep"
    m2.normalized_sample_id += "rep"
    pd.testing.assert_frame_equal(run(c, v, m)[0],
        run(pd.concat([c,c2]), pd.concat([v,v2]), pd.concat([m,m2]))[0])


def test_missing_voc_uses_cognate_samples_before_averaging():
    c, v, m = fixture()
    v.loc["s0", "v"] = np.nan
    r = run(c, v, m)[0].query("asv == 'a' and normalization == 'relative_abundance'").iloc[0]
    assert np.isclose(r.rho, spearmanr([.3,.2,.4,.5,.6,.8,.9], [3,2,4,5,6,8,9]).statistic)


def test_small_and_constant_inputs_report_unavailable():
    c, v, m = fixture()
    r = run(c.iloc[:3], v.iloc[:3], m.iloc[:3])[0]
    assert r.status.eq("too_few_patients").all()
    assert r.p_value.isna().all()
    v[:] = 1
    r = run(c, v, m)[0]
    assert r.query("asv == 'a'").status.eq("constant_voc").all()


def test_exact_mann_whitney_with_ties_matches_enumeration():
    a, b = np.array([1,1,3]), np.array([2,3,4])
    u, p, draws, method = voc.permutation_mann_whitney(a,b,999)
    ranks = rankdata(np.r_[a,b]); center = 3*7/2
    obs = abs(ranks[:3].sum()-center)
    expected = np.mean([abs(ranks[list(i)].sum()-center) >= obs-1e-12
                        for i in itertools.combinations(range(6),3)])
    assert p == expected and draws == 20 and method == "exact_permutation"
    assert u == ranks[:3].sum()-6


def test_unknown_status_not_silently_control():
    assert voc.normalize_case_status(np.nan) == "Unknown"
    assert voc.normalize_case_status("not recorded") == "Unknown"


def test_invalid_patient_id_fails():
    c,v,m = fixture(); m.loc[0,"patient_id"] = np.nan
    import pytest
    with pytest.raises(ValueError, match="valid patient ID"):
        run(c,v,m)


def test_fdr_union_does_not_change_for_candidate_order_or_duplicates():
    c,v,m = fixture()
    pd.testing.assert_frame_equal(run(c,v,m,["a","b"])[0], run(c,v,m,["b","a","b"])[0])


def test_cli_end_to_end_on_synthetic_samples(tmp_path, monkeypatch):
    import sys
    import json
    c,v,m = fixture()
    c.T.to_csv(tmp_path / "counts.tsv", sep="\t")
    v.rename_axis("sample").to_csv(tmp_path / "voc.tsv", sep="\t")
    meta = m.rename(columns={"normalized_sample_id":"Sample", "patient_id":"Participant_ID",
                             "case_status":"Case"})
    meta.loc[meta.Participant_ID.isin(["p5","p6","p7"]), "Case"] = "Cancer"
    meta["Type_Group"] = "Bronchial Brush"
    meta = pd.concat([meta.assign(ASV_ID=a) for a in c.columns], ignore_index=True)
    meta.to_csv(tmp_path / "meta.tsv", sep="\t", index=False)
    monkeypatch.setattr(voc, "save_clustermap", lambda *a, **k: None)
    monkeypatch.setattr(voc, "save_case_voc_barplots", lambda *a, **k: None)
    monkeypatch.setattr(sys, "argv", [str(SCRIPT), "--asv-counts",str(tmp_path/"counts.tsv"),
        "--asv-meta",str(tmp_path/"meta.tsv"),"--voc",str(tmp_path/"voc.tsv"),
        "--outdir",str(tmp_path/"out"),"--sample-id-mode","exact","--patient-permutations","19",
        "--indicspecies-glob",str(tmp_path/"no_isa_*.tsv")])
    voc.main()
    summary = json.loads((tmp_path/"out/patient_inference_summary.json").read_text())
    assert summary["patients"] == 7
    assert summary["matched_samples"] == 8
    result = pd.read_csv(tmp_path/"out/patient_asv_voc_permutation_long.tsv", sep="\t")
    assert set(result.normalization) == {"relative_abundance","clr"}
    assert result.n_patients.max() == 7
    cases = pd.read_csv(tmp_path/"out/patient_voc_case_tests_brush.tsv", sep="\t")
    assert cases.n_cancer_patients.iloc[0] == 3
    assert cases.n_control_patients.iloc[0] == 4
