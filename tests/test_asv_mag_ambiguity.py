import importlib.util
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
import pytest


SCRIPT = Path(__file__).parents[1] / "processes" / "asv_mag_network" / "asv_mag_network.py"
SPEC = importlib.util.spec_from_file_location("asv_mag_network", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_paper_outputs_preserve_recomputed_multimag_multiplicity():
    graph = nx.Graph()
    graph.add_node("ASV1")
    validation = pd.DataFrame([
        {
            "ASV_ID": "ASV1", "genome_id": "MAG1", "mag_match_id": "MAG1",
            "analysis_eligible_pair": True, "eligible_mag_count": 2,
        },
        {
            "ASV_ID": "ASV1", "genome_id": "MAG2", "mag_match_id": "MAG2",
            "analysis_eligible_pair": True, "eligible_mag_count": 2,
        },
    ])

    nodes, _edges, _paper_graph = MODULE.make_paper_outputs(
        graph,
        pd.DataFrame({"ASV_ID": ["ASV1"]}),
        validation,
        pd.DataFrame(),
    )

    row = nodes.set_index("ASV_ID").loc["ASV1"]
    assert row["eligible_mag_count"] == 2
    assert row["eligible_mag_ids"] == "MAG1|MAG2"


def test_plate_suffix_harmonization_recovers_shared_samples():
    asv = pd.DataFrame(
        [[1, 2]], index=["ASV1"], columns=["SI001_10m_plate8", "SI001_20m_plate8"],
    )
    mag = pd.DataFrame(
        [[3, 4]], index=["MAG1"], columns=["SI001_10m", "SI001_20m"],
    )
    asv_out, _mag_out, audit = MODULE.harmonize_abundance_samples(asv, mag)

    assert list(asv_out.columns) == ["SI001_10m", "SI001_20m"]
    assert audit["sample_harmonization_mode"] == "strip_terminal_plate"
    assert audit["shared_sample_count_used"] == 2


def test_heterogeneous_network_keeps_full_asv_context_and_correlation():
    graph = nx.Graph()
    graph.add_edge("ASV1", "ASV2", weight=0.4)
    paper_nodes = pd.DataFrame([
        {"ASV_ID": "ASV1", "node_type": "ASV"},
        {"ASV_ID": "ASV2", "node_type": "ASV"},
    ])
    paper_edges = pd.DataFrame([
        {"source": "ASV1", "target": "ASV2", "edge_type": "asv_association"},
    ])
    validation = pd.DataFrame([{
        "ASV_ID": "ASV1", "genome_id": "MAG1", "mag_match_id": "MAG1",
        "analysis_eligible_pair": True, "link_pident": 100.0, "link_qcov": 100.0,
    }])
    agreement = pd.DataFrame([{
        "ASV_ID": "ASV1", "genome_id": "MAG1", "mag_match_id": "MAG1",
        "shared_samples": 12, "abundance_metric": "relative_abundance",
        "transform": "log1p", "spearman_rho": -0.5,
    }])

    nodes, edges, hetero, mag_pairwise, excluded = MODULE.make_heterogeneous_outputs(
        graph, paper_nodes, paper_edges, validation, agreement, {},
    )

    assert set(nodes["id"]) == {"ASV1", "ASV2", "MAG1"}
    assert bool(nodes.set_index("id").loc["ASV1", "has_accepted_mag"])
    assert not bool(nodes.set_index("id").loc["ASV2", "has_accepted_mag"])
    link = edges.loc[edges["edge_type"] == "sequence_match"].iloc[0]
    assert link["spearman_rho"] == -0.5
    assert hetero.nodes["MAG1"]["node_type"] == "MAG"
    assert mag_pairwise.empty
    assert excluded.empty


def test_mag_abundance_neighbors_subcluster_candidates_for_same_asv():
    graph = nx.Graph()
    graph.add_edge("ASV1", "ASV2", weight=0.4)
    paper_nodes = pd.DataFrame([
        {"ASV_ID": "ASV1", "node_type": "ASV"},
        {"ASV_ID": "ASV2", "node_type": "ASV"},
    ])
    paper_edges = pd.DataFrame([
        {"source": "ASV1", "target": "ASV2", "edge_type": "asv_association"},
    ])
    validation = pd.DataFrame([
        {"ASV_ID": "ASV1", "genome_id": "MAG1", "mag_match_id": "MAG1", "analysis_eligible_pair": True},
        {"ASV_ID": "ASV1", "genome_id": "MAG2", "mag_match_id": "MAG2", "analysis_eligible_pair": True},
        {"ASV_ID": "ASV1", "genome_id": "MAG3", "mag_match_id": "MAG3", "analysis_eligible_pair": True},
    ])
    agreement = pd.DataFrame([
        {"ASV_ID": "ASV1", "genome_id": mag, "mag_match_id": mag, "spearman_rho": rho}
        for mag, rho in (("MAG1", 0.8), ("MAG2", 0.7), ("MAG3", -0.2))
    ])
    samples = ["S1", "S2", "S3", "S4"]
    matrices = {"samples": samples, "mag": pd.DataFrame(
        [[1, 2, 3, 4], [1, 2, 4, 3], [4, 3, 2, 1]],
        index=["MAG1", "MAG2", "MAG3"], columns=samples,
    )}

    _nodes, edges, hetero, pairwise, _excluded = MODULE.make_heterogeneous_outputs(
        graph, paper_nodes, paper_edges, validation, agreement, matrices, mag_knn=1,
    )

    assert len(pairwise) == 3
    mag_edges = edges.loc[edges["edge_type"] == "mag_abundance_association"]
    assert not mag_edges.empty
    assert any(
        data.get("edge_type") == "mag_abundance_association"
        for _left, _right, data in hetero.edges(data=True)
    )


def test_mags_linked_only_to_filtered_asvs_are_audited_not_orphaned():
    graph = nx.Graph()
    graph.add_edge("ASV1", "ASV2", weight=0.4)
    paper_nodes = pd.DataFrame([
        {"ASV_ID": "ASV1", "node_type": "ASV"},
        {"ASV_ID": "ASV2", "node_type": "ASV"},
    ])
    paper_edges = pd.DataFrame([
        {"source": "ASV1", "target": "ASV2", "edge_type": "asv_association"},
    ])
    validation = pd.DataFrame([
        {"ASV_ID": "ASV1", "genome_id": "MAG1", "mag_match_id": "MAG1", "analysis_eligible_pair": True},
        {"ASV_ID": "ASV_FILTERED", "genome_id": "MAG_ORPHAN", "mag_match_id": "MAG_ORPHAN", "analysis_eligible_pair": True},
    ])
    agreement = pd.DataFrame([
        {"ASV_ID": "ASV1", "genome_id": "MAG1", "mag_match_id": "MAG1", "spearman_rho": 0.8},
    ])

    nodes, _edges, hetero, _pairwise, excluded = MODULE.make_heterogeneous_outputs(
        graph, paper_nodes, paper_edges, validation, agreement, {},
    )

    assert "MAG_ORPHAN" not in set(nodes["id"])
    assert "MAG_ORPHAN" not in hetero
    assert excluded["genome_id"].tolist() == ["MAG_ORPHAN"]
    assert excluded["heterogeneous_exclusion_reason"].iloc[0] == "ASV_not_present_in_SPIEC_EASI_network"


def test_median_ratio_normalization_does_not_use_column_proportions():
    counts = pd.DataFrame(
        [[10, 20, 40], [5, 10, 20], [0, 4, 8]],
        index=["MAG1", "MAG2", "MAG3"],
        columns=["S1", "S2", "S3"],
    )
    normalized, factors = MODULE.median_ratio_normalize(counts)
    assert np.all(np.isfinite(factors))
    assert np.all(factors > 0)
    assert not np.allclose(normalized.sum(axis=0).to_numpy(), 1.0)
    standardized = MODULE.standardize_profiles(np.log1p(normalized))
    assert np.allclose(standardized.mean(axis=1).to_numpy(), 0.0, atol=1e-10)


def test_seqkit_paired_library_fpm_normalization(tmp_path):
    seqkit_path = tmp_path / "seqkit.tsv"
    pd.DataFrame({
        "file": [
            "fastq_renamed_metag/SI034_100m/SI034_100m_pe.1.fq.gz",
            "fastq_renamed_metag/SI034_100m/SI034_100m_pe.2.fq.gz",
            "fastq_renamed_metag/SI034_10m/SI034_10m_R1_001.fastq.gz",
            "fastq_renamed_metag/SI034_10m/SI034_10m_R2_001.fastq.gz",
        ],
        "num_seqs": [40_000_000, 40_000_000, 20_000_000, 20_000_000],
    }).to_csv(seqkit_path, sep="\t", index=False)
    counts = pd.DataFrame(
        {"SI034_100m": [4_000], "SI034_10m": [1_000]}, index=["MAG1"]
    )

    normalized, denominators, method, audit = (
        MODULE.normalize_mag_recruitment(
            counts,
            method="input_fragment_fpm",
            seqkit_path=seqkit_path,
        )
    )

    assert method == "input_fragment_fpm"
    assert denominators["SI034_100m"] == 40_000_000
    assert normalized.loc["MAG1", "SI034_100m"] == 100.0
    assert normalized.loc["MAG1", "SI034_10m"] == 50.0
    assert len(audit) == 4


def test_seqkit_missing_mate_is_rejected():
    seqkit = pd.DataFrame({
        "file": ["SI034_100m_pe.1.fq.gz"],
        "num_seqs": [100],
    })
    with pytest.raises(ValueError, match="lack R1 or R2"):
        MODULE.parse_seqkit_paired_libraries(seqkit)


def test_identical_duplicate_seqkit_references_are_collapsed():
    seqkit = pd.DataFrame({
        "file": [
            "old/SI060_200m_pe.1.fq.gz",
            "old/SI060_200m_pe.2.fq.gz",
            "current/SI060_200m_pe.1.fq.gz",
            "current/SI060_200m_pe.2.fq.gz",
        ],
        "num_seqs": [100, 100, 100, 100],
    })
    fragments, audit = MODULE.parse_seqkit_paired_libraries(seqkit)

    assert fragments.to_dict() == {"SI060_200m": 100.0}
    assert len(audit) == 2
    assert set(audit["duplicate_source_rows"]) == {2}
    assert set(audit["pair_validation"]) == {
        "matched_duplicate_reference_collapsed"
    }


def test_auto_mag_normalization_requires_seqkit():
    counts = pd.DataFrame({"S1": [10]}, index=["MAG1"])
    with pytest.raises(ValueError, match="requires a paired-end SeqKit"):
        MODULE.normalize_mag_recruitment(counts, method="auto")


def test_provided_fpkm_is_preserved_without_renormalization():
    values = pd.DataFrame(
        {"S1": [12.5, 0.0], "S2": [3.25, 8.0]}, index=["MAG1", "MAG2"]
    )
    normalized, factors, method, audit = MODULE.normalize_mag_recruitment(
        values, method="provided_fpkm"
    )
    pd.testing.assert_frame_equal(normalized, values)
    assert method == "provided_fpkm"
    assert factors.empty
    assert audit.empty


def test_provided_tpm_rejects_negative_values():
    values = pd.DataFrame({"S1": [10.0, -1.0]}, index=["MAG1", "MAG2"])
    with pytest.raises(ValueError, match="negative"):
        MODULE.normalize_mag_recruitment(values, method="provided_tpm")


def test_ambiguity_outputs_plot_every_eligible_multimag_asv(tmp_path):
    for subdir in ("mapping", "abundance", "functional", "qc"):
        (tmp_path / subdir).mkdir()
    pairing = pd.DataFrame([
        {"ASV_ID": "ASV1", "genome_id": "MAG1", "mag_match_id": "MAG1", "mapping_class": "ambiguous_match", "passes_paper_thresholds": True, "taxonomy_validation_status": "taxonomy_accepted", "analysis_eligible_pair": True},
        {"ASV_ID": "ASV1", "genome_id": "MAG2", "mag_match_id": "MAG2", "mapping_class": "ambiguous_match", "passes_paper_thresholds": True, "taxonomy_validation_status": "taxonomy_accepted", "analysis_eligible_pair": True},
        {"ASV_ID": "ASV2", "genome_id": "MAG3", "mag_match_id": "MAG3", "mapping_class": "unique_match", "passes_paper_thresholds": True, "taxonomy_validation_status": "taxonomy_accepted", "analysis_eligible_pair": True},
    ])
    samples = ["S1", "S2", "S3"]
    matrices = {
        "samples": samples,
        "asv": pd.DataFrame([[0.1, 0.2, 0.4]], index=["ASV1"], columns=samples),
        "mag": pd.DataFrame([[0.1, 0.2, 0.4], [0.4, 0.1, 0.0]], index=["MAG1", "MAG2"], columns=samples),
    }
    functions = pd.DataFrame([
        {"mag_match_id": "MAG1", "module_id": "M00529", "fraction_covered": 1.0},
        {"mag_match_id": "MAG2", "module_id": "M00529", "fraction_covered": 0.0},
    ])
    agreement = pd.DataFrame([
        {"ASV_ID": "ASV1", "genome_id": "MAG1", "mag_match_id": "MAG1", "spearman_rho": 1.0},
        {"ASV_ID": "ASV1", "genome_id": "MAG2", "mag_match_id": "MAG2", "spearman_rho": -0.5},
    ])

    MODULE.ambiguity_outputs(
        pairing, agreement, matrices, functions, tmp_path, "test",
        ["M00529", "M00973", "M00595", "M00984"], ["png"],
    )

    assert (tmp_path / "mapping" / "test_ambiguous_ASV1.png").is_file()
    assert not (tmp_path / "mapping" / "test_ambiguous_ASV2.png").exists()
    assert not (tmp_path / "mapping" / "test_ambiguity_maintext.png").exists()
    ranking = pd.read_csv(tmp_path / "qc" / "test_ambiguous_asv_ranking.tsv", sep="\t")
    assert ranking["ASV_ID"].tolist() == ["ASV1"]
