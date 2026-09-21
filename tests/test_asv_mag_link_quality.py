import importlib.util
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).parents[1]
SCRIPT = ROOT / "processes" / "asv_mag_link" / "asv_mag_barrnap_linker.py"
SPEC = importlib.util.spec_from_file_location("asv_mag_barrnap_linker", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_hq_quality_gate_is_inclusive_for_completeness_and_exclusive_for_contamination(tmp_path):
    metadata = pd.DataFrame(
        {
            "mag_native_genome_id": ["PASS", "LOW_COMPLETENESS", "FIVE_CONTAM"],
            "mag_completeness": [90.0, 89.99, 95.0],
            "mag_contamination": [4.99, 1.0, 5.0],
            "mag_eligible_for_linking": [True, True, True],
        }
    )
    source = MODULE.GenomeQcSource(
        source_label="test",
        genome_qc_dir=tmp_path,
        barrnap_dir=tmp_path,
        genome_fasta_dir=tmp_path,
        allowed_genomes={"PASS", "LOW_COMPLETENESS", "FIVE_CONTAM"},
        exact_barrnap_map={},
        exact_genome_fasta_map={},
        mag_metadata=metadata,
        multi_source=False,
        id_token_index=None,
    )

    MODULE.apply_mag_quality_gates([source], 90.0, 5.0)

    assert source.allowed_genomes == {"PASS"}
    assert metadata["mag_eligible_for_linking"].tolist() == [True, False, False]
    assert metadata["mag_pass_quality_gates"].tolist() == [True, False, False]
    assert source.mag_metadata["mag_native_genome_id"].tolist() == ["PASS"]


def test_exact_barrnap_matches_cannot_bypass_quality_gate(tmp_path):
    passed = tmp_path / "pass_barrnap.fasta"
    failed = tmp_path / "failed_barrnap.fasta"
    passed.write_text(">16S_rRNA\n" + "A" * 100 + "\n")
    failed.write_text(">16S_rRNA\n" + "C" * 100 + "\n")

    source = MODULE.GenomeQcSource(
        source_label="test",
        genome_qc_dir=tmp_path,
        barrnap_dir=tmp_path,
        genome_fasta_dir=None,
        allowed_genomes={"PASS"},
        exact_barrnap_map={
            MODULE.exact_match_key(passed): "PASS",
            MODULE.exact_match_key(failed): "FAILED",
        },
        exact_genome_fasta_map={},
        mag_metadata=pd.DataFrame(),
        multi_source=False,
        id_token_index=None,
    )

    refs = MODULE.collect_rrna_fastas(source)

    assert {ref.native_genome_id for ref in refs} == {"PASS"}


def test_dual_gunc_gate_requires_both_assessments(tmp_path):
    metadata = pd.DataFrame(
        {
            "mag_native_genome_id": ["BOTH", "STANDARD_ONLY", "STRICT_ONLY"],
            "mag_gunc_assessment": [
                "likely_not_chimeric",
                "likely_not_chimeric",
                "credible_chimeric_signal",
            ],
            "mag_gunc_strict_assessment": [
                "likely_not_chimeric",
                "credible_chimeric_signal",
                "likely_not_chimeric",
            ],
            "mag_eligible_for_linking": [True, True, True],
        }
    )
    source = MODULE.GenomeQcSource(
        source_label="test",
        genome_qc_dir=tmp_path,
        barrnap_dir=tmp_path,
        genome_fasta_dir=tmp_path,
        allowed_genomes={"BOTH", "STANDARD_ONLY", "STRICT_ONLY"},
        exact_barrnap_map={},
        exact_genome_fasta_map={},
        mag_metadata=metadata,
        multi_source=False,
        id_token_index=None,
    )

    MODULE.apply_mag_quality_gates(
        [source],
        min_completeness=None,
        max_contamination=None,
        gunc_assessment_value="likely_not_chimeric",
    )

    assert source.allowed_genomes == {"BOTH"}
    assert source.mag_metadata["mag_native_genome_id"].tolist() == ["BOTH"]


def test_species_assignment_gate_rejects_missing_or_unclassified_species(tmp_path):
    metadata = pd.DataFrame(
        {
            "mag_native_genome_id": ["NAMED", "EMPTY", "UNCLASSIFIED"],
            "mag_species": ["Example species", "", "unclassified"],
            "mag_eligible_for_linking": [True, True, True],
        }
    )
    source = MODULE.GenomeQcSource(
        source_label="test",
        genome_qc_dir=tmp_path,
        barrnap_dir=tmp_path,
        genome_fasta_dir=tmp_path,
        allowed_genomes={"NAMED", "EMPTY", "UNCLASSIFIED"},
        exact_barrnap_map={},
        exact_genome_fasta_map={},
        mag_metadata=metadata,
        multi_source=False,
        id_token_index=None,
    )

    MODULE.apply_mag_quality_gates(
        [source],
        min_completeness=None,
        max_contamination=None,
        require_species_assignment=True,
    )

    assert source.allowed_genomes == {"NAMED"}
    assert source.mag_metadata["mag_native_genome_id"].tolist() == ["NAMED"]


def test_asv_taxonomy_confidence_gate_is_inclusive_and_audited(tmp_path):
    taxonomy = tmp_path / "taxonomy.tsv"
    pd.DataFrame(
        {
            "Feature ID": ["ASV1;size=10", "ASV2;size=20", "ASV3;size=30"],
            "Taxon": ["g__One", "g__Two", "g__Three"],
            "Consensus": [0.90, 0.899, 1.0],
        }
    ).to_csv(taxonomy, sep="\t", index=False)

    retained, retained_canonical = MODULE.filter_asvs_by_taxonomy_confidence(
        ["ASV1;size=10", "ASV2;size=20", "ASV3;size=30"],
        taxonomy,
        0.90,
        tmp_path,
    )

    assert retained == ["ASV1;size=10", "ASV3;size=30"]
    assert retained_canonical == {"ASV1", "ASV3"}
    audit = pd.read_csv(tmp_path / "tables" / "asv_taxonomy_confidence_filter.tsv", sep="\t")
    assert audit["passes_taxonomy_confidence"].tolist() == [True, False, True]


def test_rrna_marker_gate_requires_16s_23s_and_5s(tmp_path):
    metadata = pd.DataFrame(
        {
            "mag_native_genome_id": ["ALL", "NO_23S", "NO_5S"],
            "mag_16s_rrna": [1, 1, 1],
            "mag_23s_rrna": [1, 0, 1],
            "mag_5s_rrna": [1, 1, 0],
            "mag_eligible_for_linking": [True, True, True],
        }
    )
    source = MODULE.GenomeQcSource(
        source_label="test",
        genome_qc_dir=tmp_path,
        barrnap_dir=tmp_path,
        genome_fasta_dir=tmp_path,
        allowed_genomes={"ALL", "NO_23S", "NO_5S"},
        exact_barrnap_map={},
        exact_genome_fasta_map={},
        mag_metadata=metadata,
        multi_source=False,
        id_token_index=None,
    )

    MODULE.apply_mag_quality_gates(
        [source],
        min_completeness=None,
        max_contamination=None,
        min_rrna_marker_count=1,
    )

    assert source.allowed_genomes == {"ALL"}
    assert source.mag_metadata["mag_native_genome_id"].tolist() == ["ALL"]
