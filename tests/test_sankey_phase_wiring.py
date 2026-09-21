from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_sankey_reads_explicit_staged_inputs() -> None:
    module = (
        ROOT / "workflow/modules/metadata_ecology.nf"
    ).read_text(encoding="utf-8")
    for variable in (
        "fastq_stats",
        "filtered_stats",
        "asv_counts",
        "asv_decon_counts",
        "asv_micro_counts",
    ):
        assert module.count(f'"\\${{PWD}}/${{{variable}}}"') == 2
    assert 'stats/"${fastq_stats}"' not in module
    assert 'ASVs/"${asv_micro_counts}"' not in module
