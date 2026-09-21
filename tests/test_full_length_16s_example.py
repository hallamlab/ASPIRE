from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]


def test_full_length_16s_example_is_single_end_and_retains_v1_v9():
    config = yaml.safe_load((ROOT / "examples" / "full_length_16s.local.yml").read_text())
    assert config["resources"]["single_end"] is True
    assert config["fastp"]["trim_front_r1"] == 0
    assert config["fastp"]["trim_tail_r1"] == 0
    assert config["filter"]["min_len"] >= 1000
    assert config["filter"]["max_len"] >= 1600
    assert config["sina"]["regions"] == ["V1-V9"]
    assert config["sina"]["trim_to"] == "V1-V9"
    assert config["mito"]["enabled"] is False
    assert config["filter_counts"]["enabled"] is False
    assert config["metadata_plots"]["enabled"] is False


def test_full_length_manifest_has_blank_r2_column():
    lines = (ROOT / "examples" / "full_length_16s_manifest.template.tsv").read_text().splitlines()
    assert lines[0].split("\t") == ["sample_id", "fastq_r1", "fastq_r2"]
    assert all(line.endswith("\t") for line in lines[1:] if line.strip())
