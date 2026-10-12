import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BUILD = ROOT / "processes/preprocessing_dataset/build_preprocessing_dataset.py"
VALIDATE = ROOT / "processes/preprocessing_dataset/validate_preprocessing_dataset.py"


def _inputs(tmp_path: Path) -> dict[str, Path]:
    values = {}
    for name in (
        "asv_counts",
        "asv_filtered",
        "asv_target",
        "asv_target_micro",
        "asv_target_mito",
        "asv_target_decon",
        "taxonomy",
        "fastq_stats",
        "fastp_stats",
        "filtered_stats",
        "concat_stats",
        "sample_manifest",
    ):
        path = tmp_path / f"{name}.tsv"
        path.write_text("id\tvalue\nrow1\t1\n", encoding="utf-8")
        values[name] = path
    fasta = tmp_path / "asv_sequences.fasta.gz"
    fasta.write_bytes(b"synthetic fasta payload")
    values["asv_sequences"] = fasta
    return values


def test_build_and_validate_preprocessing_dataset(tmp_path: Path) -> None:
    inputs = _inputs(tmp_path)
    dataset = tmp_path / "canonical"
    command = [
        sys.executable,
        str(BUILD),
        "--outdir",
        str(dataset),
        "--preprocessing-parameters",
        json.dumps({"filter": {"max_ee": 1.0}}),
    ]
    for name, path in inputs.items():
        command.extend([f"--{name.replace('_', '-')}", str(path)])
    subprocess.run(command, check=True)
    subprocess.run([sys.executable, str(VALIDATE), str(dataset)], check=True)

    manifest = json.loads((dataset / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["schema_version"] == "1.2"
    assert len(manifest["files"]) == 13
    assert (dataset / "sample_manifest.tsv").is_file()
    assert (dataset / "files.tsv").is_file()


def test_validation_rejects_modified_artifact(tmp_path: Path) -> None:
    inputs = _inputs(tmp_path)
    dataset = tmp_path / "canonical"
    command = [
        sys.executable,
        str(BUILD),
        "--outdir",
        str(dataset),
        "--preprocessing-parameters",
        "{}",
    ]
    for name, path in inputs.items():
        command.extend([f"--{name.replace('_', '-')}", str(path)])
    subprocess.run(command, check=True)
    (dataset / "asv_counts.tsv").write_text("changed\n", encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(VALIDATE), str(dataset)],
        text=True,
        capture_output=True,
    )
    assert result.returncode != 0
    assert "Checksum mismatch" in result.stderr


def test_control_and_primer_audits_survive_handoff(tmp_path):
    inputs = _inputs(tmp_path)
    primer = tmp_path / "primer.json"
    primer.write_text('{"family": "515F_806R"}')
    controls = tmp_path / "controls"
    controls.mkdir()
    (controls / "sample_qc.tsv").write_text("sample_id\tdecontam_class\nblank\ttechnical\n")
    cleaned = tmp_path / "cleaned.tsv"
    cleaned.write_text("ASV\tbio1\nASV1\t6000\n")
    dataset = tmp_path / "canonical"
    command = [sys.executable, str(BUILD), "--outdir", str(dataset),
               "--preprocessing-parameters", json.dumps({"primer_trimming": {"enabled": True},
                                                        "control_decontam": {"enabled": True}}),
               "--primer-audit", str(primer), "--control-audit", str(controls),
               "--control-cleaned", str(cleaned)]
    for name, path in inputs.items():
        command.extend([f"--{name.replace('_', '-')}", str(path)])
    subprocess.run(command, check=True)
    subprocess.run([sys.executable, str(VALIDATE), str(dataset)], check=True)
    assert (dataset / "control_audit/sample_qc.tsv").read_text() == (controls / "sample_qc.tsv").read_text()
    (dataset / "control_audit/sample_qc.tsv").write_text("changed")
    result = subprocess.run([sys.executable, str(VALIDATE), str(dataset)], capture_output=True, text=True)
    assert result.returncode != 0 and "Checksum mismatch" in result.stderr
