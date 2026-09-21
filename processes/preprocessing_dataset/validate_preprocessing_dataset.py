#!/usr/bin/env python3
"""Validate an ASPIRE canonical preprocessing dataset before analysis."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


SCHEMA_VERSION = "1.1"
REQUIRED_ROLES = {
    "asv_counts",
    "asv_filtered",
    "asv_sequences",
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
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("dataset")
    args = parser.parse_args()
    dataset = Path(args.dataset)
    manifest_path = dataset / "manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Missing preprocessing manifest: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(f"Unsupported preprocessing schema: {manifest.get('schema_version')!r}")
    records = manifest.get("files", [])
    roles = {record.get("role") for record in records}
    missing_roles = sorted(REQUIRED_ROLES - roles)
    if missing_roles:
        raise ValueError(f"Preprocessing manifest lacks required roles: {', '.join(missing_roles)}")
    for record in records:
        path = dataset / record["file"]
        if not path.is_file():
            raise FileNotFoundError(f"Manifest artifact is missing: {path}")
        observed = sha256(path)
        if observed != record["sha256"]:
            raise ValueError(f"Checksum mismatch for {path}: {observed} != {record['sha256']}")
    parameter_path = dataset / manifest["preprocessing_parameters_file"]
    if sha256(parameter_path) != manifest["preprocessing_parameters_sha256"]:
        raise ValueError(f"Checksum mismatch for {parameter_path}")
    print(f"Validated {len(records)} preprocessing artifacts in {dataset}")


if __name__ == "__main__":
    main()
