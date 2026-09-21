#!/usr/bin/env python3
"""Build and validate ASPIRE's canonical preprocessing handoff dataset."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path


SCHEMA_VERSION = "1.1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def table_shape(path: Path) -> tuple[int | None, int | None]:
    if path.suffix not in {".tsv", ".csv"}:
        return None, None
    delimiter = "\t" if path.suffix == ".tsv" else ","
    rows = 0
    columns = None
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        for row in reader:
            if columns is None:
                columns = len(row)
            else:
                rows += 1
    return rows, columns


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--preprocessing-parameters", required=True)
    for name in (
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
    ):
        parser.add_argument(f"--{name.replace('_', '-')}", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=False)

    names = {
        "asv_counts": "asv_counts.tsv",
        "asv_filtered": "asv_filtered.tsv",
        "asv_sequences": "asv_sequences.filtered.fasta.gz",
        "asv_target": "asv_target.tsv",
        "asv_target_micro": "asv_target.micro.tsv",
        "asv_target_mito": "asv_target.mito.tsv",
        "asv_target_decon": "asv_target.decon.tsv",
        "taxonomy": "taxonomy.tsv",
        "fastq_stats": "fastq_stats.tsv",
        "fastp_stats": "fastp_fastqs.tsv",
        "filtered_stats": "filtered_fastas.tsv",
        "concat_stats": "concat_fastas.tsv",
        "sample_manifest": "sample_manifest.tsv",
    }
    records = []
    for role, destination_name in names.items():
        source = Path(getattr(args, role))
        if not source.is_file():
            raise FileNotFoundError(f"Required preprocessing artifact is missing: {source}")
        destination = outdir / destination_name
        shutil.copy2(source, destination)
        rows, columns = table_shape(destination)
        records.append(
            {
                "role": role,
                "file": destination_name,
                "bytes": destination.stat().st_size,
                "sha256": sha256(destination),
                "data_rows": rows,
                "columns": columns,
            }
        )

    parameters = json.loads(args.preprocessing_parameters)
    with (outdir / "preprocessing_parameters.json").open("w", encoding="utf-8") as handle:
        json.dump(parameters, handle, indent=2, sort_keys=True)
        handle.write("\n")

    with (outdir / "files.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("role", "file", "bytes", "sha256", "data_rows", "columns"),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(records)

    manifest = {
        "schema": "ASPIRE canonical preprocessing dataset",
        "schema_version": SCHEMA_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "files": records,
        "preprocessing_parameters_file": "preprocessing_parameters.json",
        "preprocessing_parameters_sha256": sha256(outdir / "preprocessing_parameters.json"),
    }
    with (outdir / "manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")


if __name__ == "__main__":
    main()
