#!/usr/bin/env python3
"""Collect variable-isolated TITAN2 outputs into stable run-wide tables."""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import pandas as pd


SCHEMAS = {
    "titan_taxon_results.tsv": [
        "environmental_variable", "ASV_ID", "taxonomy", "change_point",
        "response_direction", "indicator_score", "z_score", "purity",
        "reliability", "bootstrap_cp_q05", "bootstrap_cp_q10",
        "bootstrap_cp_q50", "bootstrap_cp_q90", "bootstrap_cp_q95",
        "number_of_samples_used", "prevalence", "mean_relative_abundance",
        "occurrence_frequency", "response_group_code", "indval_change_point",
        "zscore_change_point", "observed_indval_probability",
        "median_bootstrap_z_score", "titan_filter_code", "passes_purity",
        "passes_reliability", "passes_purity_and_reliability",
    ],
    "titan_community_thresholds.tsv": [
        "environmental_variable", "z_minus_threshold", "z_plus_threshold",
        "z_minus_bootstrap_q05", "z_minus_bootstrap_q10", "z_minus_bootstrap_q50",
        "z_minus_bootstrap_q90", "z_minus_bootstrap_q95", "z_plus_bootstrap_q05",
        "z_plus_bootstrap_q10", "z_plus_bootstrap_q50", "z_plus_bootstrap_q90",
        "z_plus_bootstrap_q95", "filtered_z_minus_threshold",
        "filtered_z_plus_threshold", "filtered_z_minus_bootstrap_q05",
        "filtered_z_minus_bootstrap_q95", "filtered_z_plus_bootstrap_q05",
        "filtered_z_plus_bootstrap_q95", "number_of_samples",
        "number_of_ASVs_tested", "number_of_pure_reliable_z_minus",
        "number_of_pure_reliable_z_plus",
    ],
    "titan_community_response_curve.tsv": [
        "environmental_variable", "environmental_value", "sum_z_minus",
        "sum_z_plus", "filtered_sum_z_minus", "filtered_sum_z_plus",
    ],
    "titan_community_bootstrap_thresholds.tsv": [
        "environmental_variable", "bootstrap_replicate", "z_minus_threshold",
        "z_plus_threshold", "filtered_z_minus_threshold", "filtered_z_plus_threshold",
    ],
}


def collect(raw_dir: Path, filename: str, columns: list[str]) -> pd.DataFrame:
    frames = []
    for path in sorted((raw_dir / "variables").glob(f"*/{filename}")):
        frame = pd.read_csv(path, sep="\t", low_memory=False)
        missing = [column for column in columns if column not in frame.columns]
        if missing:
            raise ValueError(f"{path} lacks stable-schema columns: {missing}")
        frames.append(frame.loc[:, columns])
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=columns)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepared-dir", type=Path, required=True)
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument("--installation-manifest", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()

    tables = args.outdir / "tables"
    audit = args.outdir / "audit"
    raw = args.outdir / "raw"
    for directory in (tables, audit, raw):
        directory.mkdir(parents=True, exist_ok=True)

    for filename, columns in SCHEMAS.items():
        collect(args.raw_dir, filename, columns).to_csv(
            tables / filename, sep="\t", index=False
        )

    status = pd.read_csv(args.raw_dir / "titan_variable_status.tsv", sep="\t")
    status.to_csv(tables / "titan_variable_status.tsv", sep="\t", index=False)
    shutil.copy2(args.prepared_dir / "titan_input_manifest.tsv", audit / "titan_input_manifest.tsv")
    shutil.copy2(args.prepared_dir / "titan_filtered_asv_audit.tsv", audit / "titan_filtered_asv_audit.tsv")
    shutil.copy2(args.prepared_dir / "titan_spieceasi_cohort.tsv", audit / "titan_spieceasi_cohort.tsv")
    shutil.copy2(args.prepared_dir / "titan_spieceasi_filter_audit.tsv", audit / "titan_spieceasi_filter_audit.tsv")
    shutil.copy2(args.prepared_dir / "titan_sample_matching_audit.tsv", audit / "titan_sample_matching_audit.tsv")
    shutil.copy2(args.prepared_dir / "titan_prepare_config.json", audit / "titan_prepare_config.json")
    shutil.copy2(args.raw_dir / "titan_analysis_config.json", audit / "titan_analysis_config.json")
    shutil.copy2(args.raw_dir / "titan_session_info.txt", audit / "titan_session_info.txt")
    shutil.copy2(args.installation_manifest, audit / "titan2_installation.tsv")
    for error_path in sorted((args.raw_dir / "variables").glob("*/titan_error.txt")):
        shutil.copy2(error_path, audit / f"{error_path.parent.name}_titan_error.txt")
    for result_path in sorted((args.raw_dir / "variables").glob("*/titan_result.rds")):
        destination = raw / result_path.parent.name
        destination.mkdir(parents=True, exist_ok=True)
        shutil.copy2(result_path, destination / result_path.name)

    summary = {
        "variables_requested": int(len(status)),
        "variables_completed": int(status["status"].eq("completed").sum()),
        "variables_skipped": int(status["status"].eq("skipped").sum()),
        "variables_failed": int(status["status"].eq("failed").sum()),
        "taxon_result_rows": int(len(pd.read_csv(tables / "titan_taxon_results.tsv", sep="\t"))),
    }
    (audit / "titan_collection_summary.json").write_text(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
