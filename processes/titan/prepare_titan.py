#!/usr/bin/env python3
"""Prepare SPIEC-EASI-cohort relative abundances for variable-wise TITAN2."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


def read_table(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix.lower() in {".tsv", ".tab", ".txt"} else ","
    return pd.read_csv(path, sep=sep, low_memory=False)


def parse_list(value: str | None) -> list[str]:
    if not value:
        return []
    return [item.strip() for item in re.split(r"[,|]", value) if item.strip()]


def truthy(values: pd.Series) -> pd.Series:
    return values.astype(str).str.strip().str.lower().isin({"true", "t", "1", "yes", "y"})


def safe_slug(value: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip()).strip("._-")
    return slug or "measurement"


def load_counts(path: Path, transpose: bool, asv_id_col: str) -> pd.DataFrame:
    table = read_table(path)
    if table.empty or table.shape[1] < 2:
        raise ValueError("ASV count table must contain identifiers and at least one data column")
    identifier = asv_id_col if asv_id_col in table.columns else table.columns[0]
    ids = table[identifier].astype(str).str.strip().str.replace(r";.*$", "", regex=True)
    values = table.drop(columns=[identifier]).apply(pd.to_numeric, errors="coerce").fillna(0.0)
    if (values < 0).any().any():
        raise ValueError("ASV count table contains negative values")
    if transpose:
        values.index = ids
        if values.index.duplicated().any():
            duplicates = values.index[values.index.duplicated()].unique().tolist()
            raise ValueError(f"ASV identifiers are not unique: {duplicates[:10]}")
        counts = values.T
    else:
        values.index = ids
        counts = values
    counts.index = counts.index.astype(str).str.strip()
    counts.columns = counts.columns.astype(str).str.strip().str.replace(r";.*$", "", regex=True)
    if counts.index.duplicated().any():
        duplicates = counts.index[counts.index.duplicated()].unique().tolist()
        raise ValueError(f"Sample identifiers are not unique in ASV counts: {duplicates[:10]}")
    if counts.columns.duplicated().any():
        duplicates = counts.columns[counts.columns.duplicated()].unique().tolist()
        raise ValueError(f"ASV identifiers are not unique after normalization: {duplicates[:10]}")
    return counts


def load_taxonomy(path: Path | None) -> pd.DataFrame:
    if path is None:
        return pd.DataFrame(columns=["ASV_ID", "taxonomy"])
    table = read_table(path)
    id_col = next(
        (column for column in table.columns if column.strip().casefold() in
         {"feature id", "feature_id", "asv_id", "asv", "id"}),
        table.columns[0],
    )
    tax_col = next(
        (column for column in table.columns if column.strip().casefold() in {"taxon", "taxonomy"}),
        None,
    )
    out = pd.DataFrame({"ASV_ID": table[id_col].astype(str).str.strip().str.replace(r";.*$", "", regex=True)})
    out["taxonomy"] = table[tax_col].fillna("").astype(str) if tax_col else ""
    return out.drop_duplicates("ASV_ID", keep="first")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--counts", type=Path, required=True)
    parser.add_argument("--filter-audit", type=Path, required=True)
    parser.add_argument("--measurement-matrix", type=Path, required=True)
    parser.add_argument("--taxonomy", type=Path)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--sample-col", default="sampleID")
    parser.add_argument("--asv-id-col", default="ASV_ID")
    parser.add_argument("--variables", default="")
    parser.add_argument("--transpose", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--min-split", type=int, default=5)
    parser.add_argument("--minimum-samples", type=int, default=10)
    parser.add_argument("--minimum-occurrence", type=int, default=3)
    parser.add_argument("--minimum-prevalence", type=float, default=0.0)
    parser.add_argument("--minimum-mean-relative-abundance", type=float, default=0.0)
    args = parser.parse_args()

    if args.min_split < 3:
        raise ValueError("--min-split must be at least 3")
    if args.minimum_samples < 10:
        raise ValueError("--minimum-samples cannot be below TITAN2's hard minimum of 10")
    if args.minimum_occurrence < 3:
        raise ValueError("--minimum-occurrence cannot be below TITAN2's hard minimum of 3")
    if not 0 <= args.minimum_prevalence <= 1:
        raise ValueError("--minimum-prevalence must be between 0 and 1")
    if not 0 <= args.minimum_mean_relative_abundance <= 1:
        raise ValueError("--minimum-mean-relative-abundance must be between 0 and 1")

    args.outdir.mkdir(parents=True, exist_ok=True)
    variable_root = args.outdir / "variables"
    variable_root.mkdir(exist_ok=True)

    counts = load_counts(args.counts, args.transpose, args.asv_id_col)
    filter_audit = read_table(args.filter_audit)
    required_audit = {"ASV_ID", "retained_final"}
    if not required_audit.issubset(filter_audit.columns):
        raise ValueError(
            "SPIEC-EASI filtering audit lacks: "
            + ", ".join(sorted(required_audit.difference(filter_audit.columns)))
        )
    retained = set(
        filter_audit.loc[truthy(filter_audit["retained_final"]), "ASV_ID"]
        .astype(str).str.strip().str.replace(r";.*$", "", regex=True)
    )
    retained_order = [asv for asv in counts.columns if asv in retained]
    missing_retained = sorted(retained.difference(counts.columns))
    if missing_retained:
        raise ValueError(
            "SPIEC-EASI retained_final ASVs are absent from the supplied count table: "
            + ", ".join(missing_retained[:10])
        )
    counts = counts.loc[:, retained_order]
    if counts.empty:
        raise ValueError("No retained_final SPIEC-EASI ASVs remain")
    filter_audit = filter_audit.copy()
    filter_audit["ASV_ID"] = (
        filter_audit["ASV_ID"].astype(str).str.strip().str.replace(r";.*$", "", regex=True)
    )
    filter_audit.to_csv(args.outdir / "titan_spieceasi_filter_audit.tsv", sep="\t", index=False)

    measurements = read_table(args.measurement_matrix)
    measurement_id_col = args.sample_col if args.sample_col in measurements.columns else measurements.columns[0]
    measurements[measurement_id_col] = measurements[measurement_id_col].astype(str).str.strip()
    if measurements[measurement_id_col].duplicated().any():
        duplicates = measurements.loc[measurements[measurement_id_col].duplicated(), measurement_id_col].unique().tolist()
        raise ValueError(f"Measurement sample identifiers are not unique: {duplicates[:10]}")
    measurements = measurements.set_index(measurement_id_col)
    requested = parse_list(args.variables)
    variables = requested or list(measurements.columns)
    if len(set(map(safe_slug, variables))) != len(variables):
        raise ValueError("Environmental-variable names do not produce unique filesystem-safe names")

    taxonomy = load_taxonomy(args.taxonomy).set_index("ASV_ID")
    common_samples = [sample for sample in counts.index if sample in measurements.index]
    if not common_samples:
        raise ValueError("No explicit sample-ID overlap between ASV counts and measurement matrix")
    all_sample_ids = list(dict.fromkeys([*counts.index.tolist(), *measurements.index.tolist()]))
    sample_matching = pd.DataFrame({args.sample_col: all_sample_ids})
    sample_matching["present_in_asv_counts"] = sample_matching[args.sample_col].isin(counts.index)
    sample_matching["present_in_measurement_matrix"] = sample_matching[args.sample_col].isin(measurements.index)
    sample_matching["matched_by_explicit_identifier"] = (
        sample_matching["present_in_asv_counts"] & sample_matching["present_in_measurement_matrix"]
    )
    sample_matching.to_csv(args.outdir / "titan_sample_matching_audit.tsv", sep="\t", index=False)

    manifest_rows: list[dict[str, object]] = []
    audit_rows: list[dict[str, object]] = []
    for variable in variables:
        slug = safe_slug(variable)
        variable_dir = variable_root / slug
        variable_dir.mkdir(exist_ok=True)
        status = "ready"
        reason = ""
        negative_values_treated_missing = 0
        missing_values = len(common_samples)
        if variable not in measurements.columns:
            status, reason = "skipped", "measurement_missing"
            measured_samples: list[str] = []
            env = pd.Series(dtype=float)
        else:
            env = pd.to_numeric(measurements.loc[common_samples, variable], errors="coerce")
            negative_values_treated_missing = int(env.lt(0).fillna(False).sum())
            env = env.where(env.ge(0))
            missing_values = int(env.isna().sum())
            measured_samples = env.index[env.notna()].tolist()
            env = env.loc[measured_samples]
            if len(measured_samples) < args.minimum_samples:
                status, reason = "skipped", f"fewer_than_{args.minimum_samples}_measured_samples"
            elif env.nunique(dropna=True) < 2:
                status, reason = "skipped", "invariant_environmental_variable"
            elif len(measured_samples) < 2 * args.min_split:
                status, reason = "skipped", "insufficient_samples_for_min_split"

        tested_asvs: list[str] = []
        if measured_samples:
            subset = counts.loc[measured_samples].copy()
            library_totals = subset.sum(axis=1)
            valid_library = library_totals.gt(0)
            if not valid_library.all():
                removed_samples = subset.index[~valid_library].tolist()
                for sample in removed_samples:
                    audit_rows.append({
                        "environmental_variable": variable, "record_type": "sample",
                        "identifier": sample, "retained_for_titan": False,
                        "exclusion_reason": "zero_retained_cohort_library",
                        "prevalence": np.nan, "mean_relative_abundance": np.nan,
                    })
                subset = subset.loc[valid_library]
                env = env.loc[subset.index]
            if status == "ready" and len(env) < args.minimum_samples:
                status, reason = "skipped", "insufficient_samples_after_zero_library_screen"
            elif status == "ready" and len(env) < 2 * args.min_split:
                status, reason = "skipped", "insufficient_samples_for_min_split_after_zero_library_screen"
            relative = subset.div(subset.sum(axis=1), axis=0)
            occurrence = relative.gt(0).sum(axis=0)
            prevalence = relative.gt(0).mean(axis=0)
            mean_abundance = relative.mean(axis=0)
            variance = relative.var(axis=0)
            keep = (
                occurrence.ge(args.minimum_occurrence)
                & prevalence.ge(args.minimum_prevalence)
                & mean_abundance.ge(args.minimum_mean_relative_abundance)
                & variance.gt(0)
                & mean_abundance.gt(0)
            )
            tested_asvs = relative.columns[keep].tolist()
            for asv in relative.columns:
                reasons = []
                if occurrence[asv] < args.minimum_occurrence:
                    reasons.append("occurrence_below_titan_minimum")
                if prevalence[asv] < args.minimum_prevalence:
                    reasons.append("prevalence_below_configured_minimum")
                if mean_abundance[asv] < args.minimum_mean_relative_abundance:
                    reasons.append("mean_relative_abundance_below_configured_minimum")
                if not variance[asv] > 0:
                    reasons.append("invariant_within_measurement_cohort")
                if not mean_abundance[asv] > 0:
                    reasons.append("all_zero_within_measurement_cohort")
                audit_rows.append({
                    "environmental_variable": variable, "record_type": "ASV",
                    "identifier": asv, "retained_for_titan": bool(keep[asv]),
                    "exclusion_reason": ";".join(reasons) if reasons else "retained",
                    "prevalence": float(prevalence[asv]),
                    "mean_relative_abundance": float(mean_abundance[asv]),
                })
            if status == "ready" and not tested_asvs:
                status, reason = "skipped", "no_asvs_pass_titan_requirements"
            if status == "ready":
                pd.DataFrame({args.sample_col: env.index, "environmental_value": env.values}).to_csv(
                    variable_dir / "environment.tsv", sep="\t", index=False
                )
                taxa = relative.loc[env.index, tested_asvs]
                taxa.index.name = args.sample_col
                taxa.reset_index().to_csv(variable_dir / "taxa_relative_abundance.tsv", sep="\t", index=False)
                metadata_rows = []
                audit_index = {
                    row["identifier"]: row for row in audit_rows
                    if row["environmental_variable"] == variable and row["record_type"] == "ASV"
                }
                for asv in tested_asvs:
                    metadata_rows.append({
                        "ASV_ID": asv,
                        "taxonomy": taxonomy.at[asv, "taxonomy"] if asv in taxonomy.index else "",
                        "prevalence": audit_index[asv]["prevalence"],
                        "mean_relative_abundance": audit_index[asv]["mean_relative_abundance"],
                    })
                pd.DataFrame(metadata_rows).to_csv(variable_dir / "taxon_metadata.tsv", sep="\t", index=False)

        manifest_rows.append({
            "environmental_variable": variable,
            "variable_slug": slug,
            "status": status,
            "status_reason": reason,
            "number_of_samples": int(len(env)) if measured_samples else 0,
            "number_of_unique_values": int(env.nunique()) if measured_samples else 0,
            "sample_count_warning": "fewer_than_20_samples" if 0 < len(env) < 20 else "",
            "number_of_missing_or_invalid_values": missing_values,
            "number_of_negative_values_treated_missing": negative_values_treated_missing,
            "number_of_ASVs_tested": int(len(tested_asvs)),
            "minimum_environmental_value": float(env.min()) if measured_samples else np.nan,
            "maximum_environmental_value": float(env.max()) if measured_samples else np.nan,
        })

    pd.DataFrame(manifest_rows).to_csv(args.outdir / "titan_input_manifest.tsv", sep="\t", index=False)
    pd.DataFrame(audit_rows, columns=[
        "environmental_variable", "record_type", "identifier", "retained_for_titan",
        "exclusion_reason", "prevalence", "mean_relative_abundance",
    ]).to_csv(args.outdir / "titan_filtered_asv_audit.tsv", sep="\t", index=False)
    pd.DataFrame({"ASV_ID": retained_order}).to_csv(
        args.outdir / "titan_spieceasi_cohort.tsv", sep="\t", index=False
    )
    config = {
        key: str(value) if isinstance(value, Path) else value
        for key, value in vars(args).items()
    }
    config.update({
        "counts": str(args.counts), "filter_audit": str(args.filter_audit),
        "measurement_matrix": str(args.measurement_matrix),
        "taxonomy": str(args.taxonomy) if args.taxonomy else None,
        "variables_resolved": variables, "overlapping_samples": len(common_samples),
        "asv_count_samples": len(counts.index),
        "measurement_matrix_samples": len(measurements.index),
        "retained_final_asvs": len(retained_order),
    })
    (args.outdir / "titan_prepare_config.json").write_text(json.dumps(config, indent=2))


if __name__ == "__main__":
    main()
