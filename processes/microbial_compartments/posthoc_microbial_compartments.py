#!/usr/bin/env python3
"""Summarize microbial compartments and compare them with external labels post hoc."""

from __future__ import annotations

import argparse
import json
import platform
import re
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import sklearn
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


def read_table(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix.lower() in {".tsv", ".tab", ".txt"} else ","
    return pd.read_csv(path, sep=sep, low_memory=False)


def parse_list(value: str) -> list[str]:
    return [item.strip() for item in re.split(r"[,|]", value or "") if item.strip()]


def taxonomy_map(path: Path | None) -> dict[str, str]:
    if path is None or not path.exists():
        return {}
    table = read_table(path)
    id_col = next(
        (column for column in table.columns if column.strip().casefold() in {"feature id", "feature_id", "asv_id", "asv", "id"}),
        table.columns[0],
    )
    tax_col = next((column for column in table.columns if column.strip().casefold() in {"taxon", "taxonomy"}), None)
    if tax_col is None:
        return {}
    ids = table[id_col].astype(str).str.strip().str.replace(r";.*$", "", regex=True)
    return dict(zip(ids, table[tax_col].fillna("").astype(str)))


def copy_file(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepared-dir", type=Path, required=True)
    parser.add_argument("--inference-dir", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--taxonomy", type=Path)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--sample-col", default="sampleID")
    parser.add_argument("--environmental-compartment-cols", default="")
    parser.add_argument("--depth-col", default="Depth")
    parser.add_argument("--date-col", default="date")
    parser.add_argument("--dominant-top-n", type=int, default=10)
    args = parser.parse_args()
    if args.dominant_top_n < 1:
        raise ValueError("--dominant-top-n must be positive")

    tables = args.outdir / "tables"
    audit = args.outdir / "audit"
    tables.mkdir(parents=True, exist_ok=True)
    audit.mkdir(parents=True, exist_ok=True)

    inference_table_names = [
        "microbial_compartments.tsv", "microbial_cluster_summary_base.tsv",
        "cluster_validation.tsv", "microbial_candidate_assignments.tsv",
        "microbial_best_candidate_assignments.tsv", "pam_medoids.tsv",
        "microbial_cluster_selection_decision.tsv",
        "stability_resampling_replicates.tsv",
        "stability_cluster_recovery_replicates.tsv",
        "stability_cluster_recovery_summary.tsv",
        "season_balanced_block_replicates.tsv",
        "season_balanced_block_strata_audit.tsv",
    ]
    for name in inference_table_names:
        copy_file(args.inference_dir / name, tables / name)
    for name in ["microbial_clustering_inference_config.json", "microbial_clustering_session_info.txt"]:
        copy_file(args.inference_dir / name, audit / name)
    for name in [
        "microbial_clustering_sample_audit.tsv", "microbial_clustering_asv_filter_audit.tsv",
        "zero_replacement_audit.tsv", "microbial_clustering_preparation_config.json",
        "asv_relative_abundance.tsv", "asv_zero_replaced_composition.tsv",
        "asv_clr_matrix.tsv", "aitchison_distance_matrix.tsv",
        "microbial_clustering_stability_blocks.tsv",
    ]:
        copy_file(args.prepared_dir / name, audit / name)

    primary = read_table(args.inference_dir / "microbial_compartments.tsv")
    decision = read_table(args.inference_dir / "microbial_cluster_selection_decision.tsv").iloc[0]
    supported = str(decision["support_status"]) == "supported"
    relative = read_table(args.prepared_dir / "asv_relative_abundance.tsv")
    relative_id = relative.columns[0]
    relative[relative_id] = relative[relative_id].astype(str).str.strip()
    relative = relative.set_index(relative_id)
    primary[args.sample_col] = primary[args.sample_col].astype(str).str.strip()
    labels = primary.set_index(args.sample_col)["microbial_compartment"]
    relative = relative.loc[[sample for sample in relative.index if sample in labels.index]]
    labels = labels.loc[relative.index]
    taxa = taxonomy_map(args.taxonomy)

    dominant_rows: list[dict[str, object]] = []
    summary = read_table(args.inference_dir / "microbial_cluster_summary_base.tsv")
    if supported:
        epsilon = max(np.finfo(float).eps, float(relative[relative.gt(0)].min().min()) / 2.0)
        for compartment in sorted(labels.unique()):
            inside = labels.eq(compartment)
            mean_inside = relative.loc[inside].mean(axis=0)
            mean_outside = relative.loc[~inside].mean(axis=0)
            log2_fold = np.log2((mean_inside + epsilon) / (mean_outside + epsilon))
            ordering = pd.DataFrame({
                "ASV_ID": relative.columns,
                "mean_relative_abundance_within": mean_inside.values,
                "mean_relative_abundance_outside": mean_outside.values,
                "log2_fold_vs_other_compartments": log2_fold.values,
            }).sort_values(
                ["log2_fold_vs_other_compartments", "mean_relative_abundance_within", "ASV_ID"],
                ascending=[False, False, True],
            ).head(args.dominant_top_n)
            for rank, row in enumerate(ordering.itertuples(index=False), start=1):
                dominant_rows.append({
                    "microbial_compartment": compartment, "rank": rank,
                    "ASV_ID": row.ASV_ID, "taxonomy": taxa.get(row.ASV_ID, ""),
                    "mean_relative_abundance_within": row.mean_relative_abundance_within,
                    "mean_relative_abundance_outside": row.mean_relative_abundance_outside,
                    "log2_fold_vs_other_compartments": row.log2_fold_vs_other_compartments,
                })
        dominant = pd.DataFrame(dominant_rows)
        dominant.to_csv(tables / "microbial_compartment_dominant_asvs.tsv", sep="\t", index=False)
        dominant_asv = dominant.groupby("microbial_compartment").ASV_ID.apply(lambda values: ";".join(values))
        dominant_tax = dominant.groupby("microbial_compartment").taxonomy.apply(lambda values: ";".join(values))
        summary["dominant_ASVs"] = summary.microbial_compartment.map(dominant_asv).fillna("")
        summary["dominant_taxa"] = summary.microbial_compartment.map(dominant_tax).fillna("")
    else:
        pd.DataFrame(columns=[
            "microbial_compartment", "rank", "ASV_ID", "taxonomy",
            "mean_relative_abundance_within", "mean_relative_abundance_outside",
            "log2_fold_vs_other_compartments",
        ]).to_csv(tables / "microbial_compartment_dominant_asvs.tsv", sep="\t", index=False)
        summary["dominant_ASVs"] = pd.Series(dtype=str)
        summary["dominant_taxa"] = pd.Series(dtype=str)
    summary.to_csv(tables / "microbial_cluster_summary.tsv", sep="\t", index=False)

    metadata = read_table(args.metadata)
    if args.sample_col not in metadata.columns:
        raise ValueError(f"Metadata lacks sample column {args.sample_col}")
    metadata[args.sample_col] = metadata[args.sample_col].astype(str).str.strip()
    if metadata[args.sample_col].duplicated().any():
        raise ValueError("Metadata sample identifiers are not unique")
    joined = primary.merge(metadata, on=args.sample_col, how="left", validate="one_to_one")
    environmental_columns = parse_list(args.environmental_compartment_cols)
    missing_columns = [column for column in environmental_columns if column not in joined.columns]
    comparison_columns = [column for column in environmental_columns if column in joined.columns]
    crosswalk_rows: list[pd.DataFrame] = []
    contingency_rows: list[dict[str, object]] = []
    metric_rows: list[dict[str, object]] = []
    for column in comparison_columns:
        valid = joined[column].notna() & joined[column].astype(str).str.strip().ne("")
        subset = joined.loc[valid, [args.sample_col, "microbial_compartment", column]].copy()
        subset = subset.rename(columns={column: "environmental_compartment"})
        subset["environmental_variable"] = column
        crosswalk_rows.append(subset[[
            args.sample_col, "microbial_compartment", "environmental_variable", "environmental_compartment"
        ]])
        if not supported:
            metric_rows.append({
                "environmental_variable": column, "number_of_samples": len(subset),
                "number_of_microbial_compartments": 0,
                "number_of_environmental_compartments": subset.environmental_compartment.nunique(),
                "adjusted_rand_index": np.nan, "normalized_mutual_information": np.nan,
                "comparison_status": "not_tested_weak_microbiome_clustering",
            })
            continue
        contingency = pd.crosstab(subset.microbial_compartment, subset.environmental_compartment)
        for microbial, row in contingency.iterrows():
            for environmental, count in row.items():
                contingency_rows.append({
                    "environmental_variable": column,
                    "microbial_compartment": microbial,
                    "environmental_compartment": environmental,
                    "number_of_samples": int(count),
                })
        metric_rows.append({
            "environmental_variable": column, "number_of_samples": len(subset),
            "number_of_microbial_compartments": subset.microbial_compartment.nunique(),
            "number_of_environmental_compartments": subset.environmental_compartment.nunique(),
            "adjusted_rand_index": adjusted_rand_score(
                subset.microbial_compartment.astype(str), subset.environmental_compartment.astype(str)
            ),
            "normalized_mutual_information": normalized_mutual_info_score(
                subset.microbial_compartment.astype(str), subset.environmental_compartment.astype(str)
            ),
            "comparison_status": "completed",
        })

    crosswalk = pd.concat(crosswalk_rows, ignore_index=True) if crosswalk_rows else pd.DataFrame(columns=[
        args.sample_col, "microbial_compartment", "environmental_variable", "environmental_compartment"
    ])
    crosswalk.to_csv(tables / "microbial_environmental_compartment_crosswalk.tsv", sep="\t", index=False)
    pd.DataFrame(contingency_rows, columns=[
        "environmental_variable", "microbial_compartment", "environmental_compartment", "number_of_samples"
    ]).to_csv(tables / "microbial_environmental_contingency.tsv", sep="\t", index=False)
    pd.DataFrame(metric_rows, columns=[
        "environmental_variable", "number_of_samples", "number_of_microbial_compartments",
        "number_of_environmental_compartments", "adjusted_rand_index",
        "normalized_mutual_information", "comparison_status",
    ]).to_csv(tables / "microbial_environmental_comparison_metrics.tsv", sep="\t", index=False)

    # Retain the complete source metadata with the inferred microbial state.
    # Downstream analyses may require blocking or stratification columns (for
    # example Cruise) that are intentionally not clustering inputs and are not
    # necessarily listed among environmental comparison columns.
    posthoc_columns = [args.sample_col, "microbial_compartment", "clustering_support_status"]
    posthoc_columns.extend(
        column for column in joined.columns
        if column not in posthoc_columns
    )
    joined.loc[:, posthoc_columns].to_csv(
        tables / "microbial_compartment_posthoc_metadata.tsv", sep="\t", index=False
    )
    pd.DataFrame([
        {"variable": "ASV relative abundance", "used_for_inference": True, "use_stage": "inference"},
        {"variable": "sample identifier", "used_for_inference": True, "use_stage": "matching only"},
        {"variable": "depth", "used_for_inference": False, "use_stage": "post hoc only"},
        {"variable": "environmental measurements", "used_for_inference": False, "use_stage": "not used"},
        {"variable": "environmental compartment labels", "used_for_inference": False, "use_stage": "post hoc only"},
        {"variable": "sampling date", "used_for_inference": False, "use_stage": "post hoc only"},
    ]).to_csv(audit / "microbial_clustering_input_separation_audit.tsv", sep="\t", index=False)
    posthoc_config = {
        **{key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
        "environmental_columns_resolved": comparison_columns,
        "environmental_columns_missing": missing_columns,
        "clustering_supported": supported,
        "comparison_statistics": ["adjusted_rand_index", "normalized_mutual_information"],
        "inferential_contingency_test": "not_run_due_to_nonindependent_repeated_cruise_sampling",
        "software_versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scikit_learn": sklearn.__version__,
        },
    }
    (audit / "microbial_clustering_posthoc_config.json").write_text(json.dumps(posthoc_config, indent=2))


if __name__ == "__main__":
    main()
