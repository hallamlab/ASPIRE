#!/usr/bin/env python3
"""Calculate vertical and temporal turnover from prepared compositions."""

from __future__ import annotations

import argparse
import json
import platform
import re
from pathlib import Path

import numpy as np
import pandas as pd


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", low_memory=False)


def parse_list(value: str) -> list[str]:
    return [item.strip() for item in re.split(r"[,|]", value or "") if item.strip()]


def bootstrap_mean_ci(values: np.ndarray, replicates: int, rng: np.random.Generator) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return np.nan, np.nan
    estimates = np.empty(replicates, dtype=float)
    for index in range(replicates):
        estimates[index] = np.mean(rng.choice(values, size=len(values), replace=True))
    return float(np.quantile(estimates, 0.025)), float(np.quantile(estimates, 0.975))


def sign_flip_pvalue(values: np.ndarray, permutations: int, rng: np.random.Generator) -> float:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return np.nan
    observed = abs(float(np.mean(values)))
    exceed = 0
    for _ in range(permutations):
        signs = rng.choice(np.array([-1.0, 1.0]), size=len(values), replace=True)
        exceed += abs(float(np.mean(values * signs))) >= observed
    return float((exceed + 1) / (permutations + 1))


def benjamini_hochberg(values: pd.Series) -> pd.Series:
    result = pd.Series(np.nan, index=values.index, dtype=float)
    valid = values.dropna().astype(float).sort_values()
    if valid.empty:
        return result
    adjusted = valid.to_numpy() * len(valid) / np.arange(1, len(valid) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result.loc[valid.index] = np.minimum(adjusted, 1.0)
    return result


def consensus(values: pd.Series) -> object:
    unique = sorted({str(value).strip() for value in values.dropna() if str(value).strip()})
    if not unique:
        return np.nan
    return unique[0] if len(unique) == 1 else ";".join(unique)


def compositional_centers(
    metadata: pd.DataFrame,
    clr: pd.DataFrame,
    group_columns: list[str],
    sample_col: str,
    environmental_columns: list[str],
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for keys, group in metadata.groupby(group_columns, sort=True, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        sample_ids = group[sample_col].tolist()
        centroid = clr.loc[sample_ids].mean(axis=0).to_numpy(float)
        composition = np.exp(centroid - np.max(centroid))
        composition /= composition.sum()
        row = dict(zip(group_columns, keys))
        row.update({
            "member_sample_IDs": ";".join(sorted(sample_ids)),
            "number_of_member_samples": len(sample_ids),
            "_clr": centroid,
            "_composition": composition,
        })
        for column in environmental_columns:
            row[column] = consensus(group[column])
        rows.append(row)
    return pd.DataFrame(rows)


def distances(a: pd.Series, b: pd.Series, metrics: list[str]) -> list[tuple[str, float]]:
    output: list[tuple[str, float]] = []
    for metric in metrics:
        if metric == "aitchison":
            value = float(np.linalg.norm(a["_clr"] - b["_clr"]))
        elif metric == "braycurtis":
            denominator = float(np.sum(a["_composition"] + b["_composition"]))
            value = float(np.sum(np.abs(a["_composition"] - b["_composition"])) / denominator)
        else:
            raise ValueError(f"Unsupported distance metric: {metric}")
        output.append((metric, value))
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepared-dir", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--sample-col", default="sampleID")
    parser.add_argument("--profile-col", default="Cruise")
    parser.add_argument("--date-col", default="date")
    parser.add_argument("--depth-col", default="Depth")
    parser.add_argument("--environmental-compartment-cols", default="")
    parser.add_argument("--primary-environmental-compartment-col", default="")
    parser.add_argument("--distance-metrics", default="aitchison,braycurtis")
    parser.add_argument("--fixed-depth-min-profile-fraction", type=float, default=0.50)
    parser.add_argument("--minimum-time-difference-days", type=float, default=0.0)
    parser.add_argument("--boundary-permutations", type=int, default=9999)
    parser.add_argument("--boundary-bootstrap-replicates", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    if not 0 < args.fixed_depth_min_profile_fraction <= 1:
        raise ValueError("--fixed-depth-min-profile-fraction must be in (0, 1]")
    if args.minimum_time_difference_days < 0:
        raise ValueError("--minimum-time-difference-days cannot be negative")
    if args.boundary_permutations < 1 or args.boundary_bootstrap_replicates < 100:
        raise ValueError("Boundary permutations must be positive and bootstrap replicates at least 100")
    metrics = parse_list(args.distance_metrics)
    if not metrics or set(metrics).difference({"aitchison", "braycurtis"}):
        raise ValueError("Distance metrics must be aitchison and/or braycurtis")
    environmental_columns = parse_list(args.environmental_compartment_cols)

    args.outdir.mkdir(parents=True, exist_ok=True)
    tables = args.outdir / "tables"
    audit = args.outdir / "audit"
    tables.mkdir(exist_ok=True)
    audit.mkdir(exist_ok=True)

    clr_table = read_tsv(args.prepared_dir / "asv_clr_matrix.tsv")
    clr_id = clr_table.columns[0]
    clr_table[clr_id] = clr_table[clr_id].astype(str).str.strip()
    clr = clr_table.set_index(clr_id)
    metadata = read_tsv(args.metadata)
    required = [args.sample_col, args.profile_col, args.date_col, args.depth_col]
    missing = [column for column in required if column not in metadata.columns]
    if missing:
        raise ValueError(f"Metadata lacks required columns: {missing}")
    if metadata[args.sample_col].astype(str).str.strip().duplicated().any():
        raise ValueError("Metadata sample identifiers are not unique")
    missing_environmental = [column for column in environmental_columns if column not in metadata.columns]
    environmental_columns = [column for column in environmental_columns if column in metadata.columns]
    primary_environmental = args.primary_environmental_compartment_col.strip()
    if primary_environmental and primary_environmental not in environmental_columns:
        raise ValueError("Primary environmental compartment must be among available configured columns")

    metadata = metadata.copy()
    metadata[args.sample_col] = metadata[args.sample_col].astype(str).str.strip()
    metadata["_date"] = pd.to_datetime(metadata[args.date_col], errors="coerce")
    metadata["_depth"] = pd.to_numeric(metadata[args.depth_col], errors="coerce")
    metadata["_profile"] = metadata[args.profile_col].astype("string").str.strip()
    metadata["present_in_prepared_composition"] = metadata[args.sample_col].isin(clr.index)
    metadata["valid_profile"] = metadata["_profile"].notna() & metadata["_profile"].ne("")
    metadata["valid_date"] = metadata["_date"].notna()
    metadata["valid_depth"] = metadata["_depth"].notna() & np.isfinite(metadata["_depth"])
    metadata["retained_for_turnover"] = (
        metadata.present_in_prepared_composition & metadata.valid_profile &
        metadata.valid_date & metadata.valid_depth
    )
    metadata["turnover_exclusion_reason"] = np.select(
        [
            ~metadata.present_in_prepared_composition,
            ~metadata.valid_profile,
            ~metadata.valid_date,
            ~metadata.valid_depth,
        ],
        ["absent_from_prepared_composition", "missing_profile", "missing_or_invalid_date", "missing_or_invalid_depth"],
        default="retained",
    )
    metadata.to_csv(audit / "community_turnover_sample_audit.tsv", sep="\t", index=False)
    eligible = metadata.loc[metadata.retained_for_turnover].copy()
    if eligible.empty:
        raise ValueError("No samples have valid composition, profile, date, and depth")
    eligible = eligible.sort_values(["_date", "_profile", "_depth", args.sample_col])

    depth_centers = compositional_centers(
        eligible, clr, ["_profile", "_date", "_depth"], args.sample_col, environmental_columns
    )
    duplicate_audit = depth_centers.loc[depth_centers.number_of_member_samples.gt(1), [
        "_profile", "_date", "_depth", "number_of_member_samples", "member_sample_IDs"
    ]].rename(columns={"_profile": args.profile_col, "_date": args.date_col, "_depth": args.depth_col})
    duplicate_audit.to_csv(audit / "duplicate_profile_depth_centers.tsv", sep="\t", index=False)

    vertical_rows: list[dict[str, object]] = []
    for (profile, date), group in depth_centers.groupby(["_profile", "_date"], sort=True):
        group = group.sort_values("_depth").reset_index(drop=True)
        for index in range(len(group) - 1):
            shallow, deep = group.iloc[index], group.iloc[index + 1]
            depth_difference = float(deep["_depth"] - shallow["_depth"])
            if depth_difference <= 0:
                continue
            for metric, value in distances(shallow, deep, metrics):
                row: dict[str, object] = {
                    "profile_ID": profile, "date": date,
                    "shallow_sample_ID": shallow.member_sample_IDs,
                    "deep_sample_ID": deep.member_sample_IDs,
                    "shallow_sample_n": shallow.number_of_member_samples,
                    "deep_sample_n": deep.number_of_member_samples,
                    "shallow_depth": shallow["_depth"], "deep_depth": deep["_depth"],
                    "midpoint_depth": (shallow["_depth"] + deep["_depth"]) / 2,
                    "depth_difference": depth_difference,
                    "distance_metric": metric, "community_distance": value,
                    "turnover_per_meter": value / depth_difference,
                }
                for column in environmental_columns:
                    row[f"{column}_shallow"] = shallow[column]
                    row[f"{column}_deep"] = deep[column]
                    row[f"crosses_{column}_boundary"] = bool(shallow[column] != deep[column]) if pd.notna(shallow[column]) and pd.notna(deep[column]) else np.nan
                if primary_environmental:
                    row["environmental_compartment_shallow"] = shallow[primary_environmental]
                    row["environmental_compartment_deep"] = deep[primary_environmental]
                    row["crosses_compartment_boundary"] = row[f"crosses_{primary_environmental}_boundary"]
                vertical_rows.append(row)
    vertical_columns = [
        "profile_ID", "date", "shallow_sample_ID", "deep_sample_ID",
        "shallow_sample_n", "deep_sample_n", "shallow_depth", "deep_depth",
        "midpoint_depth", "depth_difference", "distance_metric",
        "community_distance", "turnover_per_meter",
    ]
    for column in environmental_columns:
        vertical_columns.extend([f"{column}_shallow", f"{column}_deep", f"crosses_{column}_boundary"])
    if primary_environmental:
        vertical_columns.extend([
            "environmental_compartment_shallow", "environmental_compartment_deep",
            "crosses_compartment_boundary",
        ])
    vertical = pd.DataFrame(vertical_rows, columns=vertical_columns)
    vertical.to_csv(tables / "vertical_turnover.tsv", sep="\t", index=False)

    # Compare microbial turnover across versus within each environmental
    # compartment system. Cruise-level medians make the independently sampled
    # water-column profile the replicate unit and prevent cruises with more
    # depth intervals from receiving disproportionate weight.
    boundary_contrast_rows: list[dict[str, object]] = []
    boundary_summary_rows: list[dict[str, object]] = []
    rng = np.random.default_rng(args.seed)
    display_names = {
        "o2_compartment": "O2",
        "gmm_component": "GMM",
        "o2_subcompartment_final": "Hybrid",
    }
    for metric in metrics:
        metric_data = vertical.loc[vertical.distance_metric.eq(metric)]
        for column in environmental_columns:
            boundary_column = f"crosses_{column}_boundary"
            if boundary_column not in metric_data.columns:
                continue
            usable = metric_data.loc[metric_data[boundary_column].notna()].copy()
            usable["_crosses"] = (
                usable[boundary_column].astype(str).str.strip().str.lower().map(
                    {"true": True, "false": False}
                )
            )
            usable = usable.dropna(subset=["_crosses"])
            strategy_rows: list[dict[str, object]] = []
            for profile, group in usable.groupby("profile_ID", sort=True):
                across = group.loc[group._crosses, "community_distance"]
                within = group.loc[~group._crosses, "community_distance"]
                if across.empty or within.empty:
                    continue
                across_median = float(across.median())
                within_median = float(within.median())
                row = {
                    "environmental_variable": column,
                    "compartment_strategy": display_names.get(column, column),
                    "distance_metric": metric,
                    "profile_ID": profile,
                    "crossing_transition_n": len(across),
                    "within_transition_n": len(within),
                    "median_crossing_turnover": across_median,
                    "median_within_turnover": within_median,
                    "paired_turnover_difference": across_median - within_median,
                    "paired_turnover_ratio": (
                        across_median / within_median if within_median > 0 else np.nan
                    ),
                }
                strategy_rows.append(row)
                boundary_contrast_rows.append(row)
            paired = pd.DataFrame(strategy_rows)
            differences = paired.paired_turnover_difference if not paired.empty else pd.Series(dtype=float)
            crossing_values = usable.loc[usable._crosses, "community_distance"]
            within_values = usable.loc[~usable._crosses, "community_distance"]
            difference_values = differences.to_numpy(float)
            ci_low, ci_high = bootstrap_mean_ci(
                difference_values, args.boundary_bootstrap_replicates, rng
            )
            boundary_summary_rows.append({
                "environmental_variable": column,
                "compartment_strategy": display_names.get(column, column),
                "distance_metric": metric,
                "label_complete_transition_n": len(usable),
                "crossing_transition_n": len(crossing_values),
                "within_transition_n": len(within_values),
                "paired_profile_n": len(paired),
                "pooled_median_crossing_turnover": crossing_values.median() if len(crossing_values) else np.nan,
                "pooled_median_within_turnover": within_values.median() if len(within_values) else np.nan,
                "median_paired_turnover_difference": differences.median() if len(differences) else np.nan,
                "mean_paired_turnover_difference": differences.mean() if len(differences) else np.nan,
                "mean_paired_difference_ci95_low": ci_low,
                "mean_paired_difference_ci95_high": ci_high,
                "paired_sign_flip_pvalue": sign_flip_pvalue(
                    difference_values, args.boundary_permutations, rng
                ),
                "paired_difference_q25": differences.quantile(0.25) if len(differences) else np.nan,
                "paired_difference_q75": differences.quantile(0.75) if len(differences) else np.nan,
                "profiles_with_positive_difference": int(differences.gt(0).sum()),
                "profiles_with_negative_difference": int(differences.lt(0).sum()),
                "profiles_with_zero_difference": int(differences.eq(0).sum()),
            })
    boundary_contrast_columns = [
        "environmental_variable", "compartment_strategy", "distance_metric", "profile_ID",
        "crossing_transition_n", "within_transition_n", "median_crossing_turnover",
        "median_within_turnover", "paired_turnover_difference", "paired_turnover_ratio",
    ]
    pd.DataFrame(boundary_contrast_rows, columns=boundary_contrast_columns).to_csv(
        tables / "vertical_turnover_boundary_cruise_contrasts.tsv", sep="\t", index=False
    )
    boundary_summary_columns = [
        "environmental_variable", "compartment_strategy", "distance_metric",
        "label_complete_transition_n", "crossing_transition_n", "within_transition_n",
        "paired_profile_n", "pooled_median_crossing_turnover", "pooled_median_within_turnover",
        "median_paired_turnover_difference", "mean_paired_turnover_difference",
        "mean_paired_difference_ci95_low", "mean_paired_difference_ci95_high",
        "paired_sign_flip_pvalue", "paired_sign_flip_qvalue",
        "paired_difference_q25", "paired_difference_q75",
        "profiles_with_positive_difference", "profiles_with_negative_difference",
        "profiles_with_zero_difference",
    ]
    boundary_summary = pd.DataFrame(boundary_summary_rows)
    if not boundary_summary.empty:
        boundary_summary["paired_sign_flip_qvalue"] = boundary_summary.groupby(
            "distance_metric", group_keys=False
        )["paired_sign_flip_pvalue"].apply(benjamini_hochberg)
    pd.DataFrame(boundary_summary, columns=boundary_summary_columns).to_csv(
        tables / "vertical_turnover_boundary_summary.tsv", sep="\t", index=False
    )

    profile_count = depth_centers["_profile"].nunique()
    depth_audit = depth_centers.groupby("_depth", as_index=False)["_profile"].nunique().rename(
        columns={"_depth": "depth", "_profile": "number_of_profiles"}
    )
    depth_audit["total_eligible_profiles"] = profile_count
    depth_audit["profile_fraction"] = depth_audit.number_of_profiles / profile_count
    depth_audit["retained_as_common_depth"] = depth_audit.profile_fraction.ge(args.fixed_depth_min_profile_fraction)
    depth_audit.to_csv(audit / "fixed_depth_prevalence_audit.tsv", sep="\t", index=False)
    common_depths = set(depth_audit.loc[depth_audit.retained_as_common_depth, "depth"])

    fixed_rows: list[dict[str, object]] = []
    for depth, group in depth_centers.loc[depth_centers._depth.isin(common_depths)].groupby("_depth", sort=True):
        group = group.sort_values(["_date", "_profile"]).reset_index(drop=True)
        for index in range(len(group) - 1):
            first, second = group.iloc[index], group.iloc[index + 1]
            elapsed_days = (second["_date"] - first["_date"]).total_seconds() / 86400.0
            if elapsed_days <= args.minimum_time_difference_days:
                continue
            for metric, value in distances(first, second, metrics):
                fixed_rows.append({
                    "fixed_depth": depth,
                    "profile_ID_t1": first["_profile"], "profile_ID_t2": second["_profile"],
                    "date_t1": first["_date"], "date_t2": second["_date"],
                    "sample_ID_t1": first.member_sample_IDs, "sample_ID_t2": second.member_sample_IDs,
                    "sample_n_t1": first.number_of_member_samples, "sample_n_t2": second.number_of_member_samples,
                    "depth_t1": first["_depth"], "depth_t2": second["_depth"],
                    "depth_difference": abs(second["_depth"] - first["_depth"]),
                    "depth_match_mode": "exact", "time_difference": elapsed_days,
                    "time_unit": "days", "time_difference_days": elapsed_days,
                    "distance_metric": metric, "community_distance": value,
                    "turnover_per_unit_time": value / elapsed_days,
                    "turnover_per_day": value / elapsed_days,
                })
    fixed_columns = [
        "fixed_depth", "profile_ID_t1", "profile_ID_t2", "date_t1", "date_t2",
        "sample_ID_t1", "sample_ID_t2", "sample_n_t1", "sample_n_t2",
        "depth_t1", "depth_t2", "depth_difference", "depth_match_mode",
        "time_difference", "time_unit", "time_difference_days",
        "distance_metric", "community_distance", "turnover_per_unit_time", "turnover_per_day",
    ]
    pd.DataFrame(fixed_rows, columns=fixed_columns).to_csv(
        tables / "temporal_fixed_depth_turnover.tsv", sep="\t", index=False
    )

    compartment_rows: list[dict[str, object]] = []
    compartment_audit_rows: list[dict[str, object]] = []
    for column in environmental_columns:
        usable = eligible.loc[eligible[column].notna() & eligible[column].astype(str).str.strip().ne("")]
        if usable.empty:
            compartment_audit_rows.append({"environmental_variable": column, "status": "no_observed_labels"})
            continue
        centers = compositional_centers(
            usable, clr, ["_profile", "_date", column], args.sample_col, []
        )
        for compartment, group in centers.groupby(column, sort=True):
            group = group.sort_values(["_date", "_profile"]).reset_index(drop=True)
            for index in range(len(group) - 1):
                first, second = group.iloc[index], group.iloc[index + 1]
                elapsed_days = (second["_date"] - first["_date"]).total_seconds() / 86400.0
                if elapsed_days <= args.minimum_time_difference_days:
                    continue
                for metric, value in distances(first, second, metrics):
                    compartment_rows.append({
                        "environmental_variable": column,
                        "environmental_compartment": compartment,
                        "profile_ID_t1": first["_profile"], "profile_ID_t2": second["_profile"],
                        "date_t1": first["_date"], "date_t2": second["_date"],
                        "sample_IDs_t1": first.member_sample_IDs, "sample_IDs_t2": second.member_sample_IDs,
                        "sample_n_t1": first.number_of_member_samples, "sample_n_t2": second.number_of_member_samples,
                        "time_difference": elapsed_days, "time_unit": "days",
                        "time_difference_days": elapsed_days,
                        "distance_metric": metric, "community_distance": value,
                        "turnover_per_unit_time": value / elapsed_days,
                        "turnover_per_day": value / elapsed_days,
                    })
        compartment_audit_rows.append({
            "environmental_variable": column, "status": "completed",
            "number_of_labeled_samples": len(usable),
            "number_of_profile_compartment_centers": len(centers),
            "number_of_compartments": centers[column].nunique(),
        })
    compartment_columns = [
        "environmental_variable", "environmental_compartment",
        "profile_ID_t1", "profile_ID_t2", "date_t1", "date_t2",
        "sample_IDs_t1", "sample_IDs_t2", "sample_n_t1", "sample_n_t2",
        "time_difference", "time_unit", "time_difference_days", "distance_metric",
        "community_distance", "turnover_per_unit_time", "turnover_per_day",
    ]
    pd.DataFrame(compartment_rows, columns=compartment_columns).to_csv(
        tables / "temporal_environmental_compartment_turnover.tsv", sep="\t", index=False
    )
    temporal_compartment = pd.DataFrame(compartment_rows, columns=compartment_columns)
    temporal_summary_rows: list[dict[str, object]] = []
    for (variable, compartment, metric), group in temporal_compartment.groupby(
        ["environmental_variable", "environmental_compartment", "distance_metric"], sort=True
    ):
        temporal_summary_rows.append({
            "environmental_variable": variable,
            "compartment_strategy": display_names.get(variable, variable),
            "environmental_compartment": compartment,
            "distance_metric": metric,
            "transition_n": len(group),
            "median_community_distance": group.community_distance.median(),
            "community_distance_q25": group.community_distance.quantile(0.25),
            "community_distance_q75": group.community_distance.quantile(0.75),
            "median_turnover_per_day": group.turnover_per_day.median(),
            "turnover_per_day_q25": group.turnover_per_day.quantile(0.25),
            "turnover_per_day_q75": group.turnover_per_day.quantile(0.75),
        })
    pd.DataFrame(temporal_summary_rows, columns=[
        "environmental_variable", "compartment_strategy", "environmental_compartment",
        "distance_metric", "transition_n", "median_community_distance",
        "community_distance_q25", "community_distance_q75", "median_turnover_per_day",
        "turnover_per_day_q25", "turnover_per_day_q75",
    ]).to_csv(tables / "temporal_compartment_turnover_summary.tsv", sep="\t", index=False)
    pd.DataFrame(compartment_audit_rows, columns=[
        "environmental_variable", "status", "number_of_labeled_samples",
        "number_of_profile_compartment_centers", "number_of_compartments",
    ]).to_csv(
        audit / "temporal_environmental_compartment_audit.tsv", sep="\t", index=False
    )
    config = {
        **{key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
        "distance_metrics_resolved": metrics,
        "environmental_compartment_columns_resolved": environmental_columns,
        "environmental_compartment_columns_missing": missing_environmental,
        "fixed_depth_matching": "exact_numeric_equality",
        "duplicate_profile_depth_handling": "CLR compositional center",
        "compartment_aggregation": "CLR compositional center",
        "time_unit": "days",
        "taxonomy_input": "not_required_for_sample_and_transition_level_metrics",
        "software_versions": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__},
    }
    (audit / "community_turnover_config.json").write_text(json.dumps(config, indent=2))


if __name__ == "__main__":
    main()
