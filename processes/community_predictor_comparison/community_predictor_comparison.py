#!/usr/bin/env python3
"""Matched comparisons of environmental predictors and community groupings.

The cruise analysis compares continuous BASIN cruise metrics as predictors of
depth-averaged ASV composition.  The sample analysis compares categorical
bottle-level organizations on identical observed samples.  All reported model
comparisons use identical rows within an analysis branch.
"""

from __future__ import annotations

import argparse
from itertools import combinations
import json
import math
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys as _sys
_sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shared_plot_style import install_publication_style
install_publication_style()
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.spatial.distance import pdist, squareform
from scipy.stats import f_oneway


def parse_csv(value: str | None) -> list[str]:
    return [x.strip() for x in re.split(r"[,|]", value or "") if x.strip()]


def read_table(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    return pd.read_csv(path, sep=sep)


def normalize_id(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    try:
        number = float(text)
        if np.isfinite(number) and number.is_integer():
            return str(int(number))
    except ValueError:
        pass
    return text


def read_counts(path: Path, metadata_ids: set[str]) -> pd.DataFrame:
    counts = pd.read_csv(path, sep="\t", index_col=0)
    counts.index = counts.index.map(normalize_id)
    counts.columns = [normalize_id(x) for x in counts.columns]
    row_hits = len(set(counts.index) & metadata_ids)
    col_hits = len(set(counts.columns) & metadata_ids)
    if col_hits >= row_hits:
        counts = counts.T
    counts.index = counts.index.map(normalize_id)
    counts = counts.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    counts = counts.loc[:, counts.sum(axis=0) > 0]
    return counts


def relative(counts: pd.DataFrame) -> pd.DataFrame:
    return counts.div(counts.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)


def hellinger(counts: pd.DataFrame) -> pd.DataFrame:
    return np.sqrt(relative(counts))


def cruise_depth_profiles(counts: pd.DataFrame, metadata: pd.DataFrame, cruise_col: str,
                          depth_col: str, min_prevalence: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build equal-depth-weighted Hellinger water-column profiles per cruise."""
    sample_rel = relative(counts)
    frame = metadata.loc[sample_rel.index, [cruise_col, depth_col]].copy()
    frame["_cruise"] = frame[cruise_col].map(normalize_id)
    frame["_depth"] = pd.to_numeric(frame[depth_col], errors="coerce")
    valid = frame["_cruise"].ne("") & frame["_depth"].notna()
    sample_rel = sample_rel.loc[valid]
    frame = frame.loc[valid]
    indexed = sample_rel.copy()
    indexed["_cruise"] = frame["_cruise"]
    indexed["_depth"] = frame["_depth"]
    observed = indexed.groupby(["_cruise", "_depth"]).mean()
    n_cruises = observed.index.get_level_values(0).nunique()
    coverage = observed.reset_index().groupby("_depth")["_cruise"].nunique()
    retained_depths = sorted(coverage[coverage / n_cruises >= min_prevalence].index.tolist())
    if len(retained_depths) < 2:
        raise ValueError("Fewer than two anchored depths passed cruise_depth_min_prevalence")
    rows, row_names, audit = [], [], []
    for cruise in sorted(observed.index.get_level_values(0).unique()):
        profile = observed.xs(cruise, level=0).sort_index()
        n_observed = int(profile.index.isin(retained_depths).sum())
        profile = profile.reindex(retained_depths).interpolate(method="index", limit_direction="both")
        profile = profile.div(profile.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
        values = np.sqrt(profile.to_numpy()).reshape(-1)
        rows.append(values); row_names.append(cruise)
        audit.append({"cruise": cruise, "observed_retained_depths": n_observed,
                      "retained_depths": len(retained_depths),
                      "interpolated_depths": len(retained_depths) - n_observed,
                      "coverage_fraction": n_observed / len(retained_depths)})
    columns = [f"depth_{depth:g}__{asv}" for depth in retained_depths for asv in sample_rel.columns]
    return pd.DataFrame(rows, index=row_names, columns=columns), pd.DataFrame(audit)


def design_matrix(df: pd.DataFrame, terms: list[str], continuous: set[str]) -> np.ndarray:
    pieces = [np.ones((len(df), 1), dtype=float)]
    for term in terms:
        if term in continuous:
            values = pd.to_numeric(df[term], errors="coerce").to_numpy(dtype=float)
            sd = float(np.std(values, ddof=0))
            pieces.append(((values - np.mean(values)) / sd if sd > 0 else np.zeros(len(values)))[:, None])
        else:
            categorical = pd.get_dummies(df[term].astype(str), prefix=term, drop_first=True, dtype=float)
            if categorical.shape[1]:
                pieces.append(categorical.to_numpy(dtype=float))
    return np.column_stack(pieces)


def train_test_design(
    train: pd.DataFrame, test: pd.DataFrame, terms: list[str], continuous: set[str]
) -> tuple[np.ndarray, np.ndarray]:
    """Build leakage-free train/test matrices with training-derived encodings."""
    train_parts = [np.ones((len(train), 1), dtype=float)]
    test_parts = [np.ones((len(test), 1), dtype=float)]
    for term in terms:
        if term in continuous:
            train_values = pd.to_numeric(train[term], errors="coerce").to_numpy(dtype=float)
            test_values = pd.to_numeric(test[term], errors="coerce").to_numpy(dtype=float)
            mean = float(np.mean(train_values))
            sd = float(np.std(train_values, ddof=0))
            train_parts.append(((train_values - mean) / sd if sd > 0 else np.zeros(len(train)))[:, None])
            test_parts.append(((test_values - mean) / sd if sd > 0 else np.zeros(len(test)))[:, None])
        else:
            levels = sorted(train[term].astype(str).unique().tolist())
            for level in levels[1:]:
                train_parts.append(train[term].astype(str).eq(level).to_numpy(dtype=float)[:, None])
                test_parts.append(test[term].astype(str).eq(level).to_numpy(dtype=float)[:, None])
    return np.column_stack(train_parts), np.column_stack(test_parts)


def hat(x: np.ndarray) -> np.ndarray:
    return x @ np.linalg.pinv(x)


def gower(distance: np.ndarray) -> np.ndarray:
    n = distance.shape[0]
    center = np.eye(n) - np.ones((n, n)) / n
    return -0.5 * center @ (distance ** 2) @ center


def nested_distance_stats(
    distance: np.ndarray,
    frame: pd.DataFrame,
    baseline_terms: list[str],
    tested_terms: list[str],
    continuous: set[str],
) -> dict[str, float]:
    n = len(frame)
    x0 = design_matrix(frame, baseline_terms, continuous)
    x1 = design_matrix(frame, baseline_terms + tested_terms, continuous)
    h0, h1 = hat(x0), hat(x1)
    rank0, rank1 = np.linalg.matrix_rank(x0), np.linalg.matrix_rank(x1)
    gg = gower(distance)
    total = float(np.trace(gg))
    ss_base = float(np.trace(h0 @ gg))
    ss_full = float(np.trace(h1 @ gg))
    ss_increment = max(0.0, ss_full - ss_base)
    ss_residual = max(0.0, total - ss_full)
    df_increment = max(1, rank1 - rank0)
    df_residual = max(1, n - rank1)
    pseudo_f = (ss_increment / df_increment) / (ss_residual / df_residual) if ss_residual > 0 else np.nan
    r2 = ss_full / total if total > 0 else np.nan
    baseline_r2 = ss_base / total if total > 0 else np.nan
    partial_r2 = ss_increment / (ss_increment + ss_residual) if (ss_increment + ss_residual) > 0 else np.nan
    predictors = max(0, rank1 - 1)
    adjusted_r2 = 1.0 - (1.0 - r2) * (n - 1) / (n - predictors - 1) if n > predictors + 1 else np.nan
    return {
        "n": n,
        "rank": rank1,
        "df_test": df_increment,
        "pseudo_f": pseudo_f,
        "r2": r2,
        "baseline_r2": baseline_r2,
        "incremental_r2": r2 - baseline_r2,
        "partial_r2": partial_r2,
        "adjusted_r2": adjusted_r2,
    }


def permutation_pvalue(
    distance: np.ndarray,
    frame: pd.DataFrame,
    baseline_terms: list[str],
    tested_terms: list[str],
    continuous: set[str],
    permutations: int,
    seed: int,
    blocks: pd.Series | None = None,
) -> float:
    observed = nested_distance_stats(distance, frame, baseline_terms, tested_terms, continuous)["pseudo_f"]
    if not np.isfinite(observed):
        return np.nan
    rng = np.random.default_rng(seed)
    hits = 0
    completed = 0
    for _ in range(permutations):
        permuted = frame.copy()
        if blocks is None:
            order = rng.permutation(len(frame))
            permuted.loc[:, tested_terms] = frame.iloc[order][tested_terms].to_numpy()
        else:
            block_values = blocks.reset_index(drop=True).astype(str)
            for level in block_values.unique():
                idx = np.flatnonzero(block_values.to_numpy() == level)
                if len(idx) > 1:
                    permuted.loc[idx, tested_terms] = frame.iloc[rng.permutation(idx)][tested_terms].to_numpy()
        stat = nested_distance_stats(distance, permuted, baseline_terms, tested_terms, continuous)["pseudo_f"]
        if np.isfinite(stat):
            hits += int(stat >= observed)
            completed += 1
    return (hits + 1) / (completed + 1) if completed else np.nan


def multivariate_cv(
    response: np.ndarray,
    frame: pd.DataFrame,
    models: dict[str, tuple[list[str], set[str]]],
    block_col: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    blocks = frame[block_col].astype(str)
    for fold in sorted(blocks.unique()):
        test = blocks.eq(fold).to_numpy()
        train = ~test
        if test.sum() == 0 or train.sum() < 3:
            continue
        train_mean = response[train].mean(axis=0)
        baseline_sse = float(np.square(response[test] - train_mean).sum())
        for model_name, (terms, continuous) in models.items():
            train_frame = frame.loc[train].copy()
            test_frame = frame.loc[test].copy()
            x_train, x_test = train_test_design(train_frame, test_frame, terms, continuous)
            coefficients = np.linalg.pinv(x_train) @ response[train]
            predicted = x_test @ coefficients
            sse = float(np.square(response[test] - predicted).sum())
            rows.append({
                "fold": fold,
                "model": model_name,
                "n_train": int(train.sum()),
                "n_test": int(test.sum()),
                "sse": sse,
                "baseline_sse": baseline_sse,
                "cv_r2": 1.0 - sse / baseline_sse if baseline_sse > 0 else np.nan,
                "rmse": math.sqrt(sse / (test.sum() * response.shape[1])),
            })
    folds = pd.DataFrame(rows)
    summary = (
        folds.groupby("model", as_index=False)
        .agg(folds=("fold", "nunique"), cv_r2_mean=("cv_r2", "mean"), cv_r2_median=("cv_r2", "median"),
             cv_r2_sd=("cv_r2", "std"), rmse_mean=("rmse", "mean"))
    ) if not folds.empty else pd.DataFrame()
    return folds, summary


def paired_difference(folds: pd.DataFrame, model_a: str, model_b: str, seed: int) -> pd.DataFrame:
    pivot = folds.pivot(index="fold", columns="model", values="cv_r2").dropna(subset=[model_a, model_b])
    if pivot.empty:
        return pd.DataFrame()
    diffs = (pivot[model_a] - pivot[model_b]).to_numpy(dtype=float)
    rng = np.random.default_rng(seed)
    boot = np.array([rng.choice(diffs, size=len(diffs), replace=True).mean() for _ in range(5000)])
    signs = rng.choice([-1.0, 1.0], size=(10000, len(diffs)))
    p = (1 + np.sum(np.abs((signs * diffs).mean(axis=1)) >= abs(diffs.mean()))) / 10001
    return pd.DataFrame([{
        "model_a": model_a,
        "model_b": model_b,
        "metric": "cv_r2",
        "n_folds": len(diffs),
        "mean_difference_a_minus_b": float(diffs.mean()),
        "ci_lower": float(np.quantile(boot, 0.025)),
        "ci_upper": float(np.quantile(boot, 0.975)),
        "sign_flip_p_value": float(p),
    }])


def bh_adjust(values: pd.Series) -> pd.Series:
    """Benjamini-Hochberg adjusted p-values, preserving missing values."""
    result = pd.Series(np.nan, index=values.index, dtype=float)
    valid = values.dropna().astype(float)
    if valid.empty:
        return result
    ordered = valid.sort_values()
    adjusted = (ordered.to_numpy() * len(ordered) / np.arange(1, len(ordered) + 1))[::-1]
    adjusted = np.minimum.accumulate(adjusted)[::-1]
    result.loc[ordered.index] = np.minimum(adjusted, 1.0)
    return result


def pcoa(distance: np.ndarray) -> np.ndarray:
    gg = gower(distance)
    values, vectors = np.linalg.eigh(gg)
    order = np.argsort(values)[::-1]
    values, vectors = values[order], vectors[:, order]
    keep = values > 1e-10
    return vectors[:, keep] * np.sqrt(values[keep]) if keep.any() else np.zeros((len(distance), 1))


def permdisp(distance: np.ndarray, labels: pd.Series, permutations: int, seed: int, blocks: pd.Series) -> dict:
    coords = pcoa(distance)
    label_values = labels.astype(str).to_numpy()
    def distances_to_centroid(current: np.ndarray) -> np.ndarray:
        result = np.zeros(len(current))
        for group in pd.unique(current):
            idx = np.flatnonzero(current == group)
            result[idx] = np.linalg.norm(coords[idx] - coords[idx].mean(axis=0), axis=1)
        return result
    observed_dist = distances_to_centroid(label_values)
    groups = [observed_dist[label_values == g] for g in pd.unique(label_values)]
    observed_f = float(f_oneway(*groups).statistic) if len(groups) > 1 else np.nan
    rng = np.random.default_rng(seed)
    hits = 0
    block_values = blocks.astype(str).to_numpy()
    for _ in range(permutations):
        shuffled = label_values.copy()
        for block in pd.unique(block_values):
            idx = np.flatnonzero(block_values == block)
            shuffled[idx] = shuffled[rng.permutation(idx)]
        perm_dist = distances_to_centroid(shuffled)
        perm_groups = [perm_dist[shuffled == g] for g in pd.unique(shuffled)]
        stat = float(f_oneway(*perm_groups).statistic)
        hits += int(np.isfinite(stat) and stat >= observed_f)
    return {"f_statistic": observed_f, "p_value": (hits + 1) / (permutations + 1)}


def save_figure(fig, base: Path, formats: list[str]) -> None:
    for fmt in formats:
        fig.savefig(base.with_suffix(f".{fmt}"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_model_comparison(table: pd.DataFrame, base: Path, formats: list[str], title: str) -> None:
    if table.empty:
        return
    long = table.melt(id_vars="model", value_vars=["adjusted_r2", "cv_r2_mean"], var_name="metric", value_name="value")
    long["metric"] = long["metric"].map({"adjusted_r2": "Adjusted distance R²", "cv_r2_mean": "Cross-validated R²"})
    fig, ax = plt.subplots(figsize=(10, 5.5))
    sns.barplot(data=long, x="model", y="value", hue="metric", ax=ax)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("Performance")
    ax.set_xlabel("")
    ax.set_title(title)
    ax.tick_params(axis="x", rotation=25)
    ax.legend(frameon=False)
    fig.tight_layout()
    save_figure(fig, base, formats)


def plot_cv_folds(folds: pd.DataFrame, base: Path, formats: list[str], title: str) -> None:
    if folds.empty:
        return
    fig, ax = plt.subplots(figsize=(10, 5.5))
    sns.pointplot(data=folds, x="fold", y="cv_r2", hue="model", dodge=0.25, ax=ax)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("Held-out multivariate R²")
    ax.set_xlabel("Held-out block")
    ax.set_title(title)
    ax.tick_params(axis="x", rotation=45)
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    save_figure(fig, base, formats)


def plot_group_ordination(distance: np.ndarray, frame: pd.DataFrame, groups: list[str], base: Path, formats: list[str]) -> None:
    coords = pcoa(distance)
    fig, axes = plt.subplots(1, len(groups), figsize=(7 * len(groups), 6), squeeze=False)
    for ax, group in zip(axes[0], groups):
        plot = frame.copy()
        plot["PCoA1"], plot["PCoA2"] = coords[:, 0], coords[:, 1] if coords.shape[1] > 1 else 0
        sns.scatterplot(data=plot, x="PCoA1", y="PCoA2", hue=group, s=45, alpha=0.8, ax=ax)
        ax.set_title(group)
        ax.legend(frameon=False, fontsize=8)
    fig.suptitle("Matched-sample Bray–Curtis ordination", y=1.01)
    fig.tight_layout()
    save_figure(fig, base, formats)


def plot_compartments_within_cruise_groups(distance: np.ndarray, frame: pd.DataFrame,
                                            cruise_group_col: str, compartment_col: str,
                                            base: Path, formats: list[str]) -> None:
    """Show the same community ordination split by independently defined cruise group."""
    coords = pcoa(distance)
    plot = frame.copy()
    plot["PCoA1"], plot["PCoA2"] = coords[:, 0], coords[:, 1] if coords.shape[1] > 1 else 0
    levels = sorted(plot[cruise_group_col].astype(str).unique())
    fig, axes = plt.subplots(1, len(levels), figsize=(7 * len(levels), 6), squeeze=False,
                             sharex=True, sharey=True)
    for ax, level in zip(axes[0], levels):
        subset = plot[plot[cruise_group_col].astype(str).eq(level)]
        sns.scatterplot(data=subset, x="PCoA1", y="PCoA2", hue=compartment_col,
                        s=45, alpha=.8, ax=ax)
        ax.set_title(str(level)); ax.legend(frameon=False, fontsize=7)
    fig.suptitle("Bottle compartments within environmental cruise groups", y=1.01)
    fig.tight_layout(); save_figure(fig, base, formats)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--asv-counts", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--sample-col", default="sampleID")
    parser.add_argument("--cruise-col", default="Cruise")
    parser.add_argument("--year-col", default="Year")
    parser.add_argument("--date-col", default="date")
    parser.add_argument("--season-col", default="Season")
    parser.add_argument("--depth-col", default="Depth")
    parser.add_argument("--cruise-depth-min-prevalence", type=float, default=0.50)
    parser.add_argument("--pea-col", default="pea_J_m3")
    parser.add_argument("--centroid-col", default="depth_centroid_distance")
    parser.add_argument("--cruise-group-col", default="cruise_group")
    parser.add_argument("--cruise-group-probability-col", default="cruise_group_max_prob")
    parser.add_argument("--cruise-group-uncertain-col", default="cruise_group_uncertain")
    parser.add_argument("--renewal-group-col", default="renewal_phase")
    parser.add_argument("--o2-group-col", default="o2_compartment")
    parser.add_argument("--gmm-group-col", default="gmm_component")
    parser.add_argument("--hybrid-group-col", default="o2_subcompartment_final")
    parser.add_argument("--assignment-source-col", default="o2_subcompartment_final_assignment_source")
    parser.add_argument("--observed-source-label", default="observed")
    parser.add_argument("--min-group-n", type=int, default=3)
    parser.add_argument("--permutations", type=int, default=999)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--formats", default="pdf,png,svg")
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    tables = args.outdir / "tables"
    plots = args.outdir / "plots"
    tables.mkdir(exist_ok=True)
    plots.mkdir(exist_ok=True)
    formats = parse_csv(args.formats)

    metadata = read_table(args.metadata)
    metadata[args.sample_col] = metadata[args.sample_col].map(normalize_id)
    metadata = metadata.drop_duplicates(args.sample_col).set_index(args.sample_col, drop=False)
    counts = read_counts(args.asv_counts, set(metadata.index))
    common = metadata.index.intersection(counts.index)
    metadata = metadata.loc[common].copy()
    counts = counts.loc[common].copy()

    audit_rows: list[dict] = []

    # Cruise-level matched continuous-predictor comparison.
    cruise_columns = [args.cruise_col, args.year_col, args.season_col, args.pea_col, args.centroid_col]
    if args.date_col in metadata.columns:
        cruise_columns.append(args.date_col)
    missing_cruise_columns = [c for c in cruise_columns if c not in metadata.columns]
    cruise_results = pd.DataFrame()
    if not missing_cruise_columns:
        meta_cruise = metadata.dropna(subset=cruise_columns).copy()
        cruise_ids = meta_cruise[args.cruise_col].map(normalize_id)
        cruise_response, depth_audit = cruise_depth_profiles(
            counts.loc[meta_cruise.index], meta_cruise, args.cruise_col, args.depth_col,
            args.cruise_depth_min_prevalence,
        )
        depth_audit.to_csv(tables / "cruise_depth_profile_coverage.tsv", sep="\t", index=False)
        cruise_frame = (
            meta_cruise.assign(_cruise=cruise_ids)
            .drop_duplicates("_cruise")
            .set_index("_cruise")
            .loc[cruise_response.index, cruise_columns]
            .reset_index(drop=True)
        )
        distance = squareform(pdist(cruise_response.to_numpy(), metric="braycurtis"))
        baseline = [args.season_col]
        candidates = {
            "PEA": [args.pea_col],
            "Depth centroid distance": [args.centroid_col],
            "PEA + centroid": [args.pea_col, args.centroid_col],
        }
        model_rows = []
        for name, tested in candidates.items():
            continuous = set(tested)
            stats = nested_distance_stats(distance, cruise_frame, baseline, tested, continuous)
            stats["p_value"] = permutation_pvalue(
                distance, cruise_frame, baseline, tested, continuous, args.permutations, args.seed
            )
            stats["model"] = name
            model_rows.append(stats)
        cruise_results = pd.DataFrame(model_rows)
        cv_models = {
            "Season baseline": ([args.season_col], set()),
            "PEA": ([args.season_col, args.pea_col], {args.pea_col}),
            "Depth centroid distance": ([args.season_col, args.centroid_col], {args.centroid_col}),
            "PEA + centroid": ([args.season_col, args.pea_col, args.centroid_col], {args.pea_col, args.centroid_col}),
        }
        cruise_folds, cruise_cv = multivariate_cv(
            cruise_response.to_numpy(), cruise_frame, cv_models, args.year_col
        )
        cruise_results = cruise_results.merge(cruise_cv, on="model", how="left")
        cruise_results.to_csv(tables / "cruise_model_comparison.tsv", sep="\t", index=False)
        cruise_folds.to_csv(tables / "cruise_cv_folds.tsv", sep="\t", index=False)
        paired_difference(cruise_folds, "Depth centroid distance", "PEA", args.seed).to_csv(
            tables / "cruise_centroid_vs_pea_performance.tsv", sep="\t", index=False
        )
        cruise_frame.assign(cruise_asv_samples=meta_cruise.groupby(cruise_ids).size().reindex(cruise_response.index).to_numpy()).to_csv(
            tables / "cruise_analysis_cohort.tsv", sep="\t", index=False
        )
        plot_model_comparison(cruise_results, plots / "cruise_model_performance", formats,
                              "Cruise-level community prediction")
        plot_cv_folds(cruise_folds, plots / "cruise_cv_performance", formats,
                      "Leave-year-out community prediction")

        # Neutral cruise-group comparison on a common PEA/centroid/group cohort.
        if args.cruise_group_col in meta_cruise.columns:
            group_meta = (
                meta_cruise.assign(_cruise=cruise_ids).drop_duplicates("_cruise")
                .set_index("_cruise").reindex(cruise_response.index)
            )
            matched = group_meta[args.cruise_group_col].notna().to_numpy()
            group_frame = cruise_frame.loc[matched].reset_index(drop=True).copy()
            group_frame[args.cruise_group_col] = group_meta.loc[matched, args.cruise_group_col].astype(str).to_numpy()
            for optional in (args.cruise_group_probability_col, args.cruise_group_uncertain_col):
                if optional in group_meta.columns:
                    group_frame[optional] = group_meta.loc[matched, optional].to_numpy()
            group_response = cruise_response.to_numpy()[matched]
            group_distance = squareform(pdist(group_response, metric="braycurtis"))
            stats = nested_distance_stats(
                group_distance, group_frame, [args.season_col], [args.cruise_group_col], set()
            )
            stats["p_value"] = permutation_pvalue(
                group_distance, group_frame, [args.season_col], [args.cruise_group_col], set(),
                args.permutations, args.seed, blocks=group_frame[args.year_col],
            )
            dispersion = permdisp(
                group_distance, group_frame[args.cruise_group_col], args.permutations,
                args.seed, group_frame[args.year_col],
            )
            stats["permdisp_f"] = dispersion["f_statistic"]
            stats["permdisp_p_value"] = dispersion["p_value"]
            stats["model"] = "Cruise group"
            matched_models = {
                "Season baseline": ([args.season_col], set()),
                "PEA": ([args.season_col, args.pea_col], {args.pea_col}),
                "Depth centroid distance": ([args.season_col, args.centroid_col], {args.centroid_col}),
                "PEA + centroid": ([args.season_col, args.pea_col, args.centroid_col], {args.pea_col, args.centroid_col}),
                "Cruise group": ([args.season_col, args.cruise_group_col], set()),
            }
            group_folds, group_cv = multivariate_cv(group_response, group_frame, matched_models, args.year_col)
            group_result = group_cv.copy()
            for column, value in stats.items():
                if column != "model":
                    group_result[column] = np.where(group_result["model"] == "Cruise group", value, np.nan)
            group_result.to_csv(tables / "cruise_group_model_comparison.tsv", sep="\t", index=False)
            group_folds.to_csv(tables / "cruise_group_matched_cv_folds.tsv", sep="\t", index=False)
            comparisons = []
            for comparator in ("PEA", "Depth centroid distance", "PEA + centroid"):
                current = paired_difference(group_folds, "Cruise group", comparator, args.seed)
                if not current.empty:
                    comparisons.append(current)
            if comparisons:
                pd.concat(comparisons, ignore_index=True).to_csv(
                    tables / "cruise_group_paired_performance.tsv", sep="\t", index=False
                )
            group_frame.to_csv(tables / "cruise_group_analysis_cohort.tsv", sep="\t", index=False)
            plot_cv_folds(group_folds, plots / "cruise_group_matched_cv_performance", formats,
                          "Matched leave-year-out cruise community prediction")
            plot_group_ordination(group_distance, group_frame, [args.cruise_group_col],
                                  plots / "cruise_group_ordination", formats)

            if args.cruise_group_uncertain_col in group_frame.columns:
                uncertain = group_frame[args.cruise_group_uncertain_col].astype(str).str.lower().isin(["true", "1", "yes"])
                confident_frame = group_frame.loc[~uncertain].reset_index(drop=True)
                confident_response = group_response[~uncertain.to_numpy()]
                confident_folds, confident_cv = multivariate_cv(
                    confident_response, confident_frame, matched_models, args.year_col
                )
                confident_cv.to_csv(tables / "cruise_group_confident_model_comparison.tsv", sep="\t", index=False)
                confident_folds.to_csv(tables / "cruise_group_confident_cv_folds.tsv", sep="\t", index=False)
                confident_frame.to_csv(tables / "cruise_group_confident_cohort.tsv", sep="\t", index=False)

        # Test BASIN's nitrate-qualified renewal phase at the cruise level.
        # Each cruise contributes one water-column community profile;
        # inference is adjusted for season and permutations are restricted
        # within year to avoid treating bottles as replicates.
        if args.renewal_group_col in meta_cruise.columns:
            renewal_meta = (
                meta_cruise.assign(_cruise=cruise_ids).drop_duplicates("_cruise")
                .set_index("_cruise").reindex(cruise_response.index)
            )
            matched = renewal_meta[args.renewal_group_col].notna().to_numpy()
            renewal_frame = cruise_frame.loc[matched].reset_index(drop=True).copy()
            renewal_frame[args.renewal_group_col] = (
                renewal_meta.loc[matched, args.renewal_group_col].astype(str).to_numpy()
            )
            renewal_response = cruise_response.to_numpy()[matched]
            if renewal_frame[args.renewal_group_col].nunique() >= 2 and len(renewal_frame) >= 6:
                renewal_distance = squareform(pdist(renewal_response, metric="braycurtis"))
                stats = nested_distance_stats(
                    renewal_distance, renewal_frame, [args.season_col],
                    [args.renewal_group_col], set(),
                )
                stats["p_value"] = permutation_pvalue(
                    renewal_distance, renewal_frame, [args.season_col],
                    [args.renewal_group_col], set(), args.permutations, args.seed,
                    blocks=renewal_frame[args.year_col],
                )
                dispersion = permdisp(
                    renewal_distance, renewal_frame[args.renewal_group_col],
                    args.permutations, args.seed, renewal_frame[args.year_col],
                )
                stats["permdisp_f"] = dispersion["f_statistic"]
                stats["permdisp_p_value"] = dispersion["p_value"]
                stats["model"] = "Renewal phase"
                renewal_models = {
                    "Season baseline": ([args.season_col], set()),
                    "PEA": ([args.season_col, args.pea_col], {args.pea_col}),
                    "Depth centroid distance": (
                        [args.season_col, args.centroid_col], {args.centroid_col}
                    ),
                    "PEA + centroid": (
                        [args.season_col, args.pea_col, args.centroid_col],
                        {args.pea_col, args.centroid_col},
                    ),
                    "Renewal phase": (
                        [args.season_col, args.renewal_group_col], set()
                    ),
                }
                renewal_folds, renewal_cv = multivariate_cv(
                    renewal_response, renewal_frame, renewal_models, args.year_col
                )
                renewal_result = renewal_cv.copy()
                for column, value in stats.items():
                    if column != "model":
                        renewal_result[column] = np.where(
                            renewal_result["model"] == "Renewal phase", value, np.nan
                        )
                renewal_result.to_csv(
                    tables / "renewal_group_model_comparison.tsv", sep="\t", index=False
                )
                renewal_folds.to_csv(
                    tables / "renewal_group_matched_cv_folds.tsv", sep="\t", index=False
                )
                paired = []
                for comparator in ("PEA", "Depth centroid distance", "PEA + centroid"):
                    current = paired_difference(
                        renewal_folds, "Renewal phase", comparator, args.seed
                    )
                    if not current.empty:
                        paired.append(current)
                if paired:
                    pd.concat(paired, ignore_index=True).to_csv(
                        tables / "renewal_group_paired_performance.tsv", sep="\t", index=False
                    )
                renewal_frame.to_csv(
                    tables / "renewal_group_analysis_cohort.tsv", sep="\t", index=False
                )
                plot_cv_folds(
                    renewal_folds, plots / "renewal_group_matched_cv_performance", formats,
                    "Matched leave-year-out nitrate-qualified renewal-phase assessment",
                )
                plot_group_ordination(
                    renewal_distance, renewal_frame, [args.renewal_group_col],
                    plots / "renewal_group_ordination", formats,
                )

        # Consecutive-cruise turnover directly addresses prediction of shifts.
        if args.date_col in cruise_frame.columns and len(cruise_frame) >= 6:
            ordered = cruise_frame.copy()
            ordered["_date"] = pd.to_datetime(ordered[args.date_col], errors="coerce")
            ordered["_response_row"] = np.arange(len(ordered))
            ordered = ordered.dropna(subset=["_date"]).sort_values("_date").reset_index(drop=True)
            transition_rows = []
            for idx in range(1, len(ordered)):
                previous = ordered.iloc[idx - 1]
                current = ordered.iloc[idx]
                community_shift = float(
                    pdist(
                        cruise_response.to_numpy()[
                            [int(previous["_response_row"]), int(current["_response_row"])]
                        ],
                        metric="braycurtis",
                    )[0]
                )
                transition_rows.append({
                    "from_cruise": previous[args.cruise_col],
                    "to_cruise": current[args.cruise_col],
                    "from_date": previous["_date"],
                    "to_date": current["_date"],
                    "gap_days": int((current["_date"] - previous["_date"]).days),
                    args.year_col: current[args.year_col],
                    args.season_col: current[args.season_col],
                    "community_bray_shift": community_shift,
                    "abs_pea_change": abs(float(current[args.pea_col]) - float(previous[args.pea_col])),
                    "abs_centroid_change": abs(float(current[args.centroid_col]) - float(previous[args.centroid_col])),
                })
            transitions = pd.DataFrame(transition_rows)
            transitions.to_csv(tables / "cruise_temporal_turnover_cohort.tsv", sep="\t", index=False)
            turnover_response = transitions[["community_bray_shift"]].to_numpy(dtype=float)
            turnover_models = {
                "Gap + season baseline": (["gap_days", args.season_col], {"gap_days"}),
                "PEA change": (["gap_days", args.season_col, "abs_pea_change"], {"gap_days", "abs_pea_change"}),
                "Centroid change": (["gap_days", args.season_col, "abs_centroid_change"], {"gap_days", "abs_centroid_change"}),
                "PEA + centroid change": (
                    ["gap_days", args.season_col, "abs_pea_change", "abs_centroid_change"],
                    {"gap_days", "abs_pea_change", "abs_centroid_change"},
                ),
            }
            turnover_folds, turnover_cv = multivariate_cv(
                turnover_response, transitions, turnover_models, args.year_col
            )
            turnover_cv.to_csv(tables / "cruise_temporal_turnover_model_comparison.tsv", sep="\t", index=False)
            turnover_folds.to_csv(tables / "cruise_temporal_turnover_cv_folds.tsv", sep="\t", index=False)
            paired_difference(turnover_folds, "Centroid change", "PEA change", args.seed).to_csv(
                tables / "cruise_temporal_centroid_vs_pea_performance.tsv", sep="\t", index=False
            )
            plot_cv_folds(
                turnover_folds, plots / "cruise_temporal_turnover_cv_performance", formats,
                "Prediction of consecutive-cruise Bray–Curtis shifts",
            )
        audit_rows.append({"branch": "cruise", "status": "ok", "input_samples": len(common),
                           "retained_samples": len(meta_cruise), "retained_units": len(cruise_frame),
                           "note": "matched complete PEA and centroid cruises"})
    else:
        audit_rows.append({"branch": "cruise", "status": "skipped_missing_columns", "input_samples": len(common),
                           "retained_samples": 0, "retained_units": 0, "note": ",".join(missing_cruise_columns)})

    # Bottle-level matched categorical-organization comparison.
    sample_columns = [
        args.cruise_col, args.depth_col, args.season_col,
        args.o2_group_col, args.gmm_group_col, args.hybrid_group_col,
    ]
    missing_sample_columns = [c for c in sample_columns if c not in metadata.columns]
    if not missing_sample_columns:
        sample_frame = metadata.dropna(subset=sample_columns).copy()
        if args.assignment_source_col in sample_frame.columns:
            sample_frame = sample_frame[
                sample_frame[args.assignment_source_col].astype(str).eq(args.observed_source_label)
            ]
        group_counts = sample_frame[args.hybrid_group_col].astype(str).value_counts()
        retained_groups = set(group_counts[group_counts >= args.min_group_n].index)
        sample_frame = sample_frame[sample_frame[args.hybrid_group_col].astype(str).isin(retained_groups)].copy()
        sample_counts = counts.loc[sample_frame.index]
        response = hellinger(sample_counts).to_numpy()
        distance = squareform(pdist(response, metric="braycurtis"))
        sample_frame = sample_frame.reset_index(drop=True)
        baseline = [args.depth_col, args.season_col]
        categorical_models = {
            "O2 compartment": [args.o2_group_col],
            "GMM component": [args.gmm_group_col],
            "Hybrid compartment": [args.hybrid_group_col],
        }
        model_rows = []
        for name, tested in categorical_models.items():
            stats = nested_distance_stats(distance, sample_frame, baseline, tested, set())
            stats["p_value"] = permutation_pvalue(
                distance, sample_frame, baseline, tested, set(), args.permutations, args.seed,
                blocks=sample_frame[args.cruise_col],
            )
            dispersion = permdisp(
                distance, sample_frame[tested[0]], args.permutations, args.seed,
                sample_frame[args.cruise_col],
            )
            stats["permdisp_f"] = dispersion["f_statistic"]
            stats["permdisp_p_value"] = dispersion["p_value"]
            stats["model"] = name
            model_rows.append(stats)
        sample_results = pd.DataFrame(model_rows)
        cv_models = {
            "Depth + season baseline": (baseline, set()),
            "O2 compartment": (baseline + [args.o2_group_col], set()),
            "GMM component": (baseline + [args.gmm_group_col], set()),
            "Hybrid compartment": (baseline + [args.hybrid_group_col], set()),
        }
        sample_folds, sample_cv = multivariate_cv(response, sample_frame, cv_models, args.cruise_col)
        sample_results = sample_results.merge(sample_cv, on="model", how="left")
        sample_results.to_csv(tables / "sample_grouping_model_comparison.tsv", sep="\t", index=False)
        sample_folds.to_csv(tables / "sample_grouping_cv_folds.tsv", sep="\t", index=False)
        hybrid_vs_o2 = paired_difference(
            sample_folds, "Hybrid compartment", "O2 compartment", args.seed
        )
        gmm_vs_o2 = paired_difference(
            sample_folds, "GMM component", "O2 compartment", args.seed
        )
        hybrid_vs_gmm = paired_difference(
            sample_folds, "Hybrid compartment", "GMM component", args.seed
        )
        hybrid_vs_o2.to_csv(
            tables / "sample_hybrid_vs_o2_performance.tsv", sep="\t", index=False
        )
        gmm_vs_o2.to_csv(
            tables / "sample_gmm_vs_o2_performance.tsv", sep="\t", index=False
        )
        hybrid_vs_gmm.to_csv(
            tables / "sample_hybrid_vs_gmm_performance.tsv", sep="\t", index=False
        )
        paired_compartments = pd.concat(
            [hybrid_vs_o2, gmm_vs_o2, hybrid_vs_gmm], ignore_index=True
        )
        if not paired_compartments.empty:
            paired_compartments["q_value"] = bh_adjust(
                paired_compartments["sign_flip_p_value"]
            )
        paired_compartments.to_csv(
            tables / "sample_compartment_paired_performance.tsv", sep="\t", index=False
        )

        # These nested tests ask whether the hybrid intersection adds resolution
        # within either parent classification. Permutations are restricted within
        # cruise and parent level so the parent's marginal structure is preserved.
        refinement_rows = []
        paired_lookup = {
            args.o2_group_col: hybrid_vs_o2,
            args.gmm_group_col: hybrid_vs_gmm,
        }
        parent_labels = {
            args.o2_group_col: "Legacy O2",
            args.gmm_group_col: "GMM",
        }
        for parent_col in (args.o2_group_col, args.gmm_group_col):
            stats = nested_distance_stats(
                distance, sample_frame, baseline + [parent_col],
                [args.hybrid_group_col], set(),
            )
            restricted_blocks = (
                sample_frame[args.cruise_col].astype(str)
                + "||" + sample_frame[parent_col].astype(str)
            )
            stats["p_value"] = permutation_pvalue(
                distance, sample_frame, baseline + [parent_col],
                [args.hybrid_group_col], set(), args.permutations, args.seed,
                blocks=restricted_blocks,
            )
            stats["parent_classification"] = parent_labels[parent_col]
            stats["refinement"] = "Hybrid O2-GMM"
            paired = paired_lookup[parent_col]
            if not paired.empty:
                for column in (
                    "n_folds", "mean_difference_a_minus_b", "ci_lower",
                    "ci_upper", "sign_flip_p_value",
                ):
                    stats[f"cv_{column}"] = paired.iloc[0][column]
            refinement_rows.append(stats)
        refinement = pd.DataFrame(refinement_rows)
        refinement["q_value"] = bh_adjust(refinement["p_value"])
        refinement["cv_q_value"] = bh_adjust(refinement["cv_sign_flip_p_value"])
        refinement.to_csv(
            tables / "sample_hybrid_parent_refinement.tsv", sep="\t", index=False
        )
        pd.DataFrame({
            "group": group_counts.index,
            "n": group_counts.values,
            "retained": [group in retained_groups for group in group_counts.index],
            "minimum_n": args.min_group_n,
        }).to_csv(tables / "sample_hybrid_group_audit.tsv", sep="\t", index=False)
        sample_frame[sample_columns + [args.sample_col] if args.sample_col in sample_frame else sample_columns].to_csv(
            tables / "sample_analysis_cohort.tsv", sep="\t", index=False
        )
        plot_model_comparison(sample_results, plots / "sample_grouping_model_performance", formats,
                              "Bottle-level community organization")
        plot_cv_folds(sample_folds, plots / "sample_grouping_cv_performance", formats,
                      "Leave-one-cruise-out community prediction")
        plot_group_ordination(
            distance, sample_frame,
            [args.o2_group_col, args.gmm_group_col, args.hybrid_group_col],
            plots / "sample_grouping_ordination", formats,
        )

        # Independent validation for biological comparisons: within each BASIN
        # cruise-scale state, test whether hybrid compartments organize community
        # composition after depth and season adjustment. This is run separately
        # for environmental cruise groups and nitrate-qualified renewal phase.
        outer_states = [
            (args.cruise_group_col, "cruise_group"),
            (args.renewal_group_col, "renewal_group"),
        ]
        for outer_col, output_prefix in outer_states:
            if outer_col not in metadata.columns:
                continue
            outer_values = metadata.loc[sample_counts.index, outer_col]
            sample_frame[outer_col] = outer_values.to_numpy()
            valid_group = sample_frame[outer_col].notna()
            strat_frame = sample_frame.loc[valid_group].reset_index(drop=True)
            strat_response = response[valid_group.to_numpy()]
            strat_distance = squareform(pdist(strat_response, metric="braycurtis"))
            result_parts, fold_parts, paired_parts, pairwise_rows, audit = [], [], [], [], []
            for group_level in sorted(strat_frame[outer_col].astype(str).unique()):
                in_group = strat_frame[outer_col].astype(str).eq(group_level).to_numpy()
                current = strat_frame.loc[in_group].reset_index(drop=True)
                current_response = strat_response[in_group]
                counts_by_compartment = current[args.hybrid_group_col].astype(str).value_counts()
                cruises_by_compartment = current.groupby(args.hybrid_group_col)[args.cruise_col].nunique()
                retained = counts_by_compartment[counts_by_compartment >= args.min_group_n].index
                current = current[current[args.hybrid_group_col].astype(str).isin(retained)].reset_index(drop=True)
                for compartment, n_samples in counts_by_compartment.items():
                    audit.append({"outer_state_column": outer_col, "outer_state": group_level,
                                  outer_col: group_level,
                                  "compartment": compartment,
                                  "n_samples": int(n_samples),
                                  "n_cruises": int(cruises_by_compartment.get(compartment, 0)),
                                  "retained": bool(compartment in retained), "minimum_n": args.min_group_n})
                if current[args.hybrid_group_col].nunique() < 2 or current[args.cruise_col].nunique() < 3:
                    continue
                kept_rows = np.flatnonzero(in_group)[
                    strat_frame.loc[in_group, args.hybrid_group_col].astype(str).isin(retained).to_numpy()
                ]
                current_response = strat_response[kept_rows]
                current_distance = squareform(pdist(current_response, metric="braycurtis"))
                stats = nested_distance_stats(current_distance, current, baseline,
                                              [args.hybrid_group_col], set())
                stats["p_value"] = permutation_pvalue(
                    current_distance, current, baseline, [args.hybrid_group_col], set(),
                    args.permutations, args.seed, blocks=current[args.cruise_col],
                )
                dispersion = permdisp(current_distance, current[args.hybrid_group_col],
                                      args.permutations, args.seed, current[args.cruise_col])
                stats["permdisp_f"] = dispersion["f_statistic"]
                stats["permdisp_p_value"] = dispersion["p_value"]
                models = {
                    "Depth + season baseline": (baseline, set()),
                    "Hybrid compartment": (baseline + [args.hybrid_group_col], set()),
                }
                folds, summary = multivariate_cv(current_response, current, models, args.cruise_col)
                summary["outer_state_column"] = outer_col
                summary["outer_state"] = group_level
                summary[outer_col] = group_level
                for column, value in stats.items():
                    summary[column] = np.where(summary.model.eq("Hybrid compartment"), value, np.nan)
                result_parts.append(summary)
                folds["outer_state_column"] = outer_col
                folds["outer_state"] = group_level
                folds[outer_col] = group_level
                fold_parts.append(folds)
                paired = paired_difference(folds, "Hybrid compartment", "Depth + season baseline", args.seed)
                if not paired.empty:
                    paired["outer_state_column"] = outer_col
                    paired["outer_state"] = group_level
                    paired[outer_col] = group_level
                    paired_parts.append(paired)
                for level_a, level_b in combinations(sorted(current[args.hybrid_group_col].astype(str).unique()), 2):
                    pair = current[current[args.hybrid_group_col].astype(str).isin([level_a, level_b])].reset_index(drop=True)
                    per_level_cruises = pair.groupby(args.hybrid_group_col)[args.cruise_col].nunique()
                    if len(pair) < 6 or (per_level_cruises < 2).any():
                        continue
                    pair_idx = current[args.hybrid_group_col].astype(str).isin([level_a, level_b]).to_numpy()
                    pair_response = current_response[pair_idx]
                    pair_distance = squareform(pdist(pair_response, metric="braycurtis"))
                    pair_stats = nested_distance_stats(pair_distance, pair, baseline,
                                                       [args.hybrid_group_col], set())
                    pair_stats["p_value"] = permutation_pvalue(
                        pair_distance, pair, baseline, [args.hybrid_group_col], set(),
                        args.permutations, args.seed, blocks=pair[args.cruise_col],
                    )
                    pair_dispersion = permdisp(pair_distance, pair[args.hybrid_group_col],
                                               args.permutations, args.seed, pair[args.cruise_col])
                    pairwise_rows.append({"outer_state_column": outer_col, "outer_state": group_level,
                                          outer_col: group_level,
                                          "compartment_a": level_a, "compartment_b": level_b,
                                          "n_a": int(pair[args.hybrid_group_col].astype(str).eq(level_a).sum()),
                                          "n_b": int(pair[args.hybrid_group_col].astype(str).eq(level_b).sum()),
                                          "cruises_a": int(per_level_cruises.get(level_a, 0)),
                                          "cruises_b": int(per_level_cruises.get(level_b, 0)),
                                          **pair_stats,
                                          "permdisp_f": pair_dispersion["f_statistic"],
                                          "permdisp_p_value": pair_dispersion["p_value"]})
            pd.DataFrame(audit).to_csv(tables / f"{output_prefix}_compartment_audit.tsv", sep="\t", index=False)
            if result_parts:
                pd.concat(result_parts, ignore_index=True).to_csv(
                    tables / f"{output_prefix}_compartment_model_comparison.tsv", sep="\t", index=False)
                pd.concat(fold_parts, ignore_index=True).to_csv(
                    tables / f"{output_prefix}_compartment_cv_folds.tsv", sep="\t", index=False)
            if paired_parts:
                pd.concat(paired_parts, ignore_index=True).to_csv(
                    tables / f"{output_prefix}_compartment_paired_performance.tsv", sep="\t", index=False)
            if pairwise_rows:
                pairwise = pd.DataFrame(pairwise_rows)
                pairwise["q_value"] = pairwise.groupby("outer_state", group_keys=False)["p_value"].apply(bh_adjust)
                pairwise["permdisp_q_value"] = pairwise.groupby("outer_state", group_keys=False)["permdisp_p_value"].apply(bh_adjust)
                pairwise.to_csv(tables / f"{output_prefix}_compartment_pairwise.tsv", sep="\t", index=False)
            if len(strat_frame) > 2:
                plot_compartments_within_cruise_groups(
                    strat_distance, strat_frame, outer_col, args.hybrid_group_col,
                    plots / f"{output_prefix}_compartment_ordination", formats,
                )
        audit_rows.append({"branch": "sample", "status": "ok", "input_samples": len(common),
                           "retained_samples": len(sample_frame), "retained_units": sample_frame[args.cruise_col].nunique(),
                           "note": f"observed labels; hybrid groups n>={args.min_group_n}"})
    else:
        audit_rows.append({"branch": "sample", "status": "skipped_missing_columns", "input_samples": len(common),
                           "retained_samples": 0, "retained_units": 0, "note": ",".join(missing_sample_columns)})

    pd.DataFrame(audit_rows).to_csv(tables / "community_predictor_input_audit.tsv", sep="\t", index=False)
    (tables / "run_config.json").write_text(json.dumps(vars(args), default=str, indent=2) + "\n")


if __name__ == "__main__":
    main()
