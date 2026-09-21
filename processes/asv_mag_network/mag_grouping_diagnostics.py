"""Direct MAG abundance/expression diagnostics for biochemical groupings."""

from __future__ import annotations

import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shared_plot_style import install_publication_style

install_publication_style()

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def _bh(values: pd.Series) -> pd.Series:
    out = pd.Series(np.nan, index=values.index, dtype=float)
    valid = pd.to_numeric(values, errors="coerce").dropna().sort_values()
    if valid.empty:
        return out
    n = len(valid)
    adjusted = (
        valid.mul(n).div(np.arange(1, n + 1))
        .iloc[::-1].cummin().iloc[::-1].clip(upper=1)
    )
    out.loc[adjusted.index] = adjusted
    return out


def _effect(matrix: np.ndarray, labels: np.ndarray) -> float:
    grand = matrix.mean(axis=0)
    total = float(np.square(matrix - grand).sum())
    if total <= 0:
        return 0.0
    between = 0.0
    for level in np.unique(labels):
        current = matrix[labels == level]
        between += len(current) * float(np.square(current.mean(axis=0) - grand).sum())
    return between / total


def _feature_effects(matrix: np.ndarray, labels: np.ndarray) -> np.ndarray:
    grand = matrix.mean(axis=0)
    total = np.square(matrix - grand).sum(axis=0)
    between = np.zeros(matrix.shape[1], dtype=float)
    for level in np.unique(labels):
        current = matrix[labels == level]
        between += len(current) * np.square(current.mean(axis=0) - grand)
    return np.divide(between, total, out=np.zeros_like(between), where=total > 0)


def _permute(labels: np.ndarray, blocks: np.ndarray | None, rng: np.random.Generator) -> np.ndarray:
    if blocks is None:
        return rng.permutation(labels)
    shuffled = labels.copy()
    for block in np.unique(blocks):
        indexes = np.flatnonzero(blocks == block)
        shuffled[indexes] = rng.permutation(shuffled[indexes])
    return shuffled


def _balanced_accuracy(observed: np.ndarray, predicted: np.ndarray) -> float:
    recalls = [
        float(np.mean(predicted[observed == level] == level))
        for level in np.unique(observed)
        if np.any(observed == level)
    ]
    return float(np.mean(recalls)) if recalls else np.nan


def _blocked_nearest_centroid_cv(
    matrix: np.ndarray,
    labels: np.ndarray,
    folds: np.ndarray,
) -> float:
    predictions = np.full(len(labels), "", dtype=object)
    for fold in np.unique(folds):
        test = folds == fold
        train = ~test
        train_levels = np.unique(labels[train])
        if len(train_levels) < 2:
            continue
        centroids = {
            level: matrix[train & (labels == level)].mean(axis=0)
            for level in train_levels
            if np.any(train & (labels == level))
        }
        for index in np.flatnonzero(test):
            predictions[index] = min(
                centroids,
                key=lambda level: float(np.square(matrix[index] - centroids[level]).sum()),
            )
    valid = predictions != ""
    return _balanced_accuracy(labels[valid], predictions[valid]) if valid.any() else np.nan


def _prepare(
    normalized: pd.DataFrame,
    metadata: pd.DataFrame,
    sample_col: str,
    cruise_col: str,
    grouping: str,
    cruise_level_groupings: set[str],
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray, np.ndarray | None, str]:
    samples = [sample for sample in normalized.columns if sample in set(metadata[sample_col].astype(str))]
    meta = metadata.drop_duplicates(sample_col).set_index(sample_col).reindex(samples)
    keep = meta[grouping].notna() & meta[cruise_col].notna()
    meta = meta.loc[keep].copy()
    matrix = np.log1p(normalized.loc[:, meta.index].T.apply(pd.to_numeric, errors="coerce").fillna(0.0))
    if grouping in cruise_level_groupings:
        matrix[cruise_col] = meta[cruise_col].astype(str).values
        matrix = matrix.groupby(cruise_col).mean()
        labels = meta.groupby(meta[cruise_col].astype(str))[grouping].first().reindex(matrix.index).astype(str).to_numpy()
        folds = matrix.index.astype(str).to_numpy()
        blocks = None
        unit = "cruise"
    else:
        labels = meta[grouping].astype(str).to_numpy()
        folds = meta[cruise_col].astype(str).to_numpy()
        blocks = folds
        unit = "sample_within_cruise"
    values = matrix.to_numpy(float)
    scale = values.std(axis=0)
    scale[scale == 0] = 1.0
    values = (values - values.mean(axis=0)) / scale
    # Give each omics block equal total squared weight in joint models so the
    # block with the most features cannot dominate merely by dimensionality.
    block_names = np.asarray([
        str(column).split("::", 1)[0] for column in matrix.columns
    ])
    for block in np.unique(block_names):
        indexes = block_names == block
        values[:, indexes] /= math.sqrt(int(indexes.sum()))
    return matrix, values, labels, blocks, unit


def run_mag_grouping_diagnostics(
    modalities: dict[str, pd.DataFrame],
    metadata: pd.DataFrame,
    outdir: Path,
    prefix: str,
    sample_col: str,
    cruise_col: str,
    groupings: list[str],
    cruise_level_groupings: list[str],
    permutations: int,
    seed: int,
    formats: list[str],
) -> None:
    tables = outdir / "group_diagnostics"
    plots = outdir / "group_diagnostics"
    tables.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    cruise_level_grouping_set = set(cruise_level_groupings)
    overall_rows: list[dict] = []
    mag_rows: list[dict] = []
    for modality, normalized in modalities.items():
        if normalized.empty:
            continue
        for grouping in groupings:
            if grouping not in metadata or sample_col not in metadata or cruise_col not in metadata:
                continue
            raw, matrix, labels, blocks, unit = _prepare(
                normalized,
                metadata,
                sample_col,
                cruise_col,
                grouping,
                cruise_level_grouping_set,
            )
            if len(labels) < 4 or len(np.unique(labels)) < 2:
                continue
            observed = _effect(matrix, labels)
            permuted_labels = [
                _permute(labels, blocks, rng) for _ in range(permutations)
            ]
            hits = sum(
                _effect(matrix, shuffled) >= observed - 1e-15
                for shuffled in permuted_labels
            )
            folds = (
                raw.index.astype(str).to_numpy()
                if unit == "cruise"
                else metadata.drop_duplicates(sample_col).set_index(sample_col)
                .reindex(raw.index)[cruise_col].astype(str).to_numpy()
            )
            cv = _blocked_nearest_centroid_cv(matrix, labels, folds)
            cv_hits = 0
            for shuffled in permuted_labels:
                cv_hits += _blocked_nearest_centroid_cv(matrix, shuffled, folds) >= cv - 1e-15
            feature_effects = _feature_effects(matrix, labels)
            feature_hits = np.zeros(matrix.shape[1], dtype=int)
            for shuffled in permuted_labels:
                feature_hits += (
                    _feature_effects(matrix, shuffled) >= feature_effects - 1e-15
                )
            overall_rows.append({
                "modality": modality,
                "grouping": grouping,
                "analysis_unit": unit,
                "n_units": len(labels),
                "n_groups": len(np.unique(labels)),
                "chance_balanced_accuracy": 1.0 / len(np.unique(labels)),
                "multivariate_effect_r2": observed,
                "multivariate_permutation_p_value": (hits + 1) / (permutations + 1),
                "blocked_nearest_centroid_balanced_accuracy": cv,
                "classification_permutation_p_value": (cv_hits + 1) / (permutations + 1),
            })
            for column_index, feature_id in enumerate(raw.columns):
                effect = float(feature_effects[column_index])
                means = pd.Series(raw.iloc[:, column_index].to_numpy(float)).groupby(labels).mean()
                top = str(means.idxmax())
                other = means.drop(top)
                mag_rows.append({
                    "modality": modality,
                    "grouping": grouping,
                    "feature_id": str(feature_id),
                    "feature_block": str(feature_id).split("::", 1)[0],
                    "effect_eta_squared": effect,
                    "p_value": (
                        int(feature_hits[column_index]) + 1
                    ) / (permutations + 1),
                    "highest_group_level": top,
                    "highest_group_mean_normalized_recruitment": float(means.max()),
                    "highest_to_other_mean_ratio": (
                        float(means.max() / other.mean()) if len(other) and other.mean() > 0 else np.nan
                    ),
                    "n_units": len(labels),
                })
            means = raw.assign(_group=labels).groupby("_group").mean()
            variable = means.var(axis=0).sort_values(ascending=False).head(30).index
            display = means[variable].T
            display = display.sub(display.mean(axis=1), axis=0)
            display = display.div(display.std(axis=1).replace(0, 1), axis=0)
            fig, ax = plt.subplots(figsize=(max(10, 2.0 * len(display.columns)), 14))
            sns.heatmap(display, cmap="Greys", center=0, ax=ax, cbar_kws={"label": "Group-mean z-score"})
            ax.set_xlabel(grouping)
            ax.set_ylabel("Community feature")
            fig.tight_layout()
            stem = plots / f"{prefix}_{modality}_{grouping}_feature_group_means"
            for fmt in formats:
                fig.savefig(stem.with_suffix(f".{fmt}"), dpi=300, bbox_inches="tight")
            plt.close(fig)
    overall = pd.DataFrame(overall_rows)
    if not overall.empty:
        overall["multivariate_q_value"] = _bh(overall["multivariate_permutation_p_value"])
        overall["classification_q_value"] = _bh(overall["classification_permutation_p_value"])
        chance = overall["chance_balanced_accuracy"]
        overall["chance_adjusted_balanced_accuracy"] = (
            overall["blocked_nearest_centroid_balanced_accuracy"] - chance
        ) / (1.0 - chance)
    overall.to_csv(
        tables / f"{prefix}_grouping_diagnostics_overall.tsv",
        sep="\t",
        index=False,
    )
    if not overall.empty:
        display_names = {
            "Depth": "Depth",
            "Season": "Season",
            "o2_compartment": "Legacy O₂",
            "gmm_component": "GMM",
            "o2_subcompartment_final": "Hybrid O₂ + GMM",
            "cruise_group": "Cruise group",
        }
        plot = overall.copy()
        plot["grouping_label"] = plot["grouping"].map(
            lambda value: display_names.get(value, value)
        )
        plot["modality_label"] = plot["modality"].str.replace("_", " ", regex=False)
        modality_order = list(dict.fromkeys(plot["modality_label"]))
        fig, axes = plt.subplots(
            len(modality_order),
            2,
            figsize=(15, max(5, 3.2 * len(modality_order))),
            squeeze=False,
        )
        for row_index, modality in enumerate(modality_order):
            subset = plot.loc[plot["modality_label"].eq(modality)].copy()
            subset = subset.sort_values(
                ["analysis_unit", "multivariate_effect_r2"],
                ascending=[True, False],
            )
            sns.barplot(
                data=subset,
                y="grouping_label",
                x="multivariate_effect_r2",
                color="0.35",
                ax=axes[row_index, 0],
            )
            sns.barplot(
                data=subset,
                y="grouping_label",
                x="chance_adjusted_balanced_accuracy",
                color="0.35",
                ax=axes[row_index, 1],
            )
            axes[row_index, 0].set_ylabel(modality.title())
            axes[row_index, 0].set_xlabel("Explained multivariate variation (R²)")
            axes[row_index, 1].set_ylabel("")
            axes[row_index, 1].set_xlabel("Chance-adjusted balanced accuracy")
            axes[row_index, 1].axvline(0, color="black", linestyle="--", linewidth=1)
            for axis in axes[row_index]:
                axis.set_title("")
        axes[0, 0].set_title("Observed community separation")
        axes[0, 1].set_title("Held-out predictive performance")
        fig.tight_layout()
        stem = plots / f"{prefix}_grouping_diagnostics_comparison"
        for fmt in formats:
            fig.savefig(stem.with_suffix(f".{fmt}"), dpi=300, bbox_inches="tight")
        plt.close(fig)
    per_mag = pd.DataFrame(mag_rows)
    if not per_mag.empty:
        per_mag["q_value"] = per_mag.groupby(
            ["modality", "grouping"], group_keys=False
        )["p_value"].apply(_bh)
    per_mag.to_csv(
        tables / f"{prefix}_grouping_diagnostics_per_feature.tsv",
        sep="\t",
        index=False,
    )
