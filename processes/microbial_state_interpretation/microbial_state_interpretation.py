#!/usr/bin/env python3
"""Interpret ASV-inferred microbial states using external BASIN labels and ASV modules.

Microbial-state inference is deliberately complete before this script runs. Depth,
season, BASIN labels, and ecological modules are used only for post hoc comparison.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import platform
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/aspire_matplotlib")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, is_color_like
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
import scipy
from scipy.interpolate import LinearNDInterpolator
from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree
from sklearn.compose import ColumnTransformer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    adjusted_rand_score,
    f1_score,
    log_loss,
    normalized_mutual_info_score,
)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder
import sklearn
import umap

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shared_plot_style import install_publication_style

install_publication_style()


COMPACT_PLOT_STYLE = {
    "font.size": 10,
    "axes.titlesize": 14,
    "axes.labelsize": 12,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "legend.title_fontsize": 10,
    "figure.titlesize": 15,
}


def finalize_compact_figure(
    fig: plt.Figure,
    panel_axes: list[plt.Axes] | np.ndarray = (),
) -> None:
    """Keep dense multipanel figures readable under the global style contract."""
    fig._aspire_compact_publication_typography = True
    fig._aspire_disable_panel_labels = True
    for index, ax in enumerate(np.asarray(panel_axes, dtype=object).ravel()):
        if index >= 26:
            break
        ax.text(
            -0.08, 1.04, chr(ord("A") + index), transform=ax.transAxes,
            ha="right", va="bottom", fontsize=16, fontweight="bold",
            clip_on=False,
        )


def observed_class_balanced_accuracy(truth: np.ndarray, predicted: np.ndarray) -> float:
    """Average per-class recall across classes represented in the held-out fold."""
    truth = np.asarray(truth, dtype=str)
    predicted = np.asarray(predicted, dtype=str)
    represented = np.unique(truth)
    return float(np.mean([np.mean(predicted[truth == label] == label) for label in represented]))


def read_table(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix.lower() in {".tsv", ".tab", ".txt"} else ","
    return pd.read_csv(path, sep=sep, low_memory=False)


def merge_assignments_with_metadata(
    assignments: pd.DataFrame,
    metadata: pd.DataFrame,
    sample_col: str,
) -> pd.DataFrame:
    """Merge MC assignments without allowing staged duplicate columns to gain suffixes.

    The post hoc metadata channel may already contain the assignment fields because it
    is produced by the MC augmentation stage. Assignment-table values are authoritative;
    duplicated metadata fields must agree wherever both are observed and are then removed
    before the one-to-one merge.
    """
    if sample_col not in assignments.columns or sample_col not in metadata.columns:
        raise ValueError(f"Both assignment and metadata tables must contain {sample_col!r}")
    assignments = assignments.copy()
    metadata = metadata.copy()
    assignments[sample_col] = assignments[sample_col].astype(str).str.strip()
    metadata[sample_col] = metadata[sample_col].astype(str).str.strip()
    if assignments[sample_col].duplicated().any():
        raise ValueError("Microbial-compartment assignment table contains duplicate sample identifiers")
    if metadata[sample_col].duplicated().any():
        raise ValueError("Post hoc metadata contains duplicate sample identifiers")

    overlap = sorted((set(assignments.columns) & set(metadata.columns)) - {sample_col})
    if overlap:
        assignment_index = assignments.set_index(sample_col)
        metadata_index = metadata.set_index(sample_col)
        shared_samples = assignment_index.index.intersection(metadata_index.index)
        for column in overlap:
            left = assignment_index.loc[shared_samples, column]
            right = metadata_index.loc[shared_samples, column]
            observed = left.notna() & right.notna()
            if not observed.any():
                continue
            left_observed = left.loc[observed]
            right_observed = right.loc[observed]
            left_numeric = pd.to_numeric(left_observed, errors="coerce")
            right_numeric = pd.to_numeric(right_observed, errors="coerce")
            if left_numeric.notna().all() and right_numeric.notna().all():
                agrees = np.isclose(
                    left_numeric.to_numpy(float), right_numeric.to_numpy(float),
                    rtol=1e-10, atol=1e-12,
                )
            else:
                agrees = (
                    left_observed.astype(str).str.strip().to_numpy()
                    == right_observed.astype(str).str.strip().to_numpy()
                )
            if not bool(np.all(agrees)):
                bad_samples = left_observed.index[np.flatnonzero(~agrees)[:5]].tolist()
                raise ValueError(
                    f"Conflicting staged values for assignment field {column!r}; "
                    f"example samples: {bad_samples}"
                )
        metadata = metadata.drop(columns=overlap)

    return assignments.merge(metadata, on=sample_col, how="left", validate="one_to_one")


def parse_list(value: str) -> list[str]:
    return [item.strip() for item in re.split(r"[,|]", value or "") if item.strip()]


def parse_palette(value: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in str(value).split(","):
        if "=" not in item:
            continue
        label, color = item.split("=", 1)
        if label.strip() and color.strip():
            result[label.strip()] = color.strip()
    return result


def validate_palette(palette: dict[str, str], labels: list[str], name: str) -> None:
    missing = [label for label in labels if label not in palette]
    invalid = [label for label in labels if label in palette and not is_color_like(palette[label])]
    if missing or invalid:
        raise ValueError(f"{name} palette missing={missing}, invalid={invalid}")


def label_sort(label: str) -> tuple[int, int, str]:
    text = str(label)
    oxygen = text.split("__", 1)[0].lower()
    oxygen_order = {"oxic": 0, "dysoxic": 1, "suboxic": 2, "anoxic": 3}
    gmm = re.search(r"gmm(\d+)", text, flags=re.IGNORECASE)
    mc = re.fullmatch(r"MC(\d+)", text, flags=re.IGNORECASE)
    module = re.fullmatch(r"M(\d+)", text, flags=re.IGNORECASE)
    if mc:
        return 0, int(mc.group(1)), text
    if module:
        return 0, int(module.group(1)), text
    return oxygen_order.get(oxygen, 99), int(gmm.group(1)) if gmm else 999, text


def display_label(label: str) -> str:
    match = re.fullmatch(r"([^_]+)__gmm(\d+)", str(label), flags=re.IGNORECASE)
    return f"{match.group(1)}-GMM{match.group(2)}" if match else str(label).replace("__", "-")


def entropy(values: pd.Series) -> float:
    probabilities = values.value_counts(normalize=True).to_numpy(float)
    return float(-(probabilities * np.log(probabilities)).sum()) if len(probabilities) else np.nan


def conditional_entropy(outcome: pd.Series, predictor: pd.Series) -> float:
    frame = pd.DataFrame({"outcome": outcome.astype(str), "predictor": predictor.astype(str)})
    total = len(frame)
    return float(sum(len(group) / total * entropy(group.outcome) for _, group in frame.groupby("predictor")))


def bias_corrected_cramers_v(a: pd.Series, b: pd.Series) -> float:
    table = pd.crosstab(a.astype(str), b.astype(str)).to_numpy(float)
    n = table.sum()
    if n <= 1 or min(table.shape) < 2:
        return np.nan
    expected = np.outer(table.sum(axis=1), table.sum(axis=0)) / n
    valid = expected > 0
    chi2 = float(np.sum(((table - expected) ** 2)[valid] / expected[valid]))
    phi2 = chi2 / n
    rows, columns = table.shape
    corrected = max(0.0, phi2 - ((columns - 1) * (rows - 1)) / (n - 1))
    corrected_rows = rows - ((rows - 1) ** 2) / (n - 1)
    corrected_columns = columns - ((columns - 1) ** 2) / (n - 1)
    denominator = min(corrected_rows - 1, corrected_columns - 1)
    return float(math.sqrt(corrected / denominator)) if denominator > 0 else np.nan


def majority_mapping_accuracy(outcome: pd.Series, predictor: pd.Series) -> float:
    frame = pd.DataFrame({"outcome": outcome.astype(str), "predictor": predictor.astype(str)})
    correct = sum(group.outcome.value_counts().max() for _, group in frame.groupby("predictor"))
    return float(correct / len(frame)) if len(frame) else np.nan


def association_metrics(outcome: pd.Series, predictor: pd.Series) -> dict[str, float]:
    outcome = outcome.astype(str)
    predictor = predictor.astype(str)
    base_entropy = entropy(outcome)
    conditional = conditional_entropy(outcome, predictor)
    return {
        "adjusted_rand_index": float(adjusted_rand_score(outcome, predictor)),
        "normalized_mutual_information": float(normalized_mutual_info_score(outcome, predictor)),
        "bias_corrected_cramers_v": bias_corrected_cramers_v(outcome, predictor),
        "conditional_entropy_mc_given_environment": conditional,
        "normalized_conditional_entropy_mc_given_environment": (
            conditional / base_entropy if np.isfinite(base_entropy) and base_entropy > 0 else np.nan
        ),
        "majority_mapping_accuracy": majority_mapping_accuracy(outcome, predictor),
    }


def restricted_permutations(
    frame: pd.DataFrame,
    outcome_col: str,
    predictor_col: str,
    block_col: str,
    permutations: int,
    rng: np.random.Generator,
) -> tuple[dict[str, float], pd.DataFrame]:
    observed = association_metrics(frame[outcome_col], frame[predictor_col])
    null_rows: list[dict[str, float]] = []
    block_indices = [group.index.to_numpy() for _, group in frame.groupby(block_col, sort=False)]
    values = frame[predictor_col].astype(str).to_numpy()
    for replicate in range(1, permutations + 1):
        permuted = values.copy()
        for indices in block_indices:
            permuted[indices] = rng.permutation(permuted[indices])
        metrics = association_metrics(frame[outcome_col], pd.Series(permuted, index=frame.index))
        null_rows.append({"permutation": replicate, **metrics})
    null = pd.DataFrame(null_rows)
    for metric in [
        "adjusted_rand_index", "normalized_mutual_information", "bias_corrected_cramers_v",
        "majority_mapping_accuracy",
    ]:
        observed[f"{metric}_permutation_pvalue"] = (
            1 + int(null[metric].ge(observed[metric]).sum())
        ) / (permutations + 1)
    observed["normalized_conditional_entropy_permutation_pvalue"] = (
        1 + int(
            null["normalized_conditional_entropy_mc_given_environment"].le(
                observed["normalized_conditional_entropy_mc_given_environment"]
            ).sum()
        )
    ) / (permutations + 1)
    return observed, null


def normalize_key(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.replace(r"\.0$", "", regex=True)


def categorical_pipeline(columns: list[str], seed: int) -> Pipeline:
    return Pipeline([
        ("encode", ColumnTransformer([
            ("categorical", OneHotEncoder(handle_unknown="ignore"), columns),
        ], remainder="drop")),
        ("model", LogisticRegression(
            max_iter=3000, solver="lbfgs", random_state=seed,
        )),
    ])


def bootstrap_mean_interval(values: np.ndarray, replicates: int, rng: np.random.Generator) -> tuple[float, float]:
    values = values[np.isfinite(values)]
    if not len(values):
        return np.nan, np.nan
    boot = np.empty(replicates, dtype=float)
    for index in range(replicates):
        boot[index] = np.mean(rng.choice(values, size=len(values), replace=True))
    return float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975))


def sign_flip_pvalue(values: np.ndarray, permutations: int, rng: np.random.Generator) -> float:
    values = values[np.isfinite(values)]
    if not len(values):
        return np.nan
    observed = abs(float(np.mean(values)))
    exceed = 0
    for _ in range(permutations):
        signs = rng.choice(np.array([-1.0, 1.0]), size=len(values), replace=True)
        exceed += abs(float(np.mean(values * signs))) >= observed
    return float((exceed + 1) / (permutations + 1))


def held_out_prediction(
    frame: pd.DataFrame,
    outcome_col: str,
    cruise_col: str,
    depth_col: str,
    season_col: str,
    environmental_cols: list[str],
    permutations: int,
    bootstrap_replicates: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    classes = sorted(frame[outcome_col].astype(str).unique(), key=label_sort)
    models: list[tuple[str, list[str]]] = [("depth_season", [depth_col, season_col])]
    models.extend((column, [depth_col, season_col, column]) for column in environmental_cols)
    fold_rows: list[dict[str, object]] = []
    prediction_rows: list[dict[str, object]] = []
    for model_name, predictors in models:
        for cruise in sorted(frame[cruise_col].astype(str).unique()):
            test = frame[cruise_col].astype(str).eq(cruise)
            train = ~test
            train_classes = set(frame.loc[train, outcome_col].astype(str))
            if len(train_classes) < 2 or not set(frame.loc[test, outcome_col].astype(str)).issubset(train_classes):
                fold_rows.append({
                    "model": model_name, "held_out_cruise": cruise,
                    "status": "skipped_missing_training_class", "test_samples_n": int(test.sum()),
                })
                continue
            model = categorical_pipeline(predictors, seed)
            model.fit(frame.loc[train, predictors], frame.loc[train, outcome_col].astype(str))
            predicted = model.predict(frame.loc[test, predictors])
            probabilities_raw = model.predict_proba(frame.loc[test, predictors])
            fitted_classes = list(model.named_steps["model"].classes_)
            probabilities = np.zeros((int(test.sum()), len(classes)), dtype=float)
            for source_index, label in enumerate(fitted_classes):
                probabilities[:, classes.index(label)] = probabilities_raw[:, source_index]
            truth = frame.loc[test, outcome_col].astype(str).to_numpy()
            fold_rows.append({
                "model": model_name, "held_out_cruise": cruise, "status": "completed",
                "test_samples_n": len(truth),
                "balanced_accuracy": observed_class_balanced_accuracy(truth, predicted),
                "macro_f1": f1_score(truth, predicted, labels=classes, average="macro", zero_division=0),
                "log_loss": log_loss(truth, probabilities, labels=classes),
            })
            test_rows = frame.loc[test, ["_row_id", outcome_col, cruise_col]].copy()
            test_rows["model"] = model_name
            test_rows["predicted_microbial_compartment"] = predicted
            for class_index, label in enumerate(classes):
                test_rows[f"probability_{label}"] = probabilities[:, class_index]
            prediction_rows.extend(test_rows.to_dict("records"))
    folds = pd.DataFrame(fold_rows)
    predictions = pd.DataFrame(prediction_rows)
    completed = folds.loc[folds.status.eq("completed")].copy()
    baseline = completed.loc[completed.model.eq("depth_season")].set_index("held_out_cruise")
    rng = np.random.default_rng(seed + 101)
    summary_rows: list[dict[str, object]] = []
    for model_name, _ in models:
        model_folds = completed.loc[completed.model.eq(model_name)].set_index("held_out_cruise")
        shared = baseline.index.intersection(model_folds.index)
        row: dict[str, object] = {
            "model": model_name,
            "held_out_cruises_n": len(model_folds),
            "test_samples_n": int(model_folds.test_samples_n.sum()) if len(model_folds) else 0,
            "mean_balanced_accuracy": model_folds.balanced_accuracy.mean(),
            "mean_macro_f1": model_folds.macro_f1.mean(),
            "mean_log_loss": model_folds.log_loss.mean(),
        }
        for metric in ["balanced_accuracy", "macro_f1", "log_loss"]:
            if model_name == "depth_season":
                delta = np.zeros(len(shared), dtype=float)
            elif metric == "log_loss":
                delta = (
                    baseline.loc[shared, metric] - model_folds.loc[shared, metric]
                ).to_numpy(float)
            else:
                delta = (
                    model_folds.loc[shared, metric] - baseline.loc[shared, metric]
                ).to_numpy(float)
            low, high = bootstrap_mean_interval(delta, bootstrap_replicates, rng)
            row[f"mean_improvement_{metric}"] = float(np.mean(delta)) if len(delta) else np.nan
            row[f"improvement_{metric}_ci95_low"] = low
            row[f"improvement_{metric}_ci95_high"] = high
            row[f"improvement_{metric}_sign_flip_pvalue"] = (
                np.nan if model_name == "depth_season" else sign_flip_pvalue(delta, permutations, rng)
            )
        summary_rows.append(row)
    return pd.DataFrame(summary_rows), folds, predictions


def modal_value(values: pd.Series) -> object:
    cleaned = values.dropna().astype(str).str.strip()
    cleaned = cleaned.loc[cleaned.ne("")]
    if cleaned.empty:
        return np.nan
    counts = cleaned.value_counts()
    return sorted(counts.loc[counts.eq(counts.max())].index)[0]


def boundary_concordance(
    frame: pd.DataFrame,
    outcome_col: str,
    environmental_cols: list[str],
    cruise_col: str,
    depth_col: str,
    permutations: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    depth_centers = frame.groupby([cruise_col, depth_col], as_index=False).agg(
        **{outcome_col: (outcome_col, modal_value)},
        **{column: (column, modal_value) for column in environmental_cols},
        sample_n=("_row_id", "size"),
    )
    rows: list[dict[str, object]] = []
    for cruise, group in depth_centers.groupby(cruise_col, sort=True):
        group = group.sort_values(depth_col).reset_index(drop=True)
        for index in range(len(group) - 1):
            shallow, deep = group.iloc[index], group.iloc[index + 1]
            row: dict[str, object] = {
                "cruise": cruise,
                "shallow_depth": float(shallow[depth_col]),
                "deep_depth": float(deep[depth_col]),
                "midpoint_depth": (float(shallow[depth_col]) + float(deep[depth_col])) / 2,
                "depth_difference": float(deep[depth_col]) - float(shallow[depth_col]),
                "microbial_compartment_shallow": shallow[outcome_col],
                "microbial_compartment_deep": deep[outcome_col],
                "microbial_state_transition": bool(shallow[outcome_col] != deep[outcome_col]),
            }
            for column in environmental_cols:
                row[f"{column}_shallow"] = shallow[column]
                row[f"{column}_deep"] = deep[column]
                row[f"crosses_{column}_boundary"] = bool(shallow[column] != deep[column])
            rows.append(row)
    pairs = pd.DataFrame(rows)

    def summarize(crossing: pd.Series, microbial: pd.Series) -> dict[str, float]:
        crossing = crossing.astype(bool).to_numpy()
        microbial = microbial.astype(bool).to_numpy()
        tp = int(np.sum(crossing & microbial))
        fp = int(np.sum(crossing & ~microbial))
        fn = int(np.sum(~crossing & microbial))
        tn = int(np.sum(~crossing & ~microbial))
        probability_cross = tp / (tp + fp) if tp + fp else np.nan
        probability_within = fn / (fn + tn) if fn + tn else np.nan
        sensitivity = tp / (tp + fn) if tp + fn else np.nan
        precision = tp / (tp + fp) if tp + fp else np.nan
        specificity = tn / (tn + fp) if tn + fp else np.nan
        f1 = 2 * precision * sensitivity / (precision + sensitivity) if precision + sensitivity else np.nan
        jaccard = tp / (tp + fp + fn) if tp + fp + fn else np.nan
        corrected = np.array([tp, fp, fn, tn], dtype=float) + 0.5
        odds_ratio = corrected[0] * corrected[3] / (corrected[1] * corrected[2])
        risk_ratio = (
            (corrected[0] / (corrected[0] + corrected[1])) /
            (corrected[2] / (corrected[2] + corrected[3]))
        )
        return {
            "true_positive": tp, "false_positive": fp, "false_negative": fn, "true_negative": tn,
            "mc_transition_probability_at_boundary": probability_cross,
            "mc_transition_probability_within_compartment": probability_within,
            "probability_difference": probability_cross - probability_within,
            "sensitivity": sensitivity, "precision": precision, "specificity": specificity,
            "f1_score": f1, "jaccard_overlap": jaccard,
            "haldane_corrected_odds_ratio": odds_ratio,
            "haldane_corrected_risk_ratio": risk_ratio,
        }

    rng = np.random.default_rng(seed + 202)
    summary_rows: list[dict[str, object]] = []
    null_rows: list[dict[str, object]] = []
    for column in environmental_cols:
        crossing_col = f"crosses_{column}_boundary"
        observed = summarize(pairs[crossing_col], pairs.microbial_state_transition)
        for replicate in range(1, permutations + 1):
            permuted_parts = []
            for _, group in depth_centers.groupby(cruise_col, sort=False):
                shuffled = group.sort_values(depth_col).copy()
                shuffled[column] = rng.permutation(shuffled[column].astype(str).to_numpy())
                permuted_parts.append(shuffled)
            permuted_centers = pd.concat(permuted_parts, ignore_index=True)
            crossing_values: list[bool] = []
            microbial_values: list[bool] = []
            for cruise, group in permuted_centers.groupby(cruise_col, sort=False):
                group = group.sort_values(depth_col).reset_index(drop=True)
                for index in range(len(group) - 1):
                    crossing_values.append(bool(group.iloc[index][column] != group.iloc[index + 1][column]))
                    microbial_values.append(bool(
                        group.iloc[index][outcome_col] != group.iloc[index + 1][outcome_col]
                    ))
            null_metric = summarize(pd.Series(crossing_values), pd.Series(microbial_values))
            null_rows.append({
                "environmental_variable": column, "permutation": replicate,
                "probability_difference": null_metric["probability_difference"],
                "f1_score": null_metric["f1_score"],
                "jaccard_overlap": null_metric["jaccard_overlap"],
            })
        null = pd.DataFrame([row for row in null_rows if row["environmental_variable"] == column])
        summary_rows.append({
            "environmental_variable": column,
            "adjacent_depth_pairs_n": len(pairs),
            "environmental_boundary_crossings_n": int(pairs[crossing_col].sum()),
            "microbial_state_transitions_n": int(pairs.microbial_state_transition.sum()),
            **observed,
            "probability_difference_permutation_pvalue": (
                1 + int(null.probability_difference.ge(observed["probability_difference"]).sum())
            ) / (permutations + 1),
            "f1_permutation_pvalue": (
                1 + int(null.f1_score.ge(observed["f1_score"]).sum())
            ) / (permutations + 1),
        })
    return pd.DataFrame(summary_rows), pairs, pd.DataFrame(null_rows)


def eta_squared(values: pd.Series, groups: pd.Series) -> tuple[float, float]:
    frame = pd.DataFrame({"value": pd.to_numeric(values, errors="coerce"), "group": groups.astype(str)}).dropna()
    n = len(frame)
    k = frame.group.nunique()
    if n <= k or k < 2:
        return np.nan, np.nan
    overall = frame.value.mean()
    total = float(((frame.value - overall) ** 2).sum())
    between = float(sum(len(group) * (group.value.mean() - overall) ** 2 for _, group in frame.groupby("group")))
    r2 = between / total if total > 0 else 0.0
    adjusted = 1 - (1 - r2) * (n - 1) / (n - k)
    return float(r2), float(adjusted)


def benjamini_hochberg(values: pd.Series) -> pd.Series:
    """Benjamini-Hochberg adjusted p-values, preserving missing entries."""
    numeric = pd.to_numeric(values, errors="coerce")
    valid = numeric.notna()
    adjusted = pd.Series(np.nan, index=values.index, dtype=float)
    if not valid.any():
        return adjusted
    ordered = numeric.loc[valid].sort_values()
    ranks = np.arange(1, len(ordered) + 1, dtype=float)
    corrected = (ordered.to_numpy() * len(ordered) / ranks)
    corrected = np.minimum.accumulate(corrected[::-1])[::-1]
    adjusted.loc[ordered.index] = np.minimum(corrected, 1.0)
    return adjusted


def categorical_design(values: pd.Series) -> np.ndarray:
    labels = values.astype(str)
    levels = sorted(labels.unique())
    if len(levels) < 2:
        return np.empty((len(labels), 0), dtype=float)
    return np.column_stack([labels.eq(level).to_numpy(float) for level in levels[1:]])


def regression_sse(response: np.ndarray, design: np.ndarray) -> tuple[float, int]:
    response = np.asarray(response, dtype=float)
    if response.ndim == 1:
        response = response[:, None]
    coefficients, _, rank, _ = np.linalg.lstsq(design, response, rcond=None)
    residual = response - design @ coefficients
    return float(np.square(residual).sum()), int(rank)


def depth_adjusted_group_test(
    response: np.ndarray,
    frame: pd.DataFrame,
    group_col: str,
    depth_col: str,
    cruise_col: str,
    permutations: int,
    rng: np.random.Generator,
) -> dict[str, object]:
    """Test group contribution beyond categorical depth with within-cruise permutations."""
    response = np.asarray(response, dtype=float)
    if response.ndim == 1:
        response = response[:, None]
    valid = (
        frame[[group_col, depth_col, cruise_col]].notna().all(axis=1).to_numpy()
        & np.isfinite(response).all(axis=1)
    )
    current = frame.loc[valid, [group_col, depth_col, cruise_col]].reset_index(drop=True)
    y = response[valid]
    if len(current) < 3 or current[group_col].nunique() < 2:
        return {
            "samples_n": len(current), "cruises_n": current[cruise_col].nunique(),
            "groups_n": current[group_col].nunique(), "baseline_r_squared": np.nan,
            "total_r_squared": np.nan, "incremental_r_squared": np.nan,
            "partial_r_squared": np.nan, "pseudo_f": np.nan, "permutation_pvalue": np.nan,
        }
    intercept = np.ones((len(current), 1), dtype=float)
    depth = categorical_design(current[depth_col])
    group = categorical_design(current[group_col])
    baseline = np.column_stack([intercept, depth])
    full = np.column_stack([baseline, group])
    sse0, rank0 = regression_sse(y, baseline)
    sse1, rank1 = regression_sse(y, full)
    centered = y - y.mean(axis=0, keepdims=True)
    sst = float(np.square(centered).sum())
    incremental = max(0.0, (sse0 - sse1) / sst) if sst > 0 else 0.0
    baseline_r2 = 1.0 - sse0 / sst if sst > 0 else 0.0
    total_r2 = 1.0 - sse1 / sst if sst > 0 else 0.0
    partial = max(0.0, (sse0 - sse1) / sse0) if sse0 > 0 else 0.0
    numerator_df = max(1, rank1 - rank0)
    denominator_df = max(1, len(current) - rank1)
    pseudo_f = ((sse0 - sse1) / numerator_df) / (sse1 / denominator_df) if sse1 > 0 else np.inf
    null_incremental = []
    original = current[group_col].astype(str).to_numpy()
    cruises = current[cruise_col].astype(str).to_numpy()
    for _ in range(permutations):
        shuffled = original.copy()
        for cruise in np.unique(cruises):
            indexes = np.flatnonzero(cruises == cruise)
            shuffled[indexes] = rng.permutation(shuffled[indexes])
        perm_group = categorical_design(pd.Series(shuffled))
        perm_full = np.column_stack([baseline, perm_group])
        perm_sse, _ = regression_sse(y, perm_full)
        null_incremental.append(max(0.0, (sse0 - perm_sse) / sst) if sst > 0 else 0.0)
    pvalue = (1 + int(np.sum(np.asarray(null_incremental) >= incremental - 1e-15))) / (permutations + 1)
    return {
        "samples_n": len(current), "cruises_n": current[cruise_col].nunique(),
        "groups_n": current[group_col].nunique(), "baseline_r_squared": baseline_r2,
        "total_r_squared": total_r2, "incremental_r_squared": incremental,
        "partial_r_squared": partial, "pseudo_f": pseudo_f, "permutation_pvalue": pvalue,
    }


def cruise_label_depth_adjusted_test(
    response: np.ndarray,
    frame: pd.DataFrame,
    group_col: str,
    depth_col: str,
    cruise_col: str,
    year_col: str,
    permutations: int,
    rng: np.random.Generator,
) -> dict[str, object]:
    """Test a cruise-level state beyond depth using year-restricted cruise permutations."""
    response = np.asarray(response, dtype=float)
    if response.ndim == 1:
        response = response[:, None]
    valid = (
        frame[[group_col, depth_col, cruise_col, year_col]].notna().all(axis=1).to_numpy()
        & np.isfinite(response).all(axis=1)
    )
    current = frame.loc[valid, [group_col, depth_col, cruise_col, year_col]].reset_index(drop=True)
    y = response[valid]
    cruise_label_counts = current.groupby(cruise_col)[group_col].nunique()
    if cruise_label_counts.gt(1).any():
        raise ValueError(f"{group_col} is not constant within cruise for renewal-state testing")
    if len(current) < 3 or current[group_col].nunique() < 2:
        return {
            "samples_n": len(current), "cruises_n": current[cruise_col].nunique(),
            "groups_n": current[group_col].nunique(), "baseline_r_squared": np.nan,
            "total_r_squared": np.nan, "incremental_r_squared": np.nan,
            "partial_r_squared": np.nan, "pseudo_f": np.nan,
            "permutation_pvalue": np.nan, "nonidentity_permutations_n": 0,
            "permutation_scheme": "cruise_labels_within_calendar_year",
        }
    intercept = np.ones((len(current), 1), dtype=float)
    baseline = np.column_stack([intercept, categorical_design(current[depth_col])])
    full = np.column_stack([baseline, categorical_design(current[group_col])])
    sse0, rank0 = regression_sse(y, baseline)
    sse1, rank1 = regression_sse(y, full)
    centered = y - y.mean(axis=0, keepdims=True)
    sst = float(np.square(centered).sum())
    incremental = max(0.0, (sse0 - sse1) / sst) if sst > 0 else 0.0
    baseline_r2 = 1.0 - sse0 / sst if sst > 0 else 0.0
    total_r2 = 1.0 - sse1 / sst if sst > 0 else 0.0
    partial = max(0.0, (sse0 - sse1) / sse0) if sse0 > 0 else 0.0
    numerator_df = max(1, rank1 - rank0)
    denominator_df = max(1, len(current) - rank1)
    pseudo_f = ((sse0 - sse1) / numerator_df) / (sse1 / denominator_df) if sse1 > 0 else np.inf

    cruise_table = current[[cruise_col, year_col, group_col]].drop_duplicates().reset_index(drop=True)
    original_map = dict(zip(cruise_table[cruise_col].astype(str), cruise_table[group_col].astype(str)))
    sample_cruises = current[cruise_col].astype(str)
    null_incremental: list[float] = []
    nonidentity = 0
    for _ in range(permutations):
        permuted = cruise_table.copy()
        for _, indexes in permuted.groupby(year_col, sort=False).groups.items():
            values = permuted.loc[indexes, group_col].astype(str).to_numpy()
            permuted.loc[indexes, group_col] = rng.permutation(values)
        perm_map = dict(zip(permuted[cruise_col].astype(str), permuted[group_col].astype(str)))
        shuffled = sample_cruises.map(perm_map)
        if any(perm_map[key] != original_map[key] for key in original_map):
            nonidentity += 1
        perm_full = np.column_stack([baseline, categorical_design(shuffled)])
        perm_sse, _ = regression_sse(y, perm_full)
        null_incremental.append(max(0.0, (sse0 - perm_sse) / sst) if sst > 0 else 0.0)
    pvalue = (1 + int(np.sum(np.asarray(null_incremental) >= incremental - 1e-15))) / (permutations + 1)
    return {
        "samples_n": len(current), "cruises_n": current[cruise_col].nunique(),
        "groups_n": current[group_col].nunique(), "baseline_r_squared": baseline_r2,
        "total_r_squared": total_r2, "incremental_r_squared": incremental,
        "partial_r_squared": partial, "pseudo_f": pseudo_f,
        "permutation_pvalue": pvalue, "nonidentity_permutations_n": nonidentity,
        "permutation_scheme": "cruise_labels_within_calendar_year",
    }


def module_characterization(
    relative: pd.DataFrame,
    assignments: pd.DataFrame,
    modules: pd.DataFrame,
    sample_col: str,
    environmental_cols: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    taxon_col = next((column for column in ["Taxon", "ASV_ID", "ASV", "taxon"] if column in modules.columns), modules.columns[0])
    module_col = next((column for column in ["module_label", "module", "Module"] if column in modules.columns), None)
    if module_col is None:
        raise ValueError("Ecological-module table has no module label column")
    mapping = modules[[taxon_col, module_col]].dropna().copy()
    mapping[taxon_col] = mapping[taxon_col].astype(str).str.strip()
    mapping[module_col] = mapping[module_col].astype(str).str.strip()
    available = mapping.loc[mapping[taxon_col].isin(relative.columns)].copy()
    if available.empty:
        raise ValueError("No ecological-module ASVs were present in the relative-abundance matrix")
    module_abundance = pd.DataFrame(index=relative.index)
    for module, group in available.groupby(module_col, sort=True):
        module_abundance[module] = relative[group[taxon_col].tolist()].sum(axis=1)
    module_z = (module_abundance - module_abundance.mean()) / module_abundance.std(ddof=0).replace(0, np.nan)
    module_z = module_z.fillna(0.0)
    joined = assignments.set_index(sample_col).join(module_abundance).join(
        module_z.add_suffix("__z"), how="inner"
    )
    profile_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    for module in module_abundance.columns:
        r2, adjusted = eta_squared(joined[module], joined.microbial_compartment)
        state_means = joined.groupby("microbial_compartment")[f"{module}__z"].mean()
        summary_rows.append({
            "module_label": module,
            "module_asv_n": int((available[module_col] == module).sum()),
            "samples_n": len(joined),
            "eta_squared": r2,
            "adjusted_r_squared": adjusted,
            "highest_mean_state": state_means.idxmax(),
            "highest_mean_standardized_abundance": state_means.max(),
        })
        for state, group in joined.groupby("microbial_compartment", sort=True):
            profile_rows.append({
                "module_label": module, "microbial_compartment": state,
                "samples_n": len(group),
                "mean_module_relative_abundance": group[module].mean(),
                "median_module_relative_abundance": group[module].median(),
                "module_relative_abundance_q25": group[module].quantile(0.25),
                "module_relative_abundance_q75": group[module].quantile(0.75),
                "mean_standardized_module_abundance": group[f"{module}__z"].mean(),
                "median_standardized_module_abundance": group[f"{module}__z"].median(),
            })
    basin_rows: list[dict[str, object]] = []
    for environmental in environmental_cols:
        usable = joined.loc[joined[environmental].notna() & joined[environmental].astype(str).str.strip().ne("")]
        for module in module_abundance.columns:
            r2, adjusted = eta_squared(usable[module], usable[environmental])
            basin_rows.append({
                "environmental_variable": environmental, "module_label": module,
                "samples_n": len(usable), "eta_squared": r2, "adjusted_r_squared": adjusted,
            })
    basin_stats = pd.DataFrame(basin_rows)
    wins = []
    for module, group in basin_stats.groupby("module_label", sort=True):
        best = group.adjusted_r_squared.max()
        for _, row in group.loc[np.isclose(group.adjusted_r_squared, best, equal_nan=False)].iterrows():
            wins.append({"module_label": module, "winning_environmental_variable": row.environmental_variable})
    wins_table = pd.DataFrame(wins)
    return pd.DataFrame(profile_rows), pd.DataFrame(summary_rows), basin_stats, wins_table


def module_statistical_support(
    relative: pd.DataFrame,
    joined: pd.DataFrame,
    modules: pd.DataFrame,
    args: argparse.Namespace,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Depth-adjusted module associations and a single descriptive home MC per module."""
    taxon_col = next((column for column in ["Taxon", "ASV_ID", "ASV", "taxon"] if column in modules.columns), modules.columns[0])
    module_col = next((column for column in ["module_label", "module", "Module"] if column in modules.columns), None)
    if module_col is None:
        raise ValueError("Ecological-module table has no module label column")
    mapping = modules[[taxon_col, module_col]].dropna().copy()
    mapping[taxon_col] = mapping[taxon_col].astype(str).str.strip()
    mapping[module_col] = mapping[module_col].astype(str).str.strip()
    mapping = mapping[mapping[taxon_col].isin(relative.columns)]
    module_abundance = pd.DataFrame(index=relative.index)
    for label, group in mapping.groupby(module_col, sort=True):
        module_abundance[label] = relative[group[taxon_col].tolist()].sum(axis=1)
    sample_ids = [sample for sample in joined[args.sample_col].astype(str) if sample in module_abundance.index]
    frame = joined.set_index(args.sample_col).loc[sample_ids].reset_index()
    module_abundance = module_abundance.loc[sample_ids]
    rng = np.random.default_rng(args.seed + 211)
    rows: list[dict[str, object]] = []
    for module in module_abundance.columns:
        response = module_abundance[module].to_numpy(float)
        common_mask = frame[[args.hybrid_col, "microbial_compartment"]].notna().all(axis=1).to_numpy()
        cohorts = [
            ("all_MC_samples", frame, response, ["microbial_compartment"]),
            (
                "common_hybrid_MC_samples", frame.loc[common_mask].reset_index(drop=True),
                response[common_mask], [args.hybrid_col, "microbial_compartment"],
            ),
        ]
        for cohort, cohort_frame, cohort_response, strategies in cohorts:
            for strategy in strategies:
                rows.append({
                    "module_label": module,
                    "analysis_cohort": cohort,
                    "grouping_strategy": strategy,
                    "module_asv_n": int((mapping[module_col] == module).sum()),
                    **depth_adjusted_group_test(
                        cohort_response, cohort_frame, strategy, args.depth_col, args.cruise_col,
                        args.permutations, rng,
                    ),
                })
    statistics = pd.DataFrame(rows)
    statistics["adjusted_pvalue"] = statistics.groupby(
        ["analysis_cohort", "grouping_strategy"], group_keys=False
    )["permutation_pvalue"].apply(benjamini_hochberg)

    homes: list[dict[str, object]] = []
    mc = frame["microbial_compartment"].astype(str)
    for module in module_abundance.columns:
        values = module_abundance[module].reset_index(drop=True)
        means = values.groupby(mc.reset_index(drop=True)).mean().sort_values(ascending=False)
        medians = values.groupby(mc.reset_index(drop=True)).median()
        row = statistics.loc[
            statistics.module_label.eq(module)
            & statistics.grouping_strategy.eq("microbial_compartment")
            & statistics.analysis_cohort.eq("all_MC_samples")
        ].iloc[0]
        homes.append({
            "module_label": module,
            "home_microbial_compartment": means.index[0],
            "home_mean_module_relative_abundance": means.iloc[0],
            "home_median_module_relative_abundance": medians.loc[means.index[0]],
            "second_highest_mean_module_relative_abundance": means.iloc[1] if len(means) > 1 else np.nan,
            "home_to_second_mean_ratio": means.iloc[0] / means.iloc[1] if len(means) > 1 and means.iloc[1] > 0 else np.inf,
            "depth_adjusted_incremental_r_squared": row.incremental_r_squared,
            "depth_adjusted_partial_r_squared": row.partial_r_squared,
            "permutation_pvalue": row.permutation_pvalue,
            "adjusted_pvalue": row.adjusted_pvalue,
            "home_assignment_role": "descriptive_module_profile_within_MC",
        })
    return statistics, pd.DataFrame(homes)


def renewal_state_support(
    relative: pd.DataFrame,
    joined: pd.DataFrame,
    modules: pd.DataFrame,
    args: argparse.Namespace,
    tables: Path,
    plots: Path,
) -> None:
    """Associate cruise-level renewal phase with MC composition and module abundance."""
    if args.renewal_col not in joined.columns:
        pd.DataFrame([{
            "status": "skipped", "reason": f"metadata column not found: {args.renewal_col}",
        }]).to_csv(tables / "renewal_state_analysis_audit.tsv", sep="\t", index=False)
        return
    allowed = parse_list(args.renewal_levels)
    frame = joined.loc[joined[args.renewal_col].isin(allowed)].copy()
    frame = frame.loc[frame["microbial_compartment"].notna()].reset_index(drop=True)
    if frame.empty:
        pd.DataFrame([{"status": "skipped", "reason": "no eligible renewal-phase samples"}]).to_csv(
            tables / "renewal_state_analysis_audit.tsv", sep="\t", index=False
        )
        return
    if args.year_col not in frame.columns:
        frame[args.year_col] = pd.to_datetime(frame[args.date_col], errors="coerce").dt.year
    phase_counts = frame.groupby(args.renewal_col, dropna=False).agg(
        samples_n=(args.sample_col, "size"), cruises_n=(args.cruise_col, "nunique")
    ).reset_index().rename(columns={args.renewal_col: "renewal_phase"})
    phase_counts.to_csv(tables / "renewal_state_phase_counts.tsv", sep="\t", index=False)

    definitions = [
        ("three_phase_omnibus", {"baseline": "baseline", "renewal": "renewal", "post-renewal": "post-renewal"}),
        ("active_vs_baseline", {"baseline": "baseline", "renewal": "active", "post-renewal": "active"}),
        ("renewal_vs_baseline", {"baseline": "baseline", "renewal": "renewal"}),
        ("post_renewal_vs_baseline", {"baseline": "baseline", "post-renewal": "post-renewal"}),
    ]
    rng = np.random.default_rng(args.seed + 307)
    mc_rows: list[dict[str, object]] = []
    module_rows: list[dict[str, object]] = []
    profile_rows: list[dict[str, object]] = []
    module_col = next((column for column in ["module_label", "module", "Module"] if column in modules.columns), None)
    taxon_col = next((column for column in ["Taxon", "ASV_ID", "ASV", "taxon"] if column in modules.columns), modules.columns[0])
    if module_col is None:
        raise ValueError("Ecological-module table has no module label column")
    mapping = modules[[taxon_col, module_col]].dropna().copy()
    mapping[taxon_col] = mapping[taxon_col].astype(str).str.strip()
    mapping = mapping[mapping[taxon_col].isin(relative.columns)]
    relative_ids = [sample for sample in frame[args.sample_col].astype(str) if sample in relative.index]
    frame = frame.set_index(args.sample_col).loc[relative_ids].reset_index()
    relative_frame = relative.loc[relative_ids]
    module_abundance = pd.DataFrame(index=relative_frame.index)
    for module, group in mapping.groupby(module_col, sort=True):
        module_abundance[str(module)] = relative_frame[group[taxon_col].tolist()].sum(axis=1)

    for analysis, recode in definitions:
        current = frame.loc[frame[args.renewal_col].isin(recode)].copy()
        current["_renewal_test_group"] = current[args.renewal_col].map(recode)
        response_levels = sorted(current.microbial_compartment.astype(str).unique(), key=natural_key)
        mc_response = np.column_stack([
            current.microbial_compartment.astype(str).eq(level).to_numpy(float)
            for level in response_levels
        ])
        test = cruise_label_depth_adjusted_test(
            mc_response, current, "_renewal_test_group", args.depth_col,
            args.cruise_col, args.year_col, args.permutations, rng,
        )
        mc_rows.append({
            "analysis": analysis,
            "renewal_groups": ";".join(sorted(current._renewal_test_group.unique())),
            "adjusted_rand_index": adjusted_rand_score(
                current.microbial_compartment.astype(str), current._renewal_test_group.astype(str)
            ),
            "normalized_mutual_information": normalized_mutual_info_score(
                current.microbial_compartment.astype(str), current._renewal_test_group.astype(str)
            ),
            "bias_corrected_cramers_v": bias_corrected_cramers_v(
                current.microbial_compartment.astype(str), current._renewal_test_group.astype(str)
            ),
            "inference_note": "renewal has only two cruises in the configured ASV cohort; onset-specific results are exploratory"
            if analysis == "renewal_vs_baseline" else "cruise-level phase tested beyond categorical depth",
            **test,
        })
        current_positions = frame.index[frame.index.isin(current.index)].to_numpy()
        for module in module_abundance.columns:
            response = module_abundance.iloc[current_positions][module].to_numpy(float)
            module_rows.append({
                "analysis": analysis, "module_label": module,
                "module_asv_n": int((mapping[module_col].astype(str) == module).sum()),
                **cruise_label_depth_adjusted_test(
                    response, current, "_renewal_test_group", args.depth_col,
                    args.cruise_col, args.year_col, args.permutations, rng,
                ),
            })

    mc_statistics = pd.DataFrame(mc_rows)
    mc_statistics["adjusted_pvalue"] = benjamini_hochberg(mc_statistics.permutation_pvalue)
    mc_statistics.to_csv(tables / "microbial_compartment_renewal_association.tsv", sep="\t", index=False)
    module_statistics = pd.DataFrame(module_rows)
    module_statistics["adjusted_pvalue"] = module_statistics.groupby(
        "analysis", group_keys=False
    ).permutation_pvalue.apply(benjamini_hochberg)
    module_statistics.to_csv(tables / "ecological_module_renewal_association.tsv", sep="\t", index=False)

    module_values = module_abundance.reset_index(drop=True)
    for module in module_values.columns:
        for phase, indexes in frame.groupby(args.renewal_col, sort=True).groups.items():
            values = module_values.loc[indexes, module]
            profile_rows.append({
                "module_label": module, "renewal_phase": phase,
                "samples_n": len(values), "cruises_n": frame.loc[indexes, args.cruise_col].nunique(),
                "mean_module_relative_abundance": values.mean(),
                "median_module_relative_abundance": values.median(),
                "q25": values.quantile(0.25), "q75": values.quantile(0.75),
            })
    profiles = pd.DataFrame(profile_rows)
    profiles.to_csv(tables / "ecological_module_renewal_profiles.tsv", sep="\t", index=False)

    cruise_mc = frame.groupby([args.cruise_col, args.renewal_col, "microbial_compartment"]).size().rename("n").reset_index()
    cruise_totals = cruise_mc.groupby(args.cruise_col).n.transform("sum")
    cruise_mc["within_cruise_fraction"] = cruise_mc.n / cruise_totals
    cruise_mc.to_csv(tables / "microbial_compartment_renewal_cruise_composition.tsv", sep="\t", index=False)
    phase_mc = cruise_mc.groupby([args.renewal_col, "microbial_compartment"]).within_cruise_fraction.median().unstack(fill_value=0)
    phase_mc = phase_mc.reindex(index=allowed).dropna(how="all")
    phase_modules = profiles.pivot(index="renewal_phase", columns="module_label", values="median_module_relative_abundance").reindex(phase_mc.index)
    standardized_modules = phase_modules.sub(phase_modules.mean(axis=0), axis=1).div(
        phase_modules.std(axis=0, ddof=0).replace(0, np.nan), axis=1
    ).fillna(0.0)
    with plt.rc_context(COMPACT_PLOT_STYLE):
        fig, axes = plt.subplots(
            1, 2, figsize=(14.2, 6.2),
            gridspec_kw={"width_ratios": [1.0, 1.45]},
        )
        image0 = axes[0].imshow(
            phase_mc.to_numpy(float), aspect="auto", cmap="Blues", vmin=0, vmax=1
        )
        axes[0].set_xticks(np.arange(len(phase_mc.columns)), phase_mc.columns, rotation=0)
        axes[0].set_yticks(np.arange(len(phase_mc.index)), phase_mc.index)
        axes[0].set_title("Median within-cruise MC composition", pad=10)
        image1 = axes[1].imshow(
            standardized_modules.to_numpy(float), aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2
        )
        axes[1].set_xticks(
            np.arange(len(standardized_modules.columns)),
            standardized_modules.columns,
            rotation=90,
        )
        axes[1].set_yticks(np.arange(len(standardized_modules.index)), standardized_modules.index)
        axes[1].set_title("Ecological-module abundance by renewal phase", pad=10)
        fig.subplots_adjust(left=0.10, right=0.98, top=0.88, bottom=0.30, wspace=0.26)
        left_box = axes[0].get_position()
        right_box = axes[1].get_position()
        cax0 = fig.add_axes([left_box.x0, 0.13, left_box.width, 0.025])
        cax1 = fig.add_axes([right_box.x0, 0.13, right_box.width, 0.025])
        fig.colorbar(image0, cax=cax0, orientation="horizontal", label="Median within-cruise fraction")
        fig.colorbar(
            image1, cax=cax1, orientation="horizontal",
            label="Phase median (module-wise z score)",
        )
        finalize_compact_figure(fig, axes)
    save(fig, plots / "renewal_state_mc_and_ecological_modules", parse_list(args.formats))

    pd.DataFrame([{
        "status": "completed", "eligible_samples_n": len(frame),
        "eligible_cruises_n": frame[args.cruise_col].nunique(),
        "renewal_onset_cruises_n": frame.loc[frame[args.renewal_col].eq("renewal"), args.cruise_col].nunique(),
        "interpretation_constraint": "renewal-onset inference is exploratory when represented by fewer than three cruises",
    }]).to_csv(tables / "renewal_state_analysis_audit.tsv", sep="\t", index=False)


def normalized_plot_label(value: object) -> str:
    """Use stable categorical labels while avoiding depth labels such as 10.0."""
    if pd.isna(value):
        return ""
    text = str(value).strip()
    try:
        number = float(text)
    except ValueError:
        return text
    return str(int(number)) if number.is_integer() else f"{number:g}"


def natural_key(value: object) -> list[object]:
    return [int(item) if item.isdigit() else item.casefold() for item in re.split(r"(\d+)", str(value))]


def write_mc_evidence_summary(
    microbial_dir: Path,
    correspondence: pd.DataFrame,
    prediction: pd.DataFrame,
    boundaries: pd.DataFrame,
    tables: Path,
) -> None:
    """Collect the decision-facing MC evidence without converting it into a score."""
    selection = read_table(microbial_dir / "tables" / "microbial_cluster_selection_decision.tsv").iloc[0]
    hybrid = correspondence.loc[
        correspondence.environmental_variable.eq("o2_subcompartment_final")
    ]
    hybrid_prediction = prediction.loc[prediction.model.eq("o2_subcompartment_final")]
    hybrid_boundary = boundaries.loc[
        boundaries.environmental_variable.eq("o2_subcompartment_final")
    ]
    module_stats = read_table(tables / "microbial_state_module_association_statistics.tsv")
    measurement_stats = read_table(tables / "microbial_state_measurement_association.tsv")
    composition_stats = read_table(tables / "microbial_state_composition_association.tsv")

    rows: list[dict[str, object]] = []

    def add(category: str, metric: str, value: object, role: str, source: str) -> None:
        rows.append({
            "evidence_category": category,
            "metric": metric,
            "value": value,
            "interpretation_role": role,
            "source_table": source,
        })

    selection_source = "microbial_compartments/tables/microbial_cluster_selection_decision.tsv"
    for field in [
        "support_status", "selected_K", "selected_mean_silhouette",
        "selected_median_subsampling_ARI", "selected_primary_stability_quantile",
        "selected_minimum_cluster_median_jaccard",
        "selected_median_season_balanced_block_ARI", "selected_median_prediction_strength",
    ]:
        add("MC model support", field, selection.get(field, np.nan),
            "primary_unsupervised_selection_evidence" if field != "selected_mean_silhouette" else "descriptive_cluster_separation_only",
            selection_source)

    if not hybrid.empty:
        row = hybrid.iloc[0]
        source = "microbial_state_interpretation/tables/microbial_basin_correspondence.tsv"
        for field in [
            "adjusted_rand_index", "normalized_mutual_information", "bias_corrected_cramers_v",
            "majority_mapping_accuracy", "adjusted_rand_index_permutation_pvalue",
            "normalized_mutual_information_permutation_pvalue",
        ]:
            add("MC-hybrid correspondence", field, row.get(field, np.nan),
                "independent_environmental_alignment", source)

    if not hybrid_prediction.empty:
        row = hybrid_prediction.iloc[0]
        source = "microbial_state_interpretation/tables/microbial_state_prediction_summary.tsv"
        for field in [
            "held_out_cruises_n", "test_samples_n", "mean_balanced_accuracy",
            "mean_improvement_balanced_accuracy", "improvement_balanced_accuracy_ci95_low",
            "improvement_balanced_accuracy_ci95_high", "improvement_balanced_accuracy_sign_flip_pvalue",
            "mean_improvement_macro_f1", "mean_improvement_log_loss",
        ]:
            add("Held-out MC prediction from hybrid", field, row.get(field, np.nan),
                "leave_one_cruise_out_external_validation", source)

    if not hybrid_boundary.empty:
        row = hybrid_boundary.iloc[0]
        source = "microbial_state_interpretation/tables/microbial_boundary_concordance.tsv"
        for field in [
            "probability_difference", "sensitivity", "precision", "f1_score", "jaccard_overlap",
            "probability_difference_permutation_pvalue", "f1_permutation_pvalue",
        ]:
            add("Vertical MC-hybrid boundary concordance", field, row.get(field, np.nan),
                "independent_depth_profile_alignment", source)

    for strategy in ["microbial_compartment", "o2_subcompartment_final"]:
        subset = module_stats.loc[
            module_stats.grouping_strategy.eq(strategy)
            & module_stats.analysis_cohort.eq("common_hybrid_MC_samples")
        ]
        add("Ecological-module organization", f"{strategy}_common_cohort_tested_modules_n", len(subset),
            "module_level_statistical_support", "microbial_state_interpretation/tables/microbial_state_module_association_statistics.tsv")
        add("Ecological-module organization", f"{strategy}_common_cohort_q_le_0.05_modules_n",
            int(pd.to_numeric(subset.adjusted_pvalue, errors="coerce").le(0.05).sum()),
            "module_level_statistical_support", "microbial_state_interpretation/tables/microbial_state_module_association_statistics.tsv")
        add("Ecological-module organization", f"{strategy}_common_cohort_median_depth_adjusted_incremental_r_squared",
            pd.to_numeric(subset.incremental_r_squared, errors="coerce").median(),
            "module_level_statistical_support", "microbial_state_interpretation/tables/microbial_state_module_association_statistics.tsv")

    mc_all = module_stats.loc[
        module_stats.grouping_strategy.eq("microbial_compartment")
        & module_stats.analysis_cohort.eq("all_MC_samples")
    ]
    add("Ecological-module organization", "microbial_compartment_all_samples_q_le_0.05_modules_n",
        int(pd.to_numeric(mc_all.adjusted_pvalue, errors="coerce").le(0.05).sum()),
        "primary_module_housing_support", "microbial_state_interpretation/tables/microbial_state_module_association_statistics.tsv")
    module_comparison = read_table(tables / "ecological_module_mc_hybrid_comparison.tsv")
    add("Ecological-module organization", "modules_with_larger_MC_incremental_r_squared_n",
        int(module_comparison.larger_depth_adjusted_incremental_r_squared.eq("microbial_compartment").sum()),
        "common_cohort_descriptive_comparison", "microbial_state_interpretation/tables/ecological_module_mc_hybrid_comparison.tsv")
    add("Ecological-module organization", "modules_with_larger_hybrid_incremental_r_squared_n",
        int(module_comparison.larger_depth_adjusted_incremental_r_squared.eq("o2_subcompartment_final").sum()),
        "common_cohort_descriptive_comparison", "microbial_state_interpretation/tables/ecological_module_mc_hybrid_comparison.tsv")

    for strategy in ["microbial_compartment", "o2_subcompartment_final"]:
        subset = measurement_stats.loc[
            measurement_stats.grouping_strategy.eq(strategy)
            & measurement_stats.analysis_cohort.eq("common_hybrid_MC_samples")
        ]
        add("Physicochemical interpretation", f"{strategy}_common_cohort_tested_measurements_n", len(subset),
            "independent_environmental_interpretation", "microbial_state_interpretation/tables/microbial_state_measurement_association.tsv")
        add("Physicochemical interpretation", f"{strategy}_common_cohort_q_le_0.05_measurements_n",
            int(pd.to_numeric(subset.adjusted_pvalue, errors="coerce").le(0.05).sum()),
            "independent_environmental_interpretation", "microbial_state_interpretation/tables/microbial_state_measurement_association.tsv")

    mc_measurement_all = measurement_stats.loc[
        measurement_stats.grouping_strategy.eq("microbial_compartment")
        & measurement_stats.analysis_cohort.eq("all_MC_samples")
    ]
    add("Physicochemical interpretation", "microbial_compartment_all_samples_q_le_0.05_measurements_n",
        int(pd.to_numeric(mc_measurement_all.adjusted_pvalue, errors="coerce").le(0.05).sum()),
        "primary_MC_environmental_interpretation",
        "microbial_state_interpretation/tables/microbial_state_measurement_association.tsv")
    measurement_comparison = read_table(tables / "microbial_state_measurement_mc_hybrid_comparison.tsv")
    add("Physicochemical interpretation", "measurements_with_larger_MC_incremental_r_squared_n",
        int(measurement_comparison.larger_depth_adjusted_incremental_r_squared.eq("microbial_compartment").sum()),
        "common_cohort_descriptive_comparison", "microbial_state_interpretation/tables/microbial_state_measurement_mc_hybrid_comparison.tsv")

    renewal_mc_path = tables / "microbial_compartment_renewal_association.tsv"
    renewal_module_path = tables / "ecological_module_renewal_association.tsv"
    if renewal_mc_path.exists() and renewal_module_path.exists():
        renewal_mc = read_table(renewal_mc_path)
        renewal_modules = read_table(renewal_module_path)
        for _, row in renewal_mc.iterrows():
            add("Renewal-state association", f"{row.analysis}_MC_depth_adjusted_incremental_r_squared",
                row.incremental_r_squared, "episodic_posthoc_association",
                "microbial_state_interpretation/tables/microbial_compartment_renewal_association.tsv")
            add("Renewal-state association", f"{row.analysis}_MC_adjusted_pvalue",
                row.adjusted_pvalue, "episodic_posthoc_association",
                "microbial_state_interpretation/tables/microbial_compartment_renewal_association.tsv")
        for analysis, subset in renewal_modules.groupby("analysis", sort=False):
            add("Renewal-state association", f"{analysis}_q_le_0.05_modules_n",
                int(pd.to_numeric(subset.adjusted_pvalue, errors="coerce").le(0.05).sum()),
                "episodic_posthoc_module_association",
                "microbial_state_interpretation/tables/ecological_module_renewal_association.tsv")

    mc_composition = composition_stats.loc[composition_stats.grouping_strategy.eq("microbial_compartment")]
    if not mc_composition.empty:
        row = mc_composition.iloc[0]
        for field in ["total_r_squared", "incremental_r_squared", "partial_r_squared", "permutation_pvalue"]:
            add("In-sample ASV composition", field, row.get(field, np.nan),
                "descriptive_only_shared_data_with_MC_inference",
                "microbial_state_interpretation/tables/microbial_state_composition_association.tsv")

    pd.DataFrame(rows).to_csv(
        tables / "microbial_compartment_evidence_summary.tsv", sep="\t", index=False
    )


def taxonomy_map(path: Path | None) -> dict[str, str]:
    if path is None or not path.exists():
        return {}
    table = read_table(path)
    id_col = next((column for column in table.columns if column.strip().casefold() in {"feature id", "feature_id", "asv_id", "asv", "id"}), table.columns[0])
    tax_col = next((column for column in table.columns if column.strip().casefold() in {"taxon", "taxonomy"}), None)
    if tax_col is None:
        return {}
    ids = table[id_col].astype(str).str.strip().str.replace(r";.*$", "", regex=True)
    return dict(zip(ids, table[tax_col].fillna("").astype(str)))


def build_surface(
    data: pd.DataFrame,
    group_col: str,
    labels: list[str],
    date_col: str,
    cruise_col: str,
    depth_col: str,
    maximum_depth: float,
    depth_step: float,
    subdivisions: int,
    time_sigma_months: float,
    depth_visual_sigma_m: float,
    calendar_start: pd.Timestamp,
    calendar_end: pd.Timestamp,
) -> dict[str, object]:
    eligible = data.loc[data[group_col].notna() & data[group_col].astype(str).str.strip().ne("")].copy()
    cruises = eligible[[cruise_col, date_col]].drop_duplicates().sort_values([date_col, cruise_col]).reset_index(drop=True)
    cruises["cruise_index"] = np.arange(len(cruises), dtype=int)
    eligible = eligible.merge(cruises, on=[cruise_col, date_col], validate="many_to_one")
    rows = []
    rank = {label: index for index, label in enumerate(labels)}
    for (cruise_index, depth), group in eligible.groupby(["cruise_index", depth_col], sort=True):
        values = group[group_col].astype(str).value_counts()
        candidates = values.loc[values.eq(values.max())].index.tolist()
        selected = sorted(candidates, key=lambda value: (rank.get(value, 999), value))[0]
        rows.append({"cruise_index": cruise_index, "depth_m": float(depth), "state": selected})
    profile = pd.DataFrame(rows)
    depth_grid = np.arange(0.0, maximum_depth + depth_step * 0.5, depth_step)
    month_starts = pd.date_range(calendar_start, calendar_end, freq="MS")
    fine_dates = []
    for left, right in zip(month_starts[:-1], month_starts[1:]):
        fine_dates.extend(left + (right - left) * ((index + 0.5) / subdivisions) for index in range(subdivisions))
    fine_dates = pd.DatetimeIndex(fine_dates)
    fine_x = (fine_dates - calendar_start).total_seconds().to_numpy() / 86400.0
    cruise_x = (pd.to_datetime(cruises[date_col]) - calendar_start).dt.total_seconds().to_numpy() / 86400.0
    label_to_code = {label: index for index, label in enumerate(labels)}
    coarse = np.full((len(cruises), len(depth_grid)), np.nan)
    for cruise_index, group in profile.groupby("cruise_index", sort=True):
        group = group.sort_values("depth_m")
        depths = group.depth_m.to_numpy(float)
        codes = group.state.map(label_to_code).to_numpy(int)
        selected = np.zeros(len(depth_grid), dtype=int) if len(depths) == 1 else np.searchsorted((depths[:-1] + depths[1:]) / 2, depth_grid, side="right")
        coarse[int(cruise_index)] = codes[selected]
    support = np.zeros((len(labels), len(fine_x), len(depth_grid)))
    for code in range(len(labels)):
        for depth_index in range(len(depth_grid)):
            values = (coarse[:, depth_index] == code).astype(float)
            dated = pd.DataFrame({"x": cruise_x, "value": values}).groupby("x", as_index=False).value.mean().sort_values("x")
            support[code, :, depth_index] = np.interp(fine_x, dated.x, dated.value)
        support[code] = gaussian_filter(support[code], sigma=(time_sigma_months * subdivisions, 0), mode="nearest")
    visual = np.stack([
        gaussian_filter(item, sigma=(0, depth_visual_sigma_m / depth_step), mode="nearest")
        for item in support
    ]).argmax(axis=0)
    # Enforce observed midpoint-bounded profiles at the nearest displayed time column.
    for cruise_index, group in profile.groupby("cruise_index", sort=True):
        group = group.sort_values("depth_m")
        depths = group.depth_m.to_numpy(float)
        codes = group.state.map(label_to_code).to_numpy(int)
        selected = np.zeros(len(depth_grid), dtype=int) if len(depths) == 1 else np.searchsorted((depths[:-1] + depths[1:]) / 2, depth_grid, side="right")
        time_index = int(np.argmin(np.abs(fine_x - cruise_x[int(cruise_index)])))
        visual[time_index, :] = codes[selected]
    return {
        "eligible": eligible, "cruises": cruises, "profile": profile,
        "fine_dates": fine_dates, "fine_x": fine_x, "depth_grid": depth_grid,
        "codes": visual, "labels": labels,
    }


def interpolate_membership_surface(
    observations: pd.DataFrame,
    dates: pd.DatetimeIndex,
    depths: np.ndarray,
    max_time_support_days: float,
    max_depth_support_m: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Interpolate one state-membership field using the BASIN section geometry."""
    origin = dates.min()
    observed_days = (
        observations["date"] - origin
    ).dt.total_seconds().to_numpy() / 86400.0
    grid_days = (dates - origin).total_seconds().to_numpy() / 86400.0
    unique_days = np.unique(observed_days)
    unique_depths = np.unique(observations["depth"].to_numpy(float))
    time_scale = max(
        float(np.median(np.diff(unique_days))) if len(unique_days) > 1 else 1.0,
        1.0,
    )
    depth_scale = max(
        float(np.median(np.diff(unique_depths))) if len(unique_depths) > 1 else 1.0,
        1.0,
    )
    points = np.column_stack(
        (
            observed_days / time_scale,
            observations["depth"].to_numpy(float) / depth_scale,
        )
    )
    date_mesh, depth_mesh = np.meshgrid(grid_days, depths)
    query = np.column_stack(
        (date_mesh.ravel() / time_scale, depth_mesh.ravel() / depth_scale)
    )
    try:
        values = LinearNDInterpolator(
            points,
            observations["value"].to_numpy(float),
            fill_value=np.nan,
        )(query)
    except Exception:
        values = np.full(len(query), np.nan)
    values = np.asarray(values, dtype=float).reshape(date_mesh.shape)

    nearest_date = np.min(
        np.abs(grid_days[:, None] - unique_days[None, :]), axis=1
    )
    nearest_depth = np.min(
        np.abs(depths[:, None] - unique_depths[None, :]), axis=1
    )
    support = (
        np.isfinite(values)
        & (nearest_date[None, :] <= max_time_support_days)
        & (nearest_depth[:, None] <= max_depth_support_m)
    )
    nearest_scaled, _ = cKDTree(points).query(query, k=1)
    nearest_scaled = nearest_scaled.reshape(date_mesh.shape)
    values[~support] = np.nan
    return (
        values,
        support,
        np.broadcast_to(nearest_date, date_mesh.shape),
        np.broadcast_to(nearest_depth[:, None], date_mesh.shape),
        nearest_scaled,
    )


def build_basin_comparable_mc_surface(
    data: pd.DataFrame,
    labels: list[str],
    date_col: str,
    depth_col: str,
    maximum_depth: float,
    depth_step: float,
    time_step_days: int,
    max_time_support_days: float,
    max_depth_support_m: float,
    reference_grid: Path | None = None,
) -> dict[str, object]:
    """Render PAM labels as interpolated one-hot memberships without refitting PAM."""
    eligible = data.loc[
        data["microbial_compartment"].notna()
        & data["microbial_compartment"].astype(str).str.strip().ne("")
        & data[date_col].notna()
        & data[depth_col].notna()
    ].copy()
    eligible = eligible.loc[
        eligible["microbial_compartment"].astype(str).isin(labels)
        & eligible[depth_col].between(0.0, maximum_depth)
    ]
    if eligible.empty:
        raise ValueError("No microbial-compartment observations were available for contouring")
    domain_source = "ASPIRE microbial-compartment observations"
    if reference_grid is not None:
        reference = pd.read_csv(
            reference_grid,
            sep="\t",
            compression="infer",
            usecols=lambda column: column in {"date", "depth_m"},
        )
        if not {"date", "depth_m"}.issubset(reference.columns):
            raise ValueError(
                "The BASIN contour reference grid must contain date and depth_m columns"
            )
        reference["date"] = pd.to_datetime(reference["date"], errors="coerce")
        reference["depth_m"] = pd.to_numeric(reference["depth_m"], errors="coerce")
        reference = reference.dropna(subset=["date", "depth_m"])
        reference = reference.loc[reference["depth_m"].between(0.0, maximum_depth)]
        dates = pd.DatetimeIndex(sorted(reference["date"].unique()))
        depths = np.sort(reference["depth_m"].unique().astype(float))
        if len(dates) < 2 or len(depths) < 2:
            raise ValueError("The BASIN contour reference grid has an insufficient domain")
        domain_source = str(reference_grid)
    else:
        dates = pd.date_range(
            eligible[date_col].min().floor("D"),
            eligible[date_col].max().ceil("D"),
            freq=f"{time_step_days}D",
        ).union(pd.DatetimeIndex(eligible[date_col].dropna().unique())).sort_values()
        depths = np.arange(0.0, maximum_depth + depth_step * 0.5, depth_step)

    # Replicate samples at the same date and depth contribute fractional one-hot
    # membership rather than being resolved by arbitrary row order.
    membership = eligible[[date_col, depth_col, "microbial_compartment"]].copy()
    for label in labels:
        membership[label] = membership["microbial_compartment"].astype(str).eq(label).astype(float)
    grouped = membership.groupby([date_col, depth_col], as_index=False)[labels].mean()

    fields: list[np.ndarray] = []
    common_support: np.ndarray | None = None
    date_distance = depth_distance = scaled_distance = None
    for label in labels:
        observations = grouped[[date_col, depth_col, label]].copy()
        observations.columns = ["date", "depth", "value"]
        field, support, date_distance, depth_distance, scaled_distance = (
            interpolate_membership_surface(
                observations,
                dates,
                depths,
                max_time_support_days,
                max_depth_support_m,
            )
        )
        fields.append(field)
        common_support = support if common_support is None else common_support & support
    cube = np.clip(np.stack(fields, axis=0), 0.0, None)
    denominator = np.nansum(cube, axis=0)
    valid = common_support & np.isfinite(denominator) & (denominator > 0)
    cube[:, valid] /= denominator[valid]
    winner = np.argmax(np.where(np.isfinite(cube), cube, -np.inf), axis=0)
    winner = np.ma.array(winner, mask=~valid)
    maximum_membership = np.where(valid, np.max(cube, axis=0), np.nan)
    return {
        "eligible": eligible,
        "grouped": grouped,
        "dates": dates,
        "depths": depths,
        "labels": labels,
        "memberships": cube,
        "winner": winner,
        "valid": valid,
        "maximum_membership": maximum_membership,
        "nearest_date_days": date_distance,
        "nearest_depth_m": depth_distance,
        "nearest_scaled_distance": scaled_distance,
        "domain_source": domain_source,
    }


def write_basin_comparable_mc_outputs(
    surface: dict[str, object],
    palette: dict[str, str],
    renewal_onsets: pd.DatetimeIndex,
    args: argparse.Namespace,
    tables: Path,
    plots: Path,
    formats: list[str],
) -> None:
    """Write the directly BASIN-comparable microbial-compartment contour and audit."""
    dates = pd.DatetimeIndex(surface["dates"])
    depths = np.asarray(surface["depths"], dtype=float)
    labels = list(surface["labels"])
    winner = surface["winner"]
    valid = np.asarray(surface["valid"], dtype=bool)
    eligible = surface["eligible"]

    basin_style = {
        "font.family": "Times New Roman",
        "font.serif": ["Times New Roman"],
        "font.size": 11,
        "axes.titlesize": 13,
        "axes.labelsize": 12,
        "axes.linewidth": 0.8,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
        "legend.title_fontsize": 11,
        "legend.frameon": False,
        "figure.titlesize": 14,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    }
    figure_width = 14.5
    figure_height = 3.8
    panel_height = 3.35
    vertical_margin = (figure_height - panel_height) / 2.0
    with plt.rc_context(basin_style):
        fig = plt.figure(figsize=(figure_width, figure_height))
        ax = fig.add_axes([
            0.75 / figure_width,
            vertical_margin / figure_height,
            9.50 / figure_width,
            panel_height / figure_height,
        ])
        legend_ax = fig.add_axes([
            10.45 / figure_width,
            vertical_margin / figure_height,
            3.80 / figure_width,
            panel_height / figure_height,
        ])
        legend_ax.set_axis_off()
    fig._aspire_compact_publication_typography = True
    ax.set_facecolor("#E6E6E6")
    ax.contourf(
        dates,
        depths,
        winner,
        levels=np.arange(len(labels) + 1) - 0.5,
        cmap=ListedColormap([palette[label] for label in labels]),
        antialiased=True,
    )
    ax.scatter(
        eligible[args.date_col],
        eligible[args.depth_col],
        c=[palette[str(label)] for label in eligible["microbial_compartment"]],
        s=10,
        edgecolors="black",
        linewidths=0.25,
        zorder=5,
    )
    for event_date in renewal_onsets:
        ax.axvline(
            event_date,
            color="black",
            linestyle="--",
            linewidth=0.65,
            alpha=0.8,
            zorder=4,
        )
    ax.set_xlim(dates.min(), dates.max())
    ax.set_ylim(float(np.max(depths)), 0.0)
    ax.set_xlabel("Sampling date", fontsize=12)
    ax.set_ylabel("Depth (m)", fontsize=12)
    ax.set_title("Microbial compartments", fontsize=13)
    ax.set_yticks(np.arange(0.0, float(np.max(depths)) + 0.1, 25.0))
    ax.tick_params(axis="both", which="major", labelsize=10)
    ax.xaxis.set_major_locator(mdates.YearLocator(2))
    ax.xaxis.set_minor_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    legend_ax.legend(
        handles=[
            Patch(
                facecolor=palette[label],
                edgecolor="black",
                linewidth=0.4,
                label=label,
            )
            for label in labels
        ],
        loc="upper left",
        bbox_to_anchor=(0.0, 1.0),
        frameon=False,
        fontsize=10,
    )
    save(
        fig,
        plots / "continuous_time_depth_compartments_microbial",
        formats,
        bbox_inches=None,
    )

    flat_mask = np.ma.getmaskarray(winner).ravel()
    flat_codes = winner.data.ravel()
    grid = pd.DataFrame(
        {
            "compartment_system": "Microbial",
            "date": np.tile(dates.to_numpy(), len(depths)),
            "depth_m": np.repeat(depths, len(dates)),
            "interpolated_compartment": [
                pd.NA if masked else labels[int(code)]
                for code, masked in zip(flat_codes, flat_mask)
            ],
            "maximum_interpolated_membership": np.asarray(
                surface["maximum_membership"]
            ).ravel(),
            "supported": valid.ravel(),
            "nearest_observed_date_days": np.asarray(
                surface["nearest_date_days"]
            ).ravel(),
            "nearest_observed_depth_m": np.asarray(
                surface["nearest_depth_m"]
            ).ravel(),
            "nearest_observation_scaled_distance": np.asarray(
                surface["nearest_scaled_distance"]
            ).ravel(),
        }
    )
    for index, label in enumerate(labels):
        grid[f"membership_{label}"] = np.asarray(surface["memberships"])[index].ravel()
    grid.to_csv(
        tables / "continuous_compartment_grid_microbial.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    eligible[
        [args.sample_col, args.date_col, args.depth_col, "microbial_compartment"]
    ].to_csv(
        tables / "continuous_compartment_observations_microbial.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(
        {
            "compartment": labels,
            "color_hex": [palette[label] for label in labels],
            "palette_source": "configured ASPIRE microbial-compartment palette",
        }
    ).to_csv(
        tables / "continuous_compartment_palette_microbial.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(
        [
            {
                "rendering_method": "two-dimensional linear interpolation of one-hot PAM memberships",
                "cluster_assignments_refit": False,
                "number_of_anchor_samples": len(eligible),
                "number_of_date_depth_anchors": len(surface["grouped"]),
                "number_of_compartments": len(labels),
                "first_date": dates.min().date(),
                "last_date": dates.max().date(),
                "time_step_days": args.contour_time_step_days,
                "depth_step_m": args.depth_step,
                "maximum_depth_m": args.maximum_depth,
                "maximum_time_support_days": args.contour_max_time_support_days,
                "maximum_depth_support_m": args.contour_max_depth_support_m,
                "domain_source": surface["domain_source"],
                "figure_width_inches": 14.5,
                "figure_height_inches": 3.8,
                "panel_width_inches": 9.5,
                "panel_height_inches": 3.35,
                "supported_grid_fraction": float(valid.mean()),
            }
        ]
    ).to_csv(
        tables / "continuous_compartment_audit_microbial.tsv",
        sep="\t",
        index=False,
    )


def draw_surface(ax: plt.Axes, surface: dict[str, object], palette: dict[str, str], calendar_start: pd.Timestamp, calendar_end: pd.Timestamp) -> None:
    fine_x = np.asarray(surface["fine_x"])
    depth = np.asarray(surface["depth_grid"])
    codes = np.asarray(surface["codes"])
    labels = list(surface["labels"])
    plot_x = np.concatenate(([0.0], fine_x, [float((calendar_end - calendar_start).days)]))
    plot_codes = np.concatenate((codes[:1], codes, codes[-1:]), axis=0)
    for code, label in enumerate(labels):
        mask = (plot_codes == code).astype(float).T
        if not mask.any():
            continue
        ax.contourf(plot_x, depth, mask, levels=[0.5, 1.5], colors=[palette[label]], antialiased=True, zorder=1)
        if mask.min() < 0.5 < mask.max():
            ax.contour(plot_x, depth, mask, levels=[0.5], colors=["#4D4D4D"], linewidths=0.25, alpha=0.55, zorder=2)


def format_time_depth_axis(ax: plt.Axes, calendar_start: pd.Timestamp, calendar_end: pd.Timestamp, maximum_depth: float) -> None:
    year_starts = pd.date_range(calendar_start, calendar_end, freq="YS", inclusive="left")
    month_starts = pd.date_range(calendar_start, calendar_end, freq="MS", inclusive="left")
    ax.set_xticks((year_starts - calendar_start).total_seconds().to_numpy() / 86400.0)
    ax.set_xticklabels([str(date.year) for date in year_starts], ha="left", fontsize=8)
    ax.set_xticks((month_starts - calendar_start).total_seconds().to_numpy() / 86400.0, minor=True)
    ax.tick_params(axis="x", which="minor", length=2.5)
    ax.set_xlim(0, float((calendar_end - calendar_start).days))
    ax.set_ylim(maximum_depth, 0)
    ax.set_yticks(np.arange(0, maximum_depth + 1, 25))
    ax.set_xlabel("Calendar time (monthly intervals)")
    ax.set_ylabel("Depth (m)")


def save(
    fig: plt.Figure,
    base: Path,
    formats: list[str],
    bbox_inches: str | None = "tight",
) -> None:
    for fmt in formats:
        fig.savefig(
            base.with_suffix(f".{fmt}"),
            dpi=300 if fmt == "png" else None,
            bbox_inches=bbox_inches,
        )
    plt.close(fig)


def load_renewal_onsets(path: Path | None, date_col: str) -> pd.DatetimeIndex:
    if path is None:
        return pd.DatetimeIndex([])
    frame = read_table(path)
    if date_col not in frame.columns:
        raise ValueError(f"Renewal-event table lacks onset column {date_col!r}")
    dates = pd.to_datetime(frame[date_col], errors="coerce").dropna().drop_duplicates()
    return pd.DatetimeIndex(sorted(dates))


def draw_renewal_onsets(
    ax: plt.Axes,
    onsets: pd.DatetimeIndex,
    calendar_start: pd.Timestamp,
    calendar_end: pd.Timestamp,
) -> None:
    for onset in onsets[(onsets >= calendar_start) & (onsets <= calendar_end)]:
        x = float((onset - calendar_start).total_seconds() / 86400.0)
        ax.axvline(
            x, color="black", linestyle="--", linewidth=0.9,
            alpha=0.9, zorder=2.5,
        )


def renewal_legend_handle() -> Line2D:
    return Line2D(
        [], [], color="black", linestyle="--", linewidth=0.9,
        label="Predicted renewal onset",
    )


def curtain_outputs(
    joined: pd.DataFrame,
    relative: pd.DataFrame,
    taxonomy: dict[str, str],
    args: argparse.Namespace,
    tables: Path,
    plots: Path,
) -> None:
    formats = parse_list(args.formats)
    joined = joined.copy()
    joined[args.date_col] = pd.to_datetime(joined[args.date_col], errors="coerce")
    joined[args.depth_col] = pd.to_numeric(joined[args.depth_col], errors="coerce")
    joined = joined.dropna(subset=[args.date_col, args.depth_col, args.cruise_col])
    calendar_start = pd.Timestamp(year=int(joined[args.date_col].dt.year.min()), month=1, day=1)
    calendar_end = pd.Timestamp(year=int(joined[args.date_col].dt.year.max()) + 1, month=1, day=1)
    renewal_onsets = load_renewal_onsets(args.renewal_events, args.renewal_date_col)
    renewal_onsets = renewal_onsets[
        (renewal_onsets >= calendar_start) & (renewal_onsets <= calendar_end)
    ]
    pd.DataFrame({"renewal_onset": renewal_onsets}).to_csv(
        tables / "microbial_state_curtain_renewal_onsets.tsv", sep="\t", index=False
    )
    excluded = re.compile(args.exclude_curtain_label_pattern) if args.exclude_curtain_label_pattern else None
    hybrid_observed = set(joined[args.hybrid_col].dropna().astype(str))
    if excluded:
        hybrid_observed = {label for label in hybrid_observed if not excluded.search(label)}
    hybrid_labels = [label for label in parse_list(args.hybrid_order) if label in hybrid_observed]
    hybrid_labels.extend(sorted(hybrid_observed.difference(hybrid_labels), key=label_sort))
    mc_labels = [label for label in parse_list(args.mc_order) if label in set(joined.microbial_compartment.astype(str))]
    mc_labels.extend(sorted(set(joined.microbial_compartment.astype(str)).difference(mc_labels), key=label_sort))
    hybrid_palette = parse_palette(args.hybrid_palette)
    mc_palette = parse_palette(args.mc_palette)
    validate_palette(hybrid_palette, hybrid_labels, "hybrid")
    validate_palette(mc_palette, mc_labels, "microbial compartment")
    comparable_mc_surface = build_basin_comparable_mc_surface(
        joined,
        mc_labels,
        args.date_col,
        args.depth_col,
        args.maximum_depth,
        args.depth_step,
        args.contour_time_step_days,
        args.contour_max_time_support_days,
        args.contour_max_depth_support_m,
        args.contour_reference_grid,
    )
    write_basin_comparable_mc_outputs(
        comparable_mc_surface,
        mc_palette,
        renewal_onsets,
        args,
        tables,
        plots,
        formats,
    )
    common = joined.loc[joined[args.hybrid_col].astype(str).isin(hybrid_labels) & joined.microbial_compartment.notna()].copy()
    hybrid_surface = build_surface(
        common, args.hybrid_col, hybrid_labels, args.date_col, args.cruise_col, args.depth_col,
        args.maximum_depth, args.depth_step, args.time_subdivisions_per_month,
        args.time_sigma_months, args.depth_visual_sigma_m, calendar_start, calendar_end,
    )
    mc_all_surface = build_surface(
        joined, "microbial_compartment", mc_labels, args.date_col, args.cruise_col, args.depth_col,
        args.maximum_depth, args.depth_step, args.time_subdivisions_per_month,
        args.time_sigma_months, args.depth_visual_sigma_m, calendar_start, calendar_end,
    )
    mc_common_surface = build_surface(
        common, "microbial_compartment", mc_labels, args.date_col, args.cruise_col, args.depth_col,
        args.maximum_depth, args.depth_step, args.time_subdivisions_per_month,
        args.time_sigma_months, args.depth_visual_sigma_m, calendar_start, calendar_end,
    )
    for name, surface in [
        ("asv_only_hybrid_compartment", hybrid_surface),
        ("microbial_compartment", mc_all_surface),
        ("microbial_compartment_common_cohort", mc_common_surface),
    ]:
        grid = pd.DataFrame({
            "calendar_date": np.repeat(surface["fine_dates"], len(surface["depth_grid"])),
            "depth_m": np.tile(surface["depth_grid"], len(surface["fine_dates"])),
            "rendered_state": np.asarray(surface["labels"], dtype=object)[np.asarray(surface["codes"]).ravel()],
        })
        grid.to_csv(tables / f"{name}_curtain_grid.tsv", sep="\t", index=False)
    pd.DataFrame([
        {"surface": "hybrid", "samples_n": len(common), "cohort": "samples_with_hybrid_and_MC"},
        {"surface": "MC", "samples_n": len(joined), "cohort": "all_MC_samples"},
        {"surface": "MC_common", "samples_n": len(common), "cohort": "samples_with_hybrid_and_MC"},
    ]).to_csv(tables / "microbial_state_curtain_cohort_audit.tsv", sep="\t", index=False)

    def plot_single(surface: dict[str, object], palette: dict[str, str], title: str, output: str) -> None:
        fig, ax = plt.subplots(figsize=(15.5, 7.2))
        fig._aspire_compact_publication_typography = True
        draw_surface(ax, surface, palette, calendar_start, calendar_end)
        draw_renewal_onsets(ax, renewal_onsets, calendar_start, calendar_end)
        eligible = surface["eligible"]
        x = (pd.to_datetime(eligible[args.date_col]) - calendar_start).dt.total_seconds() / 86400.0
        ax.scatter(x, eligible[args.depth_col], s=24, facecolor=[palette[str(value)] for value in eligible[args.hybrid_col if "hybrid" in output else "microbial_compartment"]], edgecolor="black", linewidth=0.7, zorder=4)
        format_time_depth_axis(ax, calendar_start, calendar_end, args.maximum_depth)
        ax.set_title(title)
        handles = [Patch(facecolor=palette[label], edgecolor="black", linewidth=0.4, label=display_label(label)) for label in surface["labels"]]
        handles.append(renewal_legend_handle())
        ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False, fontsize=11)
        fig.subplots_adjust(right=0.82)
        save(fig, plots / output, formats)

    plot_single(hybrid_surface, hybrid_palette, r"ASV-sampled hybrid O$_2$-GMM compartments", "asv_only_hybrid_compartment_curtain")
    plot_single(mc_all_surface, mc_palette, "ASV-inferred microbial compartments", "microbial_compartment_curtain")

    fig, axes = plt.subplots(1, 2, figsize=(19.0, 7.2), sharey=True)
    fig._aspire_compact_publication_typography = True
    for ax, surface, palette, title in [
        (axes[0], hybrid_surface, hybrid_palette, r"Hybrid O$_2$-GMM"),
        (axes[1], mc_common_surface, mc_palette, "Microbial compartments"),
    ]:
        draw_surface(ax, surface, palette, calendar_start, calendar_end)
        draw_renewal_onsets(ax, renewal_onsets, calendar_start, calendar_end)
        format_time_depth_axis(ax, calendar_start, calendar_end, args.maximum_depth)
        ax.set_title(title)
    axes[1].set_ylabel("")
    hybrid_handles = [Patch(facecolor=hybrid_palette[label], label=display_label(label)) for label in hybrid_labels]
    mc_handles = [Patch(facecolor=mc_palette[label], label=label) for label in mc_labels]
    hybrid_handles.append(renewal_legend_handle())
    mc_handles.append(renewal_legend_handle())
    axes[0].legend(handles=hybrid_handles, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=3, frameon=False, fontsize=10)
    axes[1].legend(handles=mc_handles, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=4, frameon=False, fontsize=10)
    fig.subplots_adjust(bottom=0.28, wspace=0.08)
    save(fig, plots / "hybrid_vs_microbial_compartment_curtain", formats)

    common_ids = [sample for sample in common[args.sample_col].astype(str) if sample in relative.index]
    relative_common = relative.loc[common_ids]
    explicit = [item for item in parse_list(args.selected_asvs) if item in relative_common.columns]
    metrics = pd.DataFrame({
        "ASV_ID": relative_common.columns,
        "mean_relative_abundance": relative_common.mean().to_numpy(float),
        "median_relative_abundance": relative_common.median().to_numpy(float),
        "maximum_relative_abundance": relative_common.max().to_numpy(float),
        "prevalence": relative_common.gt(0).mean().to_numpy(float),
    })
    if explicit:
        selected = explicit
        selection_source = "configured_ASV_list"
    else:
        selected = metrics.sort_values(
            [args.asv_selection_metric, "prevalence", "ASV_ID"], ascending=[False, False, True]
        ).head(args.asv_top_n).ASV_ID.tolist()
        selection_source = f"top_{args.asv_top_n}_by_{args.asv_selection_metric}"
    selection = metrics.loc[metrics.ASV_ID.isin(selected)].copy()
    selection["selection_rank"] = selection.ASV_ID.map({asv: index + 1 for index, asv in enumerate(selected)})
    selection["selection_source"] = selection_source
    selection["taxonomy"] = selection.ASV_ID.map(taxonomy).fillna("")
    selection.sort_values("selection_rank").to_csv(tables / "asv_curtain_selection.tsv", sep="\t", index=False)
    overlay_dir = plots / "asv_relative_abundance_curtains"
    overlay_dir.mkdir(exist_ok=True)
    common_plot = common.set_index(args.sample_col).loc[common_ids].copy()
    x = (pd.to_datetime(common_plot[args.date_col]) - calendar_start).dt.total_seconds() / 86400.0
    global_overlay_maximum = float(relative_common[selected].to_numpy(float).max()) if selected else 0.0
    for asv in selected:
        abundance = relative_common[asv].to_numpy(float)
        present = abundance > 0
        sizes = 18.0 + 170.0 * abundance / global_overlay_maximum if global_overlay_maximum > 0 else np.zeros(len(abundance))
        fig, axes = plt.subplots(1, 2, figsize=(19.0, 7.2), sharey=True)
        fig._aspire_compact_publication_typography = True
        for ax, surface, palette, title in [
            (axes[0], hybrid_surface, hybrid_palette, r"Hybrid O$_2$-GMM"),
            (axes[1], mc_common_surface, mc_palette, "Microbial compartments"),
        ]:
            draw_surface(ax, surface, palette, calendar_start, calendar_end)
            draw_renewal_onsets(ax, renewal_onsets, calendar_start, calendar_end)
            ax.scatter(x[present], common_plot.loc[present, args.depth_col], s=sizes[present], facecolor="white", edgecolor="black", linewidth=0.75, alpha=0.92, zorder=5)
            format_time_depth_axis(ax, calendar_start, calendar_end, args.maximum_depth)
            ax.set_title(title)
        axes[1].set_ylabel("")
        tax = taxonomy.get(asv, "")
        fig.suptitle(f"{asv} relative abundance across environmental and microbial compartments\n{tax}", y=0.99)
        legend_values = [value for value in [0.0001, 0.001, 0.01, 0.1] if value <= global_overlay_maximum]
        if global_overlay_maximum > 0 and not legend_values:
            legend_values = [global_overlay_maximum]
        handles = [Line2D([], [], marker="o", linestyle="", markerfacecolor="white", markeredgecolor="black", markersize=math.sqrt(18 + 170 * value / global_overlay_maximum), label=f"{100 * value:g}%") for value in legend_values]
        handles.append(renewal_legend_handle())
        axes[1].legend(handles=handles, title="Relative abundance", loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False)
        fig.subplots_adjust(right=0.88, top=0.88, wspace=0.08)
        save(fig, overlay_dir / f"{asv}_hybrid_vs_microbial_compartment_curtain", formats)


def microbial_state_validation_outputs(
    relative: pd.DataFrame,
    joined: pd.DataFrame,
    args: argparse.Namespace,
    tables: Path,
    plots: Path,
) -> None:
    """Create the common Depth → hybrid → MC validation and interpretation set."""
    sample_ids = [sample for sample in joined[args.sample_col].astype(str) if sample in relative.index]
    frame = joined.set_index(args.sample_col).loc[sample_ids].reset_index()
    composition = relative.loc[sample_ids].astype(float)
    composition = composition.div(composition.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    positive = composition.to_numpy() > 0
    values = composition.to_numpy()
    shannon = -np.where(positive, values * np.log(np.where(positive, values, 1.0)), 0.0).sum(axis=1)
    richness = positive.sum(axis=1).astype(float)
    diversity = frame[[args.sample_col, args.cruise_col, args.depth_col, args.hybrid_col, "microbial_compartment"]].copy()
    diversity["shannon_diversity"] = shannon
    diversity["observed_asv_richness"] = richness
    diversity.to_csv(tables / "microbial_state_alpha_diversity.tsv", sep="\t", index=False)

    clr_path = args.microbial_dir / "audit" / "asv_clr_matrix.tsv"
    clr = read_table(clr_path)
    clr_id = clr.columns[0]
    clr[clr_id] = clr[clr_id].astype(str).str.strip()
    clr = clr.set_index(clr_id).loc[sample_ids]
    reducer = umap.UMAP(
        n_neighbors=min(15, max(2, len(clr) - 1)), min_dist=0.1,
        metric="euclidean", random_state=args.seed, n_components=2,
    )
    embedding = reducer.fit_transform(clr.to_numpy(float))
    ordination = diversity.copy()
    ordination["UMAP1"] = embedding[:, 0]
    ordination["UMAP2"] = embedding[:, 1]
    ordination.to_csv(tables / "microbial_state_umap_coordinates.tsv", sep="\t", index=False)

    depth_palette = parse_palette(args.depth_palette)
    hybrid_palette = parse_palette(args.hybrid_palette)
    mc_palette = parse_palette(args.mc_palette)
    formats = parse_list(args.formats)
    overlay_specs = [
        (args.depth_col, depth_palette, "Depth"),
        (args.hybrid_col, hybrid_palette, r"Hybrid O$_2$-GMM compartment"),
        ("microbial_compartment", mc_palette, "Microbial compartment"),
    ]
    with plt.rc_context(COMPACT_PLOT_STYLE):
        fig, axes = plt.subplots(1, 3, figsize=(16.0, 7.8), sharex=True, sharey=True)
        for index, (ax, (column, palette, title)) in enumerate(zip(axes, overlay_specs)):
            labels = ordination[column].map(normalized_plot_label)
            observed = sorted([label for label in labels.unique() if label], key=natural_key)
            colors = [palette.get(label, "#808080") for label in labels]
            ax.scatter(
                ordination.UMAP1, ordination.UMAP2, c=colors, s=24, alpha=0.9,
                edgecolor="white", linewidth=0.25,
            )
            ax.set_title(title, pad=9)
            ax.set_xlabel("UMAP1")
            handles = [
                Line2D(
                    [], [], marker="o", linestyle="", markersize=5,
                    markerfacecolor=palette.get(label, "#808080"),
                    markeredgecolor="0.25", markeredgewidth=0.35,
                    label=display_label(label),
                )
                for label in observed
            ]
            legend_columns = 3 if index == 0 else 2
            ax.legend(
                handles=handles, title=title, loc="upper center",
                bbox_to_anchor=(0.5, -0.18), ncol=legend_columns,
                handletextpad=0.45, columnspacing=0.9, frameon=False,
            )
        axes[0].set_ylabel("UMAP2")
        fig.subplots_adjust(left=0.07, right=0.99, top=0.90, bottom=0.34, wspace=0.18)
        finalize_compact_figure(fig, axes)
    save(fig, plots / "microbial_state_umap_depth_hybrid_mc", formats)

    rng = np.random.default_rng(args.seed + 101)
    association_rows: list[dict[str, object]] = []
    for response_name, response in (("Shannon diversity", shannon), ("Observed ASV richness", richness)):
        for strategy in (args.hybrid_col, "microbial_compartment"):
            result = depth_adjusted_group_test(
                response, frame, strategy, args.depth_col, args.cruise_col,
                args.permutations, rng,
            )
            association_rows.append({"response": response_name, "grouping_strategy": strategy, **result})
    diversity_stats = pd.DataFrame(association_rows)
    diversity_stats["adjusted_pvalue"] = diversity_stats.groupby("response", group_keys=False)["permutation_pvalue"].apply(benjamini_hochberg)
    diversity_stats.to_csv(tables / "microbial_state_alpha_diversity_association.tsv", sep="\t", index=False)

    with plt.rc_context(COMPACT_PLOT_STYLE):
        fig, axes = plt.subplots(2, 3, figsize=(16.0, 9.0), sharex="col")
        for row_index, response_name in enumerate(["shannon_diversity", "observed_asv_richness"]):
            for ax, (column, palette, title) in zip(axes[row_index], overlay_specs):
                labels = diversity[column].map(normalized_plot_label)
                observed = sorted([label for label in labels.unique() if label], key=natural_key)
                arrays = [
                    diversity.loc[labels.eq(label), response_name].dropna().to_numpy(float)
                    for label in observed
                ]
                boxes = ax.boxplot(
                    arrays, tick_labels=[display_label(label) for label in observed],
                    patch_artist=True, showfliers=False, widths=0.68,
                )
                for box, label in zip(boxes["boxes"], observed):
                    box.set_facecolor(palette.get(label, "#BDBDBD"))
                    box.set_edgecolor("0.20")
                    box.set_linewidth(0.65)
                for median in boxes["medians"]:
                    median.set_color("black")
                    median.set_linewidth(1.2)
                ax.set_title(title if row_index == 0 else "", pad=8)
                ax.tick_params(axis="x", labelrotation=90, labelsize=8)
            axes[row_index, 0].set_ylabel(
                "Shannon diversity" if response_name == "shannon_diversity" else "Observed ASV richness"
            )
        fig.subplots_adjust(left=0.07, right=0.99, top=0.92, bottom=0.25, wspace=0.20, hspace=0.22)
        finalize_compact_figure(fig, axes)
    save(fig, plots / "microbial_state_alpha_diversity_depth_hybrid_mc", formats)

    composition_rows = []
    for strategy in (args.hybrid_col, "microbial_compartment"):
        result = depth_adjusted_group_test(
            clr.to_numpy(float), frame, strategy, args.depth_col, args.cruise_col,
            args.permutations, rng,
        )
        composition_rows.append({
            "response": "CLR-transformed ASV composition",
            "distance_interpretation": "Euclidean sums of squares in CLR space (Aitchison geometry)",
            "grouping_strategy": strategy,
            "inference_role": "descriptive_for_MC_due_to_shared_ASV_input" if strategy == "microbial_compartment" else "external_BASIN_comparator",
            **result,
        })
    composition_stats = pd.DataFrame(composition_rows)
    composition_stats["adjusted_pvalue"] = benjamini_hochberg(composition_stats.permutation_pvalue)
    composition_stats.to_csv(tables / "microbial_state_composition_association.tsv", sep="\t", index=False)

    measurement = read_table(args.measurement_matrix)
    measurement_id = args.sample_col if args.sample_col in measurement.columns else measurement.columns[0]
    measurement[measurement_id] = measurement[measurement_id].astype(str).str.strip()
    measurement = measurement.rename(columns={measurement_id: args.sample_col})
    grouping_frame = frame[[
        args.sample_col, args.cruise_col, args.depth_col,
        args.hybrid_col, "microbial_compartment",
    ]].copy()
    measurement = grouping_frame.merge(
        measurement, on=args.sample_col, how="inner", validate="one_to_one"
    )
    metadata_columns = set(grouping_frame.columns)
    measurement_columns = [
        column for column in measurement.columns if column not in metadata_columns
        and pd.to_numeric(measurement[column], errors="coerce").notna().any()
    ]
    measurement_rows: list[dict[str, object]] = []
    profile_rows: list[dict[str, object]] = []
    for column in measurement_columns:
        response = pd.to_numeric(measurement[column], errors="coerce").to_numpy(float)
        common_mask = (
            measurement[[args.hybrid_col, "microbial_compartment"]].notna().all(axis=1).to_numpy()
            & np.isfinite(response)
        )
        mc_mask = measurement["microbial_compartment"].notna().to_numpy() & np.isfinite(response)
        cohorts = [
            (
                "common_hybrid_MC_samples", measurement.loc[common_mask].reset_index(drop=True),
                response[common_mask], [args.hybrid_col, "microbial_compartment"],
            ),
            (
                "all_MC_samples", measurement.loc[mc_mask].reset_index(drop=True),
                response[mc_mask], ["microbial_compartment"],
            ),
        ]
        for cohort, cohort_frame, cohort_response, strategies in cohorts:
            for strategy in strategies:
                result = depth_adjusted_group_test(
                    cohort_response, cohort_frame, strategy, args.depth_col, args.cruise_col,
                    args.permutations, rng,
                )
                measurement_rows.append({
                    "measurement": column, "analysis_cohort": cohort,
                    "grouping_strategy": strategy, **result,
                })
                cohort_frame = cohort_frame.copy()
                cohort_frame["_measurement_value"] = cohort_response
                for level, group in cohort_frame.groupby(strategy, dropna=True, sort=True):
                    numeric = group["_measurement_value"].dropna()
                    profile_rows.append({
                        "measurement": column, "analysis_cohort": cohort,
                        "grouping_strategy": strategy, "group": level,
                        "samples_n": len(numeric), "mean": numeric.mean(),
                        "standard_deviation": numeric.std(ddof=1), "median": numeric.median(),
                        "q25": numeric.quantile(0.25), "q75": numeric.quantile(0.75),
                    })
    measurement_stats = pd.DataFrame(measurement_rows)
    if not measurement_stats.empty:
        measurement_stats["adjusted_pvalue"] = measurement_stats.groupby(
            ["analysis_cohort", "grouping_strategy"], group_keys=False
        )["permutation_pvalue"].apply(benjamini_hochberg)
    measurement_stats.to_csv(tables / "microbial_state_measurement_association.tsv", sep="\t", index=False)
    measurement_profiles = pd.DataFrame(profile_rows)
    measurement_profiles.to_csv(tables / "microbial_state_measurement_profiles.tsv", sep="\t", index=False)
    if not measurement_profiles.empty:
        with plt.rc_context(COMPACT_PLOT_STYLE):
            fig, axes = plt.subplots(
                1, 2,
                figsize=(14.2, max(7.2, 0.36 * len(measurement_columns) + 3.2)),
                gridspec_kw={"width_ratios": [2.0, 1.0]},
            )
            image = None
            for ax, strategy, title in [
                (axes[0], args.hybrid_col, r"Hybrid O$_2$-GMM"),
                (axes[1], "microbial_compartment", "Microbial compartments"),
            ]:
                subset = measurement_profiles.loc[
                    measurement_profiles.grouping_strategy.eq(strategy)
                    & measurement_profiles.analysis_cohort.eq("common_hybrid_MC_samples")
                ]
                matrix = subset.pivot(index="measurement", columns="group", values="median").reindex(measurement_columns)
                row_mean = matrix.mean(axis=1)
                row_sd = matrix.std(axis=1, ddof=0).replace(0, np.nan)
                standardized = matrix.sub(row_mean, axis=0).div(row_sd, axis=0).fillna(0.0)
                image = ax.imshow(
                    standardized.to_numpy(float), aspect="auto", cmap="RdBu_r", vmin=-2.0, vmax=2.0
                )
                ax.set_xticks(
                    np.arange(len(standardized.columns)),
                    [display_label(item) for item in standardized.columns],
                    rotation=90,
                )
                ax.set_yticks(np.arange(len(standardized.index)), standardized.index)
                ax.set_title(title, pad=10)
                ax.set_xlabel("Compartment")
            axes[1].tick_params(axis="y", labelleft=False)
            axes[0].set_ylabel("Environmental measurement")
            fig.subplots_adjust(left=0.18, right=0.98, top=0.88, bottom=0.32, wspace=0.06)
            cax = fig.add_axes([0.37, 0.10, 0.42, 0.025])
            colorbar = fig.colorbar(image, cax=cax, orientation="horizontal")
            colorbar.set_label("Median within-compartment value (measurement-wise z score)")
            finalize_compact_figure(fig, axes)
        save(fig, plots / "microbial_state_measurement_profiles_hybrid_mc", formats)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--microbial-dir", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--modules", type=Path, required=True)
    parser.add_argument("--taxonomy", type=Path)
    parser.add_argument("--measurement-matrix", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--sample-col", default="sampleID")
    parser.add_argument("--cruise-col", default="Cruise")
    parser.add_argument("--depth-col", default="Depth")
    parser.add_argument("--season-col", default="Season")
    parser.add_argument("--year-col", default="Year")
    parser.add_argument("--date-col", default="Date")
    parser.add_argument("--renewal-col", default="renewal_phase")
    parser.add_argument("--renewal-levels", default="baseline,renewal,post-renewal")
    parser.add_argument("--environmental-compartment-cols", default="o2_compartment,gmm_component,o2_subcompartment_final")
    parser.add_argument("--hybrid-col", default="o2_subcompartment_final")
    parser.add_argument("--permutations", type=int, default=999)
    parser.add_argument("--bootstrap-replicates", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--hybrid-palette", required=True)
    parser.add_argument("--hybrid-order", default="")
    parser.add_argument("--mc-palette", default="MC1=#0072B2,MC2=#E69F00,MC3=#009E73,MC4=#CC79A7")
    parser.add_argument("--depth-palette", default="")
    parser.add_argument("--mc-order", default="MC1,MC2,MC3,MC4")
    parser.add_argument("--asv-top-n", type=int, default=20)
    parser.add_argument("--selected-asvs", default="")
    parser.add_argument("--asv-selection-metric", choices=["mean_relative_abundance", "median_relative_abundance", "maximum_relative_abundance", "prevalence"], default="mean_relative_abundance")
    parser.add_argument("--maximum-depth", type=float, default=210.0)
    parser.add_argument("--depth-step", type=float, default=1.0)
    parser.add_argument("--time-subdivisions-per-month", type=int, default=4)
    parser.add_argument("--time-sigma-months", type=float, default=0.75)
    parser.add_argument("--depth-visual-sigma-m", type=float, default=5.0)
    parser.add_argument("--contour-time-step-days", type=int, default=7)
    parser.add_argument("--contour-max-time-support-days", type=float, default=90.0)
    parser.add_argument("--contour-max-depth-support-m", type=float, default=30.0)
    parser.add_argument("--contour-reference-grid", type=Path)
    parser.add_argument("--exclude-curtain-label-pattern", default=r"(?i)(?:other|outlier)")
    parser.add_argument("--renewal-events", type=Path)
    parser.add_argument("--renewal-date-col", default="start_date")
    parser.add_argument("--formats", default="pdf,png,svg")
    args = parser.parse_args()
    if args.permutations < 1 or args.bootstrap_replicates < 100:
        raise ValueError("Permutation count must be positive and bootstrap replicates at least 100")
    if args.asv_top_n < 1:
        raise ValueError("--asv-top-n must be positive")
    if (
        args.contour_time_step_days < 1
        or args.contour_max_time_support_days <= 0
        or args.contour_max_depth_support_m <= 0
    ):
        raise ValueError("BASIN-comparable contour grid and support settings must be positive")
    if args.contour_reference_grid is not None and not args.contour_reference_grid.is_file():
        raise FileNotFoundError(
            f"BASIN contour reference grid not found: {args.contour_reference_grid}"
        )

    tables = args.outdir / "tables"
    plots = args.outdir / "plots"
    audit = args.outdir / "audit"
    for directory in [tables, plots, audit]:
        directory.mkdir(parents=True, exist_ok=True)
    assignments = read_table(args.microbial_dir / "tables" / "microbial_compartments.tsv")
    relative_table = read_table(args.microbial_dir / "audit" / "asv_relative_abundance.tsv")
    relative_id = relative_table.columns[0]
    relative_table[relative_id] = relative_table[relative_id].astype(str).str.strip()
    relative = relative_table.set_index(relative_id)
    metadata = read_table(args.metadata)
    environmental_cols = parse_list(args.environmental_compartment_cols)
    required = [
        args.sample_col, args.cruise_col, args.depth_col, args.season_col,
        args.year_col, args.date_col, args.renewal_col, *environmental_cols,
    ]
    missing = [column for column in required if column not in metadata.columns]
    if missing:
        raise ValueError(f"Metadata lacks required post hoc columns: {missing}")
    joined = merge_assignments_with_metadata(assignments, metadata, args.sample_col)
    joined = joined.loc[joined[args.sample_col].isin(relative.index)].copy()
    joined["_row_id"] = np.arange(len(joined), dtype=int)
    joined[args.cruise_col] = normalize_key(joined[args.cruise_col])
    joined[args.depth_col] = pd.to_numeric(joined[args.depth_col], errors="coerce")
    joined[args.season_col] = joined[args.season_col].astype("string").str.strip()
    joined[args.date_col] = pd.to_datetime(joined[args.date_col], errors="coerce")
    for column in environmental_cols:
        joined[column] = joined[column].astype("string").str.strip()
        joined.loc[joined[column].isin(["", "nan", "<NA>"]), column] = pd.NA
    joined.to_csv(tables / "microbial_state_sample_crosswalk.tsv", sep="\t", index=False)
    microbial_state_validation_outputs(relative, joined, args, tables, plots)

    common_mask = joined[["microbial_compartment", args.cruise_col, args.depth_col, args.season_col, *environmental_cols]].notna().all(axis=1)
    common = joined.loc[common_mask].reset_index(drop=True)
    if common.empty:
        raise ValueError("No common samples contained MC, depth, season, cruise, and all BASIN labels")
    rng = np.random.default_rng(args.seed)
    correspondence_rows: list[dict[str, object]] = []
    correspondence_nulls: list[pd.DataFrame] = []
    for column in environmental_cols:
        metrics, null = restricted_permutations(
            common, "microbial_compartment", column, args.cruise_col,
            args.permutations, rng,
        )
        correspondence_rows.append({
            "environmental_variable": column, "common_samples_n": len(common),
            "cruises_n": common[args.cruise_col].nunique(),
            "microbial_compartments_n": common.microbial_compartment.nunique(),
            "environmental_compartments_n": common[column].nunique(), **metrics,
        })
        null.insert(0, "environmental_variable", column)
        correspondence_nulls.append(null)
    correspondence = pd.DataFrame(correspondence_rows)
    correspondence.to_csv(tables / "microbial_basin_correspondence.tsv", sep="\t", index=False)
    pd.concat(correspondence_nulls, ignore_index=True).to_csv(tables / "microbial_basin_correspondence_permutation_null.tsv", sep="\t", index=False)

    prediction_summary, prediction_folds, prediction_samples = held_out_prediction(
        common, "microbial_compartment", args.cruise_col, args.depth_col, args.season_col,
        environmental_cols, args.permutations, args.bootstrap_replicates, args.seed,
    )
    prediction_summary.to_csv(tables / "microbial_state_prediction_summary.tsv", sep="\t", index=False)
    prediction_folds.to_csv(tables / "microbial_state_prediction_fold_metrics.tsv", sep="\t", index=False)
    prediction_samples.to_csv(tables / "microbial_state_prediction_sample_predictions.tsv", sep="\t", index=False)

    boundary_summary, boundary_pairs, boundary_null = boundary_concordance(
        common, "microbial_compartment", environmental_cols, args.cruise_col,
        args.depth_col, args.permutations, args.seed,
    )
    boundary_summary.to_csv(tables / "microbial_boundary_concordance.tsv", sep="\t", index=False)
    boundary_pairs.to_csv(tables / "microbial_boundary_pairs.tsv", sep="\t", index=False)
    boundary_null.to_csv(tables / "microbial_boundary_permutation_null.tsv", sep="\t", index=False)

    modules = read_table(args.modules)
    assignment_metadata = joined[[args.sample_col, "microbial_compartment", *environmental_cols]].copy()
    module_profiles, module_summary, basin_module, module_wins = module_characterization(
        relative, assignment_metadata, modules, args.sample_col, environmental_cols,
    )
    module_profiles.to_csv(tables / "microbial_state_module_profiles.tsv", sep="\t", index=False)
    module_summary.to_csv(tables / "microbial_state_module_association_summary.tsv", sep="\t", index=False)
    basin_module.to_csv(tables / "basin_module_association.tsv", sep="\t", index=False)
    module_wins.to_csv(tables / "basin_module_best_strategy.tsv", sep="\t", index=False)
    module_statistics, module_homes = module_statistical_support(
        relative, joined, modules, args,
    )
    module_statistics.to_csv(
        tables / "microbial_state_module_association_statistics.tsv", sep="\t", index=False
    )
    module_homes.to_csv(
        tables / "ecological_module_home_microbial_compartment.tsv", sep="\t", index=False
    )
    common_module_statistics = module_statistics.loc[
        module_statistics.analysis_cohort.eq("common_hybrid_MC_samples")
    ]
    module_comparison = common_module_statistics.pivot(
        index="module_label", columns="grouping_strategy",
        values=["incremental_r_squared", "partial_r_squared", "adjusted_pvalue"],
    )
    module_comparison.columns = [f"{metric}__{strategy}" for metric, strategy in module_comparison.columns]
    module_comparison = module_comparison.reset_index()
    mc_increment = "incremental_r_squared__microbial_compartment"
    hybrid_increment = f"incremental_r_squared__{args.hybrid_col}"
    module_comparison["MC_minus_hybrid_incremental_r_squared"] = (
        module_comparison[mc_increment] - module_comparison[hybrid_increment]
    )
    module_comparison["larger_depth_adjusted_incremental_r_squared"] = np.where(
        module_comparison["MC_minus_hybrid_incremental_r_squared"] > 0,
        "microbial_compartment",
        np.where(module_comparison["MC_minus_hybrid_incremental_r_squared"] < 0, args.hybrid_col, "tie"),
    )
    module_comparison.to_csv(
        tables / "ecological_module_mc_hybrid_comparison.tsv", sep="\t", index=False
    )
    module_taxonomy = taxonomy_map(args.taxonomy)
    module_id_col = next(
        (column for column in ["Taxon", "ASV_ID", "ASV", "taxon"] if column in modules.columns),
        modules.columns[0],
    )
    module_inventory = modules.copy()
    module_inventory[module_id_col] = module_inventory[module_id_col].astype(str).str.strip()
    module_inventory = module_inventory.merge(module_homes, on="module_label", how="left", validate="many_to_one")
    module_inventory.insert(1, "taxonomy", module_inventory[module_id_col].map(module_taxonomy).fillna(""))
    module_inventory.to_csv(
        tables / "ecological_module_asv_membership_with_mc_home.tsv", sep="\t", index=False
    )
    renewal_state_support(relative, joined, modules, args, tables, plots)

    measurement_statistics = read_table(tables / "microbial_state_measurement_association.tsv")
    common_measurement_statistics = measurement_statistics.loc[
        measurement_statistics.analysis_cohort.eq("common_hybrid_MC_samples")
    ]
    measurement_comparison = common_measurement_statistics.pivot(
        index="measurement", columns="grouping_strategy",
        values=["incremental_r_squared", "partial_r_squared", "adjusted_pvalue"],
    )
    measurement_comparison.columns = [f"{metric}__{strategy}" for metric, strategy in measurement_comparison.columns]
    measurement_comparison = measurement_comparison.reset_index()
    measurement_comparison["MC_minus_hybrid_incremental_r_squared"] = (
        measurement_comparison[mc_increment] - measurement_comparison[hybrid_increment]
    )
    measurement_comparison["larger_depth_adjusted_incremental_r_squared"] = np.where(
        measurement_comparison["MC_minus_hybrid_incremental_r_squared"] > 0,
        "microbial_compartment",
        np.where(measurement_comparison["MC_minus_hybrid_incremental_r_squared"] < 0, args.hybrid_col, "tie"),
    )
    measurement_comparison.to_csv(
        tables / "microbial_state_measurement_mc_hybrid_comparison.tsv", sep="\t", index=False
    )

    integrated_rows = []
    for column in environmental_cols:
        corr = correspondence.loc[correspondence.environmental_variable.eq(column)].iloc[0]
        pred = prediction_summary.loc[prediction_summary.model.eq(column)].iloc[0]
        bound = boundary_summary.loc[boundary_summary.environmental_variable.eq(column)].iloc[0]
        module_values = basin_module.loc[basin_module.environmental_variable.eq(column)]
        integrated_rows.append({
            "environmental_variable": column,
            "adjusted_rand_index_with_MC": corr.adjusted_rand_index,
            "normalized_mutual_information_with_MC": corr.normalized_mutual_information,
            "bias_corrected_cramers_v_with_MC": corr.bias_corrected_cramers_v,
            "normalized_conditional_entropy_MC_given_environment": corr.normalized_conditional_entropy_mc_given_environment,
            "majority_mapping_accuracy": corr.majority_mapping_accuracy,
            "held_out_balanced_accuracy": pred.mean_balanced_accuracy,
            "held_out_balanced_accuracy_improvement_over_depth_season": pred.mean_improvement_balanced_accuracy,
            "held_out_macro_f1_improvement_over_depth_season": pred.mean_improvement_macro_f1,
            "held_out_log_loss_improvement_over_depth_season": pred.mean_improvement_log_loss,
            "boundary_sensitivity_for_MC_transition": bound.sensitivity,
            "boundary_precision_for_MC_transition": bound.precision,
            "boundary_f1": bound.f1_score,
            "boundary_jaccard": bound.jaccard_overlap,
            "boundary_probability_difference": bound.probability_difference,
            "median_module_adjusted_r_squared": module_values.adjusted_r_squared.median(),
            "modules_with_largest_adjusted_r_squared_n": int((module_wins.winning_environmental_variable == column).sum()),
        })
    integrated = pd.DataFrame(integrated_rows)
    integrated.to_csv(tables / "microbial_state_integrated_basin_comparison.tsv", sep="\t", index=False)
    write_mc_evidence_summary(
        args.microbial_dir, correspondence, prediction_summary, boundary_summary, tables
    )

    formats = parse_list(args.formats)
    metric_plot = correspondence.set_index("environmental_variable")
    correspondence_metrics = [
        ("adjusted_rand_index", "Adjusted Rand index"),
        ("normalized_mutual_information", "Normalized mutual information"),
        ("bias_corrected_cramers_v", "Bias-corrected Cramér's V"),
    ]
    with plt.rc_context(COMPACT_PLOT_STYLE):
        fig, ax = plt.subplots(figsize=(8.2, 4.8))
        metric_values = [float(metric_plot.iloc[0][column]) for column, _ in correspondence_metrics]
        metric_labels = [label for _, label in correspondence_metrics]
        y = np.arange(len(metric_labels))
        bars = ax.barh(y, metric_values, color="0.72", edgecolor="0.20", linewidth=0.8)
        ax.set_yticks(y, metric_labels)
        ax.invert_yaxis()
        ax.set_xlim(0, 1)
        ax.set_xlabel("Correspondence with microbial compartments")
        ax.set_title(r"Microbial and hybrid O$_2$-GMM compartment correspondence", pad=10)
        for bar, value in zip(bars, metric_values):
            ax.text(value + 0.018, bar.get_y() + bar.get_height() / 2, f"{value:.3f}", va="center")
        ax.spines[["top", "right"]].set_visible(False)
        fig.subplots_adjust(left=0.34, right=0.96, top=0.86, bottom=0.18)
        finalize_compact_figure(fig)
    save(fig, plots / "microbial_basin_correspondence_summary", formats)

    boundary_row = boundary_summary.iloc[0]
    with plt.rc_context(COMPACT_PLOT_STYLE):
        fig, ax = plt.subplots(figsize=(8.2, 4.8))
        probabilities = [
            float(boundary_row.mc_transition_probability_at_boundary),
            float(boundary_row.mc_transition_probability_within_compartment),
        ]
        labels = ["Hybrid boundary crossed", "Within hybrid compartment"]
        y = np.arange(len(labels))
        bars = ax.barh(y, probabilities, color=["0.30", "0.78"], edgecolor="0.15", linewidth=0.8)
        ax.set_yticks(y, labels)
        ax.invert_yaxis()
        ax.set_xlim(0, 1)
        ax.set_xlabel("Probability of an adjacent-depth MC transition")
        ax.set_title("Microbial-state transitions at hybrid compartment boundaries", pad=10)
        for bar, value in zip(bars, probabilities):
            ax.text(value + 0.018, bar.get_y() + bar.get_height() / 2, f"{value:.3f}", va="center")
        ax.text(
            0.98, 0.06,
            f"Difference = {boundary_row.probability_difference:.3f}; permutation p = "
            f"{boundary_row.probability_difference_permutation_pvalue:.3f}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=9,
        )
        ax.spines[["top", "right"]].set_visible(False)
        fig.subplots_adjust(left=0.34, right=0.96, top=0.86, bottom=0.18)
        finalize_compact_figure(fig)
    save(fig, plots / "microbial_boundary_concordance", formats)

    profile_matrix = module_profiles.pivot(index="module_label", columns="microbial_compartment", values="mean_standardized_module_abundance")
    profile_matrix = profile_matrix.reindex(sorted(profile_matrix.index, key=label_sort)).reindex(columns=sorted(profile_matrix.columns, key=label_sort))
    with plt.rc_context(COMPACT_PLOT_STYLE):
        fig, ax = plt.subplots(figsize=(8.2, max(6.4, 0.34 * len(profile_matrix) + 2)))
        limit = max(1, np.nanmax(np.abs(profile_matrix.to_numpy(float))))
        image = ax.imshow(
            profile_matrix.to_numpy(float), aspect="auto", cmap="RdBu_r",
            vmin=-limit, vmax=limit,
        )
        image._aspire_preserve_categorical_cmap = True
        ax.set_xticks(np.arange(len(profile_matrix.columns)), profile_matrix.columns)
        ax.set_yticks(np.arange(len(profile_matrix.index)), profile_matrix.index)
        ax.set_xlabel("Microbial compartment")
        ax.set_ylabel("Ecological module")
        ax.set_title("Standardized ecological-module abundance by microbial state", pad=10)
        fig.subplots_adjust(left=0.16, right=0.84, top=0.90, bottom=0.14)
        cax = fig.add_axes([0.88, 0.22, 0.025, 0.58])
        colorbar = fig.colorbar(
            image, cax=cax, label="Mean standardized module abundance"
        )
        # The publication contract converts generic heatmaps to grayscale. This
        # matrix is an explicitly diverging signed z-score display, so preserve
        # the matching diverging color scale on both the image and its colorbar.
        colorbar.solids._aspire_preserve_categorical_cmap = True
        finalize_compact_figure(fig)
    save(fig, plots / "microbial_state_module_heatmap", formats)

    curtain_outputs(joined, relative, taxonomy_map(args.taxonomy), args, tables, plots)

    audit_config = {
        **{key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
        "clustering_inputs": "not modified; all analyses in this process are post hoc",
        "common_comparison_cohort_samples": len(common),
        "common_comparison_cohort_cruises": common[args.cruise_col].nunique(),
        "module_association_interpretation": "descriptive because MCs and modules share the ASV input matrix",
        "software_versions": {
            "python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__,
            "scipy": scipy.__version__, "scikit_learn": sklearn.__version__,
        },
    }
    (audit / "microbial_state_interpretation_config.json").write_text(json.dumps(audit_config, indent=2))


if __name__ == "__main__":
    main()
