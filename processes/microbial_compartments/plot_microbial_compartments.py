#!/usr/bin/env python3
"""Render microbial compartment diagnostics and explicitly post hoc comparisons."""

from __future__ import annotations

import argparse
import json
import os
import platform
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/aspire_matplotlib")
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import sklearn
from sklearn.decomposition import PCA

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shared_plot_style import install_publication_style

install_publication_style()


def formats(value: str) -> list[str]:
    return [part.strip().lower() for part in re.split(r"[,|]", value) if part.strip()]


def slug(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._-") or "group"


def save(fig: plt.Figure, path: Path, output_formats: list[str]) -> None:
    for fmt in output_formats:
        fig.savefig(path.with_suffix(f".{fmt}"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def palette(labels: list[str]) -> dict[str, tuple[float, float, float]]:
    colors = sns.color_palette("colorblind", max(3, len(labels)))
    return {label: colors[index] for index, label in enumerate(labels)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--module-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--sample-col", default="sampleID")
    parser.add_argument("--depth-col", default="Depth")
    parser.add_argument("--date-col", default="date")
    parser.add_argument("--formats", default="pdf,png,svg")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    output_formats = formats(args.formats)
    tables = args.module_dir / "tables"
    audit = args.module_dir / "audit"

    primary = pd.read_csv(tables / "microbial_compartments.tsv", sep="\t")
    best = pd.read_csv(tables / "microbial_best_candidate_assignments.tsv", sep="\t")
    validation = pd.read_csv(tables / "cluster_validation.tsv", sep="\t")
    decision = pd.read_csv(tables / "microbial_cluster_selection_decision.tsv", sep="\t").iloc[0]
    clr = pd.read_csv(audit / "asv_clr_matrix.tsv", sep="\t")
    clr_id = clr.columns[0]
    clr = clr.set_index(clr_id)
    scores = PCA(n_components=min(2, clr.shape[0] - 1, clr.shape[1]), svd_solver="full").fit_transform(clr)
    score_table = pd.DataFrame({args.sample_col: clr.index, "PC1": scores[:, 0]})
    score_table["PC2"] = scores[:, 1] if scores.shape[1] > 1 else 0.0
    score_table = score_table.merge(primary, on=args.sample_col, how="left", validate="one_to_one")
    supported = str(decision.support_status) == "supported"

    fig, ax = plt.subplots(figsize=(7.8, 6.3))
    groups = sorted(score_table.microbial_compartment.astype(str).unique())
    colors = palette(groups) if supported else {"unresolved": (0.45, 0.45, 0.45)}
    for group in groups:
        subset = score_table.loc[score_table.microbial_compartment.astype(str).eq(group)]
        ax.scatter(subset.PC1, subset.PC2, s=30, color=colors.get(group, "0.5"),
                   edgecolor="black", linewidth=0.35, label=group)
    ax.set_xlabel("CLR-PC1")
    ax.set_ylabel("CLR-PC2")
    ax.set_title("Microbial compartment CLR-PCA" if supported else "CLR-PCA: no supported discrete microbial compartments")
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    save(fig, args.outdir / "microbial_compartment_clr_pca", output_formats)

    pam = validation.loc[validation.cluster_method.eq("PAM")].sort_values("K")
    fig, ax1 = plt.subplots(figsize=(8.1, 5.2))
    ax1.plot(pam.K, pam.primary_stability_quantile, marker="s", color="#b2182b", label="Lower-tail within-cruise ARI")
    ax1.plot(pam.K, pam.minimum_cluster_median_jaccard, marker="^", color="#ef8a62", label="Minimum cluster median Jaccard")
    if "median_prediction_strength" in pam:
        ax1.plot(pam.K, pam.median_prediction_strength, marker="D", color="#762a83", label="Season-balanced prediction strength")
    if "median_season_balanced_block_ARI" in pam:
        ax1.plot(pam.K, pam.median_season_balanced_block_ARI, marker="v", color="#2166ac", label="Median season-balanced blocked ARI")
    ax1.plot(pam.K, pam.mean_silhouette, marker="o", color="0.45", label="Mean silhouette (descriptive)")
    ax1.axhline(float(decision.required_primary_stability_quantile), color="#b2182b", linestyle="--", linewidth=1)
    ax1.axhline(float(decision.required_minimum_cluster_median_jaccard), color="#ef8a62", linestyle="--", linewidth=1)
    ax1.set_xlabel("Candidate K")
    ax1.set_ylabel("Validation statistic")
    ax1.set_ylim(min(-0.05, float(pam.mean_silhouette.min()) - 0.05), 1.02)
    selected_k = pd.to_numeric(pd.Series([decision.selected_K]), errors="coerce").iloc[0]
    best_k = pd.to_numeric(pd.Series([decision.best_candidate_K]), errors="coerce").iloc[0]
    marked_k = selected_k if np.isfinite(selected_k) else best_k
    if np.isfinite(marked_k):
        ax1.axvline(marked_k, color="0.25", linestyle=":", linewidth=1)
    ax1.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    ax1.set_title(f"PAM stability and blocked robustness ({decision.support_status})")
    save(fig, args.outdir / "microbial_cluster_validation", output_formats)

    plot_assignments = primary[[args.sample_col, "microbial_compartment", "silhouette_width_for_sample"]].copy()
    title = "Selected PAM compartment silhouettes"
    output_name = "microbial_compartment_silhouette"
    if not supported:
        plot_assignments = best.rename(columns={
            "best_candidate_compartment": "microbial_compartment",
            "best_candidate_silhouette_width": "silhouette_width_for_sample",
        })[["sample_ID", "microbial_compartment", "silhouette_width_for_sample"]]
        plot_assignments = plot_assignments.rename(columns={"sample_ID": args.sample_col})
        title = "Best PAM candidate silhouettes (weak; not selected)"
        output_name = "microbial_compartment_silhouette_weak_candidate"
    plot_assignments = plot_assignments.sort_values(["microbial_compartment", "silhouette_width_for_sample"])
    fig, ax = plt.subplots(figsize=(8.4, 5.5))
    y = np.arange(len(plot_assignments))
    group_colors = palette(sorted(plot_assignments.microbial_compartment.astype(str).unique()))
    ax.barh(y, plot_assignments.silhouette_width_for_sample,
            color=[group_colors[str(group)] for group in plot_assignments.microbial_compartment], height=0.9)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_yticks([])
    ax.set_xlabel("Sample silhouette width")
    ax.set_ylabel("Samples grouped by compartment")
    ax.set_title(title)
    save(fig, args.outdir / output_name, output_formats)

    metrics = pd.read_csv(tables / "microbial_environmental_comparison_metrics.tsv", sep="\t")
    contingency = pd.read_csv(tables / "microbial_environmental_contingency.tsv", sep="\t")
    if supported and not contingency.empty:
        for variable in contingency.environmental_variable.unique():
            subset = contingency.loc[contingency.environmental_variable.eq(variable)]
            matrix = subset.pivot(index="microbial_compartment", columns="environmental_compartment", values="number_of_samples").fillna(0)
            fig_width = max(7.0, 0.55 * matrix.shape[1] + 3.0)
            fig_height = max(4.5, 0.48 * matrix.shape[0] + 2.2)
            fig, ax = plt.subplots(figsize=(fig_width, fig_height))
            sns.heatmap(matrix, cmap="Greys", annot=True, fmt=".0f", cbar_kws={"label": "Samples"}, ax=ax)
            metric = metrics.loc[metrics.environmental_variable.eq(variable)].iloc[0]
            ax.set_title(
                f"Microbial vs {variable}\nARI={metric.adjusted_rand_index:.3f}; NMI={metric.normalized_mutual_information:.3f}"
            )
            ax.set_xlabel(variable)
            ax.set_ylabel("Microbial compartment")
            save(fig, args.outdir / f"microbial_vs_{slug(variable)}_contingency", output_formats)

    posthoc = pd.read_csv(tables / "microbial_compartment_posthoc_metadata.tsv", sep="\t")
    if supported and args.depth_col in posthoc.columns and args.date_col in posthoc.columns:
        posthoc[args.date_col] = pd.to_datetime(posthoc[args.date_col], errors="coerce")
        posthoc[args.depth_col] = pd.to_numeric(posthoc[args.depth_col], errors="coerce")
        usable = posthoc.dropna(subset=[args.date_col, args.depth_col])
        if not usable.empty:
            groups = sorted(usable.microbial_compartment.astype(str).unique())
            colors = palette(groups)
            fig, ax = plt.subplots(figsize=(11.0, 5.8))
            for group in groups:
                subset = usable.loc[usable.microbial_compartment.astype(str).eq(group)]
                ax.scatter(subset[args.date_col], subset[args.depth_col], s=28,
                           color=colors[group], edgecolor="black", linewidth=0.35, label=group)
            ax.invert_yaxis()
            ax.set_xlabel("Sampling date")
            ax.set_ylabel("Depth")
            ax.set_title("Post hoc distribution of microbial compartments")
            ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
            save(fig, args.outdir / "microbial_compartment_depth_time_posthoc", output_formats)

    plot_config = {
        "formats": output_formats,
        "sample_col": args.sample_col,
        "depth_col": args.depth_col,
        "date_col": args.date_col,
        "clustering_supported": supported,
        "software_versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "matplotlib": mpl.__version__,
            "seaborn": sns.__version__,
            "scikit_learn": sklearn.__version__,
        },
    }
    (audit / "microbial_clustering_plot_config.json").write_text(json.dumps(plot_config, indent=2))


if __name__ == "__main__":
    main()
