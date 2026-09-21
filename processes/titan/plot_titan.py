#!/usr/bin/env python3
"""Create publication-style plots from stable, collected TITAN2 tables."""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/aspire_matplotlib")

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shared_plot_style import install_publication_style

install_publication_style()

DIRECTION_COLORS = {"z-": "#2166ac", "z+": "#b2182b"}


def as_bool(value: object) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    return str(value).strip().lower() in {"true", "t", "1", "yes", "y"}


def parse_formats(value: str) -> list[str]:
    return [item.strip().lower() for item in re.split(r"[,|]", value) if item.strip()]


def slug(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._-") or "measurement"


def save(fig: plt.Figure, base: Path, formats: list[str]) -> None:
    for fmt in formats:
        fig.savefig(base.with_suffix(f".{fmt}"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_distribution(taxa: pd.DataFrame, variable: str, output: Path, formats: list[str]) -> None:
    fig, ax = plt.subplots(figsize=(8.2, 5.3))
    for direction in ("z-", "z+"):
        values = pd.to_numeric(
            taxa.loc[taxa["response_direction"].eq(direction), "change_point"], errors="coerce"
        ).dropna()
        if values.empty:
            continue
        bins = min(30, max(6, int(np.sqrt(len(values)) * 2)))
        ax.hist(values, bins=bins, density=True, histtype="stepfilled", alpha=0.25,
                color=DIRECTION_COLORS[direction], label=f"{direction} (n={len(values)})")
        if len(values) >= 3 and values.nunique() >= 3:
            sns.kdeplot(x=values, ax=ax, color=DIRECTION_COLORS[direction], linewidth=1.6)
        ax.plot(values, np.full(len(values), -0.006 if direction == "z-" else -0.012), "|",
                color=DIRECTION_COLORS[direction], markersize=6, alpha=0.65)
    ax.set_xlabel(variable)
    ax.set_ylabel("Change-point density")
    ax.set_title(f"TITAN taxon change points: {variable}")
    ax.legend(frameon=False)
    save(fig, output, formats)


def plot_community(curve: pd.DataFrame, threshold: pd.Series, variable: str,
                   output: Path, formats: list[str]) -> None:
    fig, ax = plt.subplots(figsize=(8.4, 5.4))
    x = pd.to_numeric(curve["environmental_value"], errors="coerce")
    ax.plot(x, curve["sum_z_minus"], color=DIRECTION_COLORS["z-"], linewidth=1.6, label="sum(z−)")
    ax.plot(x, curve["sum_z_plus"], color=DIRECTION_COLORS["z+"], linewidth=1.6, label="sum(z+)")
    if pd.to_numeric(curve["filtered_sum_z_minus"], errors="coerce").notna().any():
        ax.plot(x, curve["filtered_sum_z_minus"], color=DIRECTION_COLORS["z-"], linewidth=1.0,
                linestyle="--", label="filtered sum(z−)")
    if pd.to_numeric(curve["filtered_sum_z_plus"], errors="coerce").notna().any():
        ax.plot(x, curve["filtered_sum_z_plus"], color=DIRECTION_COLORS["z+"], linewidth=1.0,
                linestyle="--", label="filtered sum(z+)")
    for direction, prefix in (("z-", "z_minus"), ("z+", "z_plus")):
        center = pd.to_numeric(pd.Series([threshold.get(f"{prefix}_threshold")]), errors="coerce").iloc[0]
        low = pd.to_numeric(pd.Series([threshold.get(f"{prefix}_bootstrap_q05")]), errors="coerce").iloc[0]
        high = pd.to_numeric(pd.Series([threshold.get(f"{prefix}_bootstrap_q95")]), errors="coerce").iloc[0]
        if np.isfinite(low) and np.isfinite(high):
            ax.axvspan(low, high, color=DIRECTION_COLORS[direction], alpha=0.08, linewidth=0)
        if np.isfinite(center):
            ax.axvline(center, color=DIRECTION_COLORS[direction], linestyle=":", linewidth=1.1)
    ax.set_xlabel(variable)
    ax.set_ylabel("Summed indicator z-score")
    ax.set_title(f"TITAN community response: {variable}")
    ax.legend(frameon=False, ncol=2)
    save(fig, output, formats)


def plot_ranked(taxa: pd.DataFrame, variable: str, output: Path,
                formats: list[str], top_n: int) -> None:
    data = taxa.copy()
    data["z_score"] = pd.to_numeric(data["z_score"], errors="coerce")
    data["change_point"] = pd.to_numeric(data["change_point"], errors="coerce")
    data = data.dropna(subset=["change_point", "z_score"])
    if top_n > 0 and len(data) > top_n:
        data = data.sort_values(["passes_purity_and_reliability", "z_score"], ascending=False).head(top_n)
    data = data.sort_values(["response_direction", "change_point", "ASV_ID"])
    height = max(5.5, 0.22 * len(data) + 2.0)
    fig, ax = plt.subplots(figsize=(9.0, height))
    y = np.arange(len(data))
    low = pd.to_numeric(data["bootstrap_cp_q05"], errors="coerce")
    high = pd.to_numeric(data["bootstrap_cp_q95"], errors="coerce")
    center = data["change_point"].to_numpy(float)
    for index, row in enumerate(data.itertuples(index=False)):
        color = DIRECTION_COLORS.get(row.response_direction, "0.45")
        if np.isfinite(low.iloc[index]) and np.isfinite(high.iloc[index]):
            ax.plot([low.iloc[index], high.iloc[index]], [index, index], color=color, linewidth=0.8)
        face = color if as_bool(row.passes_purity_and_reliability) else "white"
        ax.scatter(center[index], index, s=25, facecolor=face, edgecolor=color, linewidth=0.8, zorder=3)
    ax.set_yticks(y, data["ASV_ID"].astype(str), fontsize=6.5)
    ax.set_xlabel(variable)
    ax.set_ylabel("ASV")
    ax.set_title(f"Ranked TITAN taxon change points: {variable}")
    ax.invert_yaxis()
    save(fig, output, formats)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tables-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--formats", default="pdf,png,svg")
    parser.add_argument("--ranked-top-n", type=int, default=50)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    formats = parse_formats(args.formats)
    taxa = pd.read_csv(args.tables_dir / "titan_taxon_results.tsv", sep="\t", low_memory=False)
    thresholds = pd.read_csv(args.tables_dir / "titan_community_thresholds.tsv", sep="\t", low_memory=False)
    curves = pd.read_csv(args.tables_dir / "titan_community_response_curve.tsv", sep="\t", low_memory=False)
    for variable in thresholds["environmental_variable"].dropna().astype(str):
        variable_taxa = taxa.loc[taxa["environmental_variable"].astype(str).eq(variable)].copy()
        variable_curve = curves.loc[curves["environmental_variable"].astype(str).eq(variable)].copy()
        threshold = thresholds.loc[thresholds["environmental_variable"].astype(str).eq(variable)].iloc[0]
        base = args.outdir / slug(variable)
        plot_distribution(variable_taxa, variable, base.with_name(base.name + "_taxon_change_point_distribution"), formats)
        plot_community(variable_curve, threshold, variable, base.with_name(base.name + "_community_response"), formats)
        plot_ranked(variable_taxa, variable, base.with_name(base.name + "_taxon_change_points_ranked"), formats, args.ranked_top_n)


if __name__ == "__main__":
    main()
