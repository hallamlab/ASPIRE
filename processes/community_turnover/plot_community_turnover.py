#!/usr/bin/env python3
"""Render publication-style community turnover and LCBD diagnostics."""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/aspire_matplotlib")
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
import seaborn as sns

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shared_plot_style import install_publication_style

install_publication_style()


def output_formats(value: str) -> list[str]:
    return [item.strip().lower() for item in re.split(r"[,|]", value) if item.strip()]


def save(fig: plt.Figure, base: Path, formats: list[str]) -> None:
    for fmt in formats:
        fig.savefig(base.with_suffix(f".{fmt}"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def read_optional(path: Path) -> pd.DataFrame:
    try:
        return pd.read_csv(path, sep="\t", low_memory=False)
    except (pd.errors.EmptyDataError, FileNotFoundError):
        return pd.DataFrame()


def format_calendar_axis(ax: plt.Axes) -> None:
    locator = mdates.AutoDateLocator(minticks=4, maxticks=8)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--module-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--primary-metric", default="aitchison")
    parser.add_argument("--depth-col", default="Depth")
    parser.add_argument("--date-col", default="date")
    parser.add_argument("--formats", default="pdf,png,svg")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    formats = output_formats(args.formats)
    tables = args.module_dir / "tables"

    vertical = read_optional(tables / "vertical_turnover.tsv")
    if not vertical.empty:
        data = vertical.loc[vertical.distance_metric.eq(args.primary_metric)].copy()
        data["date"] = pd.to_datetime(data.date, errors="coerce")
        data["midpoint_depth"] = pd.to_numeric(data.midpoint_depth, errors="coerce")
        data["community_distance"] = pd.to_numeric(data.community_distance, errors="coerce")
        data = data.dropna(subset=["profile_ID", "midpoint_depth", "community_distance"])
        profiles = data.profile_ID.dropna().astype(str).unique().tolist()

        # Interpolate each observed cruise trace only within its sampled midpoint
        # range. This provides a common display grid for the median and IQR while
        # retaining every unmodified cruise trace in the background.
        grid = np.arange(np.ceil(data.midpoint_depth.min()), np.floor(data.midpoint_depth.max()) + 1, 1.0)
        interpolated = []
        fig, ax = plt.subplots(figsize=(8.2, 7.2))
        for profile in profiles:
            subset = data.loc[data.profile_ID.astype(str).eq(profile)].sort_values("midpoint_depth")
            collapsed = subset.groupby("midpoint_depth", as_index=False).community_distance.mean()
            depth = collapsed.midpoint_depth.to_numpy(float)
            turnover = collapsed.community_distance.to_numpy(float)
            if len(depth) < 2:
                continue
            values = np.interp(grid, depth, turnover)
            values[(grid < depth.min()) | (grid > depth.max())] = np.nan
            interpolated.append(values)
            ax.plot(
                turnover, depth, color="0.72", linewidth=0.65, alpha=0.55,
                solid_capstyle="round", zorder=1,
            )
        if not interpolated:
            raise ValueError("No cruise contained enough adjacent-depth transitions for profile plotting")
        matrix = np.vstack(interpolated)
        median = np.nanmedian(matrix, axis=0)
        q25 = np.nanquantile(matrix, 0.25, axis=0)
        q75 = np.nanquantile(matrix, 0.75, axis=0)
        supported = np.isfinite(median) & np.isfinite(q25) & np.isfinite(q75)
        ax.fill_betweenx(
            grid[supported], q25[supported], q75[supported],
            color="0.72", alpha=0.45, linewidth=0, zorder=2,
        )
        ax.plot(
            median[supported], grid[supported], color="black", linewidth=3.2,
            solid_capstyle="round", label="Cross-cruise median", zorder=4,
        )

        boundary_counts = pd.DataFrame()
        if "crosses_compartment_boundary" in data:
            boundary = data.crosses_compartment_boundary.astype(str).str.strip().str.lower().map(
                {"true": True, "false": False}
            )
            boundary_counts = (
                data.assign(_boundary=boundary)
                .dropna(subset=["_boundary"])
                .groupby("midpoint_depth", as_index=False)
                .agg(boundary_crossing_cruises=("_boundary", "sum"))
            )
            boundary_counts = boundary_counts.loc[boundary_counts.boundary_crossing_cruises.gt(0)].copy()
            if not boundary_counts.empty:
                boundary_counts["summary_turnover"] = np.interp(
                    boundary_counts.midpoint_depth, grid[supported], median[supported]
                )
                boundary_counts["boundary_crossing_cruises"] = pd.to_numeric(
                    boundary_counts.boundary_crossing_cruises, errors="raise"
                ).astype(float)
                boundary_counts["marker_area"] = 30.0 + 22.0 * boundary_counts.boundary_crossing_cruises
                ax.scatter(
                    boundary_counts.summary_turnover, boundary_counts.midpoint_depth,
                    s=boundary_counts.marker_area, facecolor="white", edgecolor="black",
                    linewidth=1.2, zorder=6,
                )

        handles = [
            Line2D([0], [0], color="black", linewidth=3.2, label="Cross-cruise median"),
            Patch(facecolor="0.72", alpha=0.45, edgecolor="none", label="Interquartile range"),
            Line2D([0], [0], color="0.72", linewidth=0.8, label="Individual cruises"),
        ]
        if not boundary_counts.empty:
            observed_counts = sorted(boundary_counts.boundary_crossing_cruises.astype(int).unique())
            legend_counts = sorted(set([observed_counts[0], observed_counts[len(observed_counts) // 2], observed_counts[-1]]))
            handles.extend([
                plt.scatter([], [], s=30 + 22 * count, facecolor="white", edgecolor="black", linewidth=1.2,
                            label=f"{count} hybrid-boundary crossing{'s' if count != 1 else ''}")
                for count in legend_counts
            ])
        ax.invert_yaxis()
        ax.set_xlabel(f"Adjacent-depth turnover ({args.primary_metric.capitalize()} distance)")
        ax.set_ylabel("Midpoint depth (m)")
        ax.set_title("Vertical microbial community turnover")
        ax.legend(handles=handles, frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
        save(fig, args.outdir / "vertical_turnover_profiles", formats)

        contrast = read_optional(tables / "vertical_turnover_boundary_cruise_contrasts.tsv")
        summary = read_optional(tables / "vertical_turnover_boundary_summary.tsv")
        if not contrast.empty and not summary.empty:
            contrast = contrast.loc[contrast.distance_metric.eq(args.primary_metric)].copy()
            summary = summary.loc[summary.distance_metric.eq(args.primary_metric)].copy()
            strategy_order = [
                strategy for strategy in ["O2", "GMM", "Hybrid"]
                if strategy in set(summary.compartment_strategy.astype(str))
            ]
            strategy_order.extend(sorted(set(summary.compartment_strategy.astype(str)).difference(strategy_order)))
            if strategy_order:
                fig, ax = plt.subplots(figsize=(8.8, 4.8))
                rng = np.random.default_rng(42)
                for position, strategy in enumerate(strategy_order):
                    values = pd.to_numeric(
                        contrast.loc[
                            contrast.compartment_strategy.astype(str).eq(strategy),
                            "paired_turnover_difference",
                        ],
                        errors="coerce",
                    ).dropna()
                    row = summary.loc[summary.compartment_strategy.astype(str).eq(strategy)].iloc[0]
                    jitter = rng.uniform(-0.10, 0.10, len(values))
                    ax.scatter(
                        values, np.full(len(values), position) + jitter,
                        s=28, color="0.70", edgecolor="0.35", linewidth=0.35, alpha=0.8, zorder=2,
                    )
                    ax.hlines(
                        position, float(row.paired_difference_q25), float(row.paired_difference_q75),
                        color="black", linewidth=3.0, zorder=4,
                    )
                    ax.scatter(
                        float(row.median_paired_turnover_difference), position,
                        marker="D", s=78, color="black", edgecolor="white", linewidth=0.7, zorder=5,
                    )
                ax.axvline(0, color="0.25", linestyle="--", linewidth=1.0, zorder=1)
                strategy_labels = {"O2": r"O$_2$", "GMM": "GMM", "Hybrid": "Hybrid"}
                tick_labels = []
                for strategy in strategy_order:
                    label = strategy_labels.get(strategy, strategy)
                    row = summary.loc[summary.compartment_strategy.astype(str).eq(strategy)].iloc[0]
                    qvalue = pd.to_numeric(pd.Series([row.get("paired_sign_flip_qvalue", np.nan)]), errors="coerce").iloc[0]
                    if np.isfinite(qvalue):
                        label = f"{label} (q={qvalue:.3g})"
                    tick_labels.append(label)
                ax.set_yticks(
                    np.arange(len(strategy_order)),
                    tick_labels,
                )
                ax.invert_yaxis()
                ax.set_xlabel(
                    f"Cross-boundary minus within-compartment turnover\n({args.primary_metric.capitalize()} distance)"
                )
                ax.set_ylabel("Compartment strategy")
                ax.set_title("Microbial turnover at environmental compartment boundaries")
                handles = [
                    Line2D([0], [0], marker="o", linestyle="none", markerfacecolor="0.70",
                           markeredgecolor="0.35", markersize=6, label="Cruise-level paired contrast"),
                    Line2D([0], [0], marker="D", linestyle="none", color="black", markersize=7,
                           label="Median; line = IQR"),
                ]
                ax.legend(handles=handles, frameon=False, bbox_to_anchor=(1.02, 0.5), loc="center left")
                save(fig, args.outdir / "vertical_turnover_compartment_boundary_contrast", formats)

        heat = data.pivot_table(index="midpoint_depth", columns="date", values="community_distance", aggfunc="mean")
        if not heat.empty:
            fig, ax = plt.subplots(figsize=(max(9, 0.20 * heat.shape[1] + 4), 6.4))
            sns.heatmap(heat.sort_index(ascending=True), cmap="viridis", cbar_kws={"label": f"{args.primary_metric} turnover"}, ax=ax)
            tick_count = min(10, heat.shape[1])
            tick_indices = np.unique(np.linspace(0, heat.shape[1] - 1, tick_count, dtype=int))
            tick_dates = pd.to_datetime(heat.columns[tick_indices], errors="coerce")
            ax.set_xticks(tick_indices + 0.5)
            ax.set_xticklabels(
                [value.strftime("%Y-%m-%d") if not pd.isna(value) else "" for value in tick_dates],
                rotation=45,
                ha="right",
            )
            ax.invert_yaxis()
            ax.set_xlabel("Sampling date")
            ax.set_ylabel("Midpoint depth")
            ax.set_title("Adjacent-depth microbial turnover through time")
            save(fig, args.outdir / "vertical_turnover_depth_time_heatmap", formats)

    fixed = read_optional(tables / "temporal_fixed_depth_turnover.tsv")
    if not fixed.empty:
        data = fixed.loc[fixed.distance_metric.eq(args.primary_metric)].copy()
        data["date_t2"] = pd.to_datetime(data.date_t2, errors="coerce")
        fig, ax = plt.subplots(figsize=(11, 6.2))
        for depth, subset in data.groupby("fixed_depth", sort=True):
            ax.plot(subset.date_t2, subset.community_distance, marker="o", linewidth=0.9, markersize=3.5, label=f"{depth:g} m")
        ax.set_xlabel("Later sampling date")
        ax.set_ylabel(f"Temporal turnover ({args.primary_metric})")
        ax.set_title("Exact-depth microbial turnover")
        format_calendar_axis(ax)
        ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left", ncol=2)
        save(fig, args.outdir / "temporal_fixed_depth_turnover", formats)

    compartment = read_optional(tables / "temporal_environmental_compartment_turnover.tsv")
    if not compartment.empty:
        data = compartment.loc[compartment.distance_metric.eq(args.primary_metric)].copy()
        data["date_t2"] = pd.to_datetime(data.date_t2, errors="coerce")
        for variable, subset_variable in data.groupby("environmental_variable", sort=True):
            fig, ax = plt.subplots(figsize=(11, 5.8))
            for state, subset in subset_variable.groupby("environmental_compartment", sort=True):
                ax.plot(subset.date_t2, subset.community_distance, marker="o", linewidth=1, markersize=4, label=str(state))
            ax.set_xlabel("Later sampling date")
            ax.set_ylabel(f"Temporal turnover ({args.primary_metric})")
            ax.set_title(f"Compartment-following microbial turnover: {variable}")
            format_calendar_axis(ax)
            ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
            slug = re.sub(r"[^A-Za-z0-9._-]+", "_", str(variable)).strip("_")
            save(fig, args.outdir / f"temporal_compartment_turnover_{slug}", formats)

    lcbd = read_optional(tables / "lcbd.tsv")
    if not lcbd.empty:
        global_lcbd = lcbd.loc[lcbd.reference_population.eq("global")].copy()
        depth_column = args.depth_col if args.depth_col in global_lcbd.columns else None
        date_column = args.date_col if args.date_col in global_lcbd.columns else None
        if depth_column:
            global_lcbd[depth_column] = pd.to_numeric(global_lcbd[depth_column], errors="coerce")
            fig, ax = plt.subplots(figsize=(7.4, 6.2))
            significant = pd.to_numeric(global_lcbd.adjusted_pvalue, errors="coerce").le(0.05)
            ax.scatter(global_lcbd.LCBD, global_lcbd[depth_column], s=30, color="0.65", edgecolor="black", linewidth=0.35)
            ax.scatter(global_lcbd.loc[significant, "LCBD"], global_lcbd.loc[significant, depth_column],
                       s=52, color="#b2182b", edgecolor="black", linewidth=0.5, label="BH q ≤ 0.05")
            ax.invert_yaxis()
            ax.set_xlabel("Global LCBD")
            ax.set_ylabel("Depth")
            ax.set_title("Sample compositional uniqueness by depth")
            if significant.any():
                ax.legend(frameon=False)
            save(fig, args.outdir / "global_lcbd_by_depth", formats)
        if depth_column and date_column:
            global_lcbd[date_column] = pd.to_datetime(global_lcbd[date_column], errors="coerce")
            fig, ax = plt.subplots(figsize=(11, 6.2))
            points = ax.scatter(global_lcbd[date_column], global_lcbd[depth_column], c=global_lcbd.LCBD,
                                cmap="magma", s=42, edgecolor="black", linewidth=0.35)
            ax.invert_yaxis()
            ax.set_xlabel("Sampling date")
            ax.set_ylabel("Depth")
            ax.set_title("Global LCBD across depth and time")
            fig.colorbar(points, ax=ax, label="Global LCBD")
            format_calendar_axis(ax)
            save(fig, args.outdir / "global_lcbd_depth_time", formats)


if __name__ == "__main__":
    main()
