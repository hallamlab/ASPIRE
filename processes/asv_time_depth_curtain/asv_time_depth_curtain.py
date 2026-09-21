#!/usr/bin/env python3
"""Render a smoothed hybrid-compartment curtain from ASPIRE ASV metadata."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import is_color_like
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.ndimage import gaussian_filter

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shared_plot_style import install_publication_style

install_publication_style()


def read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t" if path.suffix.lower() in {".tsv", ".txt"} else ",")


def parse_palette(spec: str) -> dict[str, str]:
    palette: dict[str, str] = {}
    for item in str(spec).split(","):
        if "=" not in item:
            continue
        label, color = item.split("=", 1)
        if label.strip() and color.strip():
            palette[label.strip()] = color.strip()
    return palette


def validate_palette(palette: dict[str, str], labels: list[str], raw_spec: str) -> None:
    if not palette:
        raise ValueError(
            "The configured palette contained no label=color entries. "
            f"Received: {raw_spec!r}"
        )
    missing_colors = [label for label in labels if label not in palette]
    if missing_colors:
        raise ValueError(
            "The configured palette has no color for observed curtain states: "
            + ", ".join(missing_colors)
        )
    invalid_colors = [label for label in labels if not is_color_like(palette[label])]
    if invalid_colors:
        raise ValueError(
            "The configured palette has invalid colors for curtain states: "
            + ", ".join(invalid_colors)
        )


def label_sort(label: str) -> tuple[int, int, str]:
    oxygen = str(label).split("__", 1)[0].lower()
    oxygen_order = {"oxic": 0, "dysoxic": 1, "suboxic": 2, "anoxic": 3}
    match = re.search(r"gmm(\d+)", str(label), flags=re.IGNORECASE)
    return oxygen_order.get(oxygen, 99), int(match.group(1)) if match else 999, str(label)


def display_label(label: str) -> str:
    match = re.fullmatch(r"([^_]+)__gmm(\d+)", str(label), flags=re.IGNORECASE)
    return f"{match.group(1)}-GMM{match.group(2)}" if match else str(label).replace("__", "-")


def draw_categorical_contours(
    ax: plt.Axes,
    x: np.ndarray,
    depth: np.ndarray,
    codes: np.ndarray,
    colors: list[str],
) -> None:
    """Render independently filled categorical regions with smooth-looking edges."""
    for code, color in enumerate(colors):
        mask = (codes == code).astype(float).T
        if not mask.any():
            continue
        ax.contourf(
            x, depth, mask,
            levels=[0.5, 1.5], colors=[color],
            antialiased=True, zorder=1,
        )
        if mask.min() < 0.5 < mask.max():
            ax.contour(
                x, depth, mask,
                levels=[0.5], colors=["#4D4D4D"],
                linewidths=0.25, alpha=0.55, zorder=2,
            )


def modal_profile(
    metadata: pd.DataFrame,
    cruise_col: str,
    depth_col: str,
    group_col: str,
    order: list[str],
) -> pd.DataFrame:
    rank = {label: index for index, label in enumerate(order)}
    rows = []
    for (cruise_index, depth), frame in metadata.groupby(["cruise_index", depth_col], sort=True):
        counts = frame[group_col].astype(str).value_counts()
        maximum = counts.max()
        candidates = counts.loc[counts.eq(maximum)].index.tolist()
        selected = sorted(candidates, key=lambda label: (rank.get(label, 999), label))[0]
        rows.append({
            "cruise_index": int(cruise_index),
            "depth_m": float(depth),
            "hybrid_compartment": selected,
            "samples_at_depth_n": int(len(frame)),
            "selected_label_samples_n": int(counts[selected]),
            "selection_tied": len(candidates) > 1,
        })
    return pd.DataFrame(rows)


def build_surface(
    profile: pd.DataFrame,
    cruise_dates: pd.Series,
    labels: list[str],
    maximum_depth: float,
    depth_step: float,
    time_subdivisions_per_month: int,
    time_sigma_months: float,
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray, np.ndarray, np.ndarray, pd.DatetimeIndex, pd.Timestamp, pd.Timestamp]:
    """Create a calendar-time surface with smoothing restricted to time.

    Observed states define midpoint-bounded vertical layers within each cruise.
    One-hot state support is interpolated between cruise dates and Gaussian-
    smoothed along calendar time only; no smoothing is applied across depth.
    """
    label_to_code = {label: index for index, label in enumerate(labels)}
    depth_centers = np.arange(0.0, maximum_depth + depth_step * 0.5, depth_step)
    depth_centers[-1] = min(depth_centers[-1], maximum_depth)
    dates = pd.to_datetime(cruise_dates, errors="coerce")
    if dates.isna().any():
        raise ValueError("All cruises require valid dates for calendar-time smoothing")
    n_cruises = len(dates)
    calendar_start = pd.Timestamp(year=int(dates.dt.year.min()), month=1, day=1)
    calendar_end = pd.Timestamp(year=int(dates.dt.year.max()) + 1, month=1, day=1)
    month_starts = pd.date_range(calendar_start, calendar_end, freq="MS")
    subdivisions = max(1, int(time_subdivisions_per_month))
    fine_dates_list: list[pd.Timestamp] = []
    for left, right in zip(month_starts[:-1], month_starts[1:]):
        width = right - left
        fine_dates_list.extend(
            left + width * ((index + 0.5) / subdivisions)
            for index in range(subdivisions)
        )
    fine_dates = pd.DatetimeIndex(fine_dates_list)
    fine_x = (fine_dates - calendar_start).total_seconds().to_numpy() / 86400.0
    coarse_x = (dates - calendar_start).dt.total_seconds().to_numpy() / 86400.0
    coarse = np.full((n_cruises, len(depth_centers)), np.nan)
    for cruise_index, frame in profile.groupby("cruise_index", sort=True):
        frame = frame.sort_values("depth_m")
        depths = frame["depth_m"].to_numpy(float)
        codes = frame["hybrid_compartment"].map(label_to_code).to_numpy(int)
        if len(depths) == 1:
            selected = np.zeros(len(depth_centers), dtype=int)
        else:
            selected = np.searchsorted(
                (depths[:-1] + depths[1:]) / 2.0, depth_centers, side="right"
            )
        coarse[int(cruise_index)] = codes[selected]

    support = np.zeros((len(labels), len(fine_x), len(depth_centers)))
    for code in range(len(labels)):
        for depth_index in range(len(depth_centers)):
            available = np.isfinite(coarse[:, depth_index])
            if not available.any():
                continue
            values = (coarse[available, depth_index] == code).astype(float)
            dated = pd.DataFrame({"x": coarse_x[available], "value": values})
            dated = dated.groupby("x", as_index=False)["value"].mean().sort_values("x")
            support[code, :, depth_index] = np.interp(fine_x, dated["x"], dated["value"])
        support[code] = gaussian_filter(
            support[code],
            sigma=(max(0.0, time_sigma_months) * subdivisions, 0.0),
            mode="nearest",
        )
    total = support.sum(axis=0)
    codes = support.argmax(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        confidence = np.divide(
            support.max(axis=0), total,
            out=np.zeros_like(total), where=total > 0,
        )
    grid = pd.DataFrame({
        "calendar_date": np.repeat(fine_dates.to_numpy(), len(depth_centers)),
        "calendar_month": np.repeat(fine_dates.to_period("M").astype(str), len(depth_centers)),
        "display_time_coordinate_days": np.repeat(fine_x, len(depth_centers)),
        "depth_m": np.tile(depth_centers, len(fine_x)),
        "display_state_code": codes.ravel(),
        "hybrid_compartment_display": np.asarray(labels, dtype=object)[codes.ravel()],
        "maximum_smoothed_support_fraction": confidence.ravel(),
    })
    return grid, fine_x, depth_centers, codes, support, fine_dates, calendar_start, calendar_end


def enforce_observed_anchors(
    render_codes: np.ndarray,
    fine_x: np.ndarray,
    depth_centers: np.ndarray,
    profile: pd.DataFrame,
    cruise_dates: pd.Series,
    labels: list[str],
    calendar_start: pd.Timestamp,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Impose midpoint-bounded observed profiles after visual smoothing."""
    anchored = np.asarray(render_codes, dtype=int).copy()
    anchor_mask = np.zeros_like(anchored, dtype=bool)
    anchor_owner = np.full_like(anchored, -1, dtype=int)
    owner_distance = np.full(len(fine_x), np.inf, dtype=float)
    label_to_code = {label: index for index, label in enumerate(labels)}
    dates = pd.to_datetime(cruise_dates, errors="coerce")

    for cruise_index, frame in profile.groupby("cruise_index", sort=True):
        frame = frame.sort_values("depth_m").drop_duplicates("depth_m", keep="first")
        depths = frame["depth_m"].to_numpy(float)
        observed_codes = frame["hybrid_compartment"].map(label_to_code).to_numpy(int)
        if len(depths) == 1:
            selected = np.zeros(len(depth_centers), dtype=int)
        else:
            selected = np.searchsorted(
                (depths[:-1] + depths[1:]) / 2.0,
                depth_centers,
                side="right",
            )
        depth_profile_codes = observed_codes[selected]
        sample_date = dates.iloc[int(cruise_index)]
        sample_x = float((sample_date - calendar_start).total_seconds() / 86400.0)
        insertion = int(np.searchsorted(fine_x, sample_x, side="left"))
        candidate_columns = {
            max(0, min(len(fine_x) - 1, insertion - 1)),
            max(0, min(len(fine_x) - 1, insertion)),
        }
        if insertion < len(fine_x) and np.isclose(fine_x[insertion], sample_x):
            candidate_columns.add(min(len(fine_x) - 1, insertion + 1))
        for column in sorted(candidate_columns):
            distance = abs(float(fine_x[column]) - sample_x)
            if distance > owner_distance[column]:
                continue
            anchored[column, :] = depth_profile_codes
            anchor_mask[column, :] = True
            anchor_owner[column, :] = int(cruise_index)
            owner_distance[column] = distance
    return anchored, anchor_mask, anchor_owner


def render_basin_surface(args: argparse.Namespace, metadata: pd.DataFrame) -> None:
    """Render BASIN's canonical curtain and emphasize BASIN cells with ASV data."""
    grid = read_table(args.basin_grid)
    cells = read_table(args.basin_cells)
    required_grid = {
        "calendar_date", "display_time_coordinate_days", "depth_m",
        "o2_subcompartment_rendered",
    }
    required_cells = {
        args.cruise_col,
        "date",
        "measurement_depth_m",
        "o2_subcompartment_final",
        "included_in_smoothed_display_support",
        "rendered_state_matches_observation",
    }
    missing_grid = sorted(required_grid.difference(grid.columns))
    missing_cells = sorted(required_cells.difference(cells.columns))
    if missing_grid:
        raise ValueError(
            "BASIN curtain grid lacks the canonical rendered-state columns: "
            + ", ".join(missing_grid)
            + ". Rerun BASIN from BIOCHEM_SPLIT_O2_BY_GMM first."
        )
    if missing_cells:
        raise ValueError(
            "BASIN curtain cells are missing the observed-anchor audit columns: "
            + ", ".join(missing_cells)
            + ". Rerun BASIN from BIOCHEM_SPLIT_O2_BY_GMM first."
        )

    grid = grid.copy()
    grid["calendar_date"] = pd.to_datetime(grid["calendar_date"], errors="coerce")
    grid["display_time_coordinate_days"] = pd.to_numeric(
        grid["display_time_coordinate_days"], errors="coerce"
    )
    grid["depth_m"] = pd.to_numeric(grid["depth_m"], errors="coerce")
    grid = grid.dropna(subset=[
        "calendar_date", "display_time_coordinate_days", "depth_m",
        "o2_subcompartment_rendered",
    ])
    configured_order = [item.strip() for item in args.group_order.split(",") if item.strip()]
    observed = sorted(grid["o2_subcompartment_rendered"].astype(str).unique(), key=label_sort)
    order = [label for label in configured_order if label in observed]
    order.extend(label for label in observed if label not in order)
    palette = parse_palette(args.palette)
    validate_palette(palette, order, args.palette)
    label_to_code = {label: index for index, label in enumerate(order)}

    fine_x = np.sort(grid["display_time_coordinate_days"].unique())
    depth_centers = np.sort(grid["depth_m"].unique())
    surface = grid.pivot_table(
        index="display_time_coordinate_days", columns="depth_m",
        values="o2_subcompartment_rendered", aggfunc="first",
    ).reindex(index=fine_x, columns=depth_centers)
    if surface.isna().any().any():
        raise ValueError("BASIN rendered curtain grid is not a complete time-by-depth surface")
    codes = surface.astype(str).replace(label_to_code).to_numpy(dtype=int)
    grid["hybrid_compartment_display"] = grid["o2_subcompartment_rendered"].astype(str)
    calendar_start = pd.Timestamp(year=int(grid["calendar_date"].dt.year.min()), month=1, day=1)
    calendar_end = pd.Timestamp(year=int(grid["calendar_date"].dt.year.max()) + 1, month=1, day=1)

    cells = cells.copy()
    cells["date"] = pd.to_datetime(cells["date"], errors="coerce").dt.normalize()
    cells["measurement_depth_m"] = pd.to_numeric(
        cells["measurement_depth_m"], errors="coerce"
    )
    cells = cells.dropna(subset=[args.cruise_col, "date", "measurement_depth_m"])
    cells["time_coordinate_days"] = (
        cells["date"] - calendar_start
    ).dt.total_seconds() / 86400.0

    asv = metadata.copy()
    required_asv = {
        args.sample_col, args.cruise_col, args.date_col, args.depth_col, args.group_col,
    }
    missing_asv = sorted(required_asv.difference(asv.columns))
    if missing_asv:
        raise ValueError("ASV curtain metadata is missing: " + ", ".join(missing_asv))
    asv[args.date_col] = pd.to_datetime(asv[args.date_col], errors="coerce").dt.normalize()
    asv[args.depth_col] = pd.to_numeric(asv[args.depth_col], errors="coerce")
    asv = asv.dropna(subset=[args.cruise_col, args.date_col, args.depth_col]).copy()
    asv[args.group_col] = asv[args.group_col].astype("string").str.strip()

    def normalized_cruise(series: pd.Series) -> pd.Series:
        return series.astype(str).str.strip().str.replace(r"\.0$", "", regex=True)

    cells["_match_cruise"] = normalized_cruise(cells[args.cruise_col])
    cells["_match_depth"] = cells["measurement_depth_m"].round(6)
    asv["_match_cruise"] = normalized_cruise(asv[args.cruise_col])
    asv["_match_depth"] = asv[args.depth_col].round(6)
    key_columns_cells = ["_match_cruise", "date", "_match_depth"]
    key_columns_asv = ["_match_cruise", args.date_col, "_match_depth"]
    asv_keys = {
        tuple(row) for row in asv[key_columns_asv].drop_duplicates().itertuples(index=False, name=None)
    }
    basin_keys = {
        tuple(row) for row in cells[key_columns_cells].drop_duplicates().itertuples(index=False, name=None)
    }
    cells["asv_data_available"] = [
        tuple(row) in asv_keys
        for row in cells[key_columns_cells].itertuples(index=False, name=None)
    ]
    rank = {label: index for index, label in enumerate(order)}

    def choose_asv_compartment(series: pd.Series) -> object:
        labels = series.dropna().astype(str)
        labels = labels.loc[labels.str.len().gt(0)]
        if labels.empty:
            return pd.NA
        counts = labels.value_counts()
        candidates = counts.loc[counts.eq(counts.max())].index.tolist()
        return sorted(candidates, key=lambda label: (rank.get(label, 999), label))[0]

    asv_compartments = (
        asv.groupby(key_columns_asv, dropna=False)[args.group_col]
        .agg(choose_asv_compartment)
        .rename("asv_sample_compartment")
    )
    asv_compartment_lookup = asv_compartments.to_dict()
    cells["asv_sample_compartment"] = [
        asv_compartment_lookup.get(tuple(row), pd.NA)
        for row in cells[key_columns_cells].itertuples(index=False, name=None)
    ]
    missing_asv_compartment = cells["asv_data_available"] & cells["asv_sample_compartment"].isna()
    if missing_asv_compartment.any():
        raise ValueError(
            f"{int(missing_asv_compartment.sum())} BASIN positions with ASV data lack "
            f"a {args.group_col!r} label"
        )
    observed_asv_compartments = sorted(
        cells.loc[cells["asv_data_available"], "asv_sample_compartment"].astype(str).unique(),
        key=label_sort,
    )
    validate_palette(palette, observed_asv_compartments, args.palette)
    basin_anchor_ok = cells["rendered_state_matches_observation"].astype(str).str.lower().map(
        {"true": True, "false": False}
    )
    basin_included = cells["included_in_smoothed_display_support"].astype(str).str.lower().map(
        {"true": True, "false": False}
    )
    if basin_anchor_ok.isna().any() or basin_included.isna().any() or not basin_anchor_ok[basin_included].all():
        raise ValueError(
            "The BASIN curtain failed its observed-anchor audit; rerun BASIN from "
            "BIOCHEM_SPLIT_O2_BY_GMM before rendering the ASPIRE curtain"
        )
    cells["asv_compartment_matches_basin_observation"] = (
        ~cells["asv_data_available"]
        | ~basin_included
        | (
            cells["asv_sample_compartment"].astype(str)
            == cells["o2_subcompartment_final"].astype(str)
        )
    )
    if not cells["asv_compartment_matches_basin_observation"].all():
        raise ValueError(
            "Matched ASV metadata contain hybrid-compartment labels that differ from "
            "the canonical BASIN observations"
        )

    renewal_rows = pd.DataFrame(columns=["renewal_date", "time_coordinate_days"])
    if args.renewal_events:
        renewal = read_table(args.renewal_events)
        if args.renewal_date_col not in renewal.columns:
            raise ValueError(
                f"Renewal-event table lacks configured onset column {args.renewal_date_col!r}"
            )
        renewal_dates = pd.to_datetime(
            renewal[args.renewal_date_col], errors="coerce"
        ).dropna().drop_duplicates().sort_values()
        renewal_rows = pd.DataFrame({"renewal_date": renewal_dates})
        renewal_rows["time_coordinate_days"] = (
            renewal_rows["renewal_date"] - calendar_start
        ).dt.total_seconds() / 86400.0
        renewal_rows = renewal_rows.loc[
            renewal_rows["time_coordinate_days"].between(
                0.0, float((calendar_end - calendar_start).days), inclusive="both"
            )
        ].reset_index(drop=True)

    tables = args.outdir / "tables"
    plots = args.outdir / "plots"
    tables.mkdir(parents=True, exist_ok=True)
    plots.mkdir(parents=True, exist_ok=True)
    renewal_rows = pd.DataFrame(columns=["renewal_date", "time_coordinate_days"])
    if args.renewal_events:
        renewal = read_table(args.renewal_events)
        if args.renewal_date_col not in renewal.columns:
            raise ValueError(
                f"Renewal-event table lacks configured onset column {args.renewal_date_col!r}"
            )
        renewal_dates = pd.to_datetime(
            renewal[args.renewal_date_col], errors="coerce"
        ).dropna().drop_duplicates().sort_values()
        renewal_rows = pd.DataFrame({"renewal_date": renewal_dates})
        renewal_rows["time_coordinate_days"] = (
            renewal_rows["renewal_date"] - calendar_start
        ).dt.total_seconds() / 86400.0
        renewal_rows = renewal_rows.loc[
            renewal_rows["time_coordinate_days"].between(
                0.0, float((calendar_end - calendar_start).days), inclusive="both"
            )
        ].reset_index(drop=True)
    renewal_rows.to_csv(
        tables / "asv_hybrid_compartment_time_depth_renewals.tsv", sep="\t", index=False
    )
    grid.to_csv(tables / "asv_hybrid_compartment_time_depth_grid.tsv", sep="\t", index=False)
    cells.drop(columns=["_match_cruise", "_match_depth"]).to_csv(
        tables / "asv_basin_time_depth_point_matches.tsv", sep="\t", index=False
    )
    renewal_rows.to_csv(
        tables / "asv_hybrid_compartment_time_depth_renewals.tsv", sep="\t", index=False
    )
    summary_row = {
        "surface_source": "BASIN canonical rendered curtain grid",
        "basin_grid": str(args.basin_grid),
        "basin_cells": str(args.basin_cells),
        "basin_measurement_positions_n": int(len(cells)),
        "basin_positions_with_asv_data_n": int(cells["asv_data_available"].sum()),
        "basin_positions_without_asv_data_n": int((~cells["asv_data_available"]).sum()),
        "asv_samples_matching_basin_positions_n": int(
            sum(
                tuple(row) in basin_keys
                for row in asv[key_columns_asv].itertuples(index=False, name=None)
            )
        ),
        "renewal_onsets_plotted_n": len(renewal_rows),
        "small_point_size": args.base_point_size,
        "asv_point_size": args.asv_point_size,
        "asv_point_edge_width": args.asv_point_edge_width,
        "maximum_display_depth_m": args.maximum_depth,
    }
    pd.DataFrame([summary_row]).to_csv(
        tables / "asv_hybrid_compartment_time_depth_summary.tsv", sep="\t", index=False
    )

    fig, ax = plt.subplots(figsize=(18.0, 7.5))
    plot_x = np.concatenate(([0.0], fine_x, [float((calendar_end - calendar_start).days)]))
    plot_codes = np.concatenate((codes[:1], codes, codes[-1:]), axis=0)
    draw_categorical_contours(
        ax, plot_x, depth_centers, plot_codes, [palette[label] for label in order]
    )
    ax.scatter(
        cells["time_coordinate_days"], cells["measurement_depth_m"],
        s=args.base_point_size, facecolor="black", edgecolor="none", zorder=4,
    )
    with_asv = cells["asv_data_available"]
    ax.scatter(
        cells.loc[with_asv, "time_coordinate_days"],
        cells.loc[with_asv, "measurement_depth_m"],
        s=args.asv_point_size,
        facecolor=[
            palette[label]
            for label in cells.loc[with_asv, "asv_sample_compartment"].astype(str)
        ],
        edgecolor="black", linewidth=args.asv_point_edge_width, zorder=5,
    )
    for renewal_x in renewal_rows["time_coordinate_days"]:
        ax.axvline(
            renewal_x, color="black", linestyle="--", linewidth=0.9,
            alpha=0.9, zorder=3,
        )
    year_starts = pd.date_range(calendar_start, calendar_end, freq="YS", inclusive="left")
    year_x = (year_starts - calendar_start).total_seconds().to_numpy() / 86400.0
    month_starts = pd.date_range(calendar_start, calendar_end, freq="MS", inclusive="left")
    month_x = (month_starts - calendar_start).total_seconds().to_numpy() / 86400.0
    ax.set_xticks(year_x)
    ax.set_xticklabels([str(date.year) for date in year_starts], ha="left", fontsize=8)
    ax.set_xticks(month_x, minor=True)
    ax.tick_params(axis="x", which="minor", length=2.5)
    ax.set_xlim(0.0, float((calendar_end - calendar_start).days))
    ax.set_ylim(args.maximum_depth, 0.0)
    ax.set_yticks(np.arange(0.0, args.maximum_depth + 1.0, 25.0))
    ax.set_xlabel("Calendar time (monthly intervals)")
    ax.set_ylabel("Depth (m)")
    ax.set_title(r"Hybrid O$_2$-GMM compartments through time and depth")
    handles = [
        Patch(facecolor=palette[label], edgecolor="none", label=display_label(label))
        for label in order
    ]
    handles.extend([
        Line2D([], [], marker="o", linestyle="", markersize=3.5,
               markerfacecolor="black", markeredgecolor="none", label="BASIN sampled depth"),
        Line2D([], [], marker="o", linestyle="", markersize=6.0,
               markerfacecolor="white", markeredgecolor="black",
               markeredgewidth=args.asv_point_edge_width,
               label="ASV sample (compartment-colored)"),
    ])
    if not renewal_rows.empty:
        handles.append(Line2D(
            [], [], color="black", linestyle="--", linewidth=0.9,
            label="Predicted renewal onset",
        ))
    ax.legend(
        handles=handles, title=r"Hybrid O$_2$-GMM compartment",
        loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False,
        fontsize=8, title_fontsize=9,
    )
    fig.subplots_adjust(left=0.07, right=0.82, bottom=0.20, top=0.91)
    for fmt in [item.strip().lower() for item in args.formats.split(",") if item.strip()]:
        fig.savefig(
            plots / f"asv_hybrid_compartment_time_depth_curtain.{fmt}",
            dpi=300 if fmt == "png" else None, bbox_inches="tight",
        )
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--sample-col", default="sampleID")
    parser.add_argument("--cruise-col", default="Cruise")
    parser.add_argument("--date-col", default="Date")
    parser.add_argument("--depth-col", default="Depth")
    parser.add_argument("--group-col", default="o2_subcompartment_final")
    parser.add_argument("--palette", required=True)
    parser.add_argument("--group-order", default="")
    parser.add_argument("--exclude-label-pattern", default="(?i)(?:other|outlier)")
    parser.add_argument("--maximum-depth", type=float, default=210.0)
    parser.add_argument("--depth-step", type=float, default=1.0)
    parser.add_argument("--time-subdivisions-per-month", type=int, default=4)
    parser.add_argument("--time-sigma-months", type=float, default=0.75)
    parser.add_argument("--contour-visual-depth-sigma-m", type=float, default=5.0)
    parser.add_argument("--basin-grid", type=Path)
    parser.add_argument("--basin-cells", type=Path)
    parser.add_argument("--renewal-events", type=Path)
    parser.add_argument("--renewal-date-col", default="start_date")
    parser.add_argument("--base-point-size", type=float, default=7.5)
    parser.add_argument("--asv-point-size", type=float, default=42.0)
    parser.add_argument("--asv-point-edge-width", type=float, default=1.2)
    parser.add_argument("--formats", default="pdf,png,svg")
    args = parser.parse_args()

    metadata = read_table(args.metadata)
    if bool(args.basin_grid) != bool(args.basin_cells):
        raise ValueError("--basin-grid and --basin-cells must be provided together")
    if args.basin_grid:
        render_basin_surface(args, metadata)
        return
    required = {args.sample_col, args.cruise_col, args.date_col, args.depth_col, args.group_col}
    missing = sorted(required.difference(metadata.columns))
    if missing:
        raise ValueError("ASV curtain metadata is missing: " + ", ".join(missing))
    data = metadata.copy()
    data[args.sample_col] = data[args.sample_col].astype(str)
    data[args.depth_col] = pd.to_numeric(data[args.depth_col], errors="coerce")
    data["_date"] = pd.to_datetime(data[args.date_col], errors="coerce")
    data = data.loc[
        data[args.depth_col].between(0.0, args.maximum_depth, inclusive="both")
        & data[args.cruise_col].notna() & data["_date"].notna()
    ].copy()
    if data.empty:
        raise ValueError("No dated ASV samples were available for the time-depth curtain")

    cruises = data[[args.cruise_col, "_date"]].drop_duplicates().sort_values(
        ["_date", args.cruise_col]
    ).reset_index(drop=True)
    cruises["cruise_index"] = np.arange(len(cruises), dtype=int)
    calendar_start = pd.Timestamp(year=int(cruises["_date"].dt.year.min()), month=1, day=1)
    cruises["time_coordinate_days"] = (
        cruises["_date"] - calendar_start
    ).dt.total_seconds() / 86400.0
    data = data.merge(cruises, on=[args.cruise_col, "_date"], validate="many_to_one")
    data["excluded_from_surface"] = (
        data[args.group_col].isna()
        | data[args.group_col].astype(str).str.contains(args.exclude_label_pattern, regex=True, na=True)
    )
    eligible = data.loc[~data["excluded_from_surface"]].copy()
    if eligible.empty:
        raise ValueError("No classified non-excluded ASV samples were available")

    palette = parse_palette(args.palette)
    configured_order = [item.strip() for item in args.group_order.split(",") if item.strip()]
    labels = sorted(eligible[args.group_col].astype(str).unique(), key=label_sort)
    order = [label for label in configured_order if label in labels]
    order.extend(label for label in labels if label not in order)
    validate_palette(palette, order, args.palette)
    profile = modal_profile(eligible, args.cruise_col, args.depth_col, args.group_col, order)
    grid, fine_x, depth_centers, codes, support, fine_dates, calendar_start, calendar_end = build_surface(
        profile, cruises.set_index("cruise_index")["_date"], order, args.maximum_depth,
        max(0.1, args.depth_step), max(1, args.time_subdivisions_per_month),
        max(0.0, args.time_sigma_months),
    )

    tables = args.outdir / "tables"
    plots = args.outdir / "plots"
    tables.mkdir(parents=True, exist_ok=True)
    plots.mkdir(parents=True, exist_ok=True)
    data[[
        args.sample_col, args.cruise_col, "_date", "cruise_index", args.depth_col,
        args.group_col, "excluded_from_surface",
    ]].rename(columns={"_date": "date"}).to_csv(
        tables / "asv_hybrid_compartment_time_depth_samples.tsv", sep="\t", index=False
    )
    profile.to_csv(tables / "asv_hybrid_compartment_time_depth_profile.tsv", sep="\t", index=False)
    grid.to_csv(tables / "asv_hybrid_compartment_time_depth_grid.tsv", sep="\t", index=False)
    summary_row = {
        "asv_samples_n": len(data),
        "classified_surface_samples_n": len(eligible),
        "excluded_or_unclassified_samples_n": int(data["excluded_from_surface"].sum()),
        "cruises_n": len(cruises),
        "display_states_n": len(order),
        "calendar_start": calendar_start,
        "calendar_end_exclusive": calendar_end,
        "calendar_months_n": (calendar_end.year - calendar_start.year) * 12,
        "time_subdivisions_per_calendar_month": max(1, args.time_subdivisions_per_month),
        "time_gaussian_sigma_months": max(0.0, args.time_sigma_months),
        "depth_grid_step_m": max(0.1, args.depth_step),
        "maximum_display_depth_m": args.maximum_depth,
        "depth_gaussian_sigma_m": 0.0,
        "smoothing_axes": "calendar_time_only",
        "rendering": "categorical_filled_contours",
        "contour_visual_depth_sigma_m": max(0.0, args.contour_visual_depth_sigma_m),
        "contour_visual_smoothing_changes_assignments": False,
        "excluded_label_pattern": args.exclude_label_pattern,
        "renewal_onsets_plotted_n": len(renewal_rows),
    }

    fig, ax = plt.subplots(figsize=(15.5, 7.2))
    plot_x = np.concatenate((
        [0.0], fine_x, [float((calendar_end - calendar_start).days)]
    ))
    contour_support = np.stack([
        gaussian_filter(
            state_support,
            sigma=(0.0, max(0.0, args.contour_visual_depth_sigma_m) / max(0.1, args.depth_step)),
            mode="nearest",
        )
        for state_support in support
    ])
    contour_codes = contour_support.argmax(axis=0)
    contour_codes, anchor_mask, anchor_owner = enforce_observed_anchors(
        contour_codes,
        fine_x,
        depth_centers,
        profile,
        cruises.set_index("cruise_index")["_date"],
        order,
        calendar_start,
    )
    grid["display_state_code"] = contour_codes.ravel()
    grid["hybrid_compartment_display"] = np.asarray(order, dtype=object)[
        contour_codes.ravel()
    ]
    grid["observed_anchor_enforced"] = anchor_mask.ravel()
    grid["anchor_cruise_index"] = np.where(anchor_mask, anchor_owner, np.nan).ravel()
    grid.to_csv(tables / "asv_hybrid_compartment_time_depth_grid.tsv", sep="\t", index=False)

    rendered_at_profile = []
    for row in profile.itertuples(index=False):
        sample_x = float(
            cruises.loc[cruises["cruise_index"].eq(row.cruise_index), "time_coordinate_days"].iloc[0]
        )
        time_index = int(np.argmin(np.abs(fine_x - sample_x)))
        depth_index = int(np.argmin(np.abs(depth_centers - float(row.depth_m))))
        rendered_at_profile.append(order[int(contour_codes[time_index, depth_index])])
    profile["rendered_state_at_measurement"] = rendered_at_profile
    profile["rendered_state_matches_observation"] = (
        profile["rendered_state_at_measurement"].astype(str)
        == profile["hybrid_compartment"].astype(str)
    )
    if not profile["rendered_state_matches_observation"].all():
        raise RuntimeError("Observed-anchor enforcement failed for the ASV-derived curtain")
    summary_row.update({
        "observed_anchor_constraint": True,
        "anchor_grid_cells_n": int(anchor_mask.sum()),
        "included_observation_mismatches_n": int(
            (~profile["rendered_state_matches_observation"]).sum()
        ),
    })
    pd.DataFrame([summary_row]).to_csv(
        tables / "asv_hybrid_compartment_time_depth_summary.tsv", sep="\t", index=False
    )
    profile.to_csv(
        tables / "asv_hybrid_compartment_time_depth_profile.tsv", sep="\t", index=False
    )
    plot_codes = np.concatenate(
        (contour_codes[:1], contour_codes, contour_codes[-1:]), axis=0
    )
    draw_categorical_contours(
        ax, plot_x, depth_centers, plot_codes, [palette[label] for label in order]
    )
    for renewal_x in renewal_rows["time_coordinate_days"]:
        ax.axvline(
            renewal_x, color="black", linestyle="--", linewidth=0.9,
            alpha=0.9, zorder=2.5,
        )
    ax.scatter(
        data["time_coordinate_days"], data[args.depth_col], s=args.asv_point_size,
        facecolor=[palette.get(str(label), "#BDBDBD") for label in data[args.group_col]],
        edgecolor="black", linewidth=args.asv_point_edge_width, zorder=4,
    )
    year_starts = pd.date_range(calendar_start, calendar_end, freq="YS", inclusive="left")
    year_x = (year_starts - calendar_start).total_seconds().to_numpy() / 86400.0
    month_starts = pd.date_range(calendar_start, calendar_end, freq="MS", inclusive="left")
    month_x = (month_starts - calendar_start).total_seconds().to_numpy() / 86400.0
    ax.set_xticks(year_x)
    ax.set_xticklabels([str(date.year) for date in year_starts], rotation=0, ha="left", fontsize=8)
    ax.set_xticks(month_x, minor=True)
    ax.tick_params(axis="x", which="minor", length=2.5)
    ax.set_xlim(0.0, float((calendar_end - calendar_start).days))
    ax.set_ylim(args.maximum_depth, 0.0)
    ax.set_yticks(np.arange(0.0, args.maximum_depth + 1.0, 25.0))
    ax.set_xlabel("Calendar time (monthly intervals)")
    ax.set_ylabel("Depth (m)")
    ax.set_title(r"ASV samples across hybrid O$_2$-GMM compartments through time and depth")
    handles = [Patch(facecolor=palette[label], edgecolor="none", label=display_label(label)) for label in order]
    handles.append(Line2D([], [], marker="o", linestyle="", markersize=6.0,
                          markerfacecolor="white", markeredgecolor="black",
                          markeredgewidth=args.asv_point_edge_width,
                          label="ASV sample (compartment-colored)"))
    if not renewal_rows.empty:
        handles.append(Line2D([], [], color="black", linestyle="--", linewidth=0.9,
                              label="Predicted renewal onset"))
    ax.legend(handles=handles, title=r"Hybrid O$_2$-GMM compartment",
              loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False,
              fontsize=8, title_fontsize=9)
    fig.subplots_adjust(left=0.07, right=0.81, bottom=0.20, top=0.90)
    for fmt in [item.strip().lower() for item in args.formats.split(",") if item.strip()]:
        fig.savefig(plots / f"asv_hybrid_compartment_time_depth_curtain.{fmt}",
                    dpi=300 if fmt == "png" else None, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
