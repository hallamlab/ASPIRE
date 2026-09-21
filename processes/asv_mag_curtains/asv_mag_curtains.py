#!/usr/bin/env python3
"""Overlay species-level MAG recruitment on environmental-state curtains."""

from __future__ import annotations

import argparse
import itertools
import json
import math
import platform
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd


PROFILE_PATTERN = re.compile(
    r"^(?P<prefix>.+)_(?P<modality>metagenome|metatranscriptome)_"
    r"(?P<asv>ASV[^_]+)_species_recruitment_by_sample\.tsv$"
)


def read_table(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t", low_memory=False)


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
    date_min: pd.Timestamp,
    date_max: pd.Timestamp,
) -> None:
    for onset in onsets[(onsets >= date_min) & (onsets <= date_max)]:
        ax.axvline(
            onset, color="black", linestyle="--", linewidth=0.9,
            alpha=0.9, zorder=2.5,
        )


def renewal_legend_handle() -> Line2D:
    return Line2D(
        [], [], color="black", linestyle="--", linewidth=0.9,
        label="Predicted renewal onset",
    )


def parse_list(value: str) -> list[str]:
    return [item.strip() for item in str(value).split(",") if item.strip()]


def parse_palette(value: str) -> dict[str, str]:
    palette: dict[str, str] = {}
    for item in parse_list(value):
        if "=" not in item:
            raise ValueError(f"Palette item lacks '=': {item}")
        label, color = item.split("=", 1)
        palette[label.strip()] = color.strip()
    return palette


def canonical_asv(value: object) -> str:
    return re.sub(r";.*$", "", str(value).strip())


def canonical_species(value: object) -> str:
    return re.sub(r"\s+", " ", str(value).replace("_", " ").strip())


def safe_name(value: object) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_") or "unknown"


def asv_sort_key(value: str) -> tuple[int, str]:
    match = re.fullmatch(r"ASV(\d+)", str(value), flags=re.IGNORECASE)
    return (int(match.group(1)), str(value)) if match else (10**12, str(value))


def load_eligible_pairs(path: Path) -> pd.DataFrame:
    frame = read_table(path)
    required = {"ASV_ID", "mag_match_id"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Eligible mapping table lacks columns: {sorted(missing)}")
    species_col = "mag_tax_species" if "mag_tax_species" in frame.columns else "mag_species"
    if species_col not in frame.columns:
        raise ValueError("Eligible mapping table lacks MAG species taxonomy")
    if "analysis_eligible_pair" in frame.columns:
        eligible = frame["analysis_eligible_pair"].astype(str).str.lower().isin({"true", "1", "yes"})
        frame = frame.loc[eligible].copy()
    frame["ASV_ID"] = frame["ASV_ID"].map(canonical_asv)
    frame["species"] = frame[species_col].map(canonical_species)
    frame["mag_match_id"] = frame["mag_match_id"].astype(str).str.strip()
    frame = frame.loc[
        frame["ASV_ID"].ne("") & frame["species"].ne("") & frame["mag_match_id"].ne("")
    ].copy()
    return frame[["ASV_ID", "species", "mag_match_id"]].drop_duplicates()


def load_recruitment_profiles(abundance_dir: Path, prefix: str) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    for path in sorted(abundance_dir.glob(f"{prefix}_*_ASV*_species_recruitment_by_sample.tsv")):
        match = PROFILE_PATTERN.match(path.name)
        if not match:
            continue
        frame = read_table(path)
        required = {
            "ASV_ID", "species", "mag_match_id", "sample", "date", "depth_m",
            "raw_recruitment_count", "normalized_recruitment", "normalization",
        }
        missing = required.difference(frame.columns)
        if missing:
            raise ValueError(f"{path.name} lacks columns: {sorted(missing)}")
        frame = frame[list(required)].copy()
        frame["modality"] = match.group("modality")
        frame["source_file"] = path.name
        rows.append(frame)
    if not rows:
        raise ValueError(f"No species recruitment profiles found in {abundance_dir}")
    frame = pd.concat(rows, ignore_index=True)
    frame["ASV_ID"] = frame["ASV_ID"].map(canonical_asv)
    frame["species"] = frame["species"].map(canonical_species)
    frame["mag_match_id"] = frame["mag_match_id"].astype(str).str.strip()
    frame["sample"] = frame["sample"].astype(str).str.strip()
    frame["date"] = pd.to_datetime(frame["date"], errors="coerce")
    frame["depth_m"] = pd.to_numeric(frame["depth_m"], errors="coerce")
    frame["raw_recruitment_count"] = pd.to_numeric(
        frame["raw_recruitment_count"], errors="coerce"
    )
    frame["normalized_recruitment"] = pd.to_numeric(
        frame["normalized_recruitment"], errors="coerce"
    )
    return frame.dropna(subset=["date", "depth_m", "normalized_recruitment"])


def aggregate_species_points(profiles: pd.DataFrame, eligible: pd.DataFrame) -> pd.DataFrame:
    linked = profiles.merge(
        eligible.assign(_eligible=True),
        on=["ASV_ID", "species", "mag_match_id"],
        how="inner",
        validate="many_to_one",
    )
    if linked.empty:
        raise ValueError("No recruitment rows matched analysis-eligible ASV-MAG links")
    keys = ["ASV_ID", "species", "modality", "sample", "date", "depth_m", "normalization"]
    points = (
        linked.groupby(keys, as_index=False, dropna=False)
        .agg(
            linked_genome_count=("mag_match_id", "nunique"),
            mean_normalized_recruitment_fpm=("normalized_recruitment", "mean"),
            median_normalized_recruitment_fpm=("normalized_recruitment", "median"),
            minimum_normalized_recruitment_fpm=("normalized_recruitment", "min"),
            maximum_normalized_recruitment_fpm=("normalized_recruitment", "max"),
            mean_raw_recruitment_count=("raw_recruitment_count", "mean"),
        )
    )
    return points.sort_values(["ASV_ID", "species", "modality", "date", "depth_m"])


def load_surface(path: Path, order: list[str]) -> dict[str, object]:
    frame = read_table(path)
    required = {"calendar_date", "depth_m", "rendered_state"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Curtain grid {path} lacks columns: {sorted(missing)}")
    frame["calendar_date"] = pd.to_datetime(frame["calendar_date"], errors="coerce")
    frame["depth_m"] = pd.to_numeric(frame["depth_m"], errors="coerce")
    frame["rendered_state"] = frame["rendered_state"].astype(str)
    frame = frame.dropna(subset=["calendar_date", "depth_m"])
    dates = np.array(sorted(frame["calendar_date"].unique()), dtype="datetime64[ns]")
    depths = np.array(sorted(frame["depth_m"].unique()), dtype=float)
    observed = set(frame["rendered_state"])
    labels = [label for label in order if label in observed]
    labels.extend(sorted(observed.difference(labels)))
    code = {label: index for index, label in enumerate(labels)}
    pivot = frame.pivot(index="depth_m", columns="calendar_date", values="rendered_state")
    pivot = pivot.reindex(index=depths, columns=pd.to_datetime(dates))
    if pivot.isna().any().any():
        raise ValueError(f"Curtain grid {path} is incomplete")
    codes = pivot.apply(lambda column: column.map(code)).to_numpy(dtype=int)
    return {"dates": pd.to_datetime(dates), "depths": depths, "labels": labels, "codes": codes}


def load_sample_anchors(
    path: Path,
    sample_col: str,
    date_col: str,
    depth_col: str,
    state_col: str,
    allowed_states: list[str],
    state_surface: str,
) -> pd.DataFrame:
    """Load the actual ASV sampling locations represented by a curtain surface."""
    frame = read_table(path)
    required = {sample_col, date_col, depth_col, state_col}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Sample crosswalk lacks anchor columns: {sorted(missing)}")
    anchors = frame[[sample_col, date_col, depth_col, state_col]].copy()
    anchors.columns = ["sample_ID", "date", "depth_m", "state"]
    anchors["sample_ID"] = anchors["sample_ID"].astype(str).str.strip()
    anchors["date"] = pd.to_datetime(anchors["date"], errors="coerce")
    anchors["depth_m"] = pd.to_numeric(anchors["depth_m"], errors="coerce")
    anchors["state"] = anchors["state"].astype(str).str.strip()
    anchors = anchors.dropna(subset=["date", "depth_m"])
    anchors = anchors.loc[
        anchors["sample_ID"].ne("") & anchors["state"].isin(set(allowed_states))
    ].drop_duplicates(subset=["sample_ID"])
    anchors["state_surface"] = state_surface
    return anchors.sort_values(["date", "depth_m", "sample_ID"]).reset_index(drop=True)


def draw_sample_anchors(ax: plt.Axes, anchors: pd.DataFrame, point_area: float) -> None:
    if anchors.empty:
        return
    ax.scatter(
        anchors["date"], anchors["depth_m"], s=point_area,
        marker="o", facecolor="black", edgecolor="black", linewidth=0,
        alpha=0.72, zorder=3,
    )


def draw_surface(ax: plt.Axes, surface: dict[str, object], palette: dict[str, str]) -> None:
    dates = pd.DatetimeIndex(surface["dates"])
    depths = np.asarray(surface["depths"], dtype=float)
    codes = np.asarray(surface["codes"], dtype=int)
    labels = list(surface["labels"])
    x = matplotlib.dates.date2num(dates.to_pydatetime())
    for index, label in enumerate(labels):
        mask = (codes == index).astype(float)
        if not np.any(mask):
            continue
        ax.contourf(x, depths, mask, levels=[0.5, 1.5], colors=[palette[label]], zorder=1)
        if mask.min() < 0.5 < mask.max():
            ax.contour(
                x, depths, mask, levels=[0.5], colors=["#4D4D4D"],
                linewidths=0.25, alpha=0.55, zorder=2,
            )


def display_state(value: str) -> str:
    return str(value).replace("__", "–").replace("_", " ")


def size_scale(values: pd.Series, maximum: float, minimum_area: float, maximum_area: float) -> np.ndarray:
    numeric = pd.to_numeric(values, errors="coerce").fillna(0).clip(lower=0).to_numpy(float)
    if maximum <= 0:
        return np.full(len(numeric), minimum_area)
    ceiling = nice_ceiling(maximum)
    scaled = maximum_area * numeric / ceiling
    return np.where(numeric > 0, np.maximum(minimum_area, scaled), minimum_area)


def contrast_size_scale(
    values: pd.Series, maximum_absolute_log2_ratio: float,
    minimum_area: float, maximum_area: float,
) -> np.ndarray:
    numeric = pd.to_numeric(values, errors="coerce").fillna(0).abs().to_numpy(float)
    if maximum_absolute_log2_ratio <= 0:
        return np.full(len(numeric), minimum_area)
    fold_difference = np.power(2.0, numeric)
    fold_ceiling = power_of_two_ceiling(2.0 ** maximum_absolute_log2_ratio)
    denominator = max(fold_ceiling - 1.0, 1.0)
    scaled = np.clip((fold_difference - 1.0) / denominator, 0.0, 1.0)
    return minimum_area + (maximum_area - minimum_area) * scaled


def calculate_species_contrasts(
    points: pd.DataFrame, pseudocount_fpm: float = 1.0, fold_threshold: float = 2.0,
) -> pd.DataFrame:
    if pseudocount_fpm <= 0 or fold_threshold <= 1:
        raise ValueError("Contrast pseudocount must be positive and fold threshold exceed one")
    result: list[pd.DataFrame] = []
    index_columns = ["ASV_ID", "modality", "sample", "date", "depth_m", "normalization"]
    for asv, asv_rows in points.groupby("ASV_ID", sort=False):
        species = sorted(asv_rows["species"].dropna().astype(str).unique())
        for numerator, denominator in itertools.combinations(species, 2):
            selected = asv_rows.loc[asv_rows["species"].isin([numerator, denominator])]
            wide = selected.pivot_table(
                index=index_columns,
                columns="species",
                values="mean_normalized_recruitment_fpm",
                aggfunc="first",
            ).reset_index()
            if numerator not in wide.columns or denominator not in wide.columns:
                continue
            wide = wide.dropna(subset=[numerator, denominator]).copy()
            if wide.empty:
                continue
            wide = wide.rename(columns={
                numerator: "numerator_mean_fpm",
                denominator: "denominator_mean_fpm",
            })
            wide["numerator_species"] = numerator
            wide["denominator_species"] = denominator
            wide["pseudocount_fpm"] = pseudocount_fpm
            wide["fold_threshold"] = fold_threshold
            wide["log2_fpm_ratio"] = np.log2(
                (wide["numerator_mean_fpm"] + pseudocount_fpm)
                / (wide["denominator_mean_fpm"] + pseudocount_fpm)
            )
            log2_threshold = math.log2(fold_threshold)
            wide["contrast_class"] = np.select(
                [
                    wide["log2_fpm_ratio"].ge(log2_threshold),
                    wide["log2_fpm_ratio"].le(-log2_threshold),
                ],
                ["numerator_at_least_threshold", "denominator_at_least_threshold"],
                default="within_threshold",
            )
            wide["higher_species"] = np.select(
                [
                    wide["log2_fpm_ratio"].gt(0),
                    wide["log2_fpm_ratio"].lt(0),
                ],
                [numerator, denominator],
                default="equal",
            )
            result.append(wide)
    if not result:
        return pd.DataFrame(columns=[
            *index_columns, "numerator_mean_fpm", "denominator_mean_fpm",
            "numerator_species", "denominator_species", "pseudocount_fpm",
            "fold_threshold", "log2_fpm_ratio", "contrast_class", "higher_species",
        ])
    return pd.concat(result, ignore_index=True).sort_values(
        ["ASV_ID", "numerator_species", "denominator_species", "modality", "date", "depth_m"]
    )


def summarize_species_contrasts(contrasts: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    group_columns = ["ASV_ID", "numerator_species", "denominator_species", "modality"]
    for keys, group in contrasts.groupby(group_columns, sort=False):
        ratio = group["log2_fpm_ratio"].to_numpy(float)
        numerator = group["numerator_mean_fpm"].to_numpy(float)
        denominator = group["denominator_mean_fpm"].to_numpy(float)
        rho = pd.Series(numerator).corr(pd.Series(denominator), method="spearman")
        threshold = math.log2(float(group["fold_threshold"].iloc[0]))
        rows.append({
            "ASV_ID": keys[0],
            "numerator_species": keys[1],
            "denominator_species": keys[2],
            "modality": keys[3],
            "matched_samples_n": len(group),
            "spearman_rho": rho,
            "numerator_median_fpm": float(np.median(numerator)),
            "denominator_median_fpm": float(np.median(denominator)),
            "median_log2_fpm_ratio": float(np.median(ratio)),
            "median_fpm_ratio": float(2 ** np.median(ratio)),
            "log2_ratio_q25": float(np.quantile(ratio, 0.25)),
            "log2_ratio_q75": float(np.quantile(ratio, 0.75)),
            "log2_ratio_minimum": float(np.min(ratio)),
            "log2_ratio_maximum": float(np.max(ratio)),
            "numerator_higher_n": int(np.sum(ratio > 0)),
            "denominator_higher_n": int(np.sum(ratio < 0)),
            "equal_n": int(np.sum(ratio == 0)),
            "numerator_at_least_fold_threshold_n": int(np.sum(ratio >= threshold)),
            "denominator_at_least_fold_threshold_n": int(np.sum(ratio <= -threshold)),
            "within_fold_threshold_n": int(np.sum(np.abs(ratio) < threshold)),
            "fold_threshold": float(group["fold_threshold"].iloc[0]),
            "pseudocount_fpm": float(group["pseudocount_fpm"].iloc[0]),
        })
    return pd.DataFrame(rows)


def nice_ceiling(value: float) -> float:
    if not np.isfinite(value) or value <= 0:
        return 1.0
    exponent = math.floor(math.log10(value))
    scaled = value / (10.0 ** exponent)
    multiplier = 1.0 if scaled <= 1 else 2.0 if scaled <= 2 else 5.0 if scaled <= 5 else 10.0
    return multiplier * (10.0 ** exponent)


def power_of_two_ceiling(value: float) -> float:
    if not np.isfinite(value) or value <= 1:
        return 1.0
    return float(2 ** math.ceil(math.log2(value)))


def reference_values(values: pd.Series, maximum_references: int = 4) -> list[float]:
    positive = pd.to_numeric(values, errors="coerce").dropna()
    positive = positive.loc[positive > 0]
    if positive.empty:
        return [0.0]
    maximum = float(positive.max())
    largest_exponent = math.floor(math.log10(maximum))
    smallest_exponent = largest_exponent - maximum_references + 1
    references = [10.0 ** exponent for exponent in range(smallest_exponent, largest_exponent + 1)]
    while len(references) > maximum_references:
        references.pop(0)
    return references


def format_fpm(value: float) -> str:
    if value == 0:
        return "0"
    if value >= 1000:
        return f"{value:,.0f}"
    if value >= 10:
        return f"{value:.1f}"
    return f"{value:.2f}"


def recruitment_unit_label(normalization: object) -> str:
    return {
        "provided_fpkm": "FPKM",
        "provided_tpm": "TPM",
        "input_fragment_fpm": "FPM",
    }.get(str(normalization), str(normalization).replace("_", " "))


def render_asv_figure(
    asv: str,
    species: list[str],
    points: pd.DataFrame,
    surface: dict[str, object],
    sample_anchors: pd.DataFrame,
    palette: dict[str, str],
    state_name: str,
    output_base: Path,
    formats: list[str],
    maximum_depth: float,
    minimum_area: float,
    maximum_area: float,
    sample_point_area: float,
    renewal_onsets: pd.DatetimeIndex,
) -> list[dict[str, object]]:
    modalities = ["metagenome", "metatranscriptome"]
    nrows = len(species)
    fig, axes = plt.subplots(
        nrows, 2, figsize=(18.5, max(6.4, 5.1 * nrows)),
        sharex=True, sharey=True, squeeze=False,
    )
    fig.suptitle(f"{asv}: species-level MAG recruitment across {state_name}", y=0.995)
    date_min = pd.DatetimeIndex(surface["dates"]).min()
    date_max = pd.DatetimeIndex(surface["dates"]).max()
    maxima = {
        modality: float(points.loc[points["modality"].eq(modality), "mean_normalized_recruitment_fpm"].max())
        if points["modality"].eq(modality).any() else 0.0
        for modality in modalities
    }
    panel_rows: list[dict[str, object]] = []
    panel_index = 0
    for row_index, taxon in enumerate(species):
        for col_index, modality in enumerate(modalities):
            ax = axes[row_index, col_index]
            draw_surface(ax, surface, palette)
            draw_renewal_onsets(ax, renewal_onsets, date_min, date_max)
            draw_sample_anchors(ax, sample_anchors, sample_point_area)
            panel = points.loc[
                points["species"].eq(taxon) & points["modality"].eq(modality)
            ].copy()
            panel_letter = chr(ord("A") + panel_index)
            panel_index += 1
            if panel.empty:
                ax.text(
                    0.5, 0.5, "No matched recruitment data", transform=ax.transAxes,
                    ha="center", va="center", fontsize=12,
                    bbox={"facecolor": "white", "edgecolor": "#555555", "alpha": 0.9},
                    zorder=6,
                )
            else:
                areas = size_scale(
                    panel["mean_normalized_recruitment_fpm"], maxima[modality],
                    minimum_area, maximum_area,
                )
                ax.scatter(
                    panel["date"], panel["depth_m"], s=areas,
                    facecolor="white", edgecolor="black", linewidth=0.85,
                    alpha=0.95, zorder=5,
                )
            ax.set_title(f"{panel_letter}. {taxon} — {modality}")
            ax.set_ylim(maximum_depth, 0)
            ax.set_xlim(date_min, date_max)
            ax.set_yticks(np.arange(0, maximum_depth + 1, 25))
            ax.xaxis.set_major_locator(matplotlib.dates.YearLocator())
            ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%Y"))
            ax.tick_params(axis="x", rotation=0)
            if col_index == 0:
                ax.set_ylabel("Depth (m)")
            if row_index == nrows - 1:
                ax.set_xlabel("Calendar time")
            panel_rows.append({
                "ASV_ID": asv,
                "species": taxon,
                "modality": modality,
                "state_surface": state_name,
                "asv_sample_anchors_n": int(len(sample_anchors)),
                "panel": panel_letter,
                "samples_with_recruitment_n": int(len(panel)),
                "nonzero_samples_n": int((panel["mean_normalized_recruitment_fpm"] > 0).sum()) if not panel.empty else 0,
                "linked_genomes_min": int(panel["linked_genome_count"].min()) if not panel.empty else 0,
                "linked_genomes_max": int(panel["linked_genome_count"].max()) if not panel.empty else 0,
                "mean_fpm": float(panel["mean_normalized_recruitment_fpm"].mean()) if not panel.empty else np.nan,
                "median_fpm": float(panel["mean_normalized_recruitment_fpm"].median()) if not panel.empty else np.nan,
                "maximum_fpm": float(panel["mean_normalized_recruitment_fpm"].max()) if not panel.empty else np.nan,
                "normalization": str(panel["normalization"].iloc[0]) if not panel.empty else "",
                "abundance_units": recruitment_unit_label(panel["normalization"].iloc[0]) if not panel.empty else "",
                "mean_normalized_recruitment": float(panel["mean_normalized_recruitment_fpm"].mean()) if not panel.empty else np.nan,
                "median_normalized_recruitment": float(panel["mean_normalized_recruitment_fpm"].median()) if not panel.empty else np.nan,
                "maximum_normalized_recruitment": float(panel["mean_normalized_recruitment_fpm"].max()) if not panel.empty else np.nan,
                "output_stem": output_base.name,
            })
    state_handles = [
        Patch(facecolor=palette[label], edgecolor="#333333", linewidth=0.4, label=display_state(label))
        for label in surface["labels"]
    ]
    state_handles.append(renewal_legend_handle())
    fig.legend(
        handles=state_handles, loc="lower center", bbox_to_anchor=(0.5, 0.015),
        ncol=min(5, max(1, len(state_handles))), frameon=False, title=state_name,
    )
    for col_index, modality in enumerate(modalities):
        modality_values = points.loc[
            points["modality"].eq(modality), "mean_normalized_recruitment_fpm"
        ]
        refs = reference_values(modality_values)
        handles = [
            plt.scatter(
                [], [], s=sample_point_area, facecolor="black", edgecolor="black",
                linewidth=0, label="ASV sample",
            ),
            *[
            plt.scatter(
                [], [], s=size_scale(pd.Series([value]), maxima[modality], minimum_area, maximum_area)[0],
                facecolor="white", edgecolor="black", linewidth=0.85,
                label=format_fpm(value),
            )
            for value in refs
            ],
        ]
        modality_modes = points.loc[
            points["modality"].eq(modality), "normalization"
        ].dropna().astype(str).unique()
        unit_label = recruitment_unit_label(modality_modes[0]) if len(modality_modes) == 1 else "normalized abundance"
        axes[0, col_index].legend(
            handles=handles, title=f"Mean {modality} {unit_label}", loc="upper left",
            bbox_to_anchor=(1.005 if col_index == 1 else 0.0, 1.0),
            frameon=True, facecolor="white", edgecolor="#666666", fontsize=9,
        )
    fig.subplots_adjust(bottom=0.12, top=0.94, hspace=0.20, wspace=0.08)
    for fmt in formats:
        fig.savefig(output_base.with_suffix(f".{fmt}"), dpi=300 if fmt == "png" else None, bbox_inches="tight")
    plt.close(fig)
    return panel_rows


def render_contrast_figure(
    contrast: pd.DataFrame,
    surface: dict[str, object],
    sample_anchors: pd.DataFrame,
    palette: dict[str, str],
    state_name: str,
    output_base: Path,
    formats: list[str],
    maximum_depth: float,
    minimum_area: float,
    maximum_area: float,
    sample_point_area: float,
    renewal_onsets: pd.DatetimeIndex,
) -> None:
    if contrast.empty:
        return
    asv = str(contrast["ASV_ID"].iloc[0])
    numerator = str(contrast["numerator_species"].iloc[0])
    denominator = str(contrast["denominator_species"].iloc[0])
    threshold = float(contrast["fold_threshold"].iloc[0])
    modalities = ["metagenome", "metatranscriptome"]
    fig, axes = plt.subplots(1, 2, figsize=(18.5, 6.8), sharex=True, sharey=True)
    fig.suptitle(
        f"{asv}: paired species recruitment contrast across {state_name}", y=0.995
    )
    date_min = pd.DatetimeIndex(surface["dates"]).min()
    date_max = pd.DatetimeIndex(surface["dates"]).max()
    maximum_absolute = float(contrast["log2_fpm_ratio"].abs().max())
    marker_styles = {
        "numerator_at_least_threshold": ("^", "black", f"{numerator} ≥ {threshold:g}-fold higher"),
        "denominator_at_least_threshold": ("v", "white", f"{denominator} ≥ {threshold:g}-fold higher"),
        "within_threshold": ("o", "#A6A6A6", f"Within {threshold:g}-fold"),
    }
    for panel_index, (ax, modality) in enumerate(zip(axes, modalities)):
        draw_surface(ax, surface, palette)
        draw_renewal_onsets(ax, renewal_onsets, date_min, date_max)
        draw_sample_anchors(ax, sample_anchors, sample_point_area)
        panel = contrast.loc[contrast["modality"].eq(modality)].copy()
        for contrast_class, (marker, facecolor, _) in marker_styles.items():
            subset = panel.loc[panel["contrast_class"].eq(contrast_class)].copy()
            if subset.empty:
                continue
            areas = contrast_size_scale(
                subset["log2_fpm_ratio"], maximum_absolute,
                minimum_area, maximum_area,
            )
            # A white halo separates the contrast symbols from every curtain color.
            ax.scatter(
                subset["date"], subset["depth_m"], s=areas + 30,
                marker=marker, facecolor="white", edgecolor="white",
                linewidth=1.8, alpha=0.98, zorder=5,
            )
            ax.scatter(
                subset["date"], subset["depth_m"], s=areas,
                marker=marker, facecolor=facecolor, edgecolor="black",
                linewidth=0.95, alpha=0.98, zorder=6,
            )
        ax.set_title(f"{chr(ord('A') + panel_index)}. {modality}")
        ax.set_ylim(maximum_depth, 0)
        ax.set_xlim(date_min, date_max)
        ax.set_yticks(np.arange(0, maximum_depth + 1, 25))
        ax.xaxis.set_major_locator(matplotlib.dates.YearLocator())
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%Y"))
        ax.set_xlabel("Calendar time")
    axes[0].set_ylabel("Depth (m)")
    axes[1].set_ylabel("")
    direction_handles = [
        plt.scatter(
            [], [], s=sample_point_area, marker="o", facecolor="black",
            edgecolor="black", linewidth=0, label="ASV sample",
        ),
        *[
        plt.scatter(
            [], [], s=95, marker=marker, facecolor=facecolor,
            edgecolor="black", linewidth=0.95, label=label,
        )
        for marker, facecolor, label in marker_styles.values()
        ],
    ]
    fold_ceiling = power_of_two_ceiling(2.0 ** maximum_absolute)
    size_refs = [1.0]
    while size_refs[-1] < fold_ceiling:
        size_refs.append(size_refs[-1] * 2.0)
    if len(size_refs) > 5:
        size_refs = [1.0, 2.0, *size_refs[-3:]]
    size_handles = [
        plt.scatter(
            [], [],
            s=contrast_size_scale(
                pd.Series([math.log2(value)]), maximum_absolute, minimum_area, maximum_area
            )[0],
            marker="o", facecolor="#A6A6A6", edgecolor="black", linewidth=0.95,
            label=f"{value:g}x",
        )
        for value in size_refs if np.isfinite(value)
    ]
    axes[1].legend(
        handles=[*direction_handles, *size_handles],
        title="Direction and absolute fold difference", loc="upper left",
        bbox_to_anchor=(1.01, 1.0), frameon=True,
        facecolor="white", edgecolor="#666666",
    )
    state_handles = [
        Patch(facecolor=palette[label], edgecolor="#333333", linewidth=0.4, label=display_state(label))
        for label in surface["labels"]
    ]
    state_handles.append(renewal_legend_handle())
    fig.legend(
        handles=state_handles, loc="lower center", bbox_to_anchor=(0.5, 0.01),
        ncol=min(5, max(1, len(state_handles))), frameon=False, title=state_name,
    )
    fig.subplots_adjust(bottom=0.17, top=0.91, wspace=0.08, right=0.84)
    for fmt in formats:
        fig.savefig(
            output_base.with_suffix(f".{fmt}"),
            dpi=300 if fmt == "png" else None,
            bbox_inches="tight",
        )
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--asv-mag-network-dir", type=Path, required=True)
    parser.add_argument("--microbial-state-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--prefix", default="asv_mag_network")
    parser.add_argument("--hybrid-palette", required=True)
    parser.add_argument("--hybrid-order", default="")
    parser.add_argument("--mc-palette", required=True)
    parser.add_argument("--mc-order", default="MC1,MC2,MC3,MC4")
    parser.add_argument("--sample-col", default="sampleID")
    parser.add_argument("--date-col", default="Date")
    parser.add_argument("--depth-col", default="Depth")
    parser.add_argument("--hybrid-col", default="o2_subcompartment_final")
    parser.add_argument("--maximum-depth", type=float, default=210.0)
    parser.add_argument("--minimum-point-area", type=float, default=1.0)
    parser.add_argument("--maximum-point-area", type=float, default=320.0)
    parser.add_argument("--asv-sample-point-area", type=float, default=6.0)
    parser.add_argument("--contrast-pseudocount-fpm", type=float, default=1.0)
    parser.add_argument("--contrast-fold-threshold", type=float, default=2.0)
    parser.add_argument("--renewal-events", type=Path)
    parser.add_argument("--renewal-date-col", default="start_date")
    parser.add_argument("--formats", default="pdf,png,svg")
    args = parser.parse_args()
    if args.maximum_depth <= 0 or args.minimum_point_area <= 0:
        raise ValueError("Depth and point areas must be positive")
    if args.maximum_point_area <= args.minimum_point_area:
        raise ValueError("Maximum point area must exceed minimum point area")
    if args.asv_sample_point_area <= 0:
        raise ValueError("ASV sample point area must be positive")
    if args.contrast_pseudocount_fpm <= 0 or args.contrast_fold_threshold <= 1:
        raise ValueError("Contrast pseudocount must be positive and fold threshold exceed one")

    matplotlib.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "font.size": 11,
        "axes.titlesize": 12,
        "axes.labelsize": 12,
        "legend.fontsize": 9,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    tables = args.outdir / "tables"
    plots = args.outdir / "plots"
    audit = args.outdir / "audit"
    for directory in (tables, plots, audit):
        directory.mkdir(parents=True, exist_ok=True)

    mappings = args.asv_mag_network_dir / "mapping" / f"{args.prefix}_analysis_eligible_mappings.tsv"
    eligible = load_eligible_pairs(mappings)
    profiles = load_recruitment_profiles(args.asv_mag_network_dir / "abundance", args.prefix)
    points = aggregate_species_points(profiles, eligible)
    points.to_csv(tables / "asv_mag_species_recruitment_curtain_points.tsv", sep="\t", index=False)
    contrasts = calculate_species_contrasts(
        points, args.contrast_pseudocount_fpm, args.contrast_fold_threshold
    )
    contrasts.to_csv(
        tables / "asv_mag_species_recruitment_pairwise_contrasts.tsv", sep="\t", index=False
    )
    contrast_summary = summarize_species_contrasts(contrasts)
    contrast_summary.to_csv(
        tables / "asv_mag_species_recruitment_pairwise_contrast_summary.tsv",
        sep="\t", index=False,
    )
    renewal_onsets = load_renewal_onsets(args.renewal_events, args.renewal_date_col)
    pd.DataFrame({"renewal_onset": renewal_onsets}).to_csv(
        tables / "asv_mag_curtain_renewal_onsets.tsv", sep="\t", index=False
    )

    crosswalk_path = (
        args.microbial_state_dir / "tables" / "microbial_state_sample_crosswalk.tsv"
    )
    surfaces = {
        "hybrid O2-GMM compartments": (
            load_surface(
                args.microbial_state_dir / "tables" / "asv_only_hybrid_compartment_curtain_grid.tsv",
                parse_list(args.hybrid_order),
            ),
            parse_palette(args.hybrid_palette),
            "hybrid",
        ),
        "microbial compartments": (
            load_surface(
                args.microbial_state_dir / "tables" / "microbial_compartment_curtain_grid.tsv",
                parse_list(args.mc_order),
            ),
            parse_palette(args.mc_palette),
            "microbial_compartment",
        ),
    }
    for state_name, (surface, palette, _) in surfaces.items():
        missing_colors = set(surface["labels"]).difference(palette)
        if missing_colors:
            raise ValueError(f"{state_name} palette lacks states: {sorted(missing_colors)}")

    anchor_specs = {
        "hybrid O2-GMM compartments": args.hybrid_col,
        "microbial compartments": "microbial_compartment",
    }
    anchors_by_surface = {
        state_name: load_sample_anchors(
            crosswalk_path, args.sample_col, args.date_col, args.depth_col,
            anchor_specs[state_name], list(surface["labels"]), state_name,
        )
        for state_name, (surface, _, _) in surfaces.items()
    }
    pd.concat(anchors_by_surface.values(), ignore_index=True).to_csv(
        tables / "asv_mag_curtain_asv_sample_anchors.tsv", sep="\t", index=False
    )

    all_asvs = sorted(set(eligible["ASV_ID"]), key=asv_sort_key)
    panel_rows: list[dict[str, object]] = []
    manifest_rows: list[dict[str, object]] = []
    contrast_manifest_rows: list[dict[str, object]] = []
    for asv in all_asvs:
        asv_points = points.loc[points["ASV_ID"].eq(asv)].copy()
        species = sorted(set(eligible.loc[eligible["ASV_ID"].eq(asv), "species"]))
        if not species:
            continue
        for state_name, (surface, palette, suffix) in surfaces.items():
            sample_anchors = anchors_by_surface[state_name]
            output_base = plots / f"{safe_name(asv)}_{suffix}_mag_recruitment_curtain"
            panel_rows.extend(render_asv_figure(
                asv, species, asv_points, surface, sample_anchors, palette, state_name,
                output_base, parse_list(args.formats), args.maximum_depth,
                args.minimum_point_area, args.maximum_point_area,
                args.asv_sample_point_area,
                renewal_onsets,
            ))
            manifest_rows.append({
                "ASV_ID": asv,
                "state_surface": state_name,
                "species_n": len(species),
                "panel_n": 2 * len(species),
                "species": " | ".join(species),
                "asv_sample_anchors_n": int(len(sample_anchors)),
                "output_stem": output_base.name,
                "formats": ",".join(parse_list(args.formats)),
            })
        asv_contrasts = contrasts.loc[contrasts["ASV_ID"].eq(asv)]
        pair_columns = ["numerator_species", "denominator_species"]
        pairs = asv_contrasts[pair_columns].drop_duplicates().itertuples(index=False, name=None)
        pair_list = list(pairs)
        for pair_index, (numerator, denominator) in enumerate(pair_list, start=1):
            pair = asv_contrasts.loc[
                asv_contrasts["numerator_species"].eq(numerator)
                & asv_contrasts["denominator_species"].eq(denominator)
            ].copy()
            pair_suffix = "" if len(pair_list) == 1 else f"_pair{pair_index}"
            for state_name, (surface, palette, suffix) in surfaces.items():
                sample_anchors = anchors_by_surface[state_name]
                output_base = plots / (
                    f"{safe_name(asv)}_{suffix}_mag_species_contrast_curtain{pair_suffix}"
                )
                render_contrast_figure(
                    pair, surface, sample_anchors, palette, state_name, output_base,
                    parse_list(args.formats), args.maximum_depth,
                    args.minimum_point_area, args.maximum_point_area,
                    args.asv_sample_point_area,
                    renewal_onsets,
                )
                contrast_manifest_rows.append({
                    "ASV_ID": asv,
                    "numerator_species": numerator,
                    "denominator_species": denominator,
                    "state_surface": state_name,
                    "asv_sample_anchors_n": int(len(sample_anchors)),
                    "output_stem": output_base.name,
                    "formats": ",".join(parse_list(args.formats)),
                })

    panel_summary = pd.DataFrame(panel_rows)
    panel_summary.to_csv(tables / "asv_mag_species_recruitment_curtain_panel_summary.tsv", sep="\t", index=False)
    pd.DataFrame(manifest_rows).to_csv(
        tables / "asv_mag_species_recruitment_curtain_manifest.tsv", sep="\t", index=False
    )
    pd.DataFrame(contrast_manifest_rows).to_csv(
        tables / "asv_mag_species_recruitment_pairwise_contrast_manifest.tsv",
        sep="\t", index=False,
    )
    eligible.groupby(["ASV_ID", "species"], as_index=False).agg(
        eligible_genome_count=("mag_match_id", "nunique")
    ).to_csv(audit / "asv_mag_species_curtain_eligibility.tsv", sep="\t", index=False)
    configuration = {
        "species_aggregation": "arithmetic_mean_across_analysis_eligible_linked_genomes",
        "abundance_unit": "recorded_per_modality_from_normalization_column",
        "point_size_transform": "linear_marker_area_on_mean_normalized_abundance_to_1-2-5_ceiling_within_ASV_and_modality",
        "modality_scales": "separate_metagenome_and_metatranscriptome",
        "contrast_definition": "log2((numerator_mean_normalized_abundance+pseudocount)/(denominator_mean_normalized_abundance+pseudocount)); contrasts are within modality and units",
        "contrast_pseudocount_fpm": args.contrast_pseudocount_fpm,
        "contrast_fold_threshold": args.contrast_fold_threshold,
        "asv_sample_point_area": args.asv_sample_point_area,
        "asv_sample_anchor_source": str(crosswalk_path),
        "renewal_event_source": str(args.renewal_events) if args.renewal_events else None,
        "renewal_date_column": args.renewal_date_col,
        "renewal_onsets_n": len(renewal_onsets),
        "contrast_direction": "species_names_sorted_lexicographically_define_numerator_then_denominator",
        "asvs_n": len(all_asvs),
        "asv_species_pairs_n": int(eligible[["ASV_ID", "species"]].drop_duplicates().shape[0]),
        "formats": parse_list(args.formats),
        "python": platform.python_version(),
        "pandas": pd.__version__,
        "numpy": np.__version__,
        "matplotlib": matplotlib.__version__,
    }
    (audit / "asv_mag_curtain_parameters.json").write_text(json.dumps(configuration, indent=2) + "\n")


if __name__ == "__main__":
    main()
