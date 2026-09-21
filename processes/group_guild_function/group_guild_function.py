#!/usr/bin/env python3
"""Integrate environmental groups, ASV network modules, MAGs, and KEGG modules."""

from __future__ import annotations

import argparse
import colorsys
import json
import re
from pathlib import Path
import sys as _sys
_sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shared_plot_style import install_publication_style
install_publication_style()

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
from matplotlib import font_manager, colors as mcolors
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact, mannwhitneyu


def read_table(path: Path | None) -> pd.DataFrame:
    if path is None or not path.is_file() or path.stat().st_size == 0:
        return pd.DataFrame()
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    return pd.read_csv(path, sep=sep, low_memory=False)


def csv_list(value: str) -> list[str]:
    return [x.strip() for x in str(value).split(",") if x.strip()]


def bh(values: pd.Series) -> pd.Series:
    out = pd.Series(np.nan, index=values.index, dtype=float)
    valid = pd.to_numeric(values, errors="coerce").dropna().sort_values()
    if valid.empty:
        return out
    n = len(valid)
    adjusted = (valid * n / np.arange(1, n + 1)).iloc[::-1].cummin().iloc[::-1].clip(upper=1)
    out.loc[adjusted.index] = adjusted
    return out


def select_well_covered_cycle_modules(
    functions: pd.DataFrame,
    minimum_fraction: float,
) -> tuple[list[str], pd.DataFrame]:
    """Select N/S KEGG modules represented at the configured MAG coverage."""
    if functions.empty:
        return [], pd.DataFrame()
    work = functions.copy()
    work["module_name"] = work["module_name"].fillna("").astype(str)
    work["fraction_covered"] = pd.to_numeric(
        work["fraction_covered"], errors="coerce"
    ).fillna(0.0)
    nitrogen = (
        r"nitrogen fixation|nitrification|denitrification|"
        r"(?:assimilatory|dissimilatory) nitrate reduction|"
        r"nitrate assimilation|anammox"
    )
    sulfur = (
        r"(?:assimilatory|dissimilatory) sulfate reduction|"
        r"sulfur oxidation|sulfide oxidation|sulfur reduction|"
        r"dimethylsulfoniopropionate .*degradation|DMSP degradation|"
        r"sulfoquinovose degradation"
    )
    work["cycle"] = np.select(
        [
            work["module_name"].str.contains(nitrogen, case=False, regex=True),
            work["module_name"].str.contains(sulfur, case=False, regex=True),
        ],
        ["N", "S"],
        default="",
    )
    per_mag = (
        work.loc[work["cycle"].ne("")]
        .groupby(
            ["module_id", "module_name", "cycle", "mag_match_id"],
            as_index=False,
        )
        .agg(fraction_covered=("fraction_covered", "max"))
    )
    audit = (
        per_mag.groupby(["module_id", "module_name", "cycle"], as_index=False)
        .agg(
            linked_mags=("mag_match_id", "nunique"),
            median_fraction_covered=("fraction_covered", "median"),
            maximum_fraction_covered=("fraction_covered", "max"),
            mags_meeting_threshold=(
                "fraction_covered",
                lambda values: int((values >= minimum_fraction).sum()),
            ),
        )
    )
    audit["minimum_fraction_threshold"] = float(minimum_fraction)
    audit["selected"] = audit["mags_meeting_threshold"].gt(0)
    selected = sorted(
        audit.loc[audit["selected"], "module_id"].dropna().astype(str).unique()
    )
    return selected, audit.sort_values(["cycle", "module_id"])


def norm_asv(value: object) -> str:
    text = str(value).strip()
    return re.sub(r"^n(?=\d+$)", "", text)


def configure_font(font: str, strict: bool) -> None:
    available = {item.name for item in font_manager.fontManager.ttflist}
    if strict and font not in available:
        raise SystemExit(
            f"Required publication font '{font}' is not installed. "
            "Install it and refresh the Matplotlib font cache before rerunning."
        )
    plt.rcParams.update({
        "font.family": font,
        "font.size": 22,
        "axes.labelsize": 24,
        "axes.titlesize": 24,
        "xtick.labelsize": 22,
        "ytick.labelsize": 22,
        "legend.fontsize": 22,
        "figure.titlesize": 24,
        "axes.linewidth": 1.2,
    })


def load_counts(path: Path, sample_ids: set[str]) -> pd.DataFrame:
    raw = pd.read_csv(path, sep="\t", index_col=0)
    raw.index = raw.index.map(norm_asv)
    raw.columns = raw.columns.astype(str)
    if len(set(raw.columns) & sample_ids) >= len(set(raw.index) & sample_ids):
        return raw.apply(pd.to_numeric, errors="coerce").fillna(0)
    out = raw.T
    out.index = out.index.astype(str)
    out.columns = out.columns.map(norm_asv)
    return out.T.apply(pd.to_numeric, errors="coerce").fillna(0)


def load_isa(indir: Path, groupings: list[str]) -> pd.DataFrame:
    rows: list[dict] = []
    suffix = "_indicator_species_summary.tsv"
    inputs: list[tuple[str, Path]] = [(grouping, indir / f"{grouping}{suffix}") for grouping in groupings]
    stratified_suffix = "_indicator_species_summary.tsv"
    for path in sorted(indir.glob(f"stratified_*{stratified_suffix}")):
        inputs.append(("", path))
    seen_paths: set[Path] = set()
    for configured_grouping, path in inputs:
        path = path.resolve()
        if path in seen_paths:
            continue
        seen_paths.add(path)
        table = read_table(path)
        if table.empty:
            continue
        grouping = configured_grouping
        if not grouping and "stratified_group_col" in table:
            values = table["stratified_group_col"].dropna().astype(str).unique()
            grouping = values[0] if len(values) == 1 else ""
        if not grouping:
            match = re.match(r"stratified_(.+?)_within_", path.name)
            grouping = match.group(1) if match else path.stem
        asv_col = "ASV" if "ASV" in table else table.columns[0]
        q_col = "q.value" if "q.value" in table else "p.value" if "p.value" in table else None
        stat_col = next((c for c in ("stat", "index", "A", "B") if c in table), None)
        sign_cols = [c for c in table if str(c).startswith("s.")]
        for _, row in table.iterrows():
            levels = [str(c)[2:] for c in sign_cols if pd.to_numeric(row.get(c), errors="coerce") == 1]
            if not levels and "index" in table:
                levels = [str(row.get("index"))]
            scope = "single_group" if len(levels) == 1 else "multigroup" if len(levels) > 1 else "unassigned"
            rows.append({
                "grouping": grouping,
                "group_level": "+".join(levels),
                "association_group_count": len(levels),
                "association_scope": scope,
                "stratified_within_col": str(row.get("stratified_within_col", "")),
                "stratified_within_value": str(row.get("stratified_within_value", "")),
                "ASV_ID": norm_asv(row.get(asv_col)),
                "isa_statistic": pd.to_numeric(row.get(stat_col), errors="coerce") if stat_col else np.nan,
                "isa_q_value": pd.to_numeric(row.get(q_col), errors="coerce") if q_col else np.nan,
            })
    out = pd.DataFrame(rows)
    keys = ["grouping", "group_level", "association_scope", "stratified_within_col", "stratified_within_value", "ASV_ID"]
    return out.drop_duplicates(keys) if not out.empty else out


def module_scores(counts: pd.DataFrame, modules: pd.DataFrame) -> pd.DataFrame:
    modules = modules.copy()
    asv_col = "Taxon" if "Taxon" in modules else "ASV_ID"
    modules["ASV_ID"] = modules[asv_col].map(norm_asv)
    totals = counts.sum(axis=0).replace(0, np.nan)
    relative = counts.div(totals, axis=1).fillna(0)
    rows = []
    for label, part in modules.groupby("module_label", sort=True):
        members = sorted(set(part["ASV_ID"]) & set(relative.index))
        if not members:
            continue
        score = relative.loc[members].sum(axis=0)
        rows.extend({"sampleID": str(sample), "ecological_module": str(label),
                     "module_relative_abundance": float(value)} for sample, value in score.items())
    return pd.DataFrame(rows)


OXYGEN_NAMES = ("oxic", "dysoxic", "suboxic", "anoxic")
O2_COLORS = {
    "oxic": "red",
    "dysoxic": "green",
    "suboxic": "lightblue",
    "anoxic": "purple",
}


def hybrid_sort_key(value: object) -> tuple[int, int]:
    match = re.search(r"(?:hybrid_c|hyb_C)([0-3])_(?:g|G)([0-9]+)", str(value))
    return (int(match.group(1)), int(match.group(2))) if match else (99, 99)


def hybrid_palette(compartments: list[str]) -> dict[str, tuple[float, float, float]]:
    parsed = [(*hybrid_sort_key(value), value) for value in compartments]
    gmm_values = sorted({gmm for oxygen, gmm, _ in parsed if oxygen < 99})
    if len(gmm_values) == 1:
        lightness = {gmm_values[0]: 0.52}
    else:
        lightness = dict(
            zip(gmm_values, np.linspace(0.30, 0.76, len(gmm_values)))
        )
    palette: dict[str, tuple[float, float, float]] = {}
    for oxygen, gmm, raw in parsed:
        if oxygen >= len(OXYGEN_NAMES):
            continue
        hue, _, saturation = colorsys.rgb_to_hls(
            *mcolors.to_rgb(O2_COLORS[OXYGEN_NAMES[oxygen]])
        )
        palette[raw] = colorsys.hls_to_rgb(
            hue, float(lightness[gmm]), max(0.58, saturation)
        )
    return palette


def lighten(color: object, fraction: float = 0.80) -> tuple[float, float, float]:
    rgb = np.asarray(mcolors.to_rgb(color), dtype=float)
    return tuple(rgb + (1.0 - rgb) * fraction)


def merge_basin_pca_background(
    scores: pd.DataFrame,
    hybrid_assignments: pd.DataFrame,
) -> pd.DataFrame:
    """Join the exact BASIN PCA scores and soft hybrid assignments."""
    if "cruise_year_month_depth" in scores and "cruise_year_month_depth" in hybrid_assignments:
        keys = ["cruise_year_month_depth"]
    else:
        keys = ["Cruise", "date", "Depth_anchored"]
    missing = (
        set(keys + ["PC1", "PC2"]) - set(scores)
    ) | (set(keys) - set(hybrid_assignments))
    if missing:
        raise ValueError(
            "BASIN PCA/hybrid inputs lack required columns: "
            + ", ".join(sorted(missing))
        )
    merged = hybrid_assignments.merge(
        scores[keys + ["PC1", "PC2"]],
        on=keys,
        how="inner",
        validate="one_to_one",
    )
    hybrid_columns = [
        column for column in merged
        if re.fullmatch(r"hyb_C[0-3]_G[0-9]+", str(column))
    ]
    if merged.empty or not hybrid_columns:
        raise ValueError("BASIN PCA and hybrid assignments produced no usable samples")
    merged["hard_hybrid"] = (
        merged[hybrid_columns].apply(pd.to_numeric, errors="coerce")
        .fillna(0.0).idxmax(axis=1)
    )
    return merged


def build_asv_pca_positions(
    counts: pd.DataFrame,
    modules: pd.DataFrame,
    metadata: pd.DataFrame,
    pca_scores: pd.DataFrame,
    sample_col: str,
) -> pd.DataFrame:
    """Calculate observed ASV PC barycenters weighted by relative abundance.

    Multiple libraries at the same cruise-depth coordinate are first collapsed
    by their median relative abundance so replicated sampling locations receive
    one contribution per ASV.
    """
    required_metadata = {sample_col, "Cruise", "Depth"}
    if not required_metadata.issubset(metadata):
        raise ValueError(
            "ASV metadata lacks PCA join columns: "
            + ", ".join(sorted(required_metadata - set(metadata)))
        )
    depth_column = "Depth_anchored" if "Depth_anchored" in pca_scores else "Depth"
    required_scores = {"Cruise", depth_column, "PC1", "PC2"}
    if not required_scores.issubset(pca_scores):
        raise ValueError(
            "BASIN PCA scores lack ASV join columns: "
            + ", ".join(sorted(required_scores - set(pca_scores)))
        )

    meta = metadata[[sample_col, "Cruise", "Depth"]].drop_duplicates(sample_col).copy()
    meta[sample_col] = meta[sample_col].astype(str)
    meta["_cruise_key"] = pd.to_numeric(meta["Cruise"], errors="coerce")
    meta["_depth_key"] = pd.to_numeric(meta["Depth"], errors="coerce").round(6)
    score_coords = pca_scores[["Cruise", depth_column, "PC1", "PC2"]].copy()
    score_coords["_cruise_key"] = pd.to_numeric(
        score_coords["Cruise"], errors="coerce"
    )
    score_coords["_depth_key"] = pd.to_numeric(
        score_coords[depth_column], errors="coerce"
    ).round(6)
    score_coords["PC1"] = pd.to_numeric(score_coords["PC1"], errors="coerce")
    score_coords["PC2"] = pd.to_numeric(score_coords["PC2"], errors="coerce")
    score_coords = (
        score_coords.dropna(subset=["_cruise_key", "_depth_key", "PC1", "PC2"])
        .drop_duplicates(["_cruise_key", "_depth_key"])
    )
    sample_coords = meta.merge(
        score_coords[["_cruise_key", "_depth_key", "PC1", "PC2"]],
        on=["_cruise_key", "_depth_key"],
        how="inner",
        validate="many_to_one",
    )
    sample_ids = [
        sample for sample in sample_coords[sample_col].astype(str)
        if sample in counts.columns
    ]
    if not sample_ids:
        raise ValueError("No ASV samples matched BASIN cruise-depth PCA coordinates")

    relative = counts[sample_ids].div(
        counts[sample_ids].sum(axis=0).replace(0, np.nan), axis=1
    ).fillna(0.0)
    module_table = modules.copy()
    asv_column = "Taxon" if "Taxon" in module_table else "ASV_ID"
    module_table["ASV_ID"] = module_table[asv_column].map(norm_asv)
    module_table = module_table[["ASV_ID", "module_label"]].drop_duplicates("ASV_ID")
    network_asvs = [
        asv for asv in module_table["ASV_ID"] if asv in relative.index
    ]
    sample_lookup = sample_coords.set_index(sample_col).loc[sample_ids]
    rows: list[dict[str, object]] = []
    for asv_id in network_asvs:
        observations = pd.DataFrame({
            "relative_abundance": pd.to_numeric(
                relative.loc[asv_id, sample_ids], errors="coerce"
            ).to_numpy(),
            "PC1": sample_lookup["PC1"].to_numpy(),
            "PC2": sample_lookup["PC2"].to_numpy(),
            "_cruise_key": sample_lookup["_cruise_key"].to_numpy(),
            "_depth_key": sample_lookup["_depth_key"].to_numpy(),
        })
        collapsed = (
            observations.groupby(
                ["_cruise_key", "_depth_key", "PC1", "PC2"], as_index=False
            )["relative_abundance"].median()
        )
        positive = collapsed["relative_abundance"].gt(0)
        collapsed = collapsed.loc[positive].copy()
        weight_sum = float(collapsed["relative_abundance"].sum())
        if collapsed.empty or weight_sum <= 0:
            continue
        rows.append({
            "ASV_ID": asv_id,
            "module_label": module_table.set_index("ASV_ID").at[
                asv_id, "module_label"
            ],
            "observed_PC1": float(np.average(
                collapsed["PC1"], weights=collapsed["relative_abundance"]
            )),
            "observed_PC2": float(np.average(
                collapsed["PC2"], weights=collapsed["relative_abundance"]
            )),
            "matched_environmental_positions_n": int(len(collapsed)),
            "relative_abundance_weight_sum": weight_sum,
            "position_interpretation": (
                "Observed relative-abundance-weighted ASV barycenter; "
                "not a predicted niche"
            ),
        })
    return pd.DataFrame(rows).sort_values(
        ["module_label", "ASV_ID"],
        key=lambda values: (
            values.map(lambda value: hybrid_sort_key(value))
            if values.name == "module_label"
            else values
        ),
    ).reset_index(drop=True)


def build_module_phylum_pca_positions(
    asv_positions: pd.DataFrame,
    taxonomy: pd.DataFrame,
) -> pd.DataFrame:
    """Pool observed ASV positions within each ecological-module/phylum pair."""
    required = {"Feature ID", "Taxon"}
    if not required.issubset(taxonomy):
        raise ValueError(
            "ASV taxonomy lacks required columns: "
            + ", ".join(sorted(required - set(taxonomy)))
        )
    tax = taxonomy[["Feature ID", "Taxon"]].copy()
    tax["ASV_ID"] = (
        tax["Feature ID"].astype(str).str.replace(r";.*$", "", regex=True)
        .map(norm_asv)
    )
    tax["phylum"] = tax["Taxon"].astype(str).str.extract(
        r"(?:^|;\s*)p__([^;]+)", expand=False
    )
    tax["phylum"] = (
        tax["phylum"].fillna("Unclassified").replace("", "Unclassified")
    )
    tax = tax[["ASV_ID", "phylum"]].drop_duplicates("ASV_ID")
    merged = asv_positions.merge(
        tax, on="ASV_ID", how="left", validate="many_to_one"
    )
    if merged["phylum"].isna().any():
        missing = merged.loc[merged["phylum"].isna(), "ASV_ID"].tolist()
        raise ValueError(
            "Taxonomy was unavailable for network ASVs: "
            + ", ".join(map(str, missing[:10]))
        )

    rows = []
    for (module, phylum), frame in merged.groupby(
        ["module_label", "phylum"], sort=False
    ):
        weights = pd.to_numeric(
            frame["relative_abundance_weight_sum"], errors="coerce"
        ).fillna(0.0)
        if weights.sum() <= 0:
            continue
        rows.append({
            "module_label": module,
            "phylum": phylum,
            "asvs_n": int(frame["ASV_ID"].nunique()),
            "ASV_IDs": ";".join(sorted(frame["ASV_ID"].astype(str).unique())),
            "observed_PC1": float(np.average(
                frame["observed_PC1"], weights=weights
            )),
            "observed_PC2": float(np.average(
                frame["observed_PC2"], weights=weights
            )),
            "relative_abundance_weight_sum": float(weights.sum()),
            "position_interpretation": (
                "Observed abundance-weighted barycenter of ASVs sharing an "
                "ecological module and phylum; not a predicted niche"
            ),
        })
    return pd.DataFrame(rows).sort_values(
        ["module_label", "phylum"],
        key=lambda values: (
            values.map(
                lambda value: int(re.search(r"\d+", str(value)).group())
                if re.search(r"\d+", str(value)) else 999
            )
            if values.name == "module_label" else values
        ),
    ).reset_index(drop=True)


def explained_axis_labels(explained: pd.DataFrame) -> tuple[str, str]:
    ratios = {}
    if {"PC", "explained_variance_ratio"}.issubset(explained):
        ratios = dict(zip(
            explained["PC"].astype(str),
            pd.to_numeric(explained["explained_variance_ratio"], errors="coerce"),
        ))
    return (
        f"PC1 ({100 * ratios.get('PC1', np.nan):.1f}%)",
        f"PC2 ({100 * ratios.get('PC2', np.nan):.1f}%)",
    )


def resolve_centroid_label_collisions(
    fig,
    ax,
    annotations: list,
    forbidden_xy: list[tuple[float, float]],
) -> None:
    """Move opaque centroid labels away from one another and centroid markers."""
    if not annotations:
        return
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    forbidden_display = [
        ax.transData.transform(point) for point in forbidden_xy
    ]
    shifts = [(0, 0)]
    for radius in (16, 28, 42, 60, 82, 108, 136, 170, 210):
        shifts.extend([
            (
                radius * np.cos(np.deg2rad(angle)),
                radius * np.sin(np.deg2rad(angle)),
            )
            for angle in range(0, 360, 30)
        ])
    accepted = []
    for annotation in annotations:
        original_data = annotation.get_position()
        original_display = ax.transData.transform(original_data)
        initial_box = annotation.get_window_extent(renderer=renderer)
        chosen = (0, 0)
        for dx, dy in shifts:
            candidate = mtransforms.Bbox.from_extents(
                initial_box.x0 + dx, initial_box.y0 + dy,
                initial_box.x1 + dx, initial_box.y1 + dy,
            )
            padded = mtransforms.Bbox.from_extents(
                candidate.x0 - 6, candidate.y0 - 6,
                candidate.x1 + 6, candidate.y1 + 6,
            )
            covers_point = any(
                padded.x0 <= px <= padded.x1
                and padded.y0 <= py <= padded.y1
                for px, py in forbidden_display
            )
            overlaps_label = any(padded.overlaps(previous) for previous in accepted)
            if not covers_point and not overlaps_label:
                chosen = (dx, dy)
                accepted.append(padded)
                break
        annotation.set_position(
            ax.transData.inverted().transform(
                (original_display[0] + chosen[0], original_display[1] + chosen[1])
            )
        )
    fig.canvas.draw()


def select_pca_label_asvs(counts: pd.DataFrame, positions: pd.DataFrame,
                          threshold: float = 0.01) -> pd.DataFrame:
    """Select before joining positions: denominator includes every input ASV."""
    if not 0 <= threshold < 1:
        raise ValueError("PCA label abundance threshold must be in [0, 1)")
    relative = counts.div(counts.sum(axis=0).replace(0, np.nan), axis=1)
    mean = relative.mean(axis=1).rename("mean_relative_abundance")
    selected = mean[mean > threshold].rename_axis("ASV_ID").reset_index()
    selected = selected.sort_values(
        ["mean_relative_abundance", "ASV_ID"], ascending=[False, True]
    )
    selected["selection_threshold"] = threshold
    selected["selection_samples_n"] = int(relative.notna().any(axis=0).sum())
    selected = selected.merge(positions, on="ASV_ID", how="left", validate="one_to_one")
    selected["plotted"] = np.isfinite(selected[["observed_PC1", "observed_PC2"]]).all(axis=1)
    return selected


def attach_label_taxonomy(selected: pd.DataFrame, taxonomy: pd.DataFrame) -> pd.DataFrame:
    tax = taxonomy[["Feature ID", "Taxon"]].copy()
    tax["ASV_ID"] = tax["Feature ID"].astype(str).str.split(";").str[0].map(norm_asv)
    for rank, prefix in [("phylum", "p"), ("family", "f"), ("genus", "g")]:
        tax[rank] = tax.Taxon.str.extract(r"(?:^|;\s*)" + prefix + r"__([^;]+)", expand=False)
        tax[rank] = tax[rank].fillna("Unclassified").replace("", "Unclassified")
    return selected.merge(tax[["ASV_ID", "phylum", "family", "genus"]].drop_duplicates("ASV_ID"),
                          on="ASV_ID", how="left", validate="one_to_one")


def attach_label_module_peaks(selected: pd.DataFrame, module_stats: pd.DataFrame) -> pd.DataFrame:
    """Attach aggregate module maxima from the matched hybrid comparison."""
    peaks = module_stats.loc[
        module_stats.grouping.eq("o2_subcompartment_final"),
        ["ecological_module", "highest_group_level", "highest_group_mean"],
    ].rename(columns={"ecological_module": "module_label",
                      "highest_group_level": "module_peak_hybrid_compartment",
                      "highest_group_mean": "module_peak_mean_abundance"})
    return selected.merge(peaks, on="module_label", how="left", validate="many_to_one")


def leader_segments_intersect(a, b, c, d, tolerance=1e-7):
    """Include crossings, collinear overlap, and touching nonshared endpoints."""
    def cross(x, y, z):
        delta, other = y - x, z - x
        return delta[0] * other[1] - delta[1] * other[0]
    def on_segment(x, y, z):
        return (abs(cross(x, y, z)) <= tolerance
                and np.all(z >= np.minimum(x, y) - tolerance)
                and np.all(z <= np.maximum(x, y) + tolerance))
    v1, v2, v3, v4 = cross(a, b, c), cross(a, b, d), cross(c, d, a), cross(c, d, b)
    return ((v1 * v2 < 0 and v3 * v4 < 0)
            or on_segment(a, b, c) or on_segment(a, b, d)
            or on_segment(c, d, a) or on_segment(c, d, b))


def local_label_layout(anchors, half_sizes, point_centers, point_radii, bounds):
    """Search nearby point-free text boxes, then optimize a collision-free layout.

    Straight leaders stop at the text box edge. Candidate layouts prohibit
    leader crossings and any leader entering another label's box.
    """
    def endpoint(anchor, center, half):
        delta = center - anchor
        scale = min(half[0] / max(abs(delta[0]), 1e-9), half[1] / max(abs(delta[1]), 1e-9))
        return center - delta * scale

    def hits_box(a, b, center, half):
        low, high = center - half, center + half
        t0, t1 = 0., 1.
        for axis in range(2):
            delta = b[axis] - a[axis]
            if abs(delta) < 1e-10:
                if a[axis] < low[axis] or a[axis] > high[axis]:
                    return False
            else:
                entry, leave = sorted(((low[axis] - a[axis]) / delta, (high[axis] - a[axis]) / delta))
                t0, t1 = max(t0, entry), min(t1, leave)
                if t0 > t1:
                    return False
        return True

    candidates = []
    for anchor, half in zip(anchors, half_sizes):
        radii = np.arange(14., 365., 5.)
        angles = np.arange(96) * 2 * np.pi / 96
        offsets = (radii[:, None, None] * np.column_stack([np.cos(angles), np.sin(angles)])[None, :, :]).reshape(-1, 2)
        centers = anchor + offsets
        valid = np.all(centers - half >= bounds[:2] + 3, axis=1) & np.all(centers + half <= bounds[2:] - 3, axis=1)
        centers = centers[valid]
        free = []
        for chunk in np.array_split(centers, max(1, len(centers) // 250)):
            # Exact circle-vs-rectangle collision check, including marker radii.
            distances = np.maximum(np.abs(chunk[:, None, :] - point_centers[None, :, :]) - half, 0)
            collision = (distances ** 2).sum(axis=2) <= point_radii[None, :] ** 2
            free.extend(chunk[~collision.any(axis=1)])
        if not free:
            raise ValueError("No point-free local label location; increase figure size")
        centers = np.asarray(free)
        candidates.append([(center, endpoint(anchor, center, half)) for center in centers])

    def compatible(i, center, end, placed):
        for j, (other, other_end) in placed.items():
            if np.all(np.abs(center - other) <= half_sizes[i] + half_sizes[j] + 2):
                return False
            if leader_segments_intersect(anchors[i], end, anchors[j], other_end):
                return False
            if hits_box(anchors[i], end, other, half_sizes[j] + 1):
                return False
            if hits_box(anchors[j], other_end, center, half_sizes[i] + 1):
                return False
        return True

    n = len(anchors)
    rng = np.random.default_rng(42)
    difficulty = sorted(range(n), key=lambda i: -np.linalg.norm(candidates[i][0][0] - anchors[i]))
    orders = [difficulty, list(reversed(difficulty)), list(range(n))]
    orders += [rng.permutation(n).tolist() for _ in range(17)]
    best, best_cost = None, np.inf
    for order in orders:
        placed = {}
        for i in order:
            for center, end in candidates[i]:
                if compatible(i, center, end, placed):
                    placed[i] = (center, end)
                    break
            else:
                break
        if len(placed) != n:
            continue
        # Move labels inward where other placements allow it.
        for _ in range(2):
            for i in order:
                old = placed.pop(i)
                for center, end in candidates[i]:
                    if compatible(i, center, end, placed):
                        placed[i] = (center, end)
                        break
                else:
                    placed[i] = old
        cost = sum(np.linalg.norm(center - anchors[i]) for i, (center, _) in placed.items())
        if cost < best_cost:
            best, best_cost = placed, cost
    if best is None:
        raise ValueError("No collision-free local ASV label layout found; no plot saved")
    return [best[i] for i in range(n)]


def add_repelled_asv_labels(ax, selected: pd.DataFrame, colors: dict,
                            obstacles: np.ndarray) -> list[dict]:
    """Place labels near their ASVs, avoiding points, text and leader crossings."""
    selected = selected.loc[selected["plotted"]].copy()
    if selected.empty:
        return []
    xy = selected[["observed_PC1", "observed_PC2"]].to_numpy(float)
    ax.scatter(xy[:, 0], xy[:, 1], marker="D", s=46,
               c=[colors.get(m, "black") for m in selected.module_label],
               edgecolors="black", linewidths=0.7, zorder=24)
    texts = [ax.text(*point, row.ASV_ID, fontsize=10, ha="center", va="center", zorder=26)
             for point, row in zip(xy, selected.itertuples(index=False))]
    ax.figure.canvas.draw()
    renderer = ax.figure.canvas.get_renderer()
    half = np.array([[t.get_window_extent(renderer).width / 2 + 3,
                      t.get_window_extent(renderer).height / 2 + 3] for t in texts])
    anchors = ax.transData.transform(xy)
    # Include all scatter layers: environmental samples, individual ASVs,
    # selected diamonds, and numbered compartment centroids.
    points, radii = [], []
    for collection in ax.collections:
        offsets = np.asarray(collection.get_offsets(), dtype=float)
        if not len(offsets):
            continue
        points.extend(ax.transData.transform(offsets))
        sizes = collection.get_sizes()
        radius = np.sqrt(max(sizes)) / 2 * ax.figure.dpi / 72 + 3 if len(sizes) else 4
        radii.extend([radius] * len(offsets))
    bounds = ax.get_window_extent()
    layout = local_label_layout(anchors, half, np.asarray(points), np.asarray(radii),
                               np.array([bounds.x0, bounds.y0, bounds.x1, bounds.y1]))
    audit = []
    for row, anchor, text, (center, end) in zip(selected.itertuples(index=False), xy, texts, layout):
        label_point, endpoint = ax.transData.inverted().transform([center, end])
        text.set_position(label_point)
        ax.plot([anchor[0], endpoint[0]], [anchor[1], endpoint[1]],
                color=colors.get(row.module_label, "0.3"), linewidth=0.7, zorder=23)
        audit.append(dict(ASV_ID=row.ASV_ID, anchor_PC1=anchor[0], anchor_PC2=anchor[1],
                          label_PC1=label_point[0], label_PC2=label_point[1],
                          label_offset_pixels=float(np.linalg.norm(center - ax.transData.transform(anchor)))))
    return audit


def _plot_asv_module_pca_overlay(
    background: pd.DataFrame,
    centroids: pd.DataFrame,
    positions: pd.DataFrame,
    explained: pd.DataFrame,
    output: Path,
    formats: list[str],
    label_positions: pd.DataFrame | None = None,
) -> None:
    """Overlay module/phylum abundance barycenters on BASIN PC1–PC2."""
    centroid_table = centroids.copy()
    required_centroids = {
        "compartment", "display_label", "PC1_centroid", "PC2_centroid"
    }
    if not required_centroids.issubset(centroid_table):
        raise ValueError(
            "Hybrid centroid table lacks: "
            + ", ".join(sorted(required_centroids - set(centroid_table)))
        )
    centroid_table["responsibility_column"] = centroid_table["compartment"].map(
        lambda value: re.sub(
            r"hybrid_c([0-3])_g([0-9]+)", r"hyb_C\1_G\2", str(value)
        )
    )
    centroid_table = centroid_table.loc[
        centroid_table["responsibility_column"].isin(background.columns)
    ].copy()
    compartments = list(centroid_table["responsibility_column"])
    palette = hybrid_palette(compartments)

    module_levels = sorted(
        positions["module_label"].dropna().astype(str).unique(),
        key=lambda value: int(re.search(r"\d+", value).group())
        if re.search(r"\d+", value) else 999,
    )
    module_colors = {
        module: mcolors.to_hex(plt.get_cmap("tab20")(index / 20.0))
        for index, module in enumerate(module_levels)
    }
    # Manual axes avoid automatic subplot letters: this is one plot plus its key.
    labeled = label_positions is not None
    n_labels = int(label_positions.plotted.sum()) if labeled else 0
    fig = plt.figure(figsize=(27, max(13, 0.60 * n_labels + 3))) if labeled else plt.figure(figsize=(14.6, 8.4))
    # This figure deliberately follows the compact BASIN biplot contract.  The
    # ASPIRE save hook still enforces Times New Roman but does not inflate its
    # carefully matched text hierarchy to the global 22-point minimum.
    fig._aspire_compact_publication_typography = True
    ax = fig.add_axes([0.04, 0.12, 0.41, 0.78] if labeled else [0.07, 0.11, 0.68, 0.82])
    legend_ax = fig.add_axes([0.855, 0.08, 0.14, 0.86] if labeled else [0.77, 0.08, 0.22, 0.86])
    legend_ax.axis("off")
    for compartment in sorted(compartments, key=hybrid_sort_key):
        mask = background["hard_hybrid"].eq(compartment)
        ax.scatter(
            background.loc[mask, "PC1"],
            background.loc[mask, "PC2"],
            s=17, color=lighten(palette[compartment], 0.80),
            edgecolor="none", alpha=1.0, zorder=1,
        )
    ax.axhline(0, linewidth=0.8, color="0.84", zorder=0)
    ax.axvline(0, linewidth=0.8, color="0.84", zorder=0)

    individual_asvs = "ASV_ID" in positions.columns and "phylum" not in positions.columns
    module_counts = (positions.groupby("module_label")["ASV_ID"].nunique()
                     if individual_asvs else positions.groupby("module_label")["asvs_n"].sum())
    # Large modules are drawn first so small modules remain visible on top.
    draw_order = sorted(
        module_levels, key=lambda module: (-module_counts[module], module)
    )
    for draw_index, module in enumerate(draw_order):
        part = positions.loc[positions["module_label"].eq(module)]
        ax.scatter(
            part["observed_PC1"], part["observed_PC2"],
            s=82, color=module_colors[module], edgecolor="white",
            linewidth=0.55, alpha=1.0, zorder=4 + draw_index / 100,
        )

    centroid_table = centroid_table.sort_values(
        "responsibility_column",
        key=lambda values: values.map(hybrid_sort_key),
    ).reset_index(drop=True)
    centroid_table["centroid_number"] = np.arange(1, len(centroid_table) + 1)
    for row in centroid_table.itertuples(index=False):
        ax.scatter(
            row.PC1_centroid, row.PC2_centroid,
            marker="o", s=230, facecolor="white",
            edgecolor="black", linewidth=1.5, zorder=20,
        )
        ax.text(
            row.PC1_centroid, row.PC2_centroid, str(row.centroid_number),
            ha="center", va="center", fontsize=8.5, fontweight="bold",
            color="black", zorder=21,
        )

    x_values = pd.to_numeric(background["PC1"], errors="coerce").dropna()
    y_values = pd.to_numeric(background["PC2"], errors="coerce").dropna()
    x_padding = max(0.35, 0.035 * (x_values.max() - x_values.min()))
    y_padding = max(0.35, 0.035 * (y_values.max() - y_values.min()))
    ax.set_xlim(x_values.min() - x_padding, x_values.max() + x_padding)
    ax.set_ylim(y_values.min() - y_padding, y_values.max() + y_padding)
    xlabel, ylabel = explained_axis_labels(explained)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(
        (r"Observed ASV distributions across hybrid O$_2$-GMM environmental space"
         if individual_asvs else
         r"Observed module–phylum distributions across hybrid O$_2$-GMM environmental space")
    )
    ax.set_box_aspect(1.0)
    module_handles = [
        Line2D(
            [], [], marker="o", linestyle="", markersize=7,
            markerfacecolor=module_colors[module], markeredgecolor="white",
            label=f"{module} (n={int(module_counts[module])} ASVs)",
        )
        for module in module_levels
    ]
    module_legend = legend_ax.legend(
        handles=module_handles, title="Ecological module",
        bbox_to_anchor=(0.0, 1.01), loc="upper left",
        frameon=False, fontsize=6.9, title_fontsize=8.5,
        ncol=2, columnspacing=0.8, handletextpad=0.35,
        labelspacing=0.45,
    )
    legend_ax.add_artist(module_legend)
    centroid_handles = [
        Line2D(
            [], [], marker="", linestyle="",
            label=f"{row.centroid_number}. {row.display_label}",
        )
        for row in centroid_table.itertuples(index=False)
    ]
    legend_ax.legend(
        handles=centroid_handles,
        title=r"Hybrid O$_2$-GMM centroid",
        bbox_to_anchor=(0.0, 0.0), loc="lower left",
        frameon=False, fontsize=6.5, title_fontsize=8.5,
        handlelength=0.0, handletextpad=0.0,
        labelspacing=0.38,
    )
    if label_positions is not None:
        obstacles = np.vstack([
            positions[["observed_PC1", "observed_PC2"]].to_numpy(float),
            centroid_table[["PC1_centroid", "PC2_centroid"]].to_numpy(float),
        ])
        layout = add_repelled_asv_labels(ax, label_positions, module_colors, obstacles)
        tax_ax = fig.add_axes([0.53, 0.06, 0.31, 0.88])
        tax_ax.axis("off")
        selected = label_positions.loc[label_positions.plotted].copy()
        selected["_module_int"] = selected.module_label.str.extract(r"(\d+)", expand=False).astype(int)
        selected["_asv_int"] = selected.ASV_ID.str.extract(r"(\d+)", expand=False).astype(int)
        selected = selected.sort_values(["_module_int", "_asv_int", "ASV_ID"])
        tax_handles, header_indices, legend_groups = [], [], []
        for module, members in selected.groupby("module_label", sort=False):
            peak = members.module_peak_hybrid_compartment.iloc[0]
            peak = str(peak).replace("__gmm", "-GMM") if pd.notna(peak) else "Unavailable"
            header_indices.append(len(tax_handles))
            tax_handles.append(Line2D([], [], linestyle="", marker="",
                                     label=f"{module} — Module peak: {peak}"))
            legend_groups.append(dict(module=module, module_peak=peak, ASVs=members.ASV_ID.tolist()))
            for row in members.itertuples(index=False):
                label = (f"  {row.ASV_ID} — Phylum: {getattr(row, 'phylum', 'Unclassified')}\n"
                         f"  Family: {getattr(row, 'family', 'Unclassified')}\n"
                         f"  Genus: {getattr(row, 'genus', 'Unclassified')}")
                tax_handles.append(Line2D([], [], marker="D", linestyle="", markersize=6,
                                         markerfacecolor=module_colors.get(module, "black"),
                                         markeredgecolor="black", label=label))
        threshold = float(label_positions.selection_threshold.iloc[0]) * 100 if len(label_positions) else 1
        tax_legend = tax_ax.legend(handles=tax_handles, title=f"ASV showcase — mean relative abundance >{threshold:g}%",
                                  loc="upper left", frameon=False, fontsize=9.5, title_fontsize=11,
                                  handletextpad=0.7, labelspacing=0.8)
        for index in header_indices:
            tax_legend.get_texts()[index].set_fontweight("bold")
        output.with_suffix(".layout.json").write_text(json.dumps({
            "labels": layout, "legend_groups": legend_groups, "label_overlaps": 0, "leader_intersections": 0,
            "labels_inside_point_axes": True, "label_point_overlaps": 0,
            "leader_label_intersections": 0}, indent=2))
    for fmt in formats:
        fig.savefig(output.with_suffix(f".{fmt}"), dpi=300, bbox_inches=None)
    plt.close(fig)


def plot_asv_module_pca_overlay(
    background: pd.DataFrame,
    centroids: pd.DataFrame,
    positions: pd.DataFrame,
    explained: pd.DataFrame,
    output: Path,
    formats: list[str],
    label_positions: pd.DataFrame | None = None,
) -> None:
    """Render the overlay with the exact typography used by BASIN biplots."""
    basin_biplot_style = {
        "font.family": "Times New Roman",
        "font.serif": ["Times New Roman"],
        "font.size": 11,
        "axes.titlesize": 13,
        "axes.labelsize": 12,
        "axes.linewidth": 0.8,
        "axes.facecolor": "white",
        "axes.grid": False,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
        "legend.title_fontsize": 11,
        "figure.titlesize": 14,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    }
    with plt.rc_context(basin_biplot_style):
        _plot_asv_module_pca_overlay(
            background, centroids, positions, explained, output, formats, label_positions
        )


def eta_squared(values: np.ndarray, labels: np.ndarray) -> float:
    grand = float(np.mean(values))
    total = float(np.sum((values - grand) ** 2))
    if total <= 0:
        return 0.0
    between = sum(len(values[labels == level]) * (float(np.mean(values[labels == level])) - grand) ** 2
                  for level in np.unique(labels))
    return between / total


def group_module_tests(scores: pd.DataFrame, metadata: pd.DataFrame, groupings: list[str],
                       sample_col: str, cruise_col: str, permutations: int, seed: int) -> pd.DataFrame:
    meta = metadata.drop_duplicates(sample_col).copy()
    meta[sample_col] = meta[sample_col].astype(str)
    data = scores.merge(meta, left_on="sampleID", right_on=sample_col, how="inner")
    rng = np.random.default_rng(seed)
    rows: list[dict] = []
    for grouping in groupings:
        if grouping not in data:
            continue
        for module, part in data.dropna(subset=[grouping]).groupby("ecological_module"):
            part = part[["module_relative_abundance", grouping, cruise_col]].dropna().copy()
            if part[grouping].nunique() < 2:
                continue
            constant = part.groupby(cruise_col)[grouping].nunique().max() <= 1
            if constant:
                test = part.groupby(cruise_col, as_index=False).agg(
                    module_relative_abundance=("module_relative_abundance", "mean"),
                    group_level=(grouping, "first"),
                )
                values = test["module_relative_abundance"].to_numpy(float)
                labels = test["group_level"].astype(str).to_numpy()
                blocks = None
                unit = "cruise"
            else:
                values = part["module_relative_abundance"].to_numpy(float)
                labels = part[grouping].astype(str).to_numpy()
                blocks = part[cruise_col].astype(str).to_numpy()
                unit = "bottle_within_cruise"
            observed = eta_squared(values, labels)
            hits = 0
            for _ in range(permutations):
                if blocks is None:
                    shuffled = rng.permutation(labels)
                else:
                    shuffled = labels.copy()
                    for block in np.unique(blocks):
                        idx = np.flatnonzero(blocks == block)
                        shuffled[idx] = rng.permutation(shuffled[idx])
                hits += eta_squared(values, shuffled) >= observed - 1e-15
            means = pd.Series(values).groupby(labels).mean()
            rows.append({
                "grouping": grouping, "ecological_module": module,
                "effect_eta_squared": observed, "p_value": (hits + 1) / (permutations + 1),
                "sample_unit": unit, "n_units": len(values), "n_levels": len(means),
                "highest_group_level": str(means.idxmax()), "highest_group_mean": float(means.max()),
            })
    out = pd.DataFrame(rows)
    if not out.empty:
        out["q_value"] = out.groupby("grouping", group_keys=False)["p_value"].apply(bh)
    return out


def matched_compartment_module_tests(
    scores: pd.DataFrame,
    metadata: pd.DataFrame,
    sample_col: str,
    cruise_col: str,
    permutations: int,
    seed: int,
    minimum_group_n: int = 3,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Compare legacy O2, GMM, and hybrid groupings on one bottle cohort.

    The hybrid grouping has more levels than either parent classification, so
    the returned table includes adjusted R2 in addition to eta-squared. This
    prevents the direct comparison from rewarding extra levels without a
    complexity penalty.
    """
    groupings = ["o2_compartment", "gmm_component", "o2_subcompartment_final"]
    required = [sample_col, cruise_col, *groupings]
    if any(column not in metadata for column in required):
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    cohort = metadata.drop_duplicates(sample_col).dropna(subset=required).copy()
    source_col = "o2_subcompartment_final_assignment_source"
    if source_col in cohort:
        cohort = cohort.loc[cohort[source_col].astype(str).eq("observed")].copy()

    hybrid_counts = cohort["o2_subcompartment_final"].astype(str).value_counts()
    retained_hybrid = set(hybrid_counts[hybrid_counts >= minimum_group_n].index)
    cohort = cohort.loc[
        cohort["o2_subcompartment_final"].astype(str).isin(retained_hybrid)
    ].copy()
    if cohort.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    tests = group_module_tests(
        scores, cohort, groupings, sample_col, cruise_col, permutations, seed
    )
    if tests.empty:
        return tests, pd.DataFrame(), cohort[required]

    denominator = tests["n_units"] - tests["n_levels"]
    tests["adjusted_r_squared"] = np.where(
        denominator > 0,
        1.0
        - (1.0 - tests["effect_eta_squared"])
        * (tests["n_units"] - 1)
        / denominator,
        np.nan,
    )
    tests["common_cohort_samples"] = cohort[sample_col].nunique()
    tests["common_cohort_cruises"] = cohort[cruise_col].nunique()
    tests["minimum_hybrid_group_n"] = minimum_group_n

    summary = (
        tests.groupby("grouping", as_index=False)
        .agg(
            common_cohort_samples=("common_cohort_samples", "first"),
            common_cohort_cruises=("common_cohort_cruises", "first"),
            n_levels=("n_levels", "max"),
            ecological_modules_tested=("ecological_module", "nunique"),
            significant_modules_q05=("q_value", lambda values: int((values < 0.05).sum())),
            median_eta_squared=("effect_eta_squared", "median"),
            maximum_eta_squared=("effect_eta_squared", "max"),
            median_adjusted_r_squared=("adjusted_r_squared", "median"),
            maximum_adjusted_r_squared=("adjusted_r_squared", "max"),
        )
    )
    labels = {
        "o2_compartment": "Legacy O2",
        "gmm_component": "GMM",
        "o2_subcompartment_final": "Hybrid O2-GMM",
    }
    summary.insert(1, "classification", summary["grouping"].map(labels))
    comparable = tests.dropna(subset=["adjusted_r_squared"])
    if comparable.empty:
        winner_counts = pd.Series(dtype=int)
    else:
        winner_counts = comparable.loc[
            comparable.groupby("ecological_module")["adjusted_r_squared"].idxmax(),
            "grouping",
        ].value_counts()
    summary["modules_with_largest_adjusted_r_squared"] = (
        summary["grouping"].map(winner_counts).fillna(0).astype(int)
    )

    cohort_columns = required + ([source_col] if source_col in cohort else [])
    return tests, summary, cohort[cohort_columns].copy()


def group_module_level_summary(
    scores: pd.DataFrame,
    metadata: pd.DataFrame,
    groupings: list[str],
    sample_col: str,
    cruise_col: str,
) -> pd.DataFrame:
    """Report direction and magnitude behind each omnibus module/group test."""
    meta = metadata.drop_duplicates(sample_col).copy()
    meta[sample_col] = meta[sample_col].astype(str)
    data = scores.merge(meta, left_on="sampleID", right_on=sample_col, how="inner")
    rows: list[dict] = []
    for grouping in groupings:
        if grouping not in data:
            continue
        for module, part in data.dropna(subset=[grouping]).groupby("ecological_module"):
            part = part[["module_relative_abundance", grouping, cruise_col]].dropna()
            if part.empty:
                continue
            constant = part.groupby(cruise_col)[grouping].nunique().max() <= 1
            if constant:
                test = part.groupby(cruise_col, as_index=False).agg(
                    module_relative_abundance=("module_relative_abundance", "mean"),
                    group_level=(grouping, "first"),
                )
                unit = "cruise"
            else:
                test = part.rename(columns={grouping: "group_level"})
                unit = "bottle_within_cruise"
            overall = float(test["module_relative_abundance"].mean())
            for level, level_data in test.groupby("group_level"):
                values = level_data["module_relative_abundance"]
                other = test.loc[
                    ~test["group_level"].astype(str).eq(str(level)),
                    "module_relative_abundance",
                ]
                mean = float(values.mean())
                rows.append({
                    "grouping": grouping,
                    "ecological_module": module,
                    "group_level": str(level),
                    "sample_unit": unit,
                    "n_units": len(values),
                    "mean_module_relative_abundance": mean,
                    "median_module_relative_abundance": float(values.median()),
                    "overall_mean_module_relative_abundance": overall,
                    "mean_difference_from_other_levels": (
                        mean - float(other.mean()) if len(other) else np.nan
                    ),
                    "mean_ratio_to_other_levels": (
                        mean / float(other.mean())
                        if len(other) and float(other.mean()) > 0
                        else np.nan
                    ),
                })
    return pd.DataFrame(rows)


SEASON_ORDER = ["Winter", "Spring", "Summer", "Fall"]


def ecological_module_order(value: object) -> tuple[int, str]:
    text = str(value)
    match = re.fullmatch(r"M(\d+)", text)
    return (int(match.group(1)), text) if match else (10**9, text)


def seasonal_module_summaries(
    group_stats: pd.DataFrame,
    group_level_stats: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return explicit module-season inferential and abundance summaries."""
    if group_stats.empty or "grouping" not in group_stats:
        return pd.DataFrame(), pd.DataFrame()
    association = group_stats.loc[group_stats["grouping"].eq("Season")].copy()
    if association.empty:
        return association, pd.DataFrame()
    association = association.sort_values(
        "ecological_module", key=lambda values: values.map(ecological_module_order)
    ).reset_index(drop=True)
    if group_level_stats.empty or "grouping" not in group_level_stats:
        return association, pd.DataFrame()
    abundance = group_level_stats.loc[
        group_level_stats["grouping"].eq("Season")
    ].copy()
    if abundance.empty:
        return association, abundance
    abundance = abundance.merge(
        association[[
            "ecological_module", "effect_eta_squared", "p_value", "q_value",
            "highest_group_level", "n_units", "n_levels",
        ]].rename(columns={"n_units": "association_cruises"}),
        on="ecological_module", how="left", validate="many_to_one",
    )
    seasonal_total = abundance.groupby("ecological_module")[
        "mean_module_relative_abundance"
    ].transform("sum")
    abundance["within_module_season_fraction"] = np.where(
        seasonal_total > 0,
        abundance["mean_module_relative_abundance"] / seasonal_total,
        np.nan,
    )
    abundance["season_rank_within_module"] = abundance.groupby(
        "ecological_module"
    )["mean_module_relative_abundance"].rank(method="min", ascending=False)
    abundance["season_order"] = abundance["group_level"].map(
        {season: index + 1 for index, season in enumerate(SEASON_ORDER)}
    )
    abundance = abundance.assign(
        _module_order=abundance["ecological_module"].map(ecological_module_order)
    ).sort_values(["_module_order", "season_order", "group_level"]).drop(
        columns="_module_order"
    ).reset_index(drop=True)
    return association, abundance


def plot_seasonal_module_abundance(
    abundance: pd.DataFrame,
    association: pd.DataFrame,
    path: Path,
    formats: list[str],
) -> None:
    """Plot comparable within-module seasonal abundance profiles."""
    if abundance.empty:
        return
    modules = sorted(
        abundance["ecological_module"].dropna().astype(str).unique(),
        key=ecological_module_order,
    )
    observed_seasons = set(abundance["group_level"].dropna().astype(str))
    seasons = [season for season in SEASON_ORDER if season in observed_seasons]
    seasons.extend(sorted(observed_seasons - set(seasons)))
    values = abundance.pivot(
        index="ecological_module",
        columns="group_level",
        values="within_module_season_fraction",
    ).reindex(index=modules, columns=seasons)
    q_values = association.set_index("ecological_module")["q_value"]
    labels = [
        f"{module}{' *' if pd.to_numeric(q_values.get(module), errors='coerce') <= 0.05 else ''}"
        for module in modules
    ]
    fig, ax = plt.subplots(figsize=(10.5, max(8.0, 0.55 * len(modules) + 2.8)))
    image = ax.imshow(
        values.to_numpy(dtype=float), aspect="auto", interpolation="nearest",
        cmap="Greys", vmin=0.0,
        vmax=max(0.4, float(np.nanmax(values.to_numpy(dtype=float)))),
    )
    for row in range(len(modules)):
        for column in range(len(seasons)):
            value = values.iloc[row, column]
            if pd.notna(value):
                ax.text(
                    column, row, f"{value:.2f}", ha="center", va="center",
                    color="white" if value >= 0.32 else "black", fontsize=13,
                )
    ax.set_xticks(range(len(seasons)), seasons)
    ax.set_yticks(range(len(modules)), labels)
    ax.set_xlabel("Season")
    ax.set_ylabel("Ecological module")
    ax.set_title(
        "Seasonal distribution of ecological-module abundance",
        loc="left", fontweight="bold", pad=12,
    )
    colorbar = fig.colorbar(image, ax=ax, pad=0.02)
    colorbar.set_label("Within-module share of seasonal mean abundance")
    ax.text(
        0.0, -0.12,
        "* Omnibus season association q ≤ 0.05; cells sum to 1 within each module.",
        transform=ax.transAxes, ha="left", va="top", fontsize=12,
    )
    fig.tight_layout()
    for fmt in formats:
        fig.savefig(path.with_suffix(f".{fmt}"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def module_indicator_concordance(
    isa: pd.DataFrame,
    modules: pd.DataFrame,
    groupings: list[str],
    permutations: int,
    seed: int,
) -> pd.DataFrame:
    """Quantify module/group agreement among significant single-group indicators."""
    membership = modules.assign(ASV_ID=modules["Taxon"].map(norm_asv))[
        ["ASV_ID", "module_label"]
    ].drop_duplicates("ASV_ID")
    rng = np.random.default_rng(seed)
    rows: list[dict] = []
    for grouping in groupings:
        evidence = isa.loc[
            isa["grouping"].eq(grouping)
            & isa["isa_q_value"].le(0.05)
            & isa["association_scope"].eq("single_group")
            & isa["stratified_within_col"].fillna("").eq("")
        ].merge(membership, on="ASV_ID", how="inner")
        evidence = evidence.drop_duplicates(["ASV_ID", "module_label", "group_level"])
        table = pd.crosstab(evidence["module_label"], evidence["group_level"])
        if table.shape[0] < 2 or table.shape[1] < 2:
            rows.append({
                "grouping": grouping,
                "significant_single_group_indicator_asvs": evidence["ASV_ID"].nunique(),
                "ecological_modules": table.shape[0],
                "group_levels": table.shape[1],
                "cramers_v": np.nan,
                "permutation_p_value": np.nan,
            })
            continue
        observed_chi2 = float(chi2_contingency(table, correction=False)[0])
        n = int(table.to_numpy().sum())
        denominator = n * min(table.shape[0] - 1, table.shape[1] - 1)
        observed_v = float(np.sqrt(observed_chi2 / denominator))
        module_labels = evidence["module_label"].to_numpy()
        group_labels = evidence["group_level"].to_numpy()
        hits = 0
        for _ in range(permutations):
            permuted = pd.crosstab(module_labels, rng.permutation(group_labels))
            chi2 = float(chi2_contingency(permuted, correction=False)[0])
            hits += chi2 >= observed_chi2 - 1e-15
        rows.append({
            "grouping": grouping,
            "significant_single_group_indicator_asvs": evidence["ASV_ID"].nunique(),
            "ecological_modules": table.shape[0],
            "group_levels": table.shape[1],
            "cramers_v": observed_v,
            "permutation_p_value": (hits + 1) / (permutations + 1),
        })
    out = pd.DataFrame(rows)
    if not out.empty:
        out["q_value"] = bh(out["permutation_p_value"])
    return out


def stratified_group_module_tests(scores: pd.DataFrame, metadata: pd.DataFrame,
                                  sample_col: str, cruise_col: str,
                                  inner_group: str, outer_groups: list[str],
                                  permutations: int, seed: int) -> pd.DataFrame:
    """Test network-module abundance by bottle group inside each cruise-scale state."""
    parts: list[pd.DataFrame] = []
    for outer_group in outer_groups:
        if outer_group not in metadata or inner_group not in metadata:
            continue
        for outer_value, subset in metadata.dropna(subset=[outer_group]).groupby(outer_group):
            if subset[inner_group].dropna().nunique() < 2 or subset[cruise_col].dropna().nunique() < 2:
                continue
            result = group_module_tests(
                scores, subset, [inner_group], sample_col, cruise_col, permutations, seed,
            )
            if result.empty:
                continue
            result["stratified_within_col"] = outer_group
            result["stratified_within_value"] = str(outer_value)
            parts.append(result)
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def fisher_enrichment(selected: set[str], background: set[str], positives: set[str]) -> tuple[float, float, int, int]:
    selected = selected & background
    other = background - selected
    a, b = len(selected & positives), len(selected - positives)
    c, d = len(other & positives), len(other - positives)
    if not selected or not other:
        return np.nan, np.nan, a, len(selected)
    odds, p = fisher_exact([[a, b], [c, d]], alternative="greater")
    return float(odds), float(p), a, len(selected)


def indicator_module_enrichment(isa: pd.DataFrame, modules: pd.DataFrame, network_asvs: set[str]) -> pd.DataFrame:
    membership = modules.assign(ASV_ID=modules["Taxon"].map(norm_asv))
    rows = []
    strata = ["grouping", "group_level", "association_scope", "stratified_within_col", "stratified_within_value"]
    for keys, hit in isa.groupby(strata, dropna=False):
        grouping, level, association_scope, within_col, within_value = keys
        positives = set(hit.loc[hit["isa_q_value"].le(.05), "ASV_ID"])
        for module, part in membership.groupby("module_label"):
            odds, p, n_hit, n_module = fisher_enrichment(set(part["ASV_ID"]), network_asvs, positives)
            rows.append({"grouping": grouping, "group_level": level,
                         "association_scope": association_scope,
                         "stratified_within_col": within_col, "stratified_within_value": within_value,
                         "ecological_module": module,
                         "indicator_asvs_in_module": n_hit, "module_asvs": n_module,
                         "odds_ratio": odds, "p_value": p})
    out = pd.DataFrame(rows)
    if not out.empty:
        out["q_value"] = out.groupby(
            ["grouping", "group_level", "association_scope", "stratified_within_col", "stratified_within_value"],
            dropna=False, group_keys=False,
        )["p_value"].apply(bh)
    return out


def functional_enrichment(selector: pd.DataFrame, selector_cols: list[str], mapping: pd.DataFrame,
                          functions: pd.DataFrame, targets: list[str], background_mags: set[str]) -> pd.DataFrame:
    links = selector.merge(mapping[["ASV_ID", "mag_match_id"]].drop_duplicates(), on="ASV_ID", how="inner")
    funcs = functions.loc[functions["module_id"].astype(str).isin(targets)].copy()
    funcs["present"] = funcs["present"].astype(str).str.lower().isin({"true", "1", "yes"})
    rows = []
    for keys, part in links.groupby(selector_cols, dropna=False):
        keys = keys if isinstance(keys, tuple) else (keys,)
        selected = set(part["mag_match_id"].dropna().astype(str))
        record_key = dict(zip(selector_cols, keys))
        for target in targets:
            current = funcs.loc[funcs["module_id"].astype(str).eq(target)]
            positive = set(current.loc[current["present"], "mag_match_id"].astype(str))
            odds, p, n_positive, n_selected = fisher_enrichment(selected, background_mags, positive)
            selected_fraction = pd.to_numeric(
                current.loc[current["mag_match_id"].astype(str).isin(selected), "fraction_covered"], errors="coerce"
            ).dropna()
            background_fraction = pd.to_numeric(
                current.loc[current["mag_match_id"].astype(str).isin(background_mags - selected), "fraction_covered"], errors="coerce"
            ).dropna()
            completeness_p = np.nan
            if len(selected_fraction) and len(background_fraction):
                completeness_p = mannwhitneyu(selected_fraction, background_fraction, alternative="greater").pvalue
            rows.append({**record_key, "functional_module_id": target,
                         "linked_unique_mags": n_selected, "function_positive_mags": n_positive,
                         "function_prevalence": n_positive / n_selected if n_selected else np.nan,
                         "median_fraction_covered": selected_fraction.median() if len(selected_fraction) else np.nan,
                         "odds_ratio": odds, "presence_p_value": p,
                         "completeness_p_value": completeness_p})
    out = pd.DataFrame(rows)
    if not out.empty:
        out["presence_q_value"] = out.groupby(selector_cols, group_keys=False)["presence_p_value"].apply(bh)
        out["completeness_q_value"] = out.groupby(selector_cols, group_keys=False)["completeness_p_value"].apply(bh)
    return out


def summarize_compartment_function_tests(
    group_function: pd.DataFrame,
    indicator_audit: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Summarize MAG-linked target-function evidence for the three compartments."""
    compartments = ["o2_compartment", "gmm_component", "o2_subcompartment_final"]
    data = pd.DataFrame()
    if not group_function.empty and "grouping" in group_function:
        data = group_function.loc[
            group_function["grouping"].isin(compartments)
            & group_function["stratified_within_col"].fillna("").eq("")
        ].copy()
    if data.empty:
        summary = pd.DataFrame({"grouping": compartments})
    else:
        summary = (
            data.groupby("grouping", as_index=False)
            .agg(
                group_levels_tested=("group_level", "nunique"),
                target_modules_tested=("functional_module_id", "nunique"),
                grouping_level_function_tests=("functional_module_id", "size"),
                presence_enrichments_q05=(
                    "presence_q_value", lambda values: int((values < 0.05).sum())
                ),
                completeness_enrichments_q05=(
                    "completeness_q_value", lambda values: int((values < 0.05).sum())
                ),
                median_function_prevalence=("function_prevalence", "median"),
                maximum_odds_ratio=("odds_ratio", "max"),
            )
        )
        summary = pd.DataFrame({"grouping": compartments}).merge(
            summary, on="grouping", how="left"
        )
    count_columns = [
        "group_levels_tested", "target_modules_tested",
        "grouping_level_function_tests", "presence_enrichments_q05",
        "completeness_enrichments_q05",
    ]
    for column in count_columns:
        if column not in summary:
            summary[column] = 0
        summary[column] = summary[column].fillna(0).astype(int)
    for column in ("median_function_prevalence", "maximum_odds_ratio"):
        if column not in summary:
            summary[column] = np.nan
    labels = {
        "o2_compartment": "Legacy O2",
        "gmm_component": "GMM",
        "o2_subcompartment_final": "Hybrid O2-GMM",
    }
    summary.insert(1, "classification", summary["grouping"].map(labels))
    summary["presence_enrichment_fraction_q05"] = (
        summary["presence_enrichments_q05"]
        / summary["grouping_level_function_tests"].replace(0, np.nan)
    )
    summary["completeness_enrichment_fraction_q05"] = (
        summary["completeness_enrichments_q05"]
        / summary["grouping_level_function_tests"].replace(0, np.nan)
    )
    if indicator_audit is not None and not indicator_audit.empty:
        audit_columns = [
            "grouping", "fdr_significant_single_group_network_asvs",
            "mag_linked_functional_indicator_asvs", "functional_test_available",
        ]
        summary = summary.merge(
            indicator_audit[[c for c in audit_columns if c in indicator_audit]],
            on="grouping", how="left",
        )
    if "functional_test_available" not in summary:
        summary["functional_test_available"] = summary[
            "grouping_level_function_tests"
        ].gt(0)
    summary["interpretation_status"] = np.where(
        summary["functional_test_available"].fillna(False),
        "tested",
        "not_tested_no_MAG_linked_single_group_indicators",
    )
    return summary


def prepare_functional_indicator_evidence(
    isa: pd.DataFrame,
    network_asvs: set[str],
    mapping: pd.DataFrame,
    groupings: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Restrict group-function tests to specific, corrected network indicators."""
    if isa.empty:
        return isa.copy(), pd.DataFrame()
    significant_single = isa.loc[
        isa["isa_q_value"].le(0.05)
        & isa["association_scope"].eq("single_group")
        & isa["ASV_ID"].isin(network_asvs)
    ].copy()
    mapped_asvs = set(mapping["ASV_ID"].dropna().astype(str))
    rows = []
    for grouping in groupings:
        all_group = isa.loc[
            isa["grouping"].eq(grouping)
            & isa["stratified_within_col"].fillna("").eq("")
        ]
        eligible = significant_single.loc[
            significant_single["grouping"].eq(grouping)
            & significant_single["stratified_within_col"].fillna("").eq("")
        ]
        linked = eligible.loc[eligible["ASV_ID"].isin(mapped_asvs)]
        rows.append({
            "grouping": grouping,
            "unstratified_indicator_asvs_total": all_group["ASV_ID"].nunique(),
            "fdr_significant_single_group_network_asvs": eligible["ASV_ID"].nunique(),
            "mag_linked_functional_indicator_asvs": linked["ASV_ID"].nunique(),
            "represented_group_levels": eligible["group_level"].nunique(),
            "functional_test_available": bool(linked["ASV_ID"].nunique()),
            "filter": "isa_q_value<=0.05; single_group; ecological_network_member",
        })
    return significant_single, pd.DataFrame(rows)


def evidence_chain(isa: pd.DataFrame, modules: pd.DataFrame, nodes: pd.DataFrame,
                   mapping: pd.DataFrame, functions: pd.DataFrame, targets: list[str]) -> pd.DataFrame:
    membership = modules.copy()
    membership["ASV_ID"] = membership["Taxon"].map(norm_asv)
    keep_module = [c for c in ["ASV_ID", "module_label", "node_stability", "is_best"] if c in membership]
    out = isa.merge(membership[keep_module].drop_duplicates("ASV_ID"), on="ASV_ID", how="left")
    node = nodes.copy()
    node_id = "Taxon" if "Taxon" in node else "ASV_ID"
    node["ASV_ID"] = node[node_id].map(norm_asv)
    central = [c for c in ["ASV_ID", "Degree", "Betweenness", "EigenCentral", "Closeness"] if c in node]
    out = out.merge(node[central].drop_duplicates("ASV_ID"), on="ASV_ID", how="left")
    map_cols = [c for c in ["ASV_ID", "genome_id", "mag_match_id", "link_pident", "link_qcov",
                            "taxonomy_validation_status", "spearman_rho", "eligible_mag_count",
                            "asv_link_weight", "mapping_class"] if c in mapping]
    out = out.merge(mapping[map_cols].drop_duplicates(), on="ASV_ID", how="left")
    funcs = functions.loc[functions["module_id"].astype(str).isin(targets)].copy()
    funcs = funcs.rename(columns={"module_id": "functional_module_id", "module_name": "functional_module_name"})
    out = out.merge(funcs, on="mag_match_id", how="left", suffixes=("", "_function"))
    return out.rename(columns={"module_label": "ecological_module"})


def complete_membership_table(
    modules: pd.DataFrame,
    nodes: pd.DataFrame,
    isa: pd.DataFrame,
    mapping: pd.DataFrame,
    functions: pd.DataFrame,
    targets: list[str],
    grades: pd.DataFrame,
) -> pd.DataFrame:
    """One auditable row per network ASV, including all group/MAG/function links."""
    membership = modules.copy()
    membership["ASV_ID"] = membership["Taxon"].map(norm_asv)
    membership = membership.rename(columns={"module_label": "ecological_module"})
    keep = [c for c in (
        "ASV_ID", "ecological_module", "module_id", "node_stability",
        "graph_variant", "method",
    ) if c in membership]
    out = membership[keep].drop_duplicates("ASV_ID")

    node = nodes.copy()
    node_id = "Taxon" if "Taxon" in node else "ASV_ID"
    node["ASV_ID"] = node[node_id].map(norm_asv)
    node_cols = [c for c in (
        "ASV_ID", "Degree", "Betweenness", "EigenCentral", "Closeness",
        "Phylum", "Class", "Order", "Family", "Genus", "Species",
    ) if c in node]
    out = out.merge(node[node_cols].drop_duplicates("ASV_ID"), on="ASV_ID", how="left")

    if not isa.empty:
        isa_work = isa.copy()
        isa_work["membership"] = (
            isa_work["grouping"].fillna("").astype(str) + "="
            + isa_work["group_level"].fillna("").astype(str)
        )
        isa_summary = isa_work.groupby("ASV_ID").agg(
            isa_group_memberships=("membership", lambda x: "|".join(sorted(set(v for v in x if v != "=")))),
            isa_grouping_count=("grouping", "nunique"),
            isa_min_q=("isa_q_value", "min"),
            isa_max_statistic=("isa_statistic", "max"),
        ).reset_index()
        out = out.merge(isa_summary, on="ASV_ID", how="left")

    map_cols = [c for c in (
        "ASV_ID", "mag_match_id", "genome_id", "mapping_class",
        "eligible_mag_count", "taxonomy_validation_status",
    ) if c in mapping]
    mapped = mapping[map_cols].drop_duplicates() if map_cols else pd.DataFrame()
    if not mapped.empty:
        map_summary = mapped.groupby("ASV_ID").agg(
            paired_mag_count=("mag_match_id", "nunique"),
            paired_mag_ids=("mag_match_id", lambda x: "|".join(sorted(set(x.dropna().astype(str))))),
            mapping_classes=("mapping_class", lambda x: "|".join(sorted(set(x.dropna().astype(str))))),
        ).reset_index()
        out = out.merge(map_summary, on="ASV_ID", how="left")
    paired_count = out["paired_mag_count"] if "paired_mag_count" in out else pd.Series(0, index=out.index)
    out["paired_mag_count"] = pd.to_numeric(paired_count, errors="coerce").fillna(0).astype(int)
    out["has_mag_pair"] = out["paired_mag_count"].gt(0)

    target = functions.loc[functions["module_id"].astype(str).isin(targets)].copy()
    if not mapped.empty and not target.empty:
        target["present_bool"] = target["present"].astype(str).str.lower().isin({"true", "1", "yes"})
        target = target.loc[target["present_bool"]]
        asv_functions = mapped[["ASV_ID", "mag_match_id"]].merge(
            target[["mag_match_id", "module_id"]].drop_duplicates(),
            on="mag_match_id", how="inner",
        )
        function_summary = asv_functions.groupby("ASV_ID").agg(
            target_function_modules=("module_id", lambda x: "|".join(sorted(set(x.astype(str))))),
            target_function_count=("module_id", "nunique"),
        ).reset_index()
        out = out.merge(function_summary, on="ASV_ID", how="left")
    function_count = out["target_function_count"] if "target_function_count" in out else pd.Series(0, index=out.index)
    out["target_function_count"] = pd.to_numeric(function_count, errors="coerce").fillna(0).astype(int)
    out["has_target_function"] = out["target_function_count"].gt(0)

    if not grades.empty:
        priority = {
            "candidate_guild": 0,
            "functionally_coherent_module": 1,
            "environment_associated_module": 2,
            "network_module": 3,
            "insufficient_MAG_coverage": 4,
        }
        grade = grades.copy()
        grade["_priority"] = grade["evidence_grade"].map(priority).fillna(9)
        grade = grade.sort_values(
            ["ecological_module", "_priority", "group_q_value", "functional_enrichment_min_q"],
            na_position="last",
        ).drop_duplicates("ecological_module")
        grade_cols = [c for c in (
            "ecological_module", "evidence_grade", "grouping", "group_q_value",
            "indicator_enrichment_min_q", "functional_enrichment_min_q",
            "mag_link_fraction",
        ) if c in grade]
        out = out.merge(grade[grade_cols], on="ecological_module", how="left")
    return out.sort_values(["ecological_module", "Degree", "ASV_ID"], ascending=[True, False, True])


def strongest_modules(membership: pd.DataFrame, max_modules: int) -> list[str]:
    grade_priority = {
        "candidate_guild": 0,
        "functionally_coherent_module": 1,
        "environment_associated_module": 2,
        "network_module": 3,
        "insufficient_MAG_coverage": 4,
    }
    stats = membership.groupby("ecological_module", dropna=False).agg(
        evidence_grade=("evidence_grade", "first"),
        mean_node_stability=("node_stability", "mean"),
        n_asvs=("ASV_ID", "nunique"),
        mag_paired_asvs=("has_mag_pair", "sum"),
        function_linked_asvs=("has_target_function", "sum"),
    ).reset_index()
    stats["_priority"] = stats["evidence_grade"].map(grade_priority).fillna(9)
    stats = stats.sort_values(
        ["_priority", "mean_node_stability", "function_linked_asvs", "mag_paired_asvs", "n_asvs"],
        ascending=[True, False, False, False, False],
    )
    return stats["ecological_module"].dropna().astype(str).head(max(1, max_modules)).tolist()


def plot_guild_heatmap(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    membership: pd.DataFrame,
    selected_modules: list[str],
    groupings: list[str],
    sample_col: str,
    output: Path,
    formats: list[str],
    top_asvs: int,
    mag_only: bool,
) -> None:
    selected = membership.loc[membership["ecological_module"].astype(str).isin(selected_modules)].copy()
    if mag_only:
        selected = selected.loc[selected["has_mag_pair"]]
    asvs = [asv for asv in selected["ASV_ID"].astype(str) if asv in counts.index]
    samples = [sample for sample in metadata[sample_col].astype(str) if sample in counts.columns]
    if not asvs or not samples:
        return
    relative = counts[samples].div(counts[samples].sum(axis=0).replace(0, np.nan), axis=1).fillna(0)
    variance = relative.loc[asvs].var(axis=1).sort_values(ascending=False)
    asvs = variance.head(max(1, top_asvs)).index.tolist()
    selected = selected.drop_duplicates("ASV_ID").set_index("ASV_ID").loc[asvs]
    meta = metadata.drop_duplicates(sample_col).set_index(sample_col).loc[samples]
    sort_cols = [g for g in groupings if g in meta]
    if sort_cols:
        samples = meta.sort_values(sort_cols, kind="stable").index.tolist()
    matrix = np.log10(1.0 + (relative.loc[asvs, samples] * 1_000_000.0))
    n_tracks = len(sort_cols)
    fig = plt.figure(figsize=(max(16, .12 * len(samples) + 8), max(10, .28 * len(asvs) + 4)))
    grid = fig.add_gridspec(n_tracks + 1, 1, height_ratios=([.24] * n_tracks) + [8], hspace=.04)
    for idx, grouping in enumerate(sort_cols):
        labels = meta.loc[samples, grouping].fillna("NA").astype(str)
        levels = list(dict.fromkeys(labels))
        colors = [mcolors.to_hex(plt.get_cmap("tab20")(i % 20)) for i in range(len(levels))]
        values = np.array([[levels.index(value) for value in labels]])
        axis = fig.add_subplot(grid[idx, 0])
        axis.set_label("<colorbar>")  # Metadata track, not an independent figure panel.
        axis.imshow(values, aspect="auto", interpolation="nearest", cmap=mcolors.ListedColormap(colors))
        axis.set_yticks([0], [grouping])
        axis.set_xticks([])
    ax = fig.add_subplot(grid[-1, 0])
    image = ax.imshow(matrix, aspect="auto", interpolation="nearest", cmap="Greys")
    labels = []
    for asv, row in selected.iterrows():
        taxon = next((str(row.get(c, "")).strip() for c in ("Genus", "Family", "Phylum") if str(row.get(c, "")).strip() not in {"", "nan"}), "")
        labels.append(f"{row['ecological_module']} | {asv}" + (f" | {taxon}" if taxon else ""))
    ax.set_yticks(range(len(asvs)), labels)
    ax.set_xticks([])
    ax.set_xlabel(f"Samples ordered by {', '.join(sort_cols)}" if sort_cols else "Samples")
    colorbar = fig.colorbar(image, ax=ax, pad=.01)
    colorbar.set_label("log10(CPM + 1)")
    fig.tight_layout()
    for fmt in formats:
        fig.savefig(output.with_suffix(f".{fmt}"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def guild_grades(group_stats: pd.DataFrame, indicator_stats: pd.DataFrame,
                 function_stats: pd.DataFrame, modules: pd.DataFrame, mapping: pd.DataFrame) -> pd.DataFrame:
    base = modules.copy()
    base["ASV_ID"] = base["Taxon"].map(norm_asv)
    coverage = base.merge(mapping[["ASV_ID", "mag_match_id"]].drop_duplicates(), on="ASV_ID", how="left")
    coverage["linked_asv_id"] = coverage["ASV_ID"].where(coverage["mag_match_id"].notna())
    coverage = coverage.groupby("module_label").agg(
        module_asvs=("ASV_ID", "nunique"),
        linked_asvs=("linked_asv_id", "nunique"),
        linked_mags=("mag_match_id", "nunique"),
    ).reset_index().rename(columns={"module_label": "ecological_module"})
    coverage["mag_link_fraction"] = coverage["linked_asvs"] / coverage["module_asvs"].replace(0, np.nan)
    rows = []
    for _, stat in group_stats.iterrows():
        grouping, module = stat["grouping"], stat["ecological_module"]
        indicator = indicator_stats.loc[(indicator_stats["grouping"] == grouping) &
                                        (indicator_stats["ecological_module"] == module)]
        if "association_scope" in indicator:
            indicator = indicator.loc[indicator["association_scope"].eq("single_group")]
        if "stratified_within_col" in indicator:
            indicator = indicator.loc[indicator["stratified_within_col"].fillna("").eq("")]
        functional = function_stats.loc[function_stats["ecological_module"] == module]
        indicator_q = indicator["q_value"].min() if not indicator.empty else np.nan
        function_q = functional["presence_q_value"].min() if not functional.empty else np.nan
        cov = coverage.loc[coverage["ecological_module"] == module]
        mag_fraction = float(cov["mag_link_fraction"].iloc[0]) if not cov.empty else 0.0
        if stat.get("q_value", 1) <= .05 and indicator_q <= .05 and function_q <= .05 and mag_fraction >= .2:
            grade = "candidate_guild"
        elif function_q <= .05:
            grade = "functionally_coherent_module"
        elif stat.get("q_value", 1) <= .05:
            grade = "environment_associated_module"
        elif mag_fraction <= 0:
            grade = "insufficient_MAG_coverage"
        else:
            grade = "network_module"
        rows.append({"grouping": grouping, "ecological_module": module,
                     "group_q_value": stat.get("q_value"), "indicator_enrichment_min_q": indicator_q,
                     "functional_enrichment_min_q": function_q, "mag_link_fraction": mag_fraction,
                     "evidence_grade": grade})
    return pd.DataFrame(rows)


def dotplot(table: pd.DataFrame, row: str, col: str, size: str, color: str,
            path: Path, formats: list[str], panel: str, title: str) -> None:
    if table.empty:
        return
    rows = sorted(table[row].dropna().astype(str).unique())
    cols = sorted(table[col].dropna().astype(str).unique())
    fig, ax = plt.subplots(figsize=(max(12, 2.2 * len(cols) + 5), max(8, .8 * len(rows) + 3)))
    work = table.copy()
    x = work[col].astype(str).map({v: i for i, v in enumerate(cols)})
    y = work[row].astype(str).map({v: i for i, v in enumerate(rows)})
    sizes = pd.to_numeric(work[size], errors="coerce").fillna(0)
    sizes = 80 + 700 * sizes / max(float(sizes.max()), 1.0)
    colors = pd.to_numeric(work[color], errors="coerce").fillna(0)
    scatter = ax.scatter(x, y, s=sizes, c=colors, cmap="Greys", vmin=0, vmax=max(1, float(colors.max())),
                         edgecolor="black", linewidth=.8)
    ax.set_xticks(range(len(cols)), cols, rotation=30, ha="right")
    ax.set_yticks(range(len(rows)), rows)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title(f"{panel}  {title}", loc="left", fontweight="bold", pad=12)
    cbar = fig.colorbar(scatter, ax=ax, pad=.02)
    cbar.ax.tick_params(labelsize=22)
    fig.tight_layout()
    for fmt in formats:
        fig.savefig(path.with_suffix(f".{fmt}"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_matched_compartment_module_performance(
    table: pd.DataFrame,
    path: Path,
    formats: list[str],
) -> None:
    """Plot matched-cohort adjusted R² values for the three compartment systems."""
    required = {
        "grouping",
        "ecological_module",
        "adjusted_r_squared",
        "q_value",
    }
    if table.empty or not required.issubset(table.columns):
        return
    labels = {
        "o2_compartment": r"O$_2$",
        "gmm_component": "GMM",
        "o2_subcompartment_final": r"O$_2$–GMM",
    }
    work = table.loc[table["grouping"].isin(labels)].copy()
    work["classification"] = work["grouping"].map(labels)
    work["adjusted_r_squared"] = pd.to_numeric(
        work["adjusted_r_squared"], errors="coerce"
    )
    work["q_value"] = pd.to_numeric(work["q_value"], errors="coerce")
    work = work.dropna(subset=["ecological_module", "adjusted_r_squared"])
    if work.empty:
        return

    def module_order(value: object) -> tuple[int, str]:
        text = str(value)
        match = re.fullmatch(r"M(\d+)", text)
        return (int(match.group(1)), text) if match else (10**9, text)

    modules = sorted(work["ecological_module"].astype(str).unique(), key=module_order)
    classifications = [r"O$_2$", "GMM", r"O$_2$–GMM"]
    values = (
        work.pivot(
            index="ecological_module",
            columns="classification",
            values="adjusted_r_squared",
        )
        .reindex(index=modules, columns=classifications)
    )
    q_values = (
        work.pivot(
            index="ecological_module",
            columns="classification",
            values="q_value",
        )
        .reindex(index=modules, columns=classifications)
    )

    fig, ax = plt.subplots(figsize=(9.5, max(8.5, 0.58 * len(modules) + 2.5)))
    image = ax.imshow(
        values.to_numpy(dtype=float),
        aspect="auto",
        interpolation="nearest",
        cmap="Greys",
        vmin=0.0,
        vmax=max(0.75, float(np.nanmax(values.to_numpy(dtype=float)))),
    )
    for row_index, module in enumerate(modules):
        significant = q_values.loc[module].lt(0.05)
        significant_values = values.loc[module].where(significant)
        winner = (
            significant_values.idxmax()
            if significant_values.notna().any()
            else None
        )
        for column_index, classification in enumerate(classifications):
            value = values.loc[module, classification]
            if pd.isna(value):
                continue
            is_winner = classification == winner
            ax.text(
                column_index,
                row_index,
                f"{value:.3f}",
                ha="center",
                va="center",
                color="white" if value >= 0.43 else "black",
                fontweight="bold" if is_winner else "normal",
                fontsize=16,
            )
            if is_winner:
                ax.add_patch(
                    plt.Rectangle(
                        (column_index - 0.48, row_index - 0.48),
                        0.96,
                        0.96,
                        fill=False,
                        edgecolor="black" if value < 0.43 else "white",
                        linewidth=2.2,
                    )
                )
    ax.set_xticks(range(len(classifications)), classifications)
    ax.set_yticks(range(len(modules)), modules)
    ax.set_xlabel("Compartment classification")
    ax.set_ylabel("Ecological module")
    ax.set_title("Ecological-module separation across compartment classifications")
    colorbar = fig.colorbar(image, ax=ax, pad=0.03)
    colorbar.set_label("Adjusted R²")
    colorbar.ax.tick_params(labelsize=18)
    fig.tight_layout()
    for fmt in formats:
        fig.savefig(path.with_suffix(f".{fmt}"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--metadata", type=Path, required=True)
    ap.add_argument("--asv-counts", type=Path, required=True)
    ap.add_argument("--modules", type=Path, required=True)
    ap.add_argument("--node-features", type=Path, required=True)
    ap.add_argument("--isa-dir", type=Path, required=True)
    ap.add_argument("--accepted-mappings", type=Path, required=True,
                    help="Analysis-eligible ASV-MAG mappings; may include unique and ambiguous candidates.")
    ap.add_argument("--functional-modules", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--groupings", default="o2_subcompartment_final,cruise_group,renewal_phase,o2_compartment,gmm_component")
    ap.add_argument("--primary-groupings", default="")
    ap.add_argument("--target-modules", default="auto_n_s")
    ap.add_argument("--target-module-min-fraction", type=float, default=0.5)
    ap.add_argument("--sample-col", default="sampleID")
    ap.add_argument("--cruise-col", default="Cruise")
    ap.add_argument("--permutations", type=int, default=999)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--max-modules", type=int, default=8)
    ap.add_argument("--heatmap-top-asvs", type=int, default=80)
    ap.add_argument("--pca-label-min-mean-abundance", type=float, default=0.01,
                    help="Label ASVs strictly above this mean relative abundance; default 0.01.")
    ap.add_argument("--pca-scores", type=Path)
    ap.add_argument("--pca-explained", type=Path)
    ap.add_argument("--hybrid-assignments", type=Path)
    ap.add_argument("--hybrid-centroids", type=Path)
    ap.add_argument("--taxonomy", type=Path)
    ap.add_argument("--formats", default="pdf,png,svg")
    ap.add_argument("--font-family", default="Times New Roman")
    ap.add_argument("--strict-font", action="store_true")
    args = ap.parse_args()

    configure_font(args.font_family, args.strict_font)
    tables, plots, audit = args.outdir / "tables", args.outdir / "plots", args.outdir / "audit"
    for directory in (tables, plots, audit):
        directory.mkdir(parents=True, exist_ok=True)
    groupings, requested_targets, formats = csv_list(args.groupings), csv_list(args.target_modules), csv_list(args.formats)
    primary_groupings = csv_list(args.primary_groupings)
    metadata = read_table(args.metadata)
    modules = read_table(args.modules)
    nodes = read_table(args.node_features)
    mapping = read_table(args.accepted_mappings)
    functions = read_table(args.functional_modules)
    if any(value.lower() == "auto_n_s" for value in requested_targets):
        targets, cycle_audit = select_well_covered_cycle_modules(
            functions,
            args.target_module_min_fraction,
        )
        cycle_audit.to_csv(
            tables / "nitrogen_sulfur_module_selection.tsv",
            sep="\t",
            index=False,
        )
        if not targets:
            raise SystemExit(
                "Automatic N/S KEGG-module selection found no adequately covered module"
            )
    else:
        targets = requested_targets
    required = {"module_label", "Taxon"} - set(modules)
    if required:
        raise SystemExit(f"Modules table missing: {', '.join(sorted(required))}")
    if mapping.empty or functions.empty:
        raise SystemExit("Analysis-eligible ASV-MAG mappings and functional-module tables are required")
    mapping["ASV_ID"] = mapping["ASV_ID"].map(norm_asv)
    counts = load_counts(args.asv_counts, set(metadata[args.sample_col].astype(str)))
    scores = module_scores(counts, modules)
    scores.to_csv(tables / "sample_ecological_module_scores.tsv", sep="\t", index=False)
    isa = load_isa(args.isa_dir, groupings)
    isa.to_csv(tables / "indicator_asv_evidence.tsv", sep="\t", index=False)
    group_stats = group_module_tests(scores, metadata, groupings, args.sample_col, args.cruise_col,
                                     args.permutations, args.seed)
    group_stats.to_csv(tables / "group_ecological_module_association.tsv", sep="\t", index=False)
    matched_stats, matched_summary, matched_cohort = matched_compartment_module_tests(
        scores, metadata, args.sample_col, args.cruise_col, args.permutations, args.seed
    )
    matched_stats.to_csv(
        tables / "compartment_matched_ecological_module_association.tsv",
        sep="\t", index=False,
    )
    matched_summary.to_csv(
        tables / "compartment_ecological_module_comparison.tsv",
        sep="\t", index=False,
    )
    matched_cohort.to_csv(
        tables / "compartment_ecological_module_cohort.tsv",
        sep="\t", index=False,
    )
    group_level_stats = group_module_level_summary(
        scores, metadata, groupings, args.sample_col, args.cruise_col
    )
    group_level_stats.to_csv(
        tables / "group_level_ecological_module_abundance.tsv",
        sep="\t",
        index=False,
    )
    season_association, season_abundance = seasonal_module_summaries(
        group_stats, group_level_stats
    )
    season_association.to_csv(
        tables / "ecological_module_season_association.tsv",
        sep="\t", index=False,
    )
    season_abundance.to_csv(
        tables / "ecological_module_season_abundance.tsv",
        sep="\t", index=False,
    )
    inner_group = "o2_subcompartment_final" if "o2_subcompartment_final" in groupings else groupings[0]
    outer_groups = [group for group in ("cruise_group", "renewal_phase") if group in groupings]
    stratified_group_stats = stratified_group_module_tests(
        scores, metadata, args.sample_col, args.cruise_col, inner_group, outer_groups,
        args.permutations, args.seed,
    )
    stratified_group_stats.to_csv(
        tables / "stratified_group_ecological_module_association.tsv", sep="\t", index=False,
    )
    network_asvs = set(modules["Taxon"].map(norm_asv))
    indicator_stats = indicator_module_enrichment(isa, modules, network_asvs)
    indicator_stats.to_csv(tables / "indicator_asv_module_enrichment.tsv", sep="\t", index=False)
    concordance = module_indicator_concordance(
        isa, modules, primary_groupings or groupings, args.permutations, args.seed
    )
    concordance.to_csv(
        tables / "ecological_module_biochemical_group_concordance.tsv",
        sep="\t",
        index=False,
    )
    background_mags = set(mapping.loc[mapping["ASV_ID"].isin(network_asvs), "mag_match_id"].dropna().astype(str))
    functional_isa, functional_isa_audit = prepare_functional_indicator_evidence(
        isa, network_asvs, mapping, groupings
    )
    functional_isa.to_csv(
        tables / "functional_indicator_asv_evidence.tsv", sep="\t", index=False
    )
    functional_isa_audit.to_csv(
        tables / "compartment_function_indicator_audit.tsv", sep="\t", index=False
    )
    isa_selector_cols = ["grouping", "group_level", "stratified_within_col", "stratified_within_value"]
    group_function = functional_enrichment(
        functional_isa, isa_selector_cols, mapping, functions, targets, background_mags
    )
    group_function.to_csv(tables / "group_target_function_enrichment.tsv", sep="\t", index=False)
    summarize_compartment_function_tests(
        group_function, functional_isa_audit
    ).to_csv(
        tables / "compartment_target_function_comparison.tsv", sep="\t", index=False
    )
    membership = modules.assign(ASV_ID=modules["Taxon"].map(norm_asv)).rename(columns={"module_label": "ecological_module"})
    module_function = functional_enrichment(membership, ["ecological_module"], mapping, functions, targets, background_mags)
    module_function.to_csv(tables / "ecological_module_target_function_enrichment.tsv", sep="\t", index=False)
    chain = evidence_chain(isa, modules, nodes, mapping, functions, targets)
    chain.to_csv(tables / "group_guild_asv_mag_function_evidence.tsv", sep="\t", index=False)
    grades = guild_grades(group_stats, indicator_stats, module_function, modules, mapping)
    grades.to_csv(tables / "candidate_guild_evidence_grades.tsv", sep="\t", index=False)
    membership_all = complete_membership_table(
        modules, nodes, isa, mapping, functions, targets, grades,
    )
    membership_all.to_csv(
        tables / "ecological_module_asv_membership_all.tsv", sep="\t", index=False,
    )
    membership_mag = membership_all.loc[membership_all["has_mag_pair"]].copy()
    membership_mag.to_csv(
        tables / "ecological_module_asv_membership_mag_paired.tsv", sep="\t", index=False,
    )
    group_membership_cols = [
        "ASV_ID", "ecological_module", "isa_group_memberships", "isa_grouping_count",
        "evidence_grade", "has_mag_pair", "paired_mag_ids", "target_function_modules",
    ]
    membership_all[[c for c in group_membership_cols if c in membership_all]].to_csv(
        tables / "guild_asv_group_membership.tsv", sep="\t", index=False,
    )
    metadata_cols = [args.sample_col, args.cruise_col] + [g for g in groupings if g in metadata]
    metadata[list(dict.fromkeys(metadata_cols))].drop_duplicates(args.sample_col).to_csv(
        tables / "guild_sample_metadata.tsv", sep="\t", index=False,
    )
    coverage = pd.DataFrame([{
        "network_asvs": len(network_asvs),
        "analysis_eligible_linked_asvs": mapping.loc[mapping["ASV_ID"].isin(network_asvs), "ASV_ID"].nunique(),
        "analysis_eligible_unique_mags": len(background_mags),
        "ambiguous_linked_asvs": mapping.loc[mapping.get("eligible_mag_count", 1) > 1, "ASV_ID"].nunique(),
        "target_module_rows": len(functions.loc[functions["module_id"].astype(str).isin(targets)]),
    }])
    coverage.to_csv(audit / "integration_coverage.tsv", sep="\t", index=False)
    (audit / "analysis_parameters.json").write_text(json.dumps(vars(args), indent=2, default=str))
    plot_group_stats = group_stats.loc[group_stats["grouping"].isin(primary_groupings)].copy() if primary_groupings else group_stats
    dotplot(plot_group_stats, "ecological_module", "grouping", "effect_eta_squared", "effect_eta_squared",
            plots / "group_ecological_module_association", formats, "A", "Environmental association of ecological network modules")
    plot_matched_compartment_module_performance(
        matched_stats,
        plots / "matched_compartment_ecological_module_performance",
        formats,
    )
    plot_seasonal_module_abundance(
        season_abundance,
        season_association,
        plots / "ecological_module_seasonal_abundance_heatmap",
        formats,
    )
    pca_overlay_inputs = (
        args.pca_scores,
        args.pca_explained,
        args.hybrid_assignments,
        args.hybrid_centroids,
    )
    if any(pca_overlay_inputs) and not all(pca_overlay_inputs):
        raise SystemExit(
            "The ecological-module PCA overlay requires --pca-scores, "
            "--pca-explained, --hybrid-assignments, and --hybrid-centroids"
        )
    if all(pca_overlay_inputs):
        if args.taxonomy is None:
            raise SystemExit(
                "The ecological-module PCA overlay requires --taxonomy"
            )
        pca_scores = read_table(args.pca_scores)
        pca_explained = read_table(args.pca_explained)
        hybrid_assignments = read_table(args.hybrid_assignments)
        hybrid_centroids = read_table(args.hybrid_centroids)
        taxonomy = read_table(args.taxonomy)
        pca_background = merge_basin_pca_background(
            pca_scores, hybrid_assignments
        )
        asv_pca_positions = build_asv_pca_positions(
            counts, modules, metadata, pca_scores, args.sample_col
        )
        if asv_pca_positions.empty:
            raise SystemExit(
                "The ecological-module PCA overlay produced no ASV positions"
            )
        asv_pca_positions.to_csv(
            tables / "ecological_module_asv_pca_positions.tsv",
            sep="\t", index=False,
        )
        module_phylum_positions = build_module_phylum_pca_positions(
            asv_pca_positions, taxonomy
        )
        module_phylum_positions.to_csv(
            tables / "ecological_module_phylum_pca_positions.tsv",
            sep="\t", index=False,
        )
        hybrid_centroids.to_csv(
            tables / "ecological_module_pca_hybrid_centroids.tsv",
            sep="\t", index=False,
        )
        plot_asv_module_pca_overlay(
            pca_background,
            hybrid_centroids,
            module_phylum_positions,
            pca_explained,
            plots / "ecological_module_asv_pca_overlay",
            formats,
        )
        plot_asv_module_pca_overlay(
            pca_background, hybrid_centroids, asv_pca_positions, pca_explained,
            plots / "ecological_module_asv_pca_overlay_asv_level", formats,
        )
        label_positions = select_pca_label_asvs(
            counts, asv_pca_positions, args.pca_label_min_mean_abundance
        )
        label_positions = attach_label_taxonomy(label_positions, taxonomy)
        label_positions = attach_label_module_peaks(label_positions, matched_stats)
        label_positions.to_csv(tables / "ecological_module_asv_pca_overlay_labels.tsv",
                               sep="\t", index=False)
        plot_asv_module_pca_overlay(
            pca_background, hybrid_centroids, asv_pca_positions, pca_explained,
            plots / "ecological_module_asv_pca_overlay_asv_level_abundant_asvs_labeled",
            formats, label_positions=label_positions,
        )
        print(f"[done] PCA overlay: {int(label_positions.plotted.sum())} labeled ASVs; "
              f"{int((~label_positions.plotted).sum())} selected ASVs lack PCA positions", flush=True)
    dotplot(module_function, "ecological_module", "functional_module_id", "linked_unique_mags", "function_prevalence",
            plots / "ecological_module_target_functions", formats, "B", "Target functions linked to ecological network modules")
    selected_modules = strongest_modules(membership_all, args.max_modules)
    pd.DataFrame({"ecological_module": selected_modules, "visualization_rank": range(1, len(selected_modules) + 1)}).to_csv(
        tables / "guild_visualization_module_selection.tsv", sep="\t", index=False,
    )
    plot_guild_heatmap(
        counts, metadata, membership_all, selected_modules, groupings, args.sample_col,
        plots / "guild_abundance_heatmap_all_asvs", formats, args.heatmap_top_asvs, mag_only=False,
    )
    plot_guild_heatmap(
        counts, metadata, membership_all, selected_modules, groupings, args.sample_col,
        plots / "guild_abundance_heatmap_mag_paired_asvs", formats, args.heatmap_top_asvs, mag_only=True,
    )


if __name__ == "__main__":
    main()
