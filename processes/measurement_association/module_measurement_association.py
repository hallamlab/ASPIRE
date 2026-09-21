#!/usr/bin/env python3
"""Relate SPIEC-EASI ecological modules to measured environmental chemistry.

The analysis deliberately keeps missing chemistry missing.  Module-level
correlations use only samples measured for the corresponding chemical, and the
member-ASV summary joins those results to the independently inferred network
modules.  No environmental value is imputed by this script.
"""

from __future__ import annotations

import argparse
import colorsys
import json
import math
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/aspire_matplotlib")

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle
try:
    import networkx as nx
except ImportError:  # The declared process environment supplies networkx.
    nx = None
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, spearmanr
from statsmodels.stats.multitest import multipletests
from adjustText import adjust_text

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shared_plot_style import install_publication_style

install_publication_style()


OXYGEN_NAMES = ["oxic", "dysoxic", "suboxic", "anoxic"]
O2_COLORS = {
    "oxic": "#ff0000",
    "dysoxic": "#008000",
    "suboxic": "#add8e6",
    "anoxic": "#800080",
}

# Informative ASV-level properties shown on the common, locked SPIEC-EASI
# layout. ``assessed_samples_n`` is intentionally excluded because it is a
# cohort-wide denominator rather than a varying node property. ``Betweenness``
# is the normalized value, so ``Betweenness_norm`` is not plotted a second time.
NETWORK_METRIC_SPECS = {
    "prevalence": ("prevalence", "prevalence"),
    "mean_relative_abundance": ("mean_relative_abundance", "mean relative abundance"),
    "max_relative_abundance": ("max_relative_abundance", "maximum relative abundance"),
    "observed_samples": ("observed_samples_n", "observed-sample count"),
    "degree": ("Degree", "degree"),
    "degree_thresholded": ("Degree_thresholded", "thresholded degree"),
    "strength": ("Strength", "weighted strength"),
    "within_module_degree": ("within_module_degree", "within-module degree"),
    "within_module_degree_z": ("within_module_degree_z", "within-module degree z-score"),
    "cross_module_degree": ("cross_module_degree", "cross-module degree"),
    "betweenness": ("Betweenness", "normalized betweenness centrality"),
    "betweenness_raw": ("Betweenness_raw", "raw betweenness centrality"),
    "closeness": ("Closeness", "closeness centrality"),
    "eigenvector": ("EigenCentral", "eigenvector centrality"),
    "participation": ("Participation", "participation coefficient"),
}


def read_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    sep = "\t" if path.suffix.lower() in {".tsv", ".tab", ".txt"} else ","
    return pd.read_csv(path, sep=sep, low_memory=False)


def csv_list(value: str | None) -> list[str]:
    if not value:
        return []
    return [item.strip() for item in re.split(r"[,|]", value) if item.strip()]


def norm_asv(value: object) -> str:
    return re.sub(r";.*$", "", str(value).strip())


def clean_taxonomy_display(value: object) -> str:
    """Format a SILVA rank for display without changing the source assignment."""
    if value is None or pd.isna(value):
        return "unclassified"
    text = re.sub(r"^[dkpcofgs]__?", "", str(value).strip(), flags=re.I)
    text = text.replace("_", " ")
    text = re.sub(r"\s+", " ", text).strip()
    if text.casefold() in {
        "", "na", "n/a", "nan", "none", "unassigned", "unknown",
        "uncultured", "unclassified",
    }:
        return "unclassified"
    return text


def anchor_taxonomy_label(row: object) -> str:
    """Return a compact, cleaned P/F/G label for an anchor legend row."""
    return (
        f"{row.ASV_ID} — P: {clean_taxonomy_display(row.asv_phylum)}; "
        f"F: {clean_taxonomy_display(row.asv_family)}; "
        f"G: {clean_taxonomy_display(row.asv_genus)}"
    )


def load_anchor_taxonomy(path: Path) -> pd.DataFrame:
    table = read_table(path)
    id_col = next(
        (column for column in table if str(column).strip().casefold() in
         {"feature id", "feature_id", "asv_id", "asv", "id"}),
        table.columns[0],
    )
    taxonomy_col = next(
        (column for column in table if str(column).strip().casefold() in {"taxon", "taxonomy"}),
        None,
    )
    rank_prefixes = {"p": "phylum", "f": "family", "g": "genus"}
    rows = []
    for row in table.itertuples(index=False, name=None):
        record = dict(zip(table.columns, row))
        parsed = {rank: "Unclassified" for rank in rank_prefixes.values()}
        taxonomy = "" if taxonomy_col is None or pd.isna(record.get(taxonomy_col)) else str(record.get(taxonomy_col))
        for part in re.split(r";|\|", taxonomy):
            match = re.match(r"^\s*([pfg])__?\s*(.*?)\s*$", part, flags=re.I)
            if not match:
                continue
            value = match.group(2).strip()
            if value and value.casefold() not in {"na", "nan", "none", "unassigned"}:
                parsed[rank_prefixes[match.group(1).lower()]] = value
        rows.append({
            "ASV_ID": norm_asv(record.get(id_col, "")),
            "asv_taxonomy": taxonomy,
            "asv_phylum": parsed["phylum"],
            "asv_family": parsed["family"],
            "asv_genus": parsed["genus"],
        })
    return pd.DataFrame(rows).drop_duplicates("ASV_ID")


def bh(values: pd.Series) -> np.ndarray:
    numeric = pd.to_numeric(values, errors="coerce").to_numpy(float)
    valid = np.isfinite(numeric)
    adjusted = np.full(len(numeric), np.nan)
    if valid.any():
        adjusted[valid] = multipletests(numeric[valid], method="fdr_bh")[1]
    return adjusted


def load_counts(path: Path) -> pd.DataFrame:
    table = read_table(path)
    id_col = "ASV_ID" if "ASV_ID" in table else table.columns[0]
    table[id_col] = table[id_col].map(norm_asv)
    table = table.drop_duplicates(id_col).set_index(id_col)
    return table.apply(pd.to_numeric, errors="coerce").fillna(0.0)


def module_membership(modules: pd.DataFrame) -> pd.DataFrame:
    asv_col = "Taxon" if "Taxon" in modules else "ASV_ID"
    if asv_col not in modules or "module_label" not in modules:
        raise ValueError("Module table requires Taxon/ASV_ID and module_label")
    out = modules[[asv_col, "module_label"]].copy()
    out["ASV_ID"] = out[asv_col].map(norm_asv)
    return out[["ASV_ID", "module_label"]].drop_duplicates("ASV_ID")


def module_scores(counts: pd.DataFrame, membership: pd.DataFrame) -> pd.DataFrame:
    relative = counts.div(counts.sum(axis=0).replace(0, np.nan), axis=1).fillna(0.0)
    rows = []
    for module, frame in membership.groupby("module_label", sort=True):
        asvs = relative.index.intersection(frame["ASV_ID"])
        if asvs.empty:
            continue
        values = relative.loc[asvs].sum(axis=0)
        rows.extend(
            {"sampleID": str(sample), "ecological_module": str(module),
             "module_relative_abundance": float(value), "module_asvs_n": int(len(asvs))}
            for sample, value in values.items()
        )
    return pd.DataFrame(rows)


def coverage_table(
    measurements: pd.DataFrame,
    metadata: pd.DataFrame,
    sample_col: str,
    cruise_col: str,
    depth_col: str,
    sparse: set[str],
) -> pd.DataFrame:
    meta = metadata.drop_duplicates(sample_col).set_index(sample_col)
    rows = []
    for measurement in measurements.columns:
        values = pd.to_numeric(measurements[measurement], errors="coerce")
        observed = values.notna()
        sample_ids = values.index[observed].intersection(meta.index)
        cruise_values = meta.loc[sample_ids, cruise_col] if cruise_col in meta else pd.Series(dtype=object)
        depths = pd.to_numeric(meta.loc[sample_ids, depth_col], errors="coerce") if depth_col in meta else pd.Series(dtype=float)
        rows.append({
            "measurement": measurement,
            "measurement_class": "sparse" if measurement in sparse else "core",
            "samples_total": int(len(values)),
            "samples_measured": int(observed.sum()),
            "samples_missing": int((~observed).sum()),
            "measured_fraction": float(observed.mean()),
            "cruises_measured": int(cruise_values.nunique(dropna=True)),
            "depths_measured": int(depths.nunique(dropna=True)),
            "minimum_depth": float(depths.min()) if depths.notna().any() else np.nan,
            "maximum_depth": float(depths.max()) if depths.notna().any() else np.nan,
        })
    return pd.DataFrame(rows)


def add_unavailable_measurements_to_coverage(
    coverage: pd.DataFrame,
    selection_audit: pd.DataFrame,
    sparse: set[str],
) -> pd.DataFrame:
    """Retain configured measurements with no usable matched observations."""
    existing = set(coverage["measurement"].astype(str))
    rows = []
    for row in selection_audit.itertuples(index=False):
        measurement = str(row.measurement)
        if measurement in existing:
            continue
        rows.append({
            "measurement": measurement,
            "measurement_class": "sparse" if measurement in sparse else "core",
            "samples_total": int(getattr(row, "overlapping_samples", 0)),
            "samples_measured": int(getattr(row, "measured_samples", 0)),
            "samples_missing": int(getattr(row, "overlapping_samples", 0)) - int(getattr(row, "measured_samples", 0)),
            "measured_fraction": 0.0,
            "cruises_measured": 0,
            "depths_measured": 0,
            "minimum_depth": np.nan,
            "maximum_depth": np.nan,
            "availability_status": str(getattr(row, "exclusion_reason", "unavailable")),
        })
    coverage = coverage.copy()
    coverage["availability_status"] = "retained"
    return pd.concat([coverage, pd.DataFrame(rows)], ignore_index=True, sort=False)


def module_measurement_tests(
    scores: pd.DataFrame,
    measurements: pd.DataFrame,
    metadata: pd.DataFrame,
    sample_col: str,
    cruise_col: str,
    depth_col: str,
    sparse: set[str],
) -> pd.DataFrame:
    wide = scores.pivot(index="sampleID", columns="ecological_module", values="module_relative_abundance")
    meta = metadata.drop_duplicates(sample_col).set_index(sample_col)
    rows = []
    measurement_order = {measurement: index for index, measurement in enumerate(measurements.columns)}
    for module in wide.columns:
        x_all = pd.to_numeric(wide[module], errors="coerce")
        for measurement in measurements.columns:
            y_all = pd.to_numeric(measurements[measurement], errors="coerce")
            common = x_all.index.intersection(y_all.index)
            x, y = x_all.loc[common], y_all.loc[common]
            valid = x.notna() & y.notna()
            ids = common[valid]
            x, y = x.loc[ids], y.loc[ids]
            if len(ids) >= 3 and x.nunique() > 1 and y.nunique() > 1:
                rho, p_value = spearmanr(x, y)
            else:
                rho, p_value = np.nan, np.nan
            cruises = meta.loc[meta.index.intersection(ids), cruise_col] if cruise_col in meta else pd.Series(dtype=object)
            depths = pd.to_numeric(meta.loc[meta.index.intersection(ids), depth_col], errors="coerce") if depth_col in meta else pd.Series(dtype=float)
            rows.append({
                "ecological_module": str(module),
                "measurement": measurement,
                "measurement_order": measurement_order[measurement],
                "measurement_class": "sparse" if measurement in sparse else "core",
                "rho": rho,
                "p_value": p_value,
                "samples_measured": int(len(ids)),
                "cruises_measured": int(cruises.nunique(dropna=True)),
                "depths_measured": int(depths.nunique(dropna=True)),
                "minimum_depth": float(depths.min()) if depths.notna().any() else np.nan,
                "maximum_depth": float(depths.max()) if depths.notna().any() else np.nan,
            })
    result = pd.DataFrame(rows)
    result["q_value"] = bh(result["p_value"])
    result["significant_q05"] = result["q_value"].le(0.05)
    return result.sort_values(["measurement_order", "ecological_module"], na_position="last")


def member_support(
    asv_correlations: pd.DataFrame,
    membership: pd.DataFrame,
    q_threshold: float,
    sparse: set[str],
) -> pd.DataFrame:
    correlations = asv_correlations.copy()
    correlations["ASV_ID"] = correlations["ASV_ID"].map(norm_asv)
    correlations["rho"] = pd.to_numeric(correlations["rho"], errors="coerce")
    correlations["q_value"] = pd.to_numeric(correlations["q_value"], errors="coerce")
    network = correlations.merge(membership, on="ASV_ID", how="inner")
    rows = []
    for measurement, background in network.groupby("measurement", sort=True):
        tested_background = background.dropna(subset=["rho"])
        background_sig = tested_background["q_value"].le(q_threshold)
        for module, frame in background.groupby("module_label", sort=True):
            tested = frame.dropna(subset=["rho"])
            significant = tested["q_value"].le(q_threshold)
            positive = significant & tested["rho"].gt(0)
            negative = significant & tested["rho"].lt(0)
            module_sig = int(significant.sum())
            module_non = int(len(tested) - module_sig)
            outside = tested_background.loc[~tested_background["ASV_ID"].isin(set(tested["ASV_ID"]))]
            outside_sig = int(outside["q_value"].le(q_threshold).sum())
            outside_non = int(len(outside) - outside_sig)
            if len(tested) and len(outside):
                odds_ratio, enrichment_p = fisher_exact(
                    [[module_sig, module_non], [outside_sig, outside_non]],
                    alternative="greater",
                )
            else:
                odds_ratio, enrichment_p = np.nan, np.nan
            sig_rho = tested.loc[significant, "rho"]
            rows.append({
                "ecological_module": str(module),
                "measurement": measurement,
                "measurement_class": "sparse" if measurement in sparse else "core",
                "module_asvs": int(membership.loc[membership["module_label"].eq(module), "ASV_ID"].nunique()),
                "asvs_tested": int(len(tested)),
                "significant_asvs": module_sig,
                "significant_fraction": float(module_sig / len(tested)) if len(tested) else np.nan,
                "significant_positive_asvs": int(positive.sum()),
                "significant_negative_asvs": int(negative.sum()),
                "sign_concordance": float(max(positive.sum(), negative.sum()) / module_sig) if module_sig else np.nan,
                "median_rho_all_asvs": float(tested["rho"].median()) if len(tested) else np.nan,
                "median_rho_significant_asvs": float(sig_rho.median()) if len(sig_rho) else np.nan,
                "maximum_absolute_rho": float(tested["rho"].abs().max()) if len(tested) else np.nan,
                "enrichment_odds_ratio": odds_ratio,
                "enrichment_p_value": enrichment_p,
                "network_asvs_tested": int(len(tested_background)),
                "network_significant_asvs": int(background_sig.sum()),
            })
    result = pd.DataFrame(rows)
    result["enrichment_q_value"] = bh(result["enrichment_p_value"])
    return result.sort_values(["enrichment_q_value", "ecological_module", "measurement"], na_position="last")


def module_order(value: object) -> tuple[int, str]:
    match = re.search(r"(\d+)", str(value))
    return (int(match.group(1)), str(value)) if match else (10**9, str(value))


def plot_heatmap(table: pd.DataFrame, output: Path, formats: list[str]) -> None:
    modules = sorted(table["ecological_module"].unique(), key=module_order)
    measurements = (
        table[["measurement", "measurement_order"]]
        .drop_duplicates().sort_values("measurement_order")["measurement"].tolist()
    )
    rho = table.pivot(index="ecological_module", columns="measurement", values="rho").reindex(index=modules, columns=measurements)
    qval = table.pivot(index="ecological_module", columns="measurement", values="q_value").reindex(index=modules, columns=measurements)
    nvals = table.pivot(index="ecological_module", columns="measurement", values="samples_measured").reindex(index=modules, columns=measurements)
    classes = table.drop_duplicates("measurement").set_index("measurement")["measurement_class"]
    fig, ax = plt.subplots(figsize=(max(13, 0.85 * len(measurements) + 5), max(9, 0.55 * len(modules) + 3)))
    fig._aspire_compact_publication_typography = True
    image = ax.imshow(rho.to_numpy(float), cmap="coolwarm", vmin=-1, vmax=1, aspect="auto")
    for i, module in enumerate(modules):
        for j, measurement in enumerate(measurements):
            value = rho.loc[module, measurement]
            if pd.isna(value):
                continue
            star = "*" if qval.loc[module, measurement] <= 0.05 else ""
            ax.text(j, i, f"{value:.2f}{star}\n(n={int(nvals.loc[module, measurement])})",
                    ha="center", va="center", fontsize=7.2,
                    color="white" if abs(value) >= 0.45 else "black")
    labels = [f"{m}\n({'sparse' if classes.get(m) == 'sparse' else 'core'})" for m in measurements]
    ax.set_xticks(range(len(measurements)), labels, rotation=45, ha="right")
    ax.set_yticks(range(len(modules)), modules)
    ax.set_xlabel("Measured physicochemical variable")
    ax.set_ylabel("Ecological module")
    ax.set_title("Ecological-module associations with measured physicochemistry")
    colorbar = fig.colorbar(image, ax=ax, pad=0.02)
    colorbar.set_label("Spearman ρ")
    fig.tight_layout()
    for fmt in formats:
        fig.savefig(output.with_suffix(f".{fmt}"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def hybrid_sort(value: object) -> tuple[int, int, str]:
    match = re.search(r"C([0-3])_G(\d+)", str(value), flags=re.I)
    return (int(match.group(1)), int(match.group(2)), str(value)) if match else (99, 99, str(value))


def hybrid_palette(columns: list[str]) -> dict[str, tuple[float, float, float]]:
    parsed = []
    for column in columns:
        match = re.search(r"C([0-3])_G(\d+)", column, flags=re.I)
        if match:
            parsed.append((int(match.group(1)), int(match.group(2)), column))
    gmms = sorted({gmm for _, gmm, _ in parsed})
    lightness = dict(zip(gmms, np.linspace(0.30, 0.76, max(1, len(gmms)))))
    palette = {}
    for oxygen, gmm, column in parsed:
        hue, _, saturation = colorsys.rgb_to_hls(*mcolors.to_rgb(O2_COLORS[OXYGEN_NAMES[oxygen]]))
        palette[column] = colorsys.hls_to_rgb(hue, float(lightness[gmm]), max(0.58, saturation))
    return palette


def lighten(color: object, fraction: float = 0.80) -> tuple[float, float, float]:
    rgb = np.asarray(mcolors.to_rgb(color), dtype=float)
    return tuple(rgb + (1.0 - rgb) * fraction)


def pca_sample_coordinates(metadata: pd.DataFrame, scores: pd.DataFrame, sample_col: str, depth_col: str) -> pd.DataFrame:
    score_depth = "Depth_anchored" if "Depth_anchored" in scores else "Depth"
    meta = metadata[[sample_col, "Cruise", depth_col]].drop_duplicates(sample_col).copy()
    meta["_cruise"] = pd.to_numeric(meta["Cruise"], errors="coerce")
    meta["_depth"] = pd.to_numeric(meta[depth_col], errors="coerce").round(6)
    coords = scores[["Cruise", score_depth, "PC1", "PC2"]].copy()
    coords["_cruise"] = pd.to_numeric(coords["Cruise"], errors="coerce")
    coords["_depth"] = pd.to_numeric(coords[score_depth], errors="coerce").round(6)
    coords["PC1"] = pd.to_numeric(coords["PC1"], errors="coerce")
    coords["PC2"] = pd.to_numeric(coords["PC2"], errors="coerce")
    return meta.merge(coords[["_cruise", "_depth", "PC1", "PC2"]].drop_duplicates(["_cruise", "_depth"]),
                      on=["_cruise", "_depth"], how="inner")


def module_pca_positions(scores: pd.DataFrame, sample_coords: pd.DataFrame, sample_col: str) -> pd.DataFrame:
    data = scores.merge(sample_coords[[sample_col, "_cruise", "_depth", "PC1", "PC2"]],
                        left_on="sampleID", right_on=sample_col, how="inner")
    collapsed = data.groupby(["ecological_module", "_cruise", "_depth", "PC1", "PC2"], as_index=False)["module_relative_abundance"].median()
    rows = []
    for module, frame in collapsed.groupby("ecological_module"):
        weights = frame["module_relative_abundance"].clip(lower=0)
        if weights.sum() <= 0:
            continue
        rows.append({
            "ecological_module": module,
            "PC1_centroid": float(np.average(frame["PC1"], weights=weights)),
            "PC2_centroid": float(np.average(frame["PC2"], weights=weights)),
            "matched_environmental_positions_n": int(len(frame)),
            "module_abundance_weight_sum": float(weights.sum()),
        })
    return pd.DataFrame(rows)


def load_asv_graph(path: Path) -> nx.Graph:
    if nx is None:
        raise ImportError("networkx is required for the ecological-module network renderer")
    graph = nx.read_graphml(path)
    mapping = {}
    for node, attributes in graph.nodes(data=True):
        asv = attributes.get(
            "Taxon", attributes.get("ASV_ID", attributes.get("name", node))
        )
        mapping[node] = norm_asv(asv)
    return nx.relabel_nodes(graph, mapping, copy=True)


def network_anchor_table(
    node_features: pd.DataFrame,
    membership: pd.DataFrame,
    anchor_top_n: int,
) -> pd.DataFrame:
    features = node_features.copy()
    if "Taxon" in features:
        features = features.rename(columns={"Taxon": "ASV_ID"})
    elif "ASV_ID" not in features:
        features = features.rename(columns={features.columns[0]: "ASV_ID"})
    features["ASV_ID"] = features["ASV_ID"].map(norm_asv)
    nodes = membership.merge(features.drop_duplicates("ASV_ID"), on="ASV_ID", how="left")
    rank_specs = {
        "Degree": "anchor_rank_degree",
        "EigenCentral": "anchor_rank_eigencentral",
        "Betweenness": "anchor_rank_betweenness",
    }
    flags = []
    for metric, rank_column in rank_specs.items():
        nodes[metric] = pd.to_numeric(nodes.get(metric), errors="coerce")
        nodes[rank_column] = (
            nodes.groupby("module_label")[metric]
            .rank(method="dense", ascending=False, na_option="bottom")
            .astype("Int64")
        )
        flag = rank_column.replace("rank", "is")
        nodes[flag] = nodes[rank_column].le(anchor_top_n).fillna(False)
        flags.append(flag)
    nodes["is_ecological_anchor"] = nodes[flags].any(axis=1)
    return nodes


def add_global_network_metrics(
    nodes: pd.DataFrame,
    graph: nx.Graph,
    raw_counts: pd.DataFrame,
) -> pd.DataFrame:
    """Add occupancy and whole-network role metrics to the network node table."""
    out = nodes.copy()
    raw = raw_counts.reindex(out["ASV_ID"]).apply(pd.to_numeric, errors="coerce").fillna(0.0)
    sample_totals = raw.sum(axis=0).replace(0, np.nan)
    relative = raw.divide(sample_totals, axis=1)
    occupancy = pd.DataFrame({
        "ASV_ID": raw.index,
        "prevalence": (raw.gt(0).sum(axis=1) / max(raw.shape[1], 1)).to_numpy(float),
        "mean_relative_abundance": relative.mean(axis=1, skipna=True).fillna(0).to_numpy(float),
        "max_relative_abundance": relative.max(axis=1, skipna=True).fillna(0).to_numpy(float),
        "observed_samples_n": raw.gt(0).sum(axis=1).to_numpy(int),
        "assessed_samples_n": int(raw.shape[1]),
    })
    out = out.merge(occupancy, on="ASV_ID", how="left", validate="one_to_one")

    module_lookup = out.set_index("ASV_ID")["module_label"].astype(str).to_dict()
    strength: dict[str, float] = {asv: 0.0 for asv in out["ASV_ID"]}
    within_degree: dict[str, int] = {asv: 0 for asv in out["ASV_ID"]}
    cross_degree: dict[str, int] = {asv: 0 for asv in out["ASV_ID"]}
    neighbor_modules: dict[str, list[str]] = {asv: [] for asv in out["ASV_ID"]}
    for source, target, attributes in graph.edges(data=True):
        source, target = norm_asv(source), norm_asv(target)
        if source not in module_lookup or target not in module_lookup:
            continue
        try:
            weight = abs(float(attributes.get("weight", 1.0)))
        except (TypeError, ValueError):
            weight = 1.0
        strength[source] += weight
        strength[target] += weight
        neighbor_modules[source].append(module_lookup[target])
        neighbor_modules[target].append(module_lookup[source])
        if module_lookup[source] == module_lookup[target]:
            within_degree[source] += 1
            within_degree[target] += 1
        else:
            cross_degree[source] += 1
            cross_degree[target] += 1

    def participation(modules: list[str]) -> float:
        if not modules:
            return 0.0
        fractions = pd.Series(modules).value_counts().to_numpy(float) / len(modules)
        return float(1.0 - np.square(fractions).sum())

    out["Strength"] = out["ASV_ID"].map(strength).fillna(0.0)
    out["within_module_degree"] = out["ASV_ID"].map(within_degree).fillna(0).astype(int)
    out["cross_module_degree"] = out["ASV_ID"].map(cross_degree).fillna(0).astype(int)
    out["Participation"] = out["ASV_ID"].map(
        lambda asv: participation(neighbor_modules.get(asv, []))
    )
    grouped = out.groupby("module_label")["within_module_degree"]
    module_mean = grouped.transform("mean")
    module_sd = grouped.transform("std").replace(0, np.nan)
    out["within_module_degree_z"] = (
        (out["within_module_degree"] - module_mean) / module_sd
    ).fillna(0.0)

    for key, (column, _) in NETWORK_METRIC_SPECS.items():
        values = pd.to_numeric(out[column], errors="coerce")
        out[f"global_rank_{key}"] = values.rank(
            method="min", ascending=False, na_option="bottom"
        ).astype("Int64")
        out[f"global_percentile_{key}"] = values.rank(
            method="average", pct=True, na_option="bottom"
        )
    return out


def top_network_players(nodes: pd.DataFrame, top_n: int) -> pd.DataFrame:
    """Return a long, deterministic table of the leading ASVs for each metric."""
    rows = []
    detail = [
        "ASV_ID", "ecological_module", "asv_taxonomy", "asv_phylum",
        "asv_family", "asv_genus", "is_ecological_anchor",
    ]
    for metric, (column, _) in NETWORK_METRIC_SPECS.items():
        ranked = nodes.sort_values(
            [column, "ASV_ID"], ascending=[False, True], na_position="last"
        ).head(top_n)
        for rank, row in enumerate(ranked.itertuples(index=False), start=1):
            record = {name: getattr(row, name) for name in detail}
            record.update({
                "metric": metric,
                "rank": rank,
                "value": getattr(row, column),
                "taxonomy_display": anchor_taxonomy_label(row).split(" — ", 1)[1],
            })
            rows.append(record)
    return pd.DataFrame(rows)[
        ["metric", "rank", "value", "taxonomy_display", *detail]
    ]


def pack_module_centers(
    original: np.ndarray,
    radii: np.ndarray,
    fixed: np.ndarray,
    iterations: int = 500,
) -> np.ndarray:
    """Separate module islands while retaining their environmental ordering."""
    centers = np.asarray(original, dtype=float).copy()
    fixed = np.asarray(fixed, dtype=float).reshape((-1, 2)) if len(fixed) else np.empty((0, 2))
    for iteration in range(iterations):
        shift = np.zeros_like(centers)
        for i in range(len(centers)):
            for j in range(i + 1, len(centers)):
                delta = centers[i] - centers[j]
                distance = float(np.linalg.norm(delta))
                required = float(radii[i] + radii[j] + 0.10)
                if distance < required:
                    if distance < 1e-12:
                        angle = (i + j + 1) * 2.399963229728653
                        direction = np.array([math.cos(angle), math.sin(angle)])
                    else:
                        direction = delta / distance
                    displacement = direction * (required - distance) * 0.24
                    shift[i] += displacement
                    shift[j] -= displacement
            for fixed_index, obstacle in enumerate(fixed):
                delta = centers[i] - obstacle
                distance = float(np.linalg.norm(delta))
                required = float(radii[i] + 0.24)
                if distance < required:
                    if distance < 1e-12:
                        angle = (i + fixed_index + 1) * 2.399963229728653
                        direction = np.array([math.cos(angle), math.sin(angle)])
                    else:
                        direction = delta / distance
                    shift[i] += direction * (required - distance) * 0.28
        attraction = 0.012 if iteration < iterations * 0.75 else 0.025
        centers += shift + attraction * (original - centers)
    return centers


def perimeter_module_centers(
    original: np.ndarray,
    radii: np.ndarray,
    background_xy: np.ndarray,
) -> np.ndarray:
    """Place module-network callouts outside the environmental sample cloud.

    Callout order follows the angular order of the unmodified module centroids,
    which limits connector crossings.  The enclosing ellipse is deliberately
    larger than the sample bounding box so every island (including its radius)
    lies outside the cloud.  Only these display coordinates are displaced; the
    scientific centroid coordinates remain unchanged.
    """
    original = np.asarray(original, dtype=float)
    radii = np.asarray(radii, dtype=float)
    background_xy = np.asarray(background_xy, dtype=float)
    finite_background = background_xy[np.isfinite(background_xy).all(axis=1)]
    if len(original) == 0:
        return original.copy()
    if len(finite_background) == 0:
        finite_background = original[np.isfinite(original).all(axis=1)]
    if len(finite_background) == 0:
        finite_background = np.zeros((1, 2), dtype=float)

    lower = finite_background.min(axis=0)
    upper = finite_background.max(axis=0)
    center = (lower + upper) / 2.0
    half_span = np.maximum((upper - lower) / 2.0, 0.5)
    normalized = (original - center) / half_span
    original_angles = np.arctan2(normalized[:, 1], normalized[:, 0])
    angular_order = np.argsort(original_angles, kind="stable")

    # Match equally spaced callout slots to the original circular ordering.
    canonical = 2.0 * np.pi * np.arange(len(original)) / len(original)
    ordered_angles = np.unwrap(original_angles[angular_order])
    offset_vectors = np.exp(1j * (ordered_angles - canonical))
    offset = float(np.angle(offset_vectors.mean()))
    slot_angles = offset + canonical

    largest_radius = float(np.nanmax(radii)) if len(radii) else 0.0
    clearance = max(0.32, 0.45 * largest_radius)
    # sqrt(2) ensures that every point on the ellipse is beyond at least one
    # edge of the expanded sample bounding box, including diagonal positions.
    base_axes = np.sqrt(2.0) * (half_span + largest_radius + clearance)
    scale = 1.0
    ordered_radii = radii[angular_order]
    for _ in range(80):
        ordered_display = center + np.column_stack((
            base_axes[0] * scale * np.cos(slot_angles),
            base_axes[1] * scale * np.sin(slot_angles),
        ))
        separated = True
        for i in range(len(ordered_display)):
            for j in range(i + 1, len(ordered_display)):
                distance = float(np.linalg.norm(ordered_display[i] - ordered_display[j]))
                required = float(ordered_radii[i] + ordered_radii[j] + 0.14)
                if distance < required:
                    separated = False
                    break
            if not separated:
                break
        if separated:
            break
        scale *= 1.04

    displayed = np.empty_like(original)
    displayed[angular_order] = ordered_display
    return displayed


def environmental_module_network(
    graph: nx.Graph,
    node_features: pd.DataFrame,
    membership: pd.DataFrame,
    taxonomy: pd.DataFrame,
    raw_counts: pd.DataFrame,
    module_centroids: pd.DataFrame,
    hybrid_centroids: pd.DataFrame,
    background: pd.DataFrame,
    anchor_top_n: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    nodes = network_anchor_table(node_features, membership, anchor_top_n)
    nodes = nodes.merge(taxonomy, on="ASV_ID", how="left")
    for column in ("asv_phylum", "asv_family", "asv_genus"):
        nodes[column] = nodes[column].fillna("Unclassified")
    nodes = add_global_network_metrics(nodes, graph, raw_counts)
    modules = sorted(module_centroids["ecological_module"].astype(str), key=module_order)
    centers = module_centroids.set_index("ecological_module").loc[modules]
    counts = nodes.groupby("module_label")["ASV_ID"].nunique().reindex(modules).fillna(0)
    # Display-only radii are intentionally generous because each complete
    # within-module network is rendered as an exterior callout.
    radii = 0.30 + 0.065 * np.sqrt(counts.to_numpy(float))

    centroid_table = hybrid_centroids.copy()
    centroid_table["responsibility_column"] = centroid_table["compartment"].map(
        lambda value: re.sub(r"hybrid_c([0-3])_g(\d+)", r"hyb_C\1_G\2", str(value))
    )
    hybrid_cols = [column for column in background if re.fullmatch(r"hyb_C[0-3]_G\d+", str(column))]
    background_hard = background[hybrid_cols].apply(pd.to_numeric, errors="coerce").fillna(0).idxmax(axis=1)
    centroid_table = centroid_table.loc[
        centroid_table["responsibility_column"].isin(set(background_hard.astype(str)))
    ].sort_values("responsibility_column", key=lambda values: values.map(hybrid_sort)).reset_index(drop=True)
    centroid_table["centroid_number"] = np.arange(1, len(centroid_table) + 1)
    original = centers[["PC1_centroid", "PC2_centroid"]].to_numpy(float)
    background_xy = background[["PC1", "PC2"]].apply(
        pd.to_numeric, errors="coerce"
    ).to_numpy(float)
    displayed = perimeter_module_centers(original, radii, background_xy)

    layout_rows = []
    edge_rows = []
    feature_index = nodes.set_index("ASV_ID")
    module_lookup = membership.set_index("ASV_ID")["module_label"].astype(str).to_dict()
    for source, target, attributes in graph.edges(data=True):
        source, target = norm_asv(source), norm_asv(target)
        source_module, target_module = module_lookup.get(source), module_lookup.get(target)
        try:
            weight = float(attributes.get("weight", 1.0))
        except (TypeError, ValueError):
            weight = np.nan
        edge_rows.append({
            "source": source,
            "target": target,
            "source_module": source_module,
            "target_module": target_module,
            "within_module": bool(source_module and source_module == target_module),
            "displayed": bool(source_module and source_module == target_module),
            "weight": weight,
            "sign": "negative" if np.isfinite(weight) and weight < 0 else "positive",
        })

    for index, module in enumerate(modules):
        module_asvs = nodes.loc[nodes["module_label"].astype(str).eq(module), "ASV_ID"].tolist()
        subgraph = graph.subgraph([asv for asv in module_asvs if asv in graph]).copy()
        if len(subgraph) == 1:
            local = {next(iter(subgraph.nodes)): np.array([0.0, 0.0])}
        elif len(subgraph) > 1:
            local = nx.spring_layout(
                subgraph,
                seed=42 + index,
                weight="weight",
                iterations=350,
                k=max(0.12, 1.15 / math.sqrt(len(subgraph))),
            )
        else:
            local = {}
        if local:
            maximum = max(float(np.linalg.norm(value)) for value in local.values())
            maximum = max(maximum, 1e-9)
        else:
            maximum = 1.0
        for asv in module_asvs:
            local_xy = np.asarray(local.get(asv, np.zeros(2)), dtype=float) / maximum * radii[index] * 0.84
            row = feature_index.loc[asv] if asv in feature_index.index else pd.Series(dtype=object)
            layout_rows.append({
                "ASV_ID": asv,
                "ecological_module": module,
                "is_ecological_anchor": bool(row.get("is_ecological_anchor", False)),
                "anchor_rank_degree": row.get("anchor_rank_degree", pd.NA),
                "anchor_rank_eigencentral": row.get("anchor_rank_eigencentral", pd.NA),
                "anchor_rank_betweenness": row.get("anchor_rank_betweenness", pd.NA),
                "Degree": row.get("Degree", np.nan),
                "Degree_thresholded": row.get("Degree_thresholded", np.nan),
                "EigenCentral": row.get("EigenCentral", np.nan),
                "Betweenness": row.get("Betweenness", np.nan),
                "Betweenness_raw": row.get("Betweenness_raw", np.nan),
                "Betweenness_norm": row.get("Betweenness_norm", np.nan),
                "Closeness": row.get("Closeness", np.nan),
                "Strength": row.get("Strength", np.nan),
                "Participation": row.get("Participation", np.nan),
                "within_module_degree": row.get("within_module_degree", np.nan),
                "cross_module_degree": row.get("cross_module_degree", np.nan),
                "within_module_degree_z": row.get("within_module_degree_z", np.nan),
                "prevalence": row.get("prevalence", np.nan),
                "observed_samples_n": row.get("observed_samples_n", np.nan),
                "assessed_samples_n": row.get("assessed_samples_n", np.nan),
                "mean_relative_abundance": row.get("mean_relative_abundance", np.nan),
                "max_relative_abundance": row.get("max_relative_abundance", np.nan),
                **{
                    column: row.get(column, np.nan)
                    for column in nodes.columns
                    if column.startswith("global_rank_") or column.startswith("global_percentile_")
                },
                "asv_taxonomy": row.get("asv_taxonomy", ""),
                "asv_phylum": row.get("asv_phylum", "Unclassified"),
                "asv_family": row.get("asv_family", "Unclassified"),
                "asv_genus": row.get("asv_genus", "Unclassified"),
                "PC1_centroid_unmodified": original[index, 0],
                "PC2_centroid_unmodified": original[index, 1],
                "network_center_x_display": displayed[index, 0],
                "network_center_y_display": displayed[index, 1],
                "network_local_x": local_xy[0],
                "network_local_y": local_xy[1],
                "network_x_display": displayed[index, 0] + local_xy[0],
                "network_y_display": displayed[index, 1] + local_xy[1],
            })
    layout = pd.DataFrame(layout_rows)
    module_layout = pd.DataFrame({
        "ecological_module": modules,
        "asvs_n": counts.to_numpy(int),
        "island_radius": radii,
        "PC1_centroid_unmodified": original[:, 0],
        "PC2_centroid_unmodified": original[:, 1],
        "network_center_x_display": displayed[:, 0],
        "network_center_y_display": displayed[:, 1],
    })
    return layout, pd.DataFrame(edge_rows), module_layout, centroid_table


def vector_table(loadings: pd.DataFrame, measurements: pd.DataFrame,
                 sample_coords: pd.DataFrame, sample_col: str, sparse: set[str]) -> pd.DataFrame:
    feature_col = loadings.columns[0]
    load = loadings.rename(columns={feature_col: "measurement"})
    rows = []
    for row in load.itertuples(index=False):
        measurement = str(getattr(row, "measurement"))
        if measurement in measurements.columns and measurement not in sparse:
            rows.append({"measurement": measurement, "measurement_class": "core",
                         "vector_PC1": float(getattr(row, "PC1")), "vector_PC2": float(getattr(row, "PC2")),
                         "samples_measured": int(measurements[measurement].notna().sum()),
                         "vector_method": "PCA loading"})
    lookup = sample_coords.drop_duplicates(sample_col).set_index(sample_col)
    for measurement in sorted(sparse):
        if measurement not in measurements:
            continue
        ids = measurements.index.intersection(lookup.index)
        values = pd.to_numeric(measurements.loc[ids, measurement], errors="coerce")
        pc1 = pd.to_numeric(lookup.loc[ids, "PC1"], errors="coerce")
        pc2 = pd.to_numeric(lookup.loc[ids, "PC2"], errors="coerce")
        valid1, valid2 = values.notna() & pc1.notna(), values.notna() & pc2.notna()
        rho1 = spearmanr(values[valid1], pc1[valid1]).statistic if valid1.sum() >= 3 else np.nan
        rho2 = spearmanr(values[valid2], pc2[valid2]).statistic if valid2.sum() >= 3 else np.nan
        rows.append({"measurement": measurement, "measurement_class": "sparse",
                     "vector_PC1": rho1, "vector_PC2": rho2,
                     "samples_measured": int(values.notna().sum()),
                     "vector_method": "Passive Spearman projection"})
    return pd.DataFrame(rows)


def plot_biplot(background: pd.DataFrame, centroids: pd.DataFrame, modules: pd.DataFrame,
                 vectors: pd.DataFrame, explained: pd.DataFrame, output: Path,
                 formats: list[str]) -> None:
    hybrid_cols = [c for c in background if re.fullmatch(r"hyb_C[0-3]_G\d+", str(c))]
    background = background.copy()
    background["hard_hybrid"] = background[hybrid_cols].apply(pd.to_numeric, errors="coerce").fillna(0).idxmax(axis=1)
    palette = hybrid_palette(hybrid_cols)
    module_levels = sorted(modules["ecological_module"].astype(str), key=module_order)
    module_colors = {m: mcolors.to_hex(plt.get_cmap("tab20")(i / 20)) for i, m in enumerate(module_levels)}
    fig = plt.figure(figsize=(15.2, 9.0))
    fig._aspire_compact_publication_typography = True
    ax = fig.add_axes([0.07, 0.11, 0.67, 0.82])
    legend_ax = fig.add_axes([0.76, 0.08, 0.23, 0.86]); legend_ax.axis("off")
    for compartment in sorted(background["hard_hybrid"].dropna().unique(), key=hybrid_sort):
        part = background.loc[background["hard_hybrid"].eq(compartment)]
        ax.scatter(part["PC1"], part["PC2"], s=17, color=lighten(palette[compartment]), edgecolor="none", zorder=1)
    ax.axhline(0, color="0.84", lw=0.8); ax.axvline(0, color="0.84", lw=0.8)
    centroid_table = centroids.copy()
    centroid_table["responsibility_column"] = centroid_table["compartment"].map(
        lambda value: re.sub(r"hybrid_c([0-3])_g(\d+)", r"hyb_C\1_G\2", str(value)))
    observed_hybrids = set(background["hard_hybrid"].dropna().astype(str))
    centroid_table = centroid_table.loc[
        centroid_table["responsibility_column"].isin(observed_hybrids)
    ].sort_values(
        "responsibility_column", key=lambda s: s.map(hybrid_sort)).reset_index(drop=True)
    centroid_table["centroid_number"] = np.arange(1, len(centroid_table) + 1)
    xspan = float(background["PC1"].max() - background["PC1"].min())
    yspan = float(background["PC2"].max() - background["PC2"].min())
    # Module markers remain at their calculated abundance-weighted environmental
    # centroids.  Only text labels may be adjusted elsewhere; moving the marker
    # itself would make its PC-space position visually ambiguous.
    module_display = modules[["PC1_centroid", "PC2_centroid"]].to_numpy(float)
    module_xy = []
    for row, displayed in zip(modules.itertuples(index=False), module_display):
        module_xy.append(tuple(displayed))
        ax.scatter(displayed[0], displayed[1], s=205,
                   color=module_colors[str(row.ecological_module)], edgecolor="white", linewidth=0.8, zorder=10)
        ax.text(displayed[0], displayed[1], str(row.ecological_module), ha="center", va="center",
                fontsize=8, fontweight="bold", color="black", zorder=11,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 0.3})
    for row in centroid_table.itertuples(index=False):
        ax.scatter(row.PC1_centroid, row.PC2_centroid, s=250, facecolor="white", edgecolor="black", linewidth=1.4, zorder=20)
        ax.text(row.PC1_centroid, row.PC2_centroid, str(row.centroid_number), ha="center", va="center", fontsize=8.5, fontweight="bold", zorder=21)
    maximum_x = max(vectors["vector_PC1"].abs().max(), 1e-9)
    maximum_y = max(vectors["vector_PC2"].abs().max(), 1e-9)
    scale = min(0.23 * xspan / maximum_x, 0.23 * yspan / maximum_y)
    vector_labels = []
    vector_anchors = []
    for index, row in enumerate(vectors.dropna(subset=["vector_PC1", "vector_PC2"]).itertuples(index=False)):
        x, y = row.vector_PC1 * scale, row.vector_PC2 * scale
        linestyle = "--" if row.measurement_class == "sparse" else "-"
        arrow = mpl.patches.FancyArrowPatch((0, 0), (x, y), arrowstyle="-|>", mutation_scale=11,
                                             linewidth=1.25, linestyle=linestyle, color="0.18", zorder=5)
        arrow.set_path_effects([path_effects.Stroke(linewidth=2.4, foreground="white"), path_effects.Normal()])
        ax.add_patch(arrow)
        vector_anchors.append((x, y))
        vector_labels.append(ax.text(
            x * 1.08, y * 1.08, row.measurement, fontsize=8,
            ha="left" if x >= 0 else "right", va="bottom" if y >= 0 else "top", zorder=6,
            bbox={"facecolor": "white", "edgecolor": "0.8", "linewidth": 0.35, "alpha": 1.0, "pad": 1.2},
        ))
    xpad, ypad = max(0.35, xspan * 0.035), max(0.35, yspan * 0.035)
    ax.set_xlim(background["PC1"].min() - xpad, background["PC1"].max() + xpad)
    ax.set_ylim(background["PC2"].min() - ypad, background["PC2"].max() + ypad)
    forbidden = module_xy + [
        (row.PC1_centroid, row.PC2_centroid)
        for row in centroid_table.itertuples(index=False)
    ]
    if vector_labels:
        adjust_text(
            vector_labels,
            x=[point[0] for point in forbidden],
            y=[point[1] for point in forbidden],
            ax=ax,
            expand=(1.08, 1.18),
            force_text=(0.35, 0.5),
            force_static=(0.5, 0.7),
            ensure_inside_axes=True,
            arrowprops={"arrowstyle": "-", "color": "0.55", "lw": 0.45},
        )
    ratios = dict(zip(explained["PC"].astype(str), pd.to_numeric(explained["explained_variance_ratio"], errors="coerce")))
    ax.set_xlabel(f"PC1 ({100 * ratios.get('PC1', np.nan):.1f}%)")
    ax.set_ylabel(f"PC2 ({100 * ratios.get('PC2', np.nan):.1f}%)")
    ax.set_title("Ecological modules and measured physicochemistry across hybrid environmental space")
    ax.set_box_aspect(1.0)
    module_handles = [Line2D([], [], marker="o", linestyle="", markerfacecolor=module_colors[m],
                             markeredgecolor="white", label=m) for m in module_levels]
    first = legend_ax.legend(handles=module_handles, title="Ecological module", loc="upper left",
                             frameon=False, ncol=2, fontsize=7, title_fontsize=8.5)
    legend_ax.add_artist(first)
    centroid_handles = [Line2D([], [], marker="", linestyle="", label=f"{r.centroid_number}. {r.display_label}")
                        for r in centroid_table.itertuples(index=False)]
    legend_ax.legend(handles=centroid_handles, title=r"Hybrid O$_2$-GMM centroid", loc="lower left",
                     frameon=False, fontsize=6.5, title_fontsize=8.5, handlelength=0)
    for fmt in formats:
        fig.savefig(output.with_suffix(f".{fmt}"), dpi=300, bbox_inches=None)
    plt.close(fig)


def plot_environmental_module_network(
    nodes: pd.DataFrame,
    edges: pd.DataFrame,
    modules: pd.DataFrame,
    hybrid_centroids: pd.DataFrame,
    background: pd.DataFrame,
    output: Path,
    formats: list[str],
    metric_column: str | None = None,
    metric_label: str | None = None,
    top_n: int = 10,
) -> None:
    module_levels = sorted(modules["ecological_module"].astype(str), key=module_order)
    module_colors = {
        module: mcolors.to_hex(plt.get_cmap("tab20")(index / 20))
        for index, module in enumerate(module_levels)
    }
    metric_mode = metric_column is not None
    if metric_mode:
        metric_values = pd.to_numeric(nodes[metric_column], errors="coerce").fillna(0.0)
        metric_percentiles = metric_values.rank(method="average", pct=True)
        node_sizes = 13.0 + 92.0 * np.power(metric_percentiles, 1.35)
        node_size_lookup = dict(zip(nodes["ASV_ID"], node_sizes))
        top_metric_nodes = nodes.assign(_metric_value=metric_values).sort_values(
            ["_metric_value", "ASV_ID"], ascending=[False, True]
        ).head(top_n)
        top_metric_ids = set(top_metric_nodes["ASV_ID"])
    else:
        metric_values = pd.Series(dtype=float)
        node_size_lookup = {}
        top_metric_nodes = nodes.iloc[0:0].copy()
        top_metric_ids = set()
    node_xy = nodes.set_index("ASV_ID")[["network_x_display", "network_y_display"]]
    fig = plt.figure(figsize=(18.5, 11.2))
    fig._aspire_compact_publication_typography = True
    ax = fig.add_axes([0.02, 0.045, 0.67, 0.90])
    legend_ax = fig.add_axes([0.71, 0.055, 0.285, 0.88])
    legend_ax.axis("off")

    if "analysis_used" in background:
        used_mask = background["analysis_used"].fillna(False).astype(bool)
    else:
        used_mask = pd.Series(True, index=background.index)
    unused_background = background.loc[~used_mask]
    used_background = background.loc[used_mask]
    ax.scatter(
        pd.to_numeric(unused_background["PC1"], errors="coerce"),
        pd.to_numeric(unused_background["PC2"], errors="coerce"),
        s=11, facecolor="0.78", edgecolor="none", alpha=0.48, zorder=-3,
    )
    ax.scatter(
        pd.to_numeric(used_background["PC1"], errors="coerce"),
        pd.to_numeric(used_background["PC2"], errors="coerce"),
        s=20, facecolor="0.68", edgecolor="none", alpha=0.62, zorder=-2,
    )

    for row in modules.itertuples(index=False):
        center = np.array([row.network_center_x_display, row.network_center_y_display])
        original = np.array([row.PC1_centroid_unmodified, row.PC2_centroid_unmodified])
        displacement = float(np.linalg.norm(center - original))
        if displacement > 1e-12:
            direction_to_centroid = (original - center) / displacement
            island_boundary = center + direction_to_centroid * row.island_radius * 0.94
            ax.plot(
                [original[0], island_boundary[0]],
                [original[1], island_boundary[1]],
                color="0.52", linewidth=0.85, zorder=0,
            )
        ax.add_patch(Circle(
            center,
            row.island_radius * 0.94,
            facecolor=mcolors.to_rgba(module_colors[row.ecological_module], 0.055),
            edgecolor=mcolors.to_rgba(module_colors[row.ecological_module], 0.72),
            linewidth=0.9,
            zorder=1,
        ))

    displayed_edges = edges.loc[edges["displayed"]].copy()
    finite_weights = pd.to_numeric(displayed_edges["weight"], errors="coerce").abs()
    weight_scale = float(finite_weights.quantile(0.95)) if finite_weights.notna().any() else 1.0
    weight_scale = max(weight_scale, 1e-9)
    for row in displayed_edges.itertuples(index=False):
        if row.source not in node_xy.index or row.target not in node_xy.index:
            continue
        source, target = node_xy.loc[row.source], node_xy.loc[row.target]
        width = 0.28
        if pd.notna(row.weight):
            width += 0.72 * min(abs(float(row.weight)) / weight_scale, 1.0)
        linestyle = "--" if row.sign == "negative" else "-"
        ax.plot(
            [source.iloc[0], target.iloc[0]],
            [source.iloc[1], target.iloc[1]],
            color="0.38", alpha=0.36, linewidth=width,
            linestyle=linestyle, zorder=2,
        )

    nonanchors = nodes.loc[~nodes["is_ecological_anchor"]]
    anchors = nodes.loc[nodes["is_ecological_anchor"]]
    for module, frame in nonanchors.groupby("ecological_module", sort=False):
        ax.scatter(
            frame["network_x_display"], frame["network_y_display"],
            s=(frame["ASV_ID"].map(node_size_lookup).to_numpy(float) if metric_mode else 27),
            facecolor=module_colors[str(module)], edgecolor="white",
            linewidth=0.38, alpha=0.96, zorder=4,
        )
    for module, frame in anchors.groupby("ecological_module", sort=False):
        ax.scatter(
            frame["network_x_display"], frame["network_y_display"],
            s=(frame["ASV_ID"].map(node_size_lookup).clip(lower=35).to_numpy(float)
               if metric_mode else 92),
            facecolor=module_colors[str(module)], edgecolor="black",
            linewidth=1.0, alpha=1.0, zorder=6,
        )

    if metric_mode and top_metric_ids:
        top_display = nodes.loc[nodes["ASV_ID"].isin(top_metric_ids)]
        ax.scatter(
            top_display["network_x_display"], top_display["network_y_display"],
            s=top_display["ASV_ID"].map(node_size_lookup).to_numpy(float) + 38,
            facecolor="none", edgecolor="black", linewidth=1.15, zorder=7,
        )

    module_labels = []
    for row in modules.itertuples(index=False):
        module_labels.append(ax.text(
            row.network_center_x_display,
            row.network_center_y_display + row.island_radius + 0.07,
            row.ecological_module,
            ha="center", va="bottom", fontsize=10.0, fontweight="bold",
            color="black", zorder=8,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.94, "pad": 0.8},
        ))

    label_nodes = top_metric_nodes if metric_mode else anchors
    anchor_labels = []
    for row in label_nodes.itertuples(index=False):
        anchor_labels.append(ax.text(
            row.network_x_display,
            row.network_y_display,
            row.ASV_ID,
            fontsize=6.4, ha="center", va="center", zorder=9,
            bbox={"facecolor": "white", "edgecolor": "0.78", "linewidth": 0.25,
                  "alpha": 0.94, "pad": 0.55},
        ))
    if anchor_labels:
        adjust_text(
            anchor_labels,
            x=nodes["network_x_display"].to_numpy(float),
            y=nodes["network_y_display"].to_numpy(float),
            ax=ax,
            expand=(1.04, 1.08),
            force_text=(0.18, 0.24),
            force_static=(0.10, 0.14),
            ensure_inside_axes=True,
            arrowprops={"arrowstyle": "-", "color": "0.55", "lw": 0.38},
        )

    # The colored, black-bordered marker is the exact abundance-weighted
    # environmental centroid.  Its associated network island is display-only.
    for row in modules.itertuples(index=False):
        module = str(row.ecological_module)
        ax.scatter(
            row.PC1_centroid_unmodified,
            row.PC2_centroid_unmodified,
            s=190,
            facecolor=module_colors[module],
            edgecolor="black",
            linewidth=1.35,
            zorder=10,
        )
        centroid_label = ax.text(
            row.PC1_centroid_unmodified,
            row.PC2_centroid_unmodified,
            module,
            ha="center",
            va="center",
            fontsize=7.2,
            fontweight="bold",
            color="black",
            zorder=11,
        )
        centroid_label.set_path_effects([
            path_effects.Stroke(linewidth=1.8, foreground="white"),
            path_effects.Normal(),
        ])

    for row in hybrid_centroids.itertuples(index=False):
        ax.scatter(
            row.PC1_centroid, row.PC2_centroid,
            s=180, facecolor="white", edgecolor="black", linewidth=1.25, zorder=12,
        )
        ax.text(
            row.PC1_centroid, row.PC2_centroid, str(row.centroid_number),
            ha="center", va="center", fontsize=7.5, fontweight="bold", zorder=13,
        )

    all_x = pd.concat([
        nodes["network_x_display"],
        hybrid_centroids["PC1_centroid"],
        pd.to_numeric(background["PC1"], errors="coerce"),
        modules["network_center_x_display"] - modules["island_radius"],
        modules["network_center_x_display"] + modules["island_radius"],
    ]).astype(float)
    all_y = pd.concat([
        nodes["network_y_display"],
        hybrid_centroids["PC2_centroid"],
        pd.to_numeric(background["PC2"], errors="coerce"),
        modules["network_center_y_display"] - modules["island_radius"],
        modules["network_center_y_display"] + modules["island_radius"],
    ]).astype(float)
    xpad = max(0.45, (all_x.max() - all_x.min()) * 0.055)
    ypad = max(0.45, (all_y.max() - all_y.min()) * 0.065)
    ax.set_xlim(all_x.min() - xpad, all_x.max() + xpad)
    ax.set_ylim(all_y.min() - ypad, all_y.max() + ypad)
    ax.set_aspect("equal", adjustable="box")
    ax.set_axis_off()
    title = "Ecological-module topology across hybrid environmental space"
    if metric_mode:
        title = f"ASV {metric_label} across hybrid environmental space"
    ax.set_title(title, pad=10)

    symbol_handles = [
        Line2D([], [], marker="o", linestyle="", markersize=3.3,
               markerfacecolor="0.78", markeredgecolor="none",
               label="Not used in ASV-module analysis"),
        Line2D([], [], marker="o", linestyle="", markersize=4.6,
               markerfacecolor="0.68", markeredgecolor="none",
               label="Used in ASV-module analysis"),
        Line2D([], [], marker="o", linestyle="", markersize=4.2,
               markerfacecolor="0.55", markeredgecolor="white", label="ASV"),
        Line2D([], [], marker="o", linestyle="", markersize=7.0,
               markerfacecolor="0.55", markeredgecolor="black", label="Module anchor ASV"),
        Line2D([], [], marker="o", linestyle="", markersize=8.0,
               markerfacecolor="0.55", markeredgecolor="black",
               markeredgewidth=1.1, label="Exact module centroid"),
        Line2D([], [], color="0.42", linewidth=0.9, label="Within-module association"),
        Line2D([], [], marker="o", linestyle="", markersize=7.0,
               markerfacecolor="white", markeredgecolor="black", label="Hybrid centroid"),
    ]
    if metric_mode:
        symbol_handles.insert(2, Line2D(
            [], [], marker="o", linestyle="", markersize=8.0,
            markerfacecolor="none", markeredgecolor="black", markeredgewidth=1.15,
            label=f"Top {top_n} {metric_label}",
        ))
    first = legend_ax.legend(
        handles=symbol_handles, title="Network element", loc="upper left",
        bbox_to_anchor=(0.0, 1.0),
        frameon=False, fontsize=7, title_fontsize=8.5,
    )
    legend_ax.add_artist(first)
    centroid_handles = [
        Line2D([], [], marker="", linestyle="", label=f"{row.centroid_number}. {row.display_label}")
        for row in hybrid_centroids.itertuples(index=False)
    ]
    second = legend_ax.legend(
        handles=centroid_handles, title=r"Hybrid O$_2$-GMM centroid", loc="upper left",
        bbox_to_anchor=(0.0, 0.84),
        frameon=False, fontsize=6.5, title_fontsize=8.5, handlelength=0,
    )
    legend_ax.add_artist(second)
    anchor_handles = []
    if metric_mode:
        for rank, row in enumerate(top_metric_nodes.itertuples(index=False), start=1):
            anchor_handles.append(Line2D(
                [], [], marker="o", linestyle="", markersize=3.7,
                markerfacecolor=module_colors[str(row.ecological_module)],
                markeredgecolor="black", markeredgewidth=0.35,
                label=(f"{rank}. {anchor_taxonomy_label(row)} "
                       f"[{row.ecological_module}; {getattr(row, metric_column):.4g}]"),
            ))
        hierarchy_title = f"Top {metric_label} ASVs and taxonomy"
    else:
        anchor_order = anchors.assign(
            _module_order=anchors["ecological_module"].map(module_order)
        ).sort_values(["_module_order", "ASV_ID"])
        for module in module_levels:
            module_rows = anchor_order.loc[
                anchor_order["ecological_module"].astype(str) == str(module)
            ]
            if module_rows.empty:
                continue
            # A module heading followed by indented anchor rows gives the legend a
            # compact tree/drop-down structure and avoids a separate module key.
            anchor_handles.append(Line2D(
                [], [], marker="o", linestyle="", markersize=5.2,
                markerfacecolor=module_colors[str(module)],
                markeredgecolor="black", markeredgewidth=0.5,
                label=str(module),
            ))
            for row in module_rows.itertuples(index=False):
                anchor_handles.append(Line2D(
                    [], [], marker="o", linestyle="", markersize=3.7,
                    markerfacecolor=module_colors[str(row.ecological_module)],
                    markeredgecolor="black", markeredgewidth=0.35,
                    label=f"    - {anchor_taxonomy_label(row)}",
                ))
        hierarchy_title = "Ecological modules and anchor taxonomy"
    hierarchy_legend = legend_ax.legend(
        handles=anchor_handles,
        title=hierarchy_title,
        loc="lower left",
        bbox_to_anchor=(0.0, -0.01),
        frameon=False,
        fontsize=4.35,
        title_fontsize=7.0,
        labelspacing=0.24,
        handletextpad=0.45,
        borderaxespad=0.0,
    )
    if not metric_mode:
        module_names = set(map(str, module_levels))
        for text in hierarchy_legend.get_texts():
            if text.get_text() in module_names:
                text.set_fontweight("bold")
                text.set_fontsize(5.2)
    for fmt in formats:
        fig.savefig(output.with_suffix(f".{fmt}"), dpi=300, bbox_inches=None)
    plt.close(fig)


def locked_spieceasi_network_layout(
    nodes: pd.DataFrame,
    layout_path: Path,
) -> pd.DataFrame:
    """Join the exact coordinates emitted by the network-analysis renderer."""
    layout = read_table(layout_path)
    required = {"ASV_ID", "spieceasi_x", "spieceasi_y"}
    missing_columns = sorted(required.difference(layout.columns))
    if missing_columns:
        raise ValueError(
            "SPIEC-EASI layout is missing required columns: "
            + ", ".join(missing_columns)
        )
    layout = layout.copy()
    layout["ASV_ID"] = layout["ASV_ID"].map(norm_asv)
    if layout["ASV_ID"].duplicated().any():
        raise ValueError("SPIEC-EASI layout contains duplicate ASV identifiers")
    coordinates = layout[["ASV_ID", "spieceasi_x", "spieceasi_y"]].copy()
    coordinates["spieceasi_x"] = pd.to_numeric(
        coordinates["spieceasi_x"], errors="coerce"
    )
    coordinates["spieceasi_y"] = pd.to_numeric(
        coordinates["spieceasi_y"], errors="coerce"
    )
    out = nodes.merge(coordinates, on="ASV_ID", how="left", validate="one_to_one")
    missing_nodes = out.loc[
        out[["spieceasi_x", "spieceasi_y"]].isna().any(axis=1), "ASV_ID"
    ].tolist()
    if missing_nodes:
        raise ValueError(
            "Network-analysis layout lacks coordinates for ASVs: "
            + ", ".join(missing_nodes[:10])
        )
    return out


def plot_spieceasi_metric_network(
    nodes: pd.DataFrame,
    edges: pd.DataFrame,
    output: Path,
    formats: list[str],
    metric_column: str,
    metric_label: str,
    top_n: int,
) -> None:
    """Plot a metric on the complete SPIEC-EASI topology without PCA islands."""
    module_levels = sorted(nodes["ecological_module"].astype(str).unique(), key=module_order)
    module_colors = {
        module: mcolors.to_hex(plt.get_cmap("tab20")(index / 20))
        for index, module in enumerate(module_levels)
    }
    metric_values = pd.to_numeric(nodes[metric_column], errors="coerce").fillna(0.0)
    percentiles = metric_values.rank(method="average", pct=True)
    node_sizes = 17.0 + 115.0 * np.power(percentiles, 1.35)
    node_size_lookup = dict(zip(nodes["ASV_ID"], node_sizes))
    top_nodes = nodes.assign(metric_plot_value=metric_values).sort_values(
        ["metric_plot_value", "ASV_ID"], ascending=[False, True], na_position="last"
    ).head(min(top_n, len(nodes))).copy()
    top_nodes["metric_rank"] = np.arange(1, len(top_nodes) + 1)
    top_ids = set(top_nodes["ASV_ID"])
    coordinates = nodes.set_index("ASV_ID")[["spieceasi_x", "spieceasi_y"]]

    fig = plt.figure(figsize=(18.5, 11.2))
    fig._aspire_compact_publication_typography = True
    ax = fig.add_axes([0.025, 0.045, 0.69, 0.90])
    legend_ax = fig.add_axes([0.735, 0.055, 0.26, 0.88])
    legend_ax.axis("off")

    finite_weights = pd.to_numeric(edges["weight"], errors="coerce").abs()
    weight_scale = float(finite_weights.quantile(0.95)) if finite_weights.notna().any() else 1.0
    weight_scale = max(weight_scale, 1e-9)
    for row in edges.itertuples(index=False):
        if row.source not in coordinates.index or row.target not in coordinates.index:
            continue
        source, target = coordinates.loc[row.source], coordinates.loc[row.target]
        width = 0.22
        if pd.notna(row.weight):
            width += 0.68 * min(abs(float(row.weight)) / weight_scale, 1.0)
        ax.plot(
            [source.iloc[0], target.iloc[0]],
            [source.iloc[1], target.iloc[1]],
            color=("0.32" if row.sign == "negative" else "0.68"),
            alpha=(0.42 if row.sign == "negative" else 0.25),
            linewidth=width,
            linestyle=("--" if row.sign == "negative" else "-"),
            zorder=1,
        )

    anchor_mask = nodes["is_ecological_anchor"].fillna(False).astype(bool)
    ordinary_nodes = nodes.loc[~anchor_mask]
    ecological_anchors = nodes.loc[anchor_mask]

    # Draw ordinary nodes first. Ecological anchors and metric-highlighted nodes
    # are redrawn afterward so nearby points cannot obscure them.
    for module, frame in ordinary_nodes.groupby("ecological_module", sort=False):
        ax.scatter(
            frame["spieceasi_x"], frame["spieceasi_y"],
            s=frame["ASV_ID"].map(node_size_lookup).to_numpy(float),
            facecolor=module_colors[str(module)], edgecolor="white",
            linewidth=0.32, alpha=0.96, zorder=3,
        )
    for module, frame in ecological_anchors.groupby("ecological_module", sort=False):
        ax.scatter(
            frame["spieceasi_x"], frame["spieceasi_y"],
            s=frame["ASV_ID"].map(node_size_lookup).clip(lower=42).to_numpy(float),
            facecolor=module_colors[str(module)], edgecolor="black",
            linewidth=0.85, alpha=1.0, zorder=4,
        )
    highlighted = nodes.loc[nodes["ASV_ID"].isin(top_ids)]
    for module, frame in highlighted.groupby("ecological_module", sort=False):
        ax.scatter(
            frame["spieceasi_x"], frame["spieceasi_y"],
            s=frame["ASV_ID"].map(node_size_lookup).to_numpy(float),
            facecolor=module_colors[str(module)], edgecolor="white",
            linewidth=0.45, alpha=1.0, zorder=5,
        )
    ax.scatter(
        highlighted["spieceasi_x"], highlighted["spieceasi_y"],
        s=highlighted["ASV_ID"].map(node_size_lookup).to_numpy(float) + 68,
        facecolor="none", edgecolor="black", linewidth=1.45, zorder=6,
    )
    rank_lookup = top_nodes.set_index("ASV_ID")["metric_rank"].to_dict()
    rank_texts = []
    for row in highlighted.itertuples(index=False):
        rank = rank_lookup[row.ASV_ID]
        angle = rank * 2.399963229728653
        offset = 0.025 + 0.008 * (rank % 3)
        label = ax.text(
            row.spieceasi_x + offset * math.cos(angle),
            row.spieceasi_y + offset * math.sin(angle),
            str(rank),
            ha="center", va="center", fontsize=8.0, fontweight="bold",
            color="black", zorder=7,
            bbox={"facecolor": "white", "edgecolor": "0.35", "linewidth": 0.35,
                  "alpha": 0.96, "pad": 0.55},
        )
        rank_texts.append(label)
    x_range = float(nodes["spieceasi_x"].max() - nodes["spieceasi_x"].min())
    y_range = float(nodes["spieceasi_y"].max() - nodes["spieceasi_y"].min())
    ax.set_xlim(nodes["spieceasi_x"].min() - max(0.04, x_range * 0.06),
                nodes["spieceasi_x"].max() + max(0.04, x_range * 0.06))
    ax.set_ylim(nodes["spieceasi_y"].min() - max(0.04, y_range * 0.06),
                nodes["spieceasi_y"].max() + max(0.04, y_range * 0.06))
    if rank_texts:
        adjust_text(
            rank_texts,
            x=nodes["spieceasi_x"].to_numpy(float),
            y=nodes["spieceasi_y"].to_numpy(float),
            target_x=highlighted["spieceasi_x"].to_numpy(float),
            target_y=highlighted["spieceasi_y"].to_numpy(float),
            ax=ax,
            expand=(1.45, 1.60),
            force_text=(1.0, 1.2),
            force_static=(0.35, 0.42),
            force_pull=(0.004, 0.004),
            force_explode=(1.2, 1.5),
            explode_radius=40,
            max_move=35,
            iter_lim=1000,
            ensure_inside_axes=True,
            arrowprops={"arrowstyle": "-", "color": "0.45", "lw": 0.42},
        )
    ax.set_aspect("equal", adjustable="box")
    ax.set_axis_off()
    ax.set_title(f"SPIEC-EASI network: ASV {metric_label}", pad=10, fontsize=16)

    network_handles = [
        Line2D([], [], color="0.68", linewidth=0.9, label="Positive association"),
        Line2D([], [], color="0.32", linewidth=0.9, linestyle="--",
               label="Negative association"),
        Line2D([], [], marker="o", linestyle="", markersize=7.0,
               markerfacecolor="none", markeredgecolor="black",
               markeredgewidth=1.2, label=f"Top {top_n} {metric_label}"),
    ]
    first = legend_ax.legend(
        handles=network_handles, title="Network element", loc="upper left",
        frameon=False, fontsize=8.2, title_fontsize=9.5,
    )
    legend_ax.add_artist(first)
    module_handles = [
        Line2D([], [], marker="o", linestyle="", markersize=5.2,
               markerfacecolor=module_colors[module], markeredgecolor="white", label=module)
        for module in module_levels
    ]
    second = legend_ax.legend(
        handles=module_handles, title="Ecological module", loc="upper left",
        bbox_to_anchor=(0.0, 0.82), frameon=False, ncol=2,
        fontsize=7.5, title_fontsize=9.5, columnspacing=0.9,
    )
    legend_ax.add_artist(second)
    taxonomy_handles = []
    for row in top_nodes.itertuples(index=False):
        taxonomy_handles.append(Line2D(
            [], [], marker="o", linestyle="", markersize=4.2,
            markerfacecolor=module_colors[str(row.ecological_module)],
            markeredgecolor="black", markeredgewidth=0.45,
            label=(f"{row.metric_rank}. {anchor_taxonomy_label(row)} "
                   f"[{row.ecological_module}; {row.metric_plot_value:.4g}]"),
        ))
    legend_ax.legend(
        handles=taxonomy_handles,
        title=f"Top {metric_label} ASVs and taxonomy",
        loc="lower left", bbox_to_anchor=(0.0, -0.01), frameon=False,
        fontsize=6.4, title_fontsize=8.8, labelspacing=0.43,
        handletextpad=0.5, borderaxespad=0.0,
    )
    for fmt in formats:
        fig.savefig(output.with_suffix(f".{fmt}"), dpi=300, bbox_inches=None)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--asv-counts", type=Path, required=True)
    parser.add_argument("--prevalence-counts", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--modules", type=Path, required=True)
    parser.add_argument("--network-graph", type=Path, required=True)
    parser.add_argument("--spieceasi-layout", type=Path, required=True)
    parser.add_argument("--node-features", type=Path, required=True)
    parser.add_argument("--taxonomy", type=Path, required=True)
    parser.add_argument("--anchor-top-n", type=int, default=1)
    parser.add_argument("--network-metric-top-n", type=int, default=10)
    parser.add_argument("--measurement-matrix", type=Path, required=True)
    parser.add_argument("--asv-correlations", type=Path, required=True)
    parser.add_argument("--measurement-selection-audit", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--sample-col", default="sampleID")
    parser.add_argument("--cruise-col", default="Cruise")
    parser.add_argument("--depth-col", default="Depth")
    parser.add_argument("--sparse-measurements", default="")
    parser.add_argument("--q-threshold", type=float, default=0.05)
    parser.add_argument("--pca-scores", type=Path)
    parser.add_argument("--pca-loadings", type=Path)
    parser.add_argument("--pca-explained", type=Path)
    parser.add_argument("--hybrid-assignments", type=Path)
    parser.add_argument("--hybrid-centroids", type=Path)
    parser.add_argument("--formats", default="pdf,png,svg")
    args = parser.parse_args()
    if args.anchor_top_n < 1:
        raise ValueError("--anchor-top-n must be at least 1")
    if args.network_metric_top_n < 1:
        raise ValueError("--network-metric-top-n must be at least 1")
    tables, plots = args.outdir / "tables", args.outdir / "plots"
    tables.mkdir(parents=True, exist_ok=True); plots.mkdir(parents=True, exist_ok=True)
    counts = load_counts(args.asv_counts)
    prevalence_counts = load_counts(args.prevalence_counts)
    metadata = read_table(args.metadata)
    metadata[args.sample_col] = metadata[args.sample_col].astype(str)
    modules = module_membership(read_table(args.modules))
    measurements = read_table(args.measurement_matrix).set_index(args.sample_col)
    measurements.index = measurements.index.astype(str)
    measurements = measurements.apply(pd.to_numeric, errors="coerce")
    sparse = set(csv_list(args.sparse_measurements))
    scores = module_scores(counts, modules)
    scores.to_csv(tables / "sample_ecological_module_scores_measurement_cohort.tsv", sep="\t", index=False)
    coverage = coverage_table(measurements, metadata, args.sample_col, args.cruise_col, args.depth_col, sparse)
    coverage = add_unavailable_measurements_to_coverage(
        coverage, read_table(args.measurement_selection_audit), sparse
    )
    coverage.to_csv(tables / "module_measurement_coverage.tsv", sep="\t", index=False)
    direct = module_measurement_tests(scores, measurements, metadata, args.sample_col,
                                      args.cruise_col, args.depth_col, sparse)
    direct.to_csv(tables / "module_measurement_association.tsv", sep="\t", index=False)
    support = member_support(read_table(args.asv_correlations), modules, args.q_threshold, sparse)
    support.to_csv(tables / "module_measurement_member_support.tsv", sep="\t", index=False)
    plot_heatmap(direct, plots / "ecological_module_measurement_association_heatmap", csv_list(args.formats))
    pca_inputs = [args.pca_scores, args.pca_loadings, args.pca_explained,
                  args.hybrid_assignments, args.hybrid_centroids]
    if any(pca_inputs) and not all(pca_inputs):
        raise ValueError("The module measurement biplot requires all five PCA/hybrid inputs")
    if all(pca_inputs):
        pca_scores = read_table(args.pca_scores)
        assignments = read_table(args.hybrid_assignments)
        if "cruise_year_month_depth" in pca_scores and "cruise_year_month_depth" in assignments:
            background = assignments.merge(pca_scores[["cruise_year_month_depth", "PC1", "PC2"]],
                                           on="cruise_year_month_depth", how="inner", validate="one_to_one")
        else:
            background = assignments.merge(pca_scores[["Cruise", "Depth_anchored", "PC1", "PC2"]],
                                           on=["Cruise", "Depth_anchored"], how="inner")
        sample_coords = pca_sample_coordinates(metadata, pca_scores, args.sample_col, args.depth_col)
        background["_cruise"] = pd.to_numeric(
            background["Cruise"] if "Cruise" in background else background["Cruise_o2"],
            errors="coerce",
        )
        background_depth_col = (
            "Depth_anchored" if "Depth_anchored" in background
            else "Depth_anchored_o2" if "Depth_anchored_o2" in background
            else "Depth_o2"
        )
        background["_depth"] = pd.to_numeric(
            background[background_depth_col], errors="coerce"
        ).round(6)
        used_sample_ids = set(scores["sampleID"].astype(str))
        matched_keys = sample_coords.loc[
            sample_coords[args.sample_col].astype(str).isin(used_sample_ids),
            ["_cruise", "_depth"],
        ].drop_duplicates()
        matched_keys["analysis_used"] = True
        background_all = background.merge(
            matched_keys, on=["_cruise", "_depth"], how="left"
        )
        background_all["analysis_used"] = (
            background_all["analysis_used"].fillna(False).astype(bool)
        )
        background_all.to_csv(
            tables / "ecological_module_environmental_network_background.tsv",
            sep="\t", index=False,
        )
        background = background_all.loc[background_all["analysis_used"]].copy()
        positions = module_pca_positions(scores, sample_coords, args.sample_col)
        positions.to_csv(tables / "module_environmental_biplot_coordinates.tsv", sep="\t", index=False)
        vectors = vector_table(read_table(args.pca_loadings), measurements, sample_coords, args.sample_col, sparse)
        vectors.to_csv(tables / "measurement_vector_coordinates.tsv", sep="\t", index=False)
        plot_biplot(background, read_table(args.hybrid_centroids), positions, vectors,
                     read_table(args.pca_explained),
                     plots / "ecological_module_measurement_biplot", csv_list(args.formats))
        network_nodes, network_edges, network_modules, network_hybrid_centroids = environmental_module_network(
            load_asv_graph(args.network_graph),
            read_table(args.node_features),
            modules,
            load_anchor_taxonomy(args.taxonomy),
            prevalence_counts,
            positions,
            read_table(args.hybrid_centroids),
            background,
            args.anchor_top_n,
        )
        network_nodes.to_csv(
            tables / "ecological_module_environmental_network_nodes.tsv", sep="\t", index=False
        )
        top_network_players(network_nodes, args.network_metric_top_n).to_csv(
            tables / "ecological_module_environmental_network_top_players.tsv",
            sep="\t", index=False,
        )
        network_edges.to_csv(
            tables / "ecological_module_environmental_network_edges.tsv", sep="\t", index=False
        )
        network_modules.to_csv(
            tables / "ecological_module_environmental_network_layout.tsv", sep="\t", index=False
        )
        plot_environmental_module_network(
            network_nodes,
            network_edges,
            network_modules,
            network_hybrid_centroids,
            background_all,
            plots / "ecological_module_environmental_network",
            csv_list(args.formats),
        )
        metric_network_nodes = locked_spieceasi_network_layout(
            network_nodes, args.spieceasi_layout
        )
        metric_network_nodes.to_csv(
            tables / "ecological_module_spieceasi_network_layout.tsv",
            sep="\t", index=False,
        )
        for suffix, (column, label) in NETWORK_METRIC_SPECS.items():
            plot_spieceasi_metric_network(
                metric_network_nodes,
                network_edges,
                plots / f"ecological_module_environmental_network_{suffix}",
                csv_list(args.formats),
                metric_column=column,
                metric_label=label,
                top_n=args.network_metric_top_n,
            )
    (args.outdir / "module_measurement_run_config.json").write_text(
        json.dumps(vars(args), indent=2, default=str)
    )


if __name__ == "__main__":
    main()
