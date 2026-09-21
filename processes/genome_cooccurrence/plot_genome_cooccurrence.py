#!/usr/bin/env python3
"""Render native and linked-ASV overlays on genome proportionality graphs."""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path
from typing import Iterable

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns

_STYLE_ROOT = str(Path(__file__).resolve().parents[1])
if _STYLE_ROOT not in sys.path:
    sys.path.insert(0, _STYLE_ROOT)
from shared_plot_style import install_publication_style

install_publication_style()
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["svg.fonttype"] = "none"


def read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t" if path.suffix.lower() in {".tsv", ".tab"} else ",", low_memory=False)


def key(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(value).strip().lower()).strip("_")


def truthy(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"true", "t", "1", "yes"})


def natural(value: object) -> tuple:
    return tuple(int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", str(value)) if part)


def contrasting_text_color(color: object) -> str:
    """Return black or white text with adequate contrast against a node fill."""
    try:
        red, green, blue = mcolors.to_rgb(color)
    except ValueError:
        return "black"
    luminance = 0.2126 * red + 0.7152 * green + 0.0722 * blue
    return "white" if luminance < 0.48 else "black"


def palette(labels: Iterable[object]) -> dict[str, str]:
    values = sorted({str(value).strip() for value in labels if str(value).strip()}, key=natural)
    known = [value for value in values if value not in {
        "Unclassified", "No linked ASV", "No significant linked-ASV association",
        "No supported renewal-phase contrast", "No supported renewal/post-renewal contrast",
    }]
    colors = sns.color_palette("husl", max(1, len(known)))
    result = {value: mcolors.to_hex(color) for value, color in zip(known, colors)}
    for value in values:
        if value not in result:
            result[value] = "#D0D0D0"
    semantic = {
        "MC1": "#0072B2", "MC2": "#E69F00", "MC3": "#009E73", "MC4": "#CC79A7",
        "Renewal onset vs stagnation": "#0072B2",
        "Post-renewal vs stagnation": "#E69F00",
        "Renewal/post-renewal vs stagnation": "#56B4E9",
        "Renewal phase vs stagnation": "#CC79A7",
        "Multiple supported renewal-module associations": "#7A5195",
        "Higher during renewal/post-renewal": "#0072B2",
        "Higher during stagnation": "#E69F00",
        "Equal pooled abundance": "#7A5195",
        "Mixed directional module links": "#7A5195",
        "No supported renewal-phase contrast": "#D0D0D0",
        "No supported renewal/post-renewal contrast": "#D0D0D0",
    }
    result.update({label: color for label, color in semantic.items() if label in values})
    return result


def graph_node_to_genome(graph: nx.Graph) -> dict[str, str]:
    return {node: str(graph.nodes[node].get("name", node)) for node in graph.nodes()}


def component_aware_layout(graph: nx.Graph, seed: int) -> dict:
    """Keep inferred topology central and place topology-free isolates legibly."""
    connected = [node for node in graph.nodes() if graph.degree(node) > 0]
    isolates = sorted(
        [node for node in graph.nodes() if graph.degree(node) == 0],
        key=lambda node: natural(graph.nodes[node].get("display_label", node)),
    )
    positions: dict = {}
    if connected:
        core = graph.subgraph(connected)
        positions.update(nx.spring_layout(
            core, seed=seed, weight="weight", iterations=2500,
            k=max(0.34, 1.7 / math.sqrt(max(core.number_of_nodes(), 1))), scale=0.78,
        ))
    if isolates:
        left_n = int(math.ceil(len(isolates) / 2))
        sides = ((isolates[:left_n], -1.0), (isolates[left_n:], 1.0))
        for nodes, x_value in sides:
            if not nodes:
                continue
            y_values = np.linspace(0.92, -0.92, len(nodes)) if len(nodes) > 1 else np.asarray([0.0])
            for node, y_value in zip(nodes, y_values):
                positions[node] = np.asarray([x_value, float(y_value)])
    ordered = sorted(positions, key=lambda node: natural(graph.nodes[node].get("display_label", node)))
    if positions:
        values = np.vstack([positions[node] for node in ordered])
        values -= values.mean(axis=0)
        extent = float(np.abs(values).max())
        if extent > 0:
            values *= 0.90 / extent
        positions = {node: values[index] for index, node in enumerate(ordered)}
    # Signed layouts can compress nodes when many attractive and repulsive
    # constraints compete. Apply the deterministic display-only collision pass
    # after normalization so its minimum spacing is not subsequently erased.
    # Topology and all edge statistics remain unchanged.
    for _ in range(180):
        moved = False
        for left in range(len(ordered)):
            for right in range(left + 1, len(ordered)):
                a, b = ordered[left], ordered[right]
                delta = np.asarray(positions[b], dtype=float) - np.asarray(positions[a], dtype=float)
                distance = float(np.linalg.norm(delta))
                # Preserve enough display-space for a two-digit integer at the
                # publication font size plus visible fill on every side.
                minimum = 0.19
                if distance >= minimum:
                    continue
                if distance <= 1e-12:
                    angle = ((left + 1) * 104729 + (right + 1) * 13007) % 360
                    unit = np.asarray([math.cos(math.radians(angle)), math.sin(math.radians(angle))])
                else:
                    unit = delta / distance
                shift = unit * (minimum - distance) * 0.51
                positions[a] = np.asarray(positions[a], dtype=float) - shift
                positions[b] = np.asarray(positions[b], dtype=float) + shift
                moved = True
        if not moved:
            break
    return positions


def renewal_module_labels(path: Path | None, max_q: float = 0.05) -> tuple[dict[str, str], pd.DataFrame]:
    """Classify modules from existing adjusted renewal-phase tests."""
    if path is None or not path.is_file():
        return {}, pd.DataFrame()
    table = read_table(path)
    required = {"analysis", "module_label", "adjusted_pvalue"}
    if not required.issubset(table.columns):
        raise ValueError(f"renewal module-association table lacks columns: {sorted(required - set(table.columns))}")
    table["adjusted_pvalue"] = pd.to_numeric(table["adjusted_pvalue"], errors="coerce")
    significant = table.loc[table["adjusted_pvalue"].le(max_q)].copy()
    labels: dict[str, str] = {}
    for module, frame in significant.groupby("module_label"):
        analyses = set(frame["analysis"].astype(str))
        onset = "renewal_vs_baseline" in analyses
        post = "post_renewal_vs_baseline" in analyses
        active = "active_vs_baseline" in analyses
        omnibus = "three_phase_omnibus" in analyses
        if onset and post:
            label = "Renewal/post-renewal vs stagnation"
        elif onset:
            label = "Renewal onset vs stagnation"
        elif post:
            label = "Post-renewal vs stagnation"
        elif active:
            label = "Renewal/post-renewal vs stagnation"
        elif omnibus:
            label = "Renewal phase vs stagnation"
        else:
            continue
        labels[str(module)] = label
    return labels, significant


def active_renewal_module_labels(
    path: Path | None, profiles_path: Path | None, max_q: float = 0.05
) -> tuple[dict[str, str], pd.DataFrame]:
    """Return significant modules labeled by pooled active-phase direction."""
    if path is None or not path.is_file() or profiles_path is None or not profiles_path.is_file():
        return {}, pd.DataFrame()
    table = read_table(path)
    required = {"analysis", "module_label", "adjusted_pvalue"}
    if not required.issubset(table.columns):
        raise ValueError(f"renewal module-association table lacks columns: {sorted(required - set(table.columns))}")
    table["adjusted_pvalue"] = pd.to_numeric(table["adjusted_pvalue"], errors="coerce")
    active = table.loc[
        table["analysis"].astype(str).eq("active_vs_baseline")
        & table["adjusted_pvalue"].le(max_q)
    ].copy()
    profiles = read_table(profiles_path)
    profile_required = {
        "module_label", "renewal_phase", "samples_n", "mean_module_relative_abundance"
    }
    if not profile_required.issubset(profiles.columns):
        raise ValueError(
            "renewal module-profile table lacks columns: "
            f"{sorted(profile_required - set(profiles.columns))}"
        )
    profiles = profiles.copy()
    profiles["samples_n"] = pd.to_numeric(profiles["samples_n"], errors="coerce")
    profiles["mean_module_relative_abundance"] = pd.to_numeric(
        profiles["mean_module_relative_abundance"], errors="coerce"
    )
    rows = []
    labels: dict[str, str] = {}
    for module in active["module_label"].dropna().astype(str).unique():
        frame = profiles.loc[profiles["module_label"].astype(str).eq(module)].copy()
        baseline = frame.loc[frame["renewal_phase"].astype(str).isin(["baseline", "stagnation"])]
        active_phases = frame.loc[frame["renewal_phase"].astype(str).isin(["renewal", "post-renewal"])]
        if baseline.empty or active_phases.empty:
            continue
        baseline_mean = float(np.average(
            baseline["mean_module_relative_abundance"], weights=baseline["samples_n"]
        ))
        active_mean = float(np.average(
            active_phases["mean_module_relative_abundance"], weights=active_phases["samples_n"]
        ))
        delta = active_mean - baseline_mean
        label = "Higher during renewal/post-renewal" if delta > 0 else (
            "Higher during stagnation" if delta < 0 else "Equal pooled abundance"
        )
        labels[module] = label
        rows.append({
            "module_label": module,
            "stagnation_mean_module_relative_abundance": baseline_mean,
            "renewal_post_renewal_mean_module_relative_abundance": active_mean,
            "active_minus_stagnation_mean_difference": delta,
            "active_to_stagnation_mean_ratio": active_mean / baseline_mean if baseline_mean > 0 else np.nan,
            "direction_label": label,
        })
    direction = pd.DataFrame(rows, columns=[
        "module_label", "stagnation_mean_module_relative_abundance",
        "renewal_post_renewal_mean_module_relative_abundance",
        "active_minus_stagnation_mean_difference",
        "active_to_stagnation_mean_ratio", "direction_label",
    ])
    active = active.merge(direction, on="module_label", how="left")
    return labels, active


def plot_proportionality_heatmaps(
    pair_statistics: pd.DataFrame,
    graph: nx.Graph,
    annotation: pd.DataFrame,
    outdir: Path,
    modality: str,
) -> None:
    """Render full and selected signed proportionality matrices."""
    required = {"feature_1", "feature_2", "rho_p", "selected_edge"}
    if not required.issubset(pair_statistics.columns):
        raise ValueError(f"pair-statistics table lacks columns: {sorted(required - set(pair_statistics.columns))}")
    node_key = {node: key(graph.nodes[node].get("name", node)) for node in graph.nodes()}
    by_key = annotation.set_index("genome_key")
    ordered = sorted(
        graph.nodes(),
        key=lambda node: (
            natural(by_key.loc[node_key[node], "proportionality_cluster"] if node_key[node] in by_key.index else ""),
            natural(graph.nodes[node].get("display_label", node)),
        ),
    )
    index = {str(node): i for i, node in enumerate(ordered)}
    full = np.eye(len(ordered), dtype=float)
    selected = np.full((len(ordered), len(ordered)), np.nan, dtype=float)
    np.fill_diagonal(selected, 1.0)
    for row in pair_statistics.itertuples(index=False):
        left, right = str(row.feature_1), str(row.feature_2)
        if left not in index or right not in index:
            continue
        i, j = index[left], index[right]
        full[i, j] = full[j, i] = float(row.rho_p)
        if truthy(pd.Series([row.selected_edge])).iloc[0]:
            selected[i, j] = selected[j, i] = float(row.rho_p)
    ylabels = [str(graph.nodes[node].get("display_label", node)) for node in ordered]
    xlabels = [str(graph.nodes[node].get("plot_number", i + 1)) for i, node in enumerate(ordered)]
    for values, stem, subtitle in (
        (full, "proportionality_heatmap_all_pairs", "All tested species pairs"),
        (selected, "proportionality_heatmap_supported_edges", "Supported network edges only"),
    ):
        fig, ax = plt.subplots(figsize=(15, 12))
        sns.heatmap(
            values, ax=ax, cmap="vlag", center=0, vmin=-1, vmax=1, square=True,
            xticklabels=xlabels, yticklabels=ylabels,
            cbar_kws={"label": r"Proportionality ($\rho_p$)", "shrink": 0.75},
            mask=np.isnan(values), linewidths=0.15, linecolor="#EFEFEF",
        )
        ax.set_xlabel("Species key number")
        ax.set_ylabel("GTDB species")
        ax.set_title(f"{modality.title()} species proportionality\n{subtitle}")
        ax.tick_params(axis="x", rotation=0)
        ax.tick_params(axis="y", rotation=0, labelsize=8)
        fig.tight_layout()
        for extension in ("pdf", "png", "svg"):
            kwargs = {"dpi": 600} if extension == "png" else {}
            fig.savefig(outdir / f"{stem}.{extension}", bbox_inches="tight", pad_inches=0.35, **kwargs)
        plt.close(fig)


def attach_metrics(graph: nx.Graph, node_features: pd.DataFrame) -> pd.DataFrame:
    id_col = "Taxon" if "Taxon" in node_features else node_features.columns[1]
    features = node_features.copy()
    features["genome_key"] = features[id_col].map(key)
    by_key = features.set_index("genome_key")
    for node, genome in graph_node_to_genome(graph).items():
        genome_key = key(genome)
        if genome_key not in by_key.index:
            continue
        row = by_key.loc[genome_key]
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]
        for column in ("Degree", "Betweenness", "Betweenness_norm", "Closeness", "EigenCentral"):
            if column in row:
                value = pd.to_numeric(row[column], errors="coerce")
                graph.nodes[node][column] = float(value) if pd.notna(value) else 0.0
    return features


def strongest_isa(
    node_metadata: pd.DataFrame,
    crosswalk: pd.DataFrame,
    summary: pd.DataFrame,
    prefix: str,
) -> pd.DataFrame:
    result = node_metadata[["Genome_Id", "genome_key"]].copy()
    required = {"ASV", "association_groups", "stat"}
    if not required.issubset(summary):
        result[f"propagated_{prefix}_label"] = "No significant linked-ASV association"
        return result
    isa = summary.copy()
    isa["ASV_ID"] = isa["ASV"].astype(str).str.replace(r";.*$", "", regex=True)
    if "significant" in isa:
        isa = isa.loc[truthy(isa["significant"])]
    if "q.value" in isa:
        isa["q_value"] = pd.to_numeric(isa["q.value"], errors="coerce")
        isa = isa.loc[isa["q_value"].le(0.05)]
    else:
        isa["q_value"] = np.nan
    isa["stat"] = pd.to_numeric(isa["stat"], errors="coerce")
    joined = crosswalk[["genome_key", "ASV_ID"]].dropna().merge(isa, on="ASV_ID", how="inner")
    rows = []
    for genome_key, frame in joined.groupby("genome_key"):
        frame = frame.sort_values(["q_value", "stat", "ASV_ID"], ascending=[True, False, True], na_position="last")
        best = frame.iloc[0]
        rows.append({
            "genome_key": genome_key,
            f"propagated_{prefix}_label": str(best["association_groups"]),
            f"propagated_{prefix}_source_asv": str(best["ASV_ID"]),
            f"propagated_{prefix}_stat": best["stat"],
            f"propagated_{prefix}_q": best["q_value"],
            f"propagated_{prefix}_significant_linked_asvs_n": frame["ASV_ID"].nunique(),
            f"propagated_{prefix}_distinct_labels": ";".join(sorted(set(frame["association_groups"].dropna().astype(str)))),
            f"propagated_{prefix}_selection_rule": "minimum q, then maximum indicator statistic, then ASV ID",
        })
    annotations = pd.DataFrame(rows)
    result = result.merge(annotations, on="genome_key", how="left")
    result[f"propagated_{prefix}_label"] = result[f"propagated_{prefix}_label"].fillna(
        "No significant linked-ASV association"
    )
    return result


def titan_annotations(
    node_metadata: pd.DataFrame,
    crosswalk: pd.DataFrame,
    titan: pd.DataFrame,
) -> pd.DataFrame:
    work = titan.copy()
    work["ASV_ID"] = work["ASV_ID"].astype(str).str.replace(r";.*$", "", regex=True)
    if "passes_purity_and_reliability" in work:
        work = work.loc[truthy(work["passes_purity_and_reliability"])]
    work["z_score"] = pd.to_numeric(work["z_score"], errors="coerce")
    joined = crosswalk[["genome_key", "ASV_ID"]].dropna().merge(work, on="ASV_ID", how="inner")
    rows = []
    for (genome_key, variable), frame in joined.groupby(["genome_key", "environmental_variable"]):
        frame = frame.assign(_abs_z=frame["z_score"].abs()).sort_values(
            ["_abs_z", "reliability", "purity", "ASV_ID"],
            ascending=[False, False, False, True],
        )
        best = frame.iloc[0]
        rows.append({
            "genome_key": genome_key,
            "environmental_variable": variable,
            "propagated_titan_direction": best["response_direction"],
            "propagated_titan_z_score": best["z_score"],
            "propagated_titan_change_point": best["change_point"],
            "propagated_titan_source_asv": best["ASV_ID"],
            "propagated_titan_linked_asvs_n": frame["ASV_ID"].nunique(),
            "propagated_titan_distinct_directions": ";".join(sorted(set(frame["response_direction"].dropna().astype(str)))),
            "propagated_titan_selection_rule": "maximum absolute z score, then reliability, purity, and ASV ID",
        })
    return pd.DataFrame(rows)


def draw_network(
    graph: nx.Graph,
    pos: dict,
    out_stem: Path,
    *,
    title: str,
    colors: dict,
    sizes: dict,
    edge_widths: dict | None = None,
    legend_handles: list | None = None,
) -> None:
    fig, ax = plt.subplots(figsize=(19, 13))
    weights = np.asarray([abs(float(data.get("weight", 0.0))) for _, _, data in graph.edges(data=True)])
    if weights.size:
        hi = max(float(np.nanmax(weights)), 1e-12)
        widths = 0.35 + 1.8 * weights / hi
    else:
        widths = 0.5
    edge_records = list(graph.edges(data=True))
    width_lookup = {
        frozenset((u, v)): float(width)
        for (u, v, _), width in zip(edge_records, np.atleast_1d(widths))
    } if edge_records else {}
    positive_edges = [(u, v) for u, v, data in edge_records if float(data.get("weight", 0.0)) >= 0]
    negative_edges = [(u, v) for u, v, data in edge_records if float(data.get("weight", 0.0)) < 0]
    if positive_edges:
        nx.draw_networkx_edges(
            graph, pos, edgelist=positive_edges, ax=ax, edge_color="#202020", alpha=0.82,
            width=[width_lookup[frozenset(edge)] for edge in positive_edges], style="solid",
        )
    if negative_edges:
        nx.draw_networkx_edges(
            graph, pos, edgelist=negative_edges, ax=ax, edge_color="#B5B5B5", alpha=0.95,
            width=[width_lookup[frozenset(edge)] for edge in negative_edges], style="dashed",
        )
    ordinary = list(graph.nodes())
    nx.draw_networkx_nodes(
        graph, pos, nodelist=ordinary, ax=ax,
        node_color=[colors.get(node, "#D0D0D0") for node in ordinary],
        node_size=[max(760.0, sizes.get(node, 760.0)) for node in ordinary],
        edgecolors=["black" for _ in ordinary],
        linewidths=[(edge_widths or {}).get(node, 0.75) for node in ordinary],
        alpha=1.0,
    )
    texts = []
    for node in graph.nodes():
        x = float(pos[node][0])
        y = float(pos[node][1])
        label = str(graph.nodes[node].get("plot_number", ""))
        fill = colors.get(node, "#D0D0D0")
        texts.append(ax.text(
            x, y, label, fontsize=11.5, weight="bold", ha="center", va="center",
            color=contrasting_text_color(fill), zorder=5,
        ))
    edge_handles = []
    if positive_edges:
        edge_handles.append(plt.Line2D([], [], color="#202020", lw=2.0, label="Positive proportionality"))
    if negative_edges:
        edge_handles.append(plt.Line2D([], [], color="#B5B5B5", lw=2.0, linestyle="--", label="Negative proportionality"))
    combined_handles = edge_handles + list(legend_handles or [])
    if combined_handles:
        fig.legend(
            handles=combined_handles,
            frameon=False,
            loc="lower center",
            bbox_to_anchor=(0.39, 0.015),
            ncol=min(4, len(combined_handles)),
            columnspacing=1.25,
            handletextpad=0.55,
        )
    species_handles = [
        mpatches.Patch(
            facecolor="none", edgecolor="none",
            label=f"{graph.nodes[node]['plot_number']}: {graph.nodes[node]['display_label']}",
        )
        for node in sorted(graph.nodes(), key=lambda item: int(graph.nodes[item]["plot_number"]))
    ]
    fig.legend(
        handles=species_handles,
        title="Species key",
        frameon=False,
        loc="center right",
        bbox_to_anchor=(0.995, 0.5),
        handlelength=0.0,
        handletextpad=0.0,
        labelspacing=0.45,
        fontsize=8.0,
        title_fontsize=9.0,
    )
    ax.set_title(title)
    ax.set_aspect("equal", adjustable="datalim")
    ax.axis("off")
    fig.subplots_adjust(bottom=0.19, right=0.72)
    out_stem.parent.mkdir(parents=True, exist_ok=True)
    for extension in ("pdf", "png", "svg"):
        kwargs = {"dpi": 600} if extension == "png" else {}
        fig.savefig(out_stem.with_suffix(f".{extension}"), bbox_inches="tight", pad_inches=0.45, **kwargs)
    plt.close(fig)


def metric_sizes(graph: nx.Graph, metric: str) -> tuple[dict, list]:
    values = np.asarray([float(graph.nodes[node].get(metric, 0.0)) for node in graph.nodes()], dtype=float)
    finite = values[np.isfinite(values)]
    lo = float(finite.min()) if finite.size else 0.0
    hi = float(finite.max()) if finite.size else 1.0
    mapper = lambda value: 760.0 + 840.0 * ((float(value) - lo) / (hi - lo) if hi > lo else 0.5)
    sizes = {node: mapper(graph.nodes[node].get(metric, 0.0)) for node in graph.nodes()}
    handles = [
        plt.scatter([], [], s=mapper(value), facecolor="#AFAFAF", edgecolor="black", label=f"{value:.3g}")
        for value in np.linspace(lo, hi, 3)
    ]
    return sizes, handles


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--graph", required=True, type=Path)
    parser.add_argument("--node-features", required=True, type=Path)
    parser.add_argument("--network-clusters", required=True, type=Path)
    parser.add_argument("--matrix", required=True, type=Path)
    parser.add_argument("--pair-statistics", required=True, type=Path)
    parser.add_argument("--node-metadata", required=True, type=Path)
    parser.add_argument("--asv-crosswalk", required=True, type=Path)
    parser.add_argument("--hybrid-isa", required=True, type=Path)
    parser.add_argument("--mc-isa", required=True, type=Path)
    parser.add_argument("--titan", required=True, type=Path)
    parser.add_argument("--module-renewal-association", type=Path)
    parser.add_argument("--module-renewal-profiles", type=Path)
    parser.add_argument("--renewal-max-q", type=float, default=0.05)
    parser.add_argument("--modality", required=True, choices=["metagenome", "metatranscriptome"])
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    graph = nx.read_graphml(args.graph)
    node_features = attach_metrics(graph, read_table(args.node_features))
    node_metadata = read_table(args.node_metadata)
    node_metadata["genome_key"] = node_metadata["Genome_Id"].map(key)
    crosswalk = read_table(args.asv_crosswalk)
    crosswalk["genome_key"] = crosswalk["Genome_Id"].map(key)
    hybrid = strongest_isa(node_metadata, crosswalk, read_table(args.hybrid_isa), "hybrid")
    mc = strongest_isa(node_metadata, crosswalk, read_table(args.mc_isa), "microbial_compartment")
    titan = titan_annotations(node_metadata, crosswalk, read_table(args.titan))
    annotation = node_metadata.merge(hybrid.drop(columns=["Genome_Id"]), on="genome_key", how="left")
    annotation = annotation.merge(mc.drop(columns=["Genome_Id"]), on="genome_key", how="left")

    clusters = read_table(args.network_clusters)
    cluster_col = "Taxon" if "Taxon" in clusters else clusters.columns[0]
    clusters["genome_key"] = clusters[cluster_col].map(key)
    annotation = annotation.merge(
        clusters[["genome_key", "proportionality_cluster", "node_stability"]],
        on="genome_key", how="left"
    )
    matrix = read_table(args.matrix).set_index("Feature_ID")
    mean_abundance = matrix.mean(axis=1).rename("mean_raw_recruitment_count")
    mean_abundance.index = mean_abundance.index.map(key)
    annotation = annotation.merge(mean_abundance.rename_axis("genome_key").reset_index(), on="genome_key", how="left")
    annotation.to_csv(args.outdir / "genome_network_node_annotations.tsv", sep="\t", index=False)
    titan.to_csv(args.outdir / "genome_network_propagated_titan.tsv", sep="\t", index=False)

    by_key = annotation.set_index("genome_key")
    node_keys = {node: key(genome) for node, genome in graph_node_to_genome(graph).items()}
    for node in graph.nodes():
        genome_key = node_keys[node]
        if genome_key in by_key.index and pd.notna(by_key.loc[genome_key, "Species"]):
            graph.nodes[node]["display_label"] = str(by_key.loc[genome_key, "Species"])
        else:
            graph.nodes[node]["display_label"] = graph_node_to_genome(graph)[node]
    ordered_nodes = sorted(graph.nodes(), key=lambda node: natural(graph.nodes[node]["display_label"]))
    for number, node in enumerate(ordered_nodes, start=1):
        graph.nodes[node]["plot_number"] = number
    pos = component_aware_layout(graph, args.seed)
    layout = pd.DataFrame([
        {
            "graph_node_id": node,
            "genome_id": graph_node_to_genome(graph)[node],
            "species_label": graph.nodes[node]["display_label"],
            "plot_number": graph.nodes[node]["plot_number"],
            "x": xy[0],
            "y": xy[1],
            "isolate": graph.degree(node) == 0,
            "seed": args.seed,
        }
        for node, xy in pos.items()
    ])
    layout.to_csv(args.outdir / "genome_network_layout.tsv", sep="\t", index=False)
    base_sizes = {node: 900.0 for node in graph.nodes()}

    summary = pd.DataFrame([{
        "modality": args.modality,
        "selected_species_in_matrix": int(matrix.shape[0]),
        "species_in_graph": int(graph.number_of_nodes()),
        "inferred_edges": int(graph.number_of_edges()),
        "isolated_species": int(len(list(nx.isolates(graph)))),
        "connected_components": int(nx.number_connected_components(graph)),
        "assay_samples": int(matrix.shape[1]),
    }])
    summary.to_csv(args.outdir / "genome_network_summary.tsv", sep="\t", index=False)

    def categorical_plot(column: str, label: str, stem: str) -> None:
        labels = {
            node: str(by_key.loc[node_keys[node], column])
            if node_keys[node] in by_key.index and pd.notna(by_key.loc[node_keys[node], column])
            else "Unclassified"
            for node in graph.nodes()
        }
        pal = palette(labels.values())
        handles = [mpatches.Patch(facecolor=pal[value], label=value) for value in sorted(pal, key=natural)]
        draw_network(
            graph, pos, args.outdir / stem,
            title=f"{args.modality.title()} species proportionality network\nNode color: {label}",
            colors={node: pal[labels[node]] for node in graph.nodes()}, sizes=base_sizes,
            legend_handles=handles,
        )

    categorical_plot("Phylum", "GTDB phylum", "network_phylum")
    categorical_plot("Family", "GTDB family", "network_family")
    categorical_plot(
        "proportionality_cluster", "native proportionality cluster",
        "network_proportionality_clusters"
    )
    annotation["asv_ecological_module_overlay"] = np.where(
        annotation["propagated_ecological_modules_n"].fillna(0).eq(0), "No linked ASV",
        np.where(annotation["propagated_ecological_modules_n"].eq(1), annotation["propagated_ecological_modules"], "Multiple linked modules")
    )
    by_key = annotation.set_index("genome_key")
    categorical_plot("asv_ecological_module_overlay", "linked-ASV ecological module", "network_propagated_asv_ecological_module")
    categorical_plot("propagated_hybrid_label", "linked-ASV hybrid indicator association", "network_propagated_hybrid")
    categorical_plot("propagated_microbial_compartment_label", "linked-ASV microbial-compartment indicator association", "network_propagated_microbial_compartment")

    renewal_labels, renewal_significant = renewal_module_labels(
        args.module_renewal_association, args.renewal_max_q
    )
    if renewal_labels:
        module_by_asv = (
            crosswalk[["genome_key", "ASV_ID", "module_label"]]
            .dropna(subset=["genome_key", "module_label"])
            .drop_duplicates()
        )
        rows = []
        for genome_key, frame in module_by_asv.groupby("genome_key"):
            supported = frame.loc[frame["module_label"].astype(str).isin(renewal_labels)]
            labels = sorted({renewal_labels[str(value)] for value in supported["module_label"]})
            rows.append({
                "genome_key": genome_key,
                "propagated_renewal_module_label": labels[0] if len(labels) == 1 else (
                    "Multiple supported renewal-module associations" if labels else
                    "No supported renewal-phase contrast"
                ),
                "propagated_renewal_source_asvs": ";".join(sorted(set(supported["ASV_ID"].astype(str)))),
                "propagated_renewal_source_modules": ";".join(sorted(set(supported["module_label"].astype(str)), key=natural)),
            })
        renewal_annotation = pd.DataFrame(rows)
        annotation = annotation.merge(renewal_annotation, on="genome_key", how="left")
        annotation["propagated_renewal_module_label"] = annotation["propagated_renewal_module_label"].fillna(
            "No supported renewal-phase contrast"
        )
        by_key = annotation.set_index("genome_key")
        categorical_plot(
            "propagated_renewal_module_label",
            "linked-ASV ecological-module renewal association",
            "network_propagated_renewal_associated_modules",
        )
        renewal_significant.to_csv(
            args.outdir / "genome_network_source_module_renewal_associations.tsv", sep="\t", index=False
        )

    active_labels, active_significant = active_renewal_module_labels(
        args.module_renewal_association, args.module_renewal_profiles, args.renewal_max_q
    )
    if active_labels:
        module_by_asv = (
            crosswalk[["genome_key", "ASV_ID", "module_label"]]
            .dropna(subset=["genome_key", "module_label"])
            .drop_duplicates()
        )
        rows = []
        for genome_key, frame in module_by_asv.groupby("genome_key"):
            supported = frame.loc[frame["module_label"].astype(str).isin(active_labels)]
            directions = sorted({active_labels[str(value)] for value in supported["module_label"]})
            rows.append({
                "genome_key": genome_key,
                "propagated_active_renewal_module_label": (
                    directions[0] if len(directions) == 1 else (
                        "Mixed directional module links" if directions
                        else "No supported renewal/post-renewal contrast"
                    )
                ),
                "propagated_active_renewal_source_asvs": ";".join(
                    sorted(set(supported["ASV_ID"].astype(str)))
                ),
                "propagated_active_renewal_source_modules": ";".join(
                    sorted(set(supported["module_label"].astype(str)), key=natural)
                ),
            })
        active_annotation = pd.DataFrame(rows)
        annotation = annotation.merge(active_annotation, on="genome_key", how="left")
        annotation["propagated_active_renewal_module_label"] = annotation[
            "propagated_active_renewal_module_label"
        ].fillna("No supported renewal/post-renewal contrast")
        by_key = annotation.set_index("genome_key")
        categorical_plot(
            "propagated_active_renewal_module_label",
            "direction among significant linked-ASV module contrasts",
            "network_propagated_renewal_post_renewal_vs_stagnation",
        )
        active_significant.to_csv(
            args.outdir / "genome_network_source_active_renewal_module_associations.tsv",
            sep="\t", index=False,
        )

    abundance = pd.to_numeric(annotation["mean_raw_recruitment_count"], errors="coerce").fillna(0.0)
    positive = abundance[abundance.gt(0)]
    lo, hi = (float(np.log10(positive.min() + 1)), float(np.log10(positive.max() + 1))) if not positive.empty else (0.0, 1.0)
    abundance_map = dict(zip(annotation["genome_key"], abundance))
    abundance_sizes = {
        node: 760.0 + 840.0 * ((math.log10(abundance_map.get(node_keys[node], 0.0) + 1) - lo) / (hi - lo) if hi > lo else 0.5)
        for node in graph.nodes()
    }
    draw_network(
        graph, pos, args.outdir / "network_mean_abundance",
        title=(f"{args.modality.title()} species proportionality network\n"
               "Node size: mean raw recruited-read count"),
        colors={node: "#AFAFAF" for node in graph.nodes()}, sizes=abundance_sizes,
    )
    for metric, stem in (
        ("Degree", "network_degree"), ("Betweenness_norm", "network_betweenness"),
        ("EigenCentral", "network_eigenvector_centrality"), ("Closeness", "network_closeness"),
    ):
        sizes, handles = metric_sizes(graph, metric)
        draw_network(
            graph, pos, args.outdir / stem,
            title=f"{args.modality.title()} species proportionality network\nNode size: {metric}",
            colors={node: "#AFAFAF" for node in graph.nodes()}, sizes=sizes,
            legend_handles=handles,
        )

    plot_proportionality_heatmaps(
        read_table(args.pair_statistics), graph, annotation, args.outdir, args.modality
    )

    if not titan.empty:
        for variable, frame in titan.groupby("environmental_variable"):
            lookup = frame.set_index("genome_key")
            directions = {
                node: str(lookup.loc[node_keys[node], "propagated_titan_direction"])
                if node_keys[node] in lookup.index else "No reliable linked-ASV TITAN response"
                for node in graph.nodes()
            }
            direction_palette = {"z+": "#111111", "z-": "#BDBDBD", "No reliable linked-ASV TITAN response": "#F2F2F2"}
            zmax = max([abs(float(value)) for value in frame["propagated_titan_z_score"].dropna()] or [1.0])
            sizes = {
                node: 760.0 + 840.0 * abs(float(lookup.loc[node_keys[node], "propagated_titan_z_score"])) / zmax
                if node_keys[node] in lookup.index else 760.0
                for node in graph.nodes()
            }
            handles = [
                mpatches.Patch(
                    facecolor=direction_palette[label], edgecolor="black",
                    label={
                        "z+": "z+ (increases along gradient)",
                        "z-": "z- (declines along gradient)",
                        "No reliable linked-ASV TITAN response": "No reliable linked-ASV TITAN response",
                    }[label],
                )
                for label in ("z+", "z-", "No reliable linked-ASV TITAN response")
            ]
            slug = re.sub(r"[^A-Za-z0-9]+", "_", str(variable)).strip("_").lower()
            draw_network(
                graph, pos, args.outdir / f"network_propagated_titan_{slug}",
                title=(f"{args.modality.title()} species proportionality network\n"
                       f"Linked-ASV TITAN response to {variable}; node size: |z|"),
                colors={node: direction_palette.get(directions[node], "#F2F2F2") for node in graph.nodes()},
                sizes=sizes, legend_handles=handles,
            )

    annotation.to_csv(args.outdir / "genome_network_node_annotations.tsv", sep="\t", index=False)
    parameters = {
        "modality": args.modality,
        "layout_seed": args.seed,
        "layout_reused_for_every_overlay": True,
        "edge_metric": "rho_p proportionality",
        "edge_interpretation": "bootstrap-stable concordant or contrasting relative recruitment profiles; not direct interaction",
        "propagated_evidence": [
            "linked-ASV ecological module", "hybrid ISA", "microbial-compartment ISA",
            "TITAN", "ecological-module renewal-phase association",
        ],
        "propagated_evidence_interpretation": "descriptive linked-ASV annotation; not direct genome-level inference",
        "renewal_display_baseline_label": "stagnation",
        "active_renewal_direction_rule": (
            "among ecological modules with active_vs_baseline adjusted p <= threshold, "
            "compare the sample-count-weighted pooled renewal/post-renewal mean module "
            "relative abundance with the stagnation mean"
        ),
        "matrix_views": ["all tested proportionalities", "supported network edges only"],
    }
    (args.outdir / "genome_network_plot_parameters.json").write_text(json.dumps(parameters, indent=2) + "\n")


if __name__ == "__main__":
    main()
