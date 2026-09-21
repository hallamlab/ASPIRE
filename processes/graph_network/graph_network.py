#!/usr/bin/env python3
"""
Network visualization CLI for SPIEC-EASI GraphML outputs + ISA/metadata overlays.

Features
- Robust CLI (argparse)
- All paths & style values are configurable
- Safe I/O + schema checks
- Reusable spring-layout (cache to JSON)
- Multiple plot “modes” (degree, abundance, ISA overlays for two group summaries, phylum×{abundance,ISA}, labeled variants)
- Consistent aesthetics with your rcParams
"""

import argparse
import gc
import json
import os
import re
import shutil
import sys
from pathlib import Path
from typing import Dict, Iterable, Tuple, Optional, List

import networkx as nx
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
_sys_path_style = str(Path(__file__).resolve().parents[1])
if _sys_path_style not in sys.path:
    sys.path.insert(0, _sys_path_style)
from shared_plot_style import install_publication_style
install_publication_style()
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import matplotlib.patheffects as mpatheffects
from matplotlib.path import Path as MplPath
import seaborn as sns
import math
from collections import Counter

# Optional: labels w/ collision-avoidance
try:
    from adjustText import adjust_text
    _HAS_ADJUSTTEXT = True
except Exception:
    _HAS_ADJUSTTEXT = False

# ---------------------------- Global style -----------------------------------
mpl.rcParams['pdf.fonttype'] = 42                 # text as text in PDF
mpl.rcParams['svg.fonttype'] = 'none'             # text as text in SVG
plt.rcParams.update({'font.size': 22})
mpl.rcParams['savefig.dpi'] = 600
plt.rcParams['font.family'] = 'Times New Roman'
sns.set_theme()
sns.set_style("white")
NOT_FOCUS_COLOR = "#D3D3D3"
MAG_PAIRED_COLOR = "#000000"
MAG_UNPAIRED_COLOR = "#D9D9D9"
EXPLICIT_GROUP_TYPE_PALETTE = {
    "Oral Rinse": "#6A3D9A",
    "BAL": "#0072B2",
    "Bronchial Brush": "#009E73",
    "Lung Brush": "#009E73",
    "BAL+Oral Rinse": "#E78AC3",
    "BAL+Bronchial Brush": "#5CC8C8",
    "BAL+Lung Brush": "#5CC8C8",
    "Bronchial Brush+Oral Rinse": "#B8E186",
    "Lung Brush+Oral Rinse": "#B8E186",
    "BAL+Bronchial Brush+Oral Rinse": "#CBB6E9",
    "BAL+Lung Brush+Oral Rinse": "#CBB6E9",
    "not_indicator": "#D3D3D3",
}
CANONICAL_GROUP_TYPE_ORDER = [
    "Oral Rinse",
    "BAL+Oral Rinse",
    "BAL",
    "BAL+Bronchial Brush",
    "Bronchial Brush",
    "Bronchial Brush+Oral Rinse",
    "BAL+Bronchial Brush+Oral Rinse",
    "not_indicator",
]
MAG_MIXED_TAXONOMY = "Mixed MAG phyla"
MAG_UNKNOWN_TAXONOMY = "Unclassified MAG phylum"
ALLOWED_MAG_MIMAG_TIERS = {"medium", "high"}
BETWEENNESS_HIGH_COLOR = "#4a4a4a"
BETWEENNESS_LOW_COLOR = "#d0d0d0"
SIZE_CONFIG = {
    "degree_size_mode": "legacy",
    "degree_min_area": 0.0,
    "abundance_size_mode": "legacy",
    "abundance_reference": 5000.0,
    "abundance_reference_area": 80.0,
}
LABEL_CONFIG = {"max_labels": 100}


# ---------------------------- Helpers ----------------------------------------
def die(msg: str, code: int = 2):
    print(f"[ERROR] {msg}", file=sys.stderr)
    sys.exit(code)


def ok(msg: str):
    print(f"[+] {msg}")


def load_table(path: str, sep: str = None, index_col: Optional[int | str] = None) -> pd.DataFrame:
    if not os.path.exists(path):
        die(f"Missing file: {path}")
    if sep is None:
        # guess by extension
        sep = '\t' if path.endswith(('.tsv', '.tab')) else ','
    try:
        return pd.read_csv(path, sep=sep, header=0, index_col=index_col)
    except Exception as e:
        die(f"Failed to read {path}: {e}")


def ensure_cols(df: pd.DataFrame, cols: Iterable[str], name: str):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        die(f"{name} missing required columns: {missing}")


def _safe_float(x, default=0.0):
    try:
        if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))):
            return default
        return float(x)
    except Exception:
        return default


def is_not_focus_color(color: object) -> bool:
    val = str(color).strip().lower()
    return val in {"lightgray", "lightgrey", "#d3d3d3", "d3d3d3"}

def preflight_node_attr_report(
    G: nx.Graph,
    name: str,
    color_attr: str,
    size_attr: str = "AxB",
    palette_values: set | None = None
) -> None:
    """Print helpful counts before plotting."""
    n_nodes = G.number_of_nodes()
    n_edges = G.number_of_edges()
    print(f"[check:{name}] nodes={n_nodes} edges={n_edges}")

    if n_nodes == 0:
        print(f"[check:{name}] graph has no nodes")
        return

    colors = [G.nodes[n].get(color_attr, None) for n in G.nodes()]
    sizes  = [G.nodes[n].get(size_attr, None) for n in G.nodes()]

    n_missing_color = sum(c is None for c in colors)
    n_missing_size  = sum(s is None for s in sizes)
    n_nonzero_size  = sum(_safe_float(s, 0.0) > 0 for s in sizes)

    # Compact distribution of color values
    color_counts = Counter(c if c is not None else "<missing>" for c in colors)
    top_colors = ", ".join([f"{k}:{v}" for k, v in color_counts.most_common(8)])
    print(f"[check:{name}] color_attr='{color_attr}': missing={n_missing_color}, top={top_colors}")

    if palette_values is not None:
        n_outside = sum((c not in palette_values) for c in colors if c is not None)
        print(f"[check:{name}] colors not in palette: {n_outside}")

    # Basic size stats
    vals = np.array([_safe_float(s, 0.0) for s in sizes], dtype=float)
    if len(vals):
        print(f"[check:{name}] size_attr='{size_attr}': >0={n_nonzero_size}, "
              f"min={vals.min():.3g}, median={np.median(vals):.3g}, max={vals.max():.3g}")
    else:
        print(f"[check:{name}] size_attr='{size_attr}': no values found")

def save_figure(figpath: str) -> None:
    """Safe save; avoid crashing on empty collections."""
    try:
        fig = plt.gcf()
        # If nothing was added (no artists), add an invisible dot so renderer has something.
        ax = plt.gca()
        if not (ax.collections or ax.patches or ax.lines):
            ax.plot([0], [0], alpha=0)  # invisible fallback
        fig.savefig(figpath, bbox_inches='tight', pad_inches=0.6)
        if figpath.endswith(".svg"):
            fig.savefig(figpath.replace(".svg", ".pdf"), bbox_inches='tight', pad_inches=0.6)
        print(f"[+] Saved: {figpath}")
    except Exception as e:
        print(f"[!] save_figure failed for {figpath}: {e}")
    finally:
        fig = plt.gcf()
        plt.close(fig)
        # Large network figures contain many mutually-referencing Matplotlib
        # artists.  Explicit collection prevents their memory accumulating
        # across a long --modes run.
        gc.collect()


def natural_sort_key(value: object) -> Tuple:
    text = str(value).strip()
    parts = re.split(r"(\d+)", text)
    out: List[object] = []
    for part in parts:
        if not part:
            continue
        if part.isdigit():
            out.append(int(part))
        else:
            out.append(part.lower())
    return tuple(out)


def ordered_present_labels(present_labels: Iterable[str], preferred_order: Optional[List[str]] = None) -> List[str]:
    present = normalize_order_list([x for x in present_labels if str(x).strip()])
    if not present:
        return []
    if preferred_order:
        preferred = normalize_order_list(preferred_order)
        ordered = [lbl for lbl in preferred if lbl in present]
        if ordered:
            ordered.extend([lbl for lbl in present if lbl not in ordered])
            return ordered
    return sorted(present, key=natural_sort_key)


def figure_ax(figsize: Tuple[float, float] = (16, 14)) -> Tuple[plt.Figure, plt.Axes]:
    fig = plt.figure(figsize=figsize)
    ax = plt.gca()
    return fig, ax


def finalize_network_axes(
    fig: plt.Figure,
    ax: plt.Axes,
    out_svg: str,
    *,
    title: str,
    legends: Optional[List[object]] = None,
    right: float = 0.74,
) -> None:
    ax.set_aspect("equal", adjustable="datalim")
    ax.autoscale(enable=True)
    ax.set_title(title)
    ax.axis("off")
    fig.subplots_adjust(right=right)
    extra = [lg for lg in (legends or []) if lg is not None]
    try:
        if not (ax.collections or ax.patches or ax.lines):
            ax.plot([0], [0], alpha=0)
        fig.savefig(out_svg, bbox_inches="tight", pad_inches=0.5, bbox_extra_artists=extra)
        if out_svg.endswith(".svg"):
            fig.savefig(out_svg.replace(".svg", ".pdf"), bbox_inches="tight", pad_inches=0.5, bbox_extra_artists=extra)
        print(f"[+] Saved: {out_svg}")
    finally:
        plt.close(fig)


def build_isa_size_legend_values(raw_scores: Iterable[float], n_legend: int = 5) -> List[float]:
    # Keep ISA legends consistent across every graph so comparisons are direct.
    shared = [1.0, 0.75, 0.50, 0.25]
    return shared[:n_legend]

def split_taxa_string(taxa_str: str, delimiter: str = ';') -> Dict[str, Optional[str]]:
    """Split a SILVA/Greengenes-like lineage into 7 standard levels."""
    levels = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    if pd.isna(taxa_str) or taxa_str == 'Unassigned':
        return {lv: ('Unassigned' if lv == 'Domain' else None) for lv in levels}
    parts = [p.strip() for p in taxa_str.split(delimiter)]
    parts = [p.split('__', 1)[1] if '__' in p else p for p in parts]
    out = {}
    for i, lv in enumerate(levels):
        out[lv] = parts[i] if i < len(parts) else None
    return out


def write_layout_json(path: str, pos: Dict[str, Tuple[float, float]]):
    try:
        with open(path, 'w') as f:
            json.dump(pos, f)
        ok(f"Cached layout: {path}")
    except Exception as e:
        print(f"[WARN] Could not write layout JSON: {e}", file=sys.stderr)


def read_layout_json(path: str) -> Optional[Dict[str, Tuple[float, float]]]:
    if not path or not os.path.exists(path):
        return None
    try:
        with open(path, 'r') as f:
            d = json.load(f)
        # keys must be strings; coords list->tuple
        return {str(k): tuple(v) for k, v in d.items()}
    except Exception as e:
        print(f"[WARN] Failed loading layout JSON ({path}): {e}", file=sys.stderr)
        return None


def spring_layout_cached(G: nx.Graph, seed: int, scale_xy: float,
                         layout_json: Optional[str]) -> Dict[str, Tuple[float, float]]:
    pos = read_layout_json(layout_json)
    if pos is None:
        ok("Computing spring_layout (NetworkX) ...")
        pos = nx.spring_layout(G, seed=seed)
        if scale_xy != 1.0:
            pos = {n: (x * scale_xy, y * scale_xy) for n, (x, y) in pos.items()}
        if layout_json:
            write_layout_json(layout_json, {k: list(v) for k, v in pos.items()})
    else:
        ok(f"Loaded layout from: {layout_json}")
        if scale_xy != 1.0:
            pos = {n: (x * scale_xy, y * scale_xy) for n, (x, y) in pos.items()}
    return pos


def add_node_attrs_from_df(G: nx.Graph, attrs_df: pd.DataFrame, keep_cols: Iterable[str]):
    """Copy selected columns from attrs_df (index must match node IDs) into G."""
    found = list(set(G.nodes()).intersection(set(attrs_df.index)))
    if not found:
        print("[WARN] No overlapping node IDs between graph and attribute table.")
        return
    for col in keep_cols:
        if col in attrs_df.columns:
            nx.set_node_attributes(G, values=attrs_df.loc[found, col].to_dict(), name=col)


def size_legend_handles(vals: Iterable[float], label_prefix: str, scale: float,
                        edge='black', face=NOT_FOCUS_COLOR) -> List[plt.Line2D]:
    return [plt.scatter([], [], s=v * scale, edgecolors=edge, facecolors=face,
                        alpha=1, label=f"{label_prefix} {v}")
            for v in vals]


def _format_legend_num(v: float) -> str:
    """Short numeric labels for legend readability."""
    if not np.isfinite(v):
        return "NA"
    if abs(v - round(v)) < 1e-9:
        return f"{int(round(v))}"
    if v == 0:
        return "0"
    av = abs(v)
    if av >= 1000:
        return f"{v:.1e}"
    if av >= 10:
        return f"{v:.0f}"
    if av >= 1:
        return f"{v:.1f}"
    return f"{v:.2f}"


def _format_legend_float(v: float, decimals: int = 2) -> str:
    if not np.isfinite(v):
        return "NA"
    text = f"{float(v):.{int(decimals)}f}"
    if "." in text:
        text = text.rstrip("0").rstrip(".")
    return text


def _nice_integer_step(max_val: float) -> int:
    if max_val <= 20:
        return 1
    if max_val <= 100:
        return 5
    if max_val <= 500:
        return 10
    if max_val <= 2000:
        return 25
    if max_val <= 10000:
        return 50
    if max_val <= 50000:
        return 100
    return 500


def _build_nice_integer_legend_values(vals: np.ndarray, n_legend: int = 5) -> List[float]:
    """Simple log ladder legend values: 1, 50, 500, 5000, 50000, ..."""
    vals = vals[np.isfinite(vals) & (vals > 0)]
    if vals.size == 0:
        return [1.0]

    vmax = float(np.max(vals))

    ladder = [1.0]
    v = 50.0
    # Build enough ladder values to cover data and one step above it.
    while v <= (vmax * 10.0):
        ladder.append(v)
        v *= 10.0

    below = [x for x in ladder if x <= vmax]
    above = [x for x in ladder if x > vmax]
    out = below if below else [1.0]

    # Add one value above observed range for context if possible.
    if above:
        out.append(above[0])

    # Keep concise.
    if len(out) > n_legend:
        out = out[:n_legend]
    return out


def build_log_size_scaler(
    raw_vals: Iterable[float],
    min_area: float = 8.0,
    max_area: float = 420.0,
    scale_power: float = 1.6,
    n_legend: int = 5,
) -> Tuple:
    """
    Build a log1p-based size scaler plus legend values from observed data.
    Interpolate in marker-radius space rather than area so decade differences
    remain visually distinct in the rendered scatter plot.
    Returns (mapper_function, legend_values).
    """
    vals = np.asarray([_safe_float(v, 0.0) for v in raw_vals], dtype=float)
    vals = vals[np.isfinite(vals) & (vals > 0)]
    if vals.size == 0:
        return (lambda _x: min_area), [1.0]

    abundance_mode = str(SIZE_CONFIG.get("abundance_size_mode", "legacy")).strip().lower()
    if abundance_mode == "decade_doubling":
        reference = max(_safe_float(SIZE_CONFIG.get("abundance_reference", 5000.0), 5000.0), 1e-12)
        reference_area = max(_safe_float(SIZE_CONFIG.get("abundance_reference_area", 80.0), 80.0), 1e-9)

        def mapper(x: float) -> float:
            xv = max(_safe_float(x, 0.0), 0.0)
            if xv <= 0.0:
                return float(min_area)
            area = reference_area * (2.0 ** math.log10(xv / reference))
            area = max(float(min_area), area)
            if max_area is not None:
                area = min(float(max_area), area)
            return float(area)

        legend_vals = _build_nice_integer_legend_values(vals, n_legend=n_legend)
        return mapper, legend_vals

    log_vals = np.log10(vals + 1.0)
    lo = float(np.min(log_vals))
    hi = float(np.max(log_vals))

    if hi <= lo + 1e-12:
        const_area = float((min_area + max_area) / 2.0)
        uniq = sorted({float(v) for v in vals.tolist()})
        legend_vals = uniq[:1] if uniq else [1.0]
        return (lambda _x: const_area), legend_vals

    min_radius = math.sqrt(max(float(min_area), 1.0))
    max_radius = math.sqrt(max(float(max_area), min_area + 1.0))

    def mapper(x: float) -> float:
        xv = max(_safe_float(x, 0.0), 0.0)
        lx = np.log10(xv + 1.0)
        t = (lx - lo) / (hi - lo)
        t = float(np.clip(t, 0.0, 1.0))
        # Expand the high end in radius space so 5k vs 50k stays readable.
        t = t ** max(float(scale_power), 0.1)
        radius = min_radius + t * (max_radius - min_radius)
        return float(radius * radius)

    legend_vals = _build_nice_integer_legend_values(vals, n_legend=n_legend)

    return mapper, legend_vals


def degree_marker_area(degree: float, degree_scale: float) -> float:
    deg = max(_safe_float(degree, 0.0), 0.0)
    degree_mode = str(SIZE_CONFIG.get("degree_size_mode", "legacy")).strip().lower()
    if degree_mode == "linear":
        min_area = max(_safe_float(SIZE_CONFIG.get("degree_min_area", 0.0), 0.0), 0.0)
        return float(min_area + (deg * max(_safe_float(degree_scale, 1.0), 0.0)))
    scale = max(_safe_float(degree_scale, 80.0), 1.0) * 0.14
    return float(max(8.0, ((deg + 1.0) ** 1.6) * scale))


def build_degree_legend_values(observed_degrees: Iterable[object]) -> List[int]:
    vals: List[int] = []
    for value in observed_degrees:
        try:
            vals.append(max(0, int(round(_safe_float(value, 0.0), 0))))
        except Exception:
            continue
    if not vals:
        return [0]

    max_obs = max(vals)
    step = 1 if max_obs <= 5 else (5 if max_obs <= 25 else 10)
    legend_max = int(math.ceil(max_obs / float(step)) * step)
    return list(range(0, legend_max + step, step))


def betweenness_label(value: object, threshold: float) -> str:
    v = _safe_float(value, 0.0)
    return f">= {threshold:.2f}" if v >= threshold else f"< {threshold:.2f}"


def get_betweenness_value(attrs: Dict, default: float = 0.0) -> float:
    if "Betweenness_norm" in attrs:
        return _safe_float(attrs.get("Betweenness_norm", default), default)
    return _safe_float(attrs.get("Betweenness", default), default)


def isa_marker_area(raw_score: float, isa_scale: float) -> float:
    scaled = max(_safe_float(raw_score, 0.0), 0.0) * max(_safe_float(isa_scale, 500.0), 1.0)
    if scaled <= 0:
        return 6.0
    return float(8.0 + (math.sqrt(scaled) * 5.5))


def edge_widths_from_weights(
    G: nx.Graph,
    scale: float = 1.0,
    min_width: float = 0.25,
    max_width: Optional[float] = None,
) -> List[float]:
    """Convert |edge weight| to visible line widths with a floor."""
    edges = list(G.edges())
    if not edges:
        return []
    widths = []
    for e in edges:
        w = abs(_safe_float(G.edges[e].get("weight", 0.0), 0.0))
        width = min_width + (w * scale)
        if max_width is not None:
            width = min(width, max_width)
        widths.append(max(min_width, width))
    return widths


def draw_edges_light(
    G: nx.Graph,
    pos: Dict,
    alpha: float = 1.0,
    edge_widths: Optional[List[float]] = None,
):
    widths = edge_widths if edge_widths is not None else edge_widths_from_weights(G, scale=1.0, min_width=0.25)
    nx.draw_networkx_edges(
        G,
        pos,
        edgelist=list(G.edges()),
        width=widths,
        edge_color=NOT_FOCUS_COLOR,
        alpha=alpha,
    )


def draw_nodes_one_by_one(G: nx.Graph, pos: Dict, color_fn, size_fn, alpha_fn=None,
                          lw_fn=None):
    """Draw heterogeneous nodes in one collection instead of one per node.

    The old implementation created thousands of PathCollections for a large
    graph.  That made multi-mode runs retain excessive renderer memory even
    after figures were closed.
    """
    nodes = list(G.nodes())
    if not nodes:
        return
    coll = nx.draw_networkx_nodes(
        G,
        pos,
        nodelist=nodes,
        node_color=[color_fn(n) for n in nodes],
        node_size=[size_fn(n) for n in nodes],
        edgecolors="black",
        linewidths=[lw_fn(n) if lw_fn else 0.25 for n in nodes],
        alpha=[alpha_fn(n) if alpha_fn else 1.0 for n in nodes],
    )
    return coll


def set_edge_artist_zorder(artists, zorder):
    """Handle NetworkX collections, arrow-patch lists, and empty edge results."""
    if artists is None:
        return
    if isinstance(artists, (list, tuple)):
        for artist in artists:
            artist.set_zorder(zorder)
    else:
        artists.set_zorder(zorder)


def draw_selected_nodes_one_by_one(G: nx.Graph, pos: Dict, nodes: List[str], color_fn, size_fn, alpha_fn=None,
                                   lw_fn=None, zorder: Optional[float] = None):
    selected = [n for n in nodes if n in G]
    if not selected:
        return None
    coll = nx.draw_networkx_nodes(
        G,
        pos,
        nodelist=selected,
        node_color=[color_fn(n) for n in selected],
        node_size=[size_fn(n) for n in selected],
        edgecolors="black",
        linewidths=[lw_fn(n) if lw_fn else 0.25 for n in selected],
        alpha=[alpha_fn(n) if alpha_fn else 1.0 for n in selected],
    )
    if zorder is not None:
        coll.set_zorder(zorder)
    return coll


def label_selected(G: nx.Graph, pos: Dict, select_nodes: List[str], text_attr: str = 'Taxon'):
    if not select_nodes and G.number_of_nodes() > 0:
        select_nodes = sorted(
            G.nodes(),
            key=lambda node: _safe_float(G.nodes[node].get("Degree", G.degree(node)), 0.0),
            reverse=True,
        )[: min(10, G.number_of_nodes())]
    max_labels = max(0, int(LABEL_CONFIG["max_labels"]))
    if max_labels and len(select_nodes) > max_labels:
        original_count = len(select_nodes)
        select_nodes = sorted(
            select_nodes,
            key=lambda node: _safe_float(G.nodes[node].get("Degree", G.degree(node)), 0.0),
            reverse=True,
        )[:max_labels]
        print(
            f"[WARN] Requested {original_count} node labels; showing the "
            f"{max_labels} highest-degree nodes to keep the figure readable and bounded."
        )
    if not _HAS_ADJUSTTEXT:
        print("[WARN] adjustText not installed; using direct labels.")
        labels = {node: str(G.nodes[node].get(text_attr, node)) for node in select_nodes}
        nx.draw_networkx_labels(G, pos, labels=labels, font_size=7, font_weight="bold")
        return
    texts = []
    for idx, n in enumerate(select_nodes):
        x, y = pos[n]
        lbl = str(G.nodes[n].get(text_attr, ""))
        dx = 0.018 if (idx % 2 == 0) else -0.018
        dy = 0.018 if ((idx // 2) % 2 == 0) else -0.018
        texts.append(plt.text(x + dx, y + dy, lbl, fontsize=8.5, weight='bold', ha='center', va='center'))
    adjust_text(
        texts,
        arrowprops=dict(arrowstyle="->", color="gray", lw=0.45),
        expand_text=(1.45, 1.65),
        expand_points=(1.35, 1.55),
        force_text=(1.2, 1.5),
        force_points=(0.9, 1.2),
        force_objects=(0.9, 1.1),
        only_move={"points": "xy", "text": "xy", "objects": "xy"},
        lim=600,
    )


def reshape_indicspecies_summary(summary_df: pd.DataFrame) -> pd.DataFrame:
    # 1) standardize IDs
    df = summary_df.rename(columns={'ASV': 'ASV_ID'}).copy()

    # 2) infer group names from the *.B columns
    b_cols = [c for c in df.columns if c.endswith('.B')]
    groups = [c[:-2] for c in b_cols]  # strip ".B"

    # 3) build rename map to create A.<group> and B.<group>
    rename_map = {}
    # A columns are the same group names w/o ".B"
    for g in groups:
        if g in df.columns:
            rename_map[g] = f"A.{g}"
        rename_map[f"{g}.B"] = f"B.{g}"

    df = df.rename(columns=rename_map)

    # 4) reshape to long; keep both A and B
    out = (
        pd.wide_to_long(
            df,
            stubnames=['A', 'B'],
            i=['ASV_ID', 'index'],
            j='Group',
            sep='.',
            suffix='.+',          # group names can include spaces/plus signs
        )
        .reset_index()[['ASV_ID', 'index', 'Group', 'A', 'B']]
        .sort_values(['ASV_ID', 'index', 'Group'], ignore_index=True)
    )
    return out

def normalize_combo(label: str) -> str:
    if not isinstance(label, str):
        return str(label)
    parts = [p.strip() for p in label.split("+")]
    parts = [p for p in parts if p]
    parts = ["Non-Cancer" if p == "Control" else p for p in parts]
    uniq: List[str] = []
    for p in parts:
        if p not in uniq:
            uniq.append(p)
    return "+".join(uniq) if uniq else label


def combo_contains_component(label: object, component: object) -> bool:
    if label is None or component is None:
        return False
    lbl = normalize_combo(str(label).strip())
    comp = normalize_combo(str(component).strip())
    if not lbl or not comp:
        return False
    lbl_parts = {p.strip() for p in lbl.split("+") if p.strip()}
    comp_parts = {p.strip() for p in comp.split("+") if p.strip()}
    if not lbl_parts or not comp_parts:
        return False
    return comp_parts.issubset(lbl_parts)


def parse_mapping(text: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if not text:
        return out
    for token in text.split(","):
        tok = token.strip()
        if not tok or "=" not in tok:
            continue
        k, v = tok.split("=", 1)
        k = k.strip()
        v = v.strip()
        if k and v:
            out[k] = v
    return out


def parse_csv_list(text: str) -> List[str]:
    if not text:
        return []
    return [t.strip() for t in str(text).split(",") if t and t.strip()]


def normalize_order_list(vals: Iterable[str]) -> List[str]:
    out: List[str] = []
    for v in vals:
        vv = normalize_combo(str(v).strip())
        if vv and vv not in out:
            out.append(vv)
    return out


def canonicalize_group_type_order(vals: Iterable[str]) -> List[str]:
    norm = normalize_order_list(vals)
    if not norm:
        return []
    known = [lbl for lbl in CANONICAL_GROUP_TYPE_ORDER if lbl in norm]
    extras = [lbl for lbl in norm if lbl not in known]
    return known + extras


def normalize_palette_keys(palette: Dict[str, str]) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for k, v in palette.items():
        kk = normalize_combo(str(k).strip())
        if kk and kk not in out:
            out[kk] = v
    return out


def infer_all_combo_label(labels: Iterable[str]) -> Optional[str]:
    """
    Infer the 'all groups' combo label as the one with the largest number of parts.
    Returns None if no combo labels are present.
    """
    norm = [normalize_combo(str(x).strip()) for x in labels if str(x).strip()]
    norm = [x for x in norm if x and x != "not_indicator"]
    if not norm:
        return None
    scored = []
    for x in norm:
        n_parts = len([p for p in x.split("+") if p.strip()])
        if n_parts > 1:
            scored.append((n_parts, x))
    if not scored:
        return None
    scored.sort(key=lambda t: (-t[0], t[1]))
    return scored[0][1]


def build_palette_from_metadata(df: pd.DataFrame, label_col: str, color_col: str, order: Optional[List[str]] = None) -> Dict[str, str]:
    ensure_cols(df, [label_col, color_col], "metadata palette")
    x = df[[label_col, color_col]].dropna().copy()
    x[label_col] = x[label_col].astype(str).str.strip().map(normalize_combo)
    x[color_col] = x[color_col].astype(str).str.strip()
    x = x[(x[label_col] != "") & (x[color_col] != "")]
    if x.empty:
        return {}

    # If multiple colors exist for a label, pick the most frequent one.
    mode_colors = (
        x.groupby(label_col)[color_col]
        .agg(lambda s: s.value_counts().index[0])
        .to_dict()
    )
    if not order:
        return mode_colors

    ordered: Dict[str, str] = {}
    for k in order:
        if k in mode_colors:
            ordered[k] = mode_colors[k]
    for k, v in mode_colors.items():
        if k not in ordered:
            ordered[k] = v
    return ordered


def palette_get(palette: Dict[str, str], label: object, default: str = NOT_FOCUS_COLOR) -> str:
    if label is None:
        return default
    raw = str(label).strip()
    if raw in palette:
        return palette[raw]
    norm = normalize_combo(raw)
    if norm in palette:
        return palette[norm]
    compact = norm.replace(" + ", "+")
    if compact in palette:
        return palette[compact]
    return default


def augment_combo_palette(palette: Dict[str, str], labels: Iterable[str]) -> Dict[str, str]:
    out = dict(palette)
    for raw_label in labels:
        label = str(raw_label).strip()
        if not label or label in out or "+" not in label:
            continue
        explicit_color = EXPLICIT_GROUP_TYPE_PALETTE.get(label)
        if explicit_color is None:
            explicit_color = EXPLICIT_GROUP_TYPE_PALETTE.get(normalize_combo(label).replace(" + ", "+"))
        if explicit_color is not None:
            out[label] = explicit_color
            out.setdefault(normalize_combo(label), explicit_color)
            out.setdefault(normalize_combo(label).replace(" + ", "+"), explicit_color)
            continue
        parts = [p.strip() for p in re.split(r"\s*\+\s*", label) if p.strip()]
        if len(parts) < 2:
            continue
        part_colors = [palette_get(out, p, "") for p in parts]
        if any(not c for c in part_colors):
            continue
        try:
            rgb = np.mean([mcolors.to_rgb(c) for c in part_colors], axis=0)
            blend = mcolors.to_hex(rgb)
            out[label] = blend
            out.setdefault(normalize_combo(label), blend)
            out.setdefault(normalize_combo(label).replace(" + ", "+"), blend)
        except Exception:
            continue
    return out


def canonicalize_group2_palette_aliases(palette: Dict[str, str]) -> Dict[str, str]:
    """Normalize palette keys without imposing study-specific label aliases."""
    return normalize_palette_keys(palette)


def ensure_palette_has_order_entries(palette: Dict[str, str], order: Iterable[str]) -> Dict[str, str]:
    out = dict(palette)
    for raw in order:
        label = normalize_combo(str(raw).strip())
        if not label or label in out:
            continue
        if label in EXPLICIT_GROUP_TYPE_PALETTE:
            out[label] = EXPLICIT_GROUP_TYPE_PALETTE[label]
    return out


def infer_index_map_from_summary(summary_df: pd.DataFrame) -> Dict[int, str]:
    # Prefer sign columns, which are guaranteed to mirror indicspecies group order.
    s_cols = [str(c).strip() for c in summary_df.columns if str(c).strip().startswith("s.")]
    if s_cols:
        groups = [c.split("s.", 1)[1].strip() for c in s_cols]
        groups = [g for g in groups if g]
        if groups:
            mapping: Dict[int, str] = {}
            n = len(groups)
            for i in range(1, (1 << n)):
                members = [groups[b] for b in range(n) if (i >> b) & 1]
                if members:
                    mapping[i] = "+".join(members) if len(members) > 1 else members[0]
            return mapping

    # Fallback: infer from available B columns.
    b_cols = [str(c).strip() for c in summary_df.columns if str(c).strip().endswith(".B")]
    groups = [c[:-2].strip() for c in b_cols if c[:-2].strip()]
    return {i + 1: g for i, g in enumerate(groups)}


def labels_from_sign_columns(summary_df: pd.DataFrame) -> Optional[pd.Series]:
    s_cols = [str(c).strip() for c in summary_df.columns if str(c).strip().startswith("s.")]
    if not s_cols:
        return None

    def _label(row: pd.Series) -> str:
        labels: List[str] = []
        for col in s_cols:
            val = row.get(col)
            try:
                keep = float(val) > 0
            except Exception:
                keep = str(val).strip().lower() in {"1", "true", "t", "yes", "y"}
            if keep:
                labels.append(col.split("s.", 1)[1].strip())
        return normalize_combo("+".join(labels)) if labels else ""

    return summary_df.apply(_label, axis=1)


def long_AB_for_group(summary_df: pd.DataFrame, index_map: Optional[Dict[int, str]] = None) -> pd.DataFrame:
    """
    Convert indicspecies *_summary.tsv into long-form and align rows by index->group mapping.
    If strict matching yields no rows, fallback to highest-B row per ASV/index.
    """
    df = reshape_indicspecies_summary(summary_df)
    if index_map is None:
        index_map = infer_index_map_from_summary(summary_df)

    def idx_to_label(x):
        try:
            if pd.isna(x):
                return None
            return index_map.get(int(x))
        except Exception:
            return None

    df["Group_mapped"] = df["index"].apply(idx_to_label)
    df["Group_norm"] = df["Group"].map(normalize_combo)
    df["Group_mapped_norm"] = df["Group_mapped"].map(lambda x: normalize_combo(x) if x is not None else None)
    out = df.loc[df["Group_norm"] == df["Group_mapped_norm"]].copy()

    if out.empty:
        print("[WARN] No direct Group/index matches found in ISA summary; using max-B fallback.")
        tmp = df.copy()
        tmp["B_num"] = pd.to_numeric(tmp["B"], errors="coerce").fillna(0.0)
        idx = tmp.groupby(["ASV_ID", "index"])["B_num"].idxmax()
        out = tmp.loc[idx].copy()

    out["A"] = pd.to_numeric(out["A"], errors="coerce").fillna(0.0)
    out["B"] = pd.to_numeric(out["B"], errors="coerce").fillna(0.0)
    out["AxB"] = (out["A"] * out["B"]).fillna(0.0)
    return out[["ASV_ID", "index", "Group", "A", "B", "AxB"]]


# ---------------------------- Plot modes -------------------------------------
def plot_degree(G: nx.Graph, pos: Dict, out_svg: str, degree_scale: float, edge_width_scale: float):
    fig, ax = figure_ax()
    # edges weighted by |weight|
    e_w = edge_widths_from_weights(G, scale=edge_width_scale, min_width=0.25)
    draw_edges_light(G, pos, alpha=0.6, edge_widths=e_w)
    degree_levels = build_degree_legend_values(G.nodes[n].get("Degree", 0.0) for n in G.nodes())
    degree_cap = degree_levels[-1]

    def size_fn(n): return degree_marker_area(min(_safe_float(G.nodes[n].get('Degree', 0), 0.0), degree_cap), degree_scale)
    def color_fn(_): return 'black'
    draw_nodes_one_by_one(G, pos, color_fn, size_fn, alpha_fn=lambda n: 0.5)

    # legend
    handles = [ax.scatter([], [], s=degree_marker_area(s, degree_scale), edgecolors='black',
                          facecolors=NOT_FOCUS_COLOR, alpha=1, label=f'{s}') for s in degree_levels]
    legend = ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(1.01, 1),
                       title="Node Degree", frameon=False, scatterpoints=1, labelspacing=1.2,
                       borderaxespad=0.0)

    finalize_network_axes(
        fig, ax, out_svg,
        title="SPIEC-EASI Network\nNode size: Degree | Edges scaled by |weight|",
        legends=[legend],
        right=0.76,
    )


def plot_abundance(
    G: nx.Graph,
    pos: Dict,
    out_svg: str,
    edge_width_scale: float = 1.0,
    abundance_min_area: float = 8.0,
    abundance_max_area: float = 420.0,
    abundance_scale_power: float = 1.6,
):
    fig, ax = figure_ax()
    e_w = edge_widths_from_weights(G, scale=edge_width_scale, min_width=0.25)
    draw_edges_light(G, pos, alpha=1.0, edge_widths=e_w)

    raw_vals = [
        _safe_float(G.nodes[n].get('mean', 0.0), 0.0)
        for n in G.nodes()
    ]
    size_mapper, size_legend_vals = build_log_size_scaler(
        raw_vals,
        min_area=abundance_min_area,
        max_area=abundance_max_area,
        scale_power=abundance_scale_power,
        n_legend=5,
    )

    def size_fn(n):
        return size_mapper(_safe_float(G.nodes[n].get('mean', 0.0), 0.0))

    def color_fn(_): return 'black'
    draw_nodes_one_by_one(G, pos, color_fn, size_fn, alpha_fn=lambda n: 0.5)

    handles = [
        ax.scatter([], [], s=size_mapper(v), edgecolors="black", facecolors=NOT_FOCUS_COLOR,
                    alpha=1, label=f"Mean abundance: {_format_legend_num(v)}")
        for v in size_legend_vals
    ]
    legend = ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(1.01, 1),
                       title="Node Attributes", frameon=False, scatterpoints=1, labelspacing=1.2,
                       borderaxespad=0.0)

    finalize_network_axes(
        fig, ax, out_svg,
        title="SPIEC-EASI Network\nNode size: Mean ASV Abundance (log-scaled)",
        legends=[legend],
        right=0.76,
    )


def plot_mag_pairing(
    G: nx.Graph,
    pos: Dict,
    out_svg: str,
    degree_scale: float,
    edge_width_scale: float,
    label: bool = False,
    title: Optional[str] = None,
):
    fig, ax = figure_ax()
    e_w = edge_widths_from_weights(G, scale=edge_width_scale, min_width=0.25)
    draw_edges_light(G, pos, alpha=0.5, edge_widths=e_w)
    degree_levels = build_degree_legend_values(
        G.nodes[n].get("Degree", 0.0) for n in G.nodes() if bool(G.nodes[n].get("has_mag_pair", False))
    )
    degree_cap = degree_levels[-1]

    def size_fn(n):
        return degree_marker_area(min(_safe_float(G.nodes[n].get("Degree", 0), 0.0), degree_cap), degree_scale)

    def color_fn(n):
        return MAG_PAIRED_COLOR if bool(G.nodes[n].get("has_mag_pair", False)) else MAG_UNPAIRED_COLOR

    def alpha_fn(n):
        return 0.95 if bool(G.nodes[n].get("has_mag_pair", False)) else 0.12

    draw_nodes_one_by_one(G, pos, color_fn, size_fn, alpha_fn=alpha_fn)

    handles = [
        mpatches.Patch(color=MAG_PAIRED_COLOR, label="Paired MAG"),
        mpatches.Patch(color=MAG_UNPAIRED_COLOR, label="No paired MAG"),
    ]
    size_handles = [
        ax.scatter([], [], s=degree_marker_area(v, degree_scale), edgecolors="black", facecolors=NOT_FOCUS_COLOR,
                   alpha=1, label=f"Degree {v}")
        for v in degree_levels
    ]
    pairing_legend = ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(1.01, 1),
                               title="ASV-to-MAG pairing", frameon=False, borderaxespad=0.0)
    legends = [pairing_legend]
    if size_handles:
        ax.add_artist(pairing_legend)
        size_legend = ax.legend(handles=size_handles, loc='upper left', bbox_to_anchor=(1.01, 0.82),
                                title="Node size", frameon=False, borderaxespad=0.0, labelspacing=1.0)
        legends.append(size_legend)

    if label:
        to_label = [n for n in G.nodes() if bool(G.nodes[n].get("has_mag_pair", False))]
        label_selected(G, pos, to_label, text_attr='Taxon')

    finalize_network_axes(
        fig, ax, out_svg,
        title=title or "SPIEC-EASI Network\nNode color: ASV with paired MAG | Node size: Degree",
        legends=legends,
        right=0.76,
    )


def plot_mag_pairing_taxonomy(
    G: nx.Graph,
    pos: Dict,
    out_svg: str,
    degree_scale: float,
    edge_width_scale: float,
    label: bool = False,
    title: Optional[str] = None,
):
    fig, ax = figure_ax()
    e_w = edge_widths_from_weights(G, scale=edge_width_scale, min_width=0.25)
    draw_edges_light(G, pos, alpha=0.5, edge_widths=e_w)
    degree_levels = build_degree_legend_values(
        G.nodes[n].get("Degree", 0.0) for n in G.nodes() if bool(G.nodes[n].get("has_mag_pair", False))
    )
    degree_cap = degree_levels[-1]

    paired_labels = sorted({
        str(G.nodes[n].get("mag_taxonomy_label", "")).strip()
        for n in G.nodes()
        if bool(G.nodes[n].get("has_mag_pair", False))
        and str(G.nodes[n].get("mag_taxonomy_label", "")).strip()
    })
    if paired_labels:
        colors = sns.color_palette("husl", n_colors=max(3, len(paired_labels)))
        taxonomy_palette = {label: mcolors.to_hex(colors[i % len(colors)]) for i, label in enumerate(paired_labels)}
    else:
        taxonomy_palette = {}

    def size_fn(n):
        return degree_marker_area(min(_safe_float(G.nodes[n].get("Degree", 0), 0.0), degree_cap), degree_scale)

    def color_fn(n):
        if not bool(G.nodes[n].get("has_mag_pair", False)):
            return MAG_UNPAIRED_COLOR
        label = str(G.nodes[n].get("mag_taxonomy_label", "")).strip()
        return taxonomy_palette.get(label, MAG_PAIRED_COLOR)

    def alpha_fn(n):
        return 0.95 if bool(G.nodes[n].get("has_mag_pair", False)) else 0.12

    draw_nodes_one_by_one(G, pos, color_fn, size_fn, alpha_fn=alpha_fn)

    handles = [mpatches.Patch(color=MAG_UNPAIRED_COLOR, label="No paired MAG")]
    handles.extend([
        mpatches.Patch(color=taxonomy_palette[label], label=label)
        for label in paired_labels
    ])
    size_handles = [
        ax.scatter([], [], s=degree_marker_area(v, degree_scale), edgecolors="black", facecolors=NOT_FOCUS_COLOR,
                   alpha=1, label=f"Degree {v}")
        for v in degree_levels
    ]
    legend_ncol = 1 if len(handles) <= 16 else 2 if len(handles) <= 32 else 3
    tax_legend = ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(1.01, 1),
                           title="MAG phylum", frameon=False, borderaxespad=0.0,
                           labelspacing=0.8, ncol=legend_ncol)
    legends = [tax_legend]
    if size_handles:
        ax.add_artist(tax_legend)
        size_legend = ax.legend(handles=size_handles, loc='upper left', bbox_to_anchor=(1.01, 0.55),
                                title="Node size", frameon=False, borderaxespad=0.0, labelspacing=1.0)
        legends.append(size_legend)

    if label:
        to_label = [n for n in G.nodes() if bool(G.nodes[n].get("has_mag_pair", False))]
        label_selected(G, pos, to_label, text_attr='Taxon')

    finalize_network_axes(
        fig, ax, out_svg,
        title=title or "SPIEC-EASI Network\nNode color: paired MAG phylum | Node size: Degree",
        legends=legends,
        right=0.76,
    )


def plot_betweenness(
    G: nx.Graph,
    pos: Dict,
    out_svg: str,
    degree_scale: float,
    edge_width_scale: float,
    *,
    betweenness_threshold: float = 0.05,
    require_mag_pair: bool = False,
    label: bool = False,
    title: Optional[str] = None,
):
    fig, ax = figure_ax()
    e_w = edge_widths_from_weights(G, scale=edge_width_scale, min_width=0.25)
    draw_edges_light(G, pos, alpha=0.5, edge_widths=e_w)

    degree_levels = build_degree_legend_values(
        G.nodes[n].get("Degree", 0.0)
        for n in G.nodes()
        if (not require_mag_pair) or bool(G.nodes[n].get("has_mag_pair", False))
    )
    degree_cap = degree_levels[-1]

    def in_focus(n: str) -> bool:
        return (not require_mag_pair) or bool(G.nodes[n].get("has_mag_pair", False))

    def size_fn(n):
        return degree_marker_area(min(_safe_float(G.nodes[n].get("Degree", 0.0), 0.0), degree_cap), degree_scale)

    def node_class(n: str) -> str:
        if not in_focus(n):
            return "background"
        return "high" if get_betweenness_value(G.nodes[n], 0.0) >= betweenness_threshold else "low"

    # Draw low-priority nodes first and high-betweenness nodes last so the
    # connector nodes remain visible when points overlap.
    draw_order = ["background", "low"]
    class_style = {
        "background": dict(color=MAG_UNPAIRED_COLOR, alpha=0.10, zorder=2),
        "low": dict(color=BETWEENNESS_LOW_COLOR, alpha=0.95, zorder=3),
        "high": dict(color=BETWEENNESS_HIGH_COLOR, alpha=0.98, zorder=6),
    }
    for klass in draw_order:
        nodes = [n for n in G.nodes() if node_class(n) == klass]
        if not nodes:
            continue
        nodes = sorted(nodes, key=lambda n: get_betweenness_value(G.nodes[n], 0.0))
        draw_selected_nodes_one_by_one(
            G,
            pos,
            nodes,
            color_fn=lambda _n, c=class_style[klass]["color"]: c,
            size_fn=size_fn,
            alpha_fn=lambda _n, a=class_style[klass]["alpha"]: a,
            zorder=class_style[klass]["zorder"],
        )

    # Final dedicated pass for the highlighted nodes only.
    high_nodes = sorted(
        [n for n in G.nodes() if node_class(n) == "high"],
        key=lambda n: get_betweenness_value(G.nodes[n], 0.0),
    )
    if high_nodes:
        # White halo first so highlighted bridge nodes stay readable over edges
        # and lower-priority points.
        draw_selected_nodes_one_by_one(
            G,
            pos,
            high_nodes,
            color_fn=lambda _n: "white",
            size_fn=lambda n: size_fn(n) * 1.18,
            alpha_fn=lambda _n: 0.98,
            lw_fn=lambda _n: 0.0,
            zorder=5,
        )
        draw_selected_nodes_one_by_one(
            G,
            pos,
            high_nodes,
            color_fn=lambda _n: class_style["high"]["color"],
            size_fn=size_fn,
            alpha_fn=lambda _n: class_style["high"]["alpha"],
            lw_fn=lambda _n: 0.45,
            zorder=class_style["high"]["zorder"],
        )

    class_handles = [
        mpatches.Patch(color=BETWEENNESS_HIGH_COLOR, label=f"Betweenness >= {betweenness_threshold:.2f}"),
        mpatches.Patch(color=BETWEENNESS_LOW_COLOR, label=f"Betweenness < {betweenness_threshold:.2f}"),
    ]
    if require_mag_pair:
        class_handles.append(mpatches.Patch(color=MAG_UNPAIRED_COLOR, label="No paired MAG"))
    class_legend = ax.legend(
        handles=class_handles,
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        title="Normalized Betweenness",
        frameon=False,
        borderaxespad=0.0,
    )
    legends = [class_legend]

    size_handles = [
        ax.scatter([], [], s=degree_marker_area(v, degree_scale), edgecolors="black", facecolors=NOT_FOCUS_COLOR,
                   alpha=1, label=f"Degree {v}")
        for v in degree_levels
    ]
    if size_handles:
        ax.add_artist(class_legend)
        size_legend = ax.legend(
            handles=size_handles,
            loc="upper left",
            bbox_to_anchor=(1.01, 0.72),
            title="Node size",
            frameon=False,
            borderaxespad=0.0,
            labelspacing=1.0,
        )
        legends.append(size_legend)

    if label:
        to_label = [
            n for n in G.nodes()
            if node_class(n) == "high"
        ]
        label_selected(G, pos, to_label, text_attr="Taxon")

    finalize_network_axes(
        fig, ax, out_svg,
        title=title or "SPIEC-EASI Network\nNode color: Normalized betweenness class | Node size: Degree",
        legends=legends,
        right=0.76,
    )

def plot_group_isa(
    G: nx.Graph,
    pos: Dict,
    out_svg: str,
    palette: Dict[str, str],
    *,
    color_attr: str,
    size_attr: str,
    label_attr: Optional[str] = None,
    isa_scale: float = 500,
    edge_width_scale: float = 1.0,
    label: bool = False,
    title: Optional[str] = None,
    legend_title: str = "Group",
    legend_order: Optional[List[str]] = None,
    focus_label: str = "",
    all_combo_label: str = "",
    require_mag_pair: bool = False,
):
    nodes = list(G.nodes())
    if not nodes:
        print("[!] plot_group_isa: graph has no nodes; skipping.")
        return

    node_colors = []
    node_components: List[List[str]] = []
    node_sizes  = []
    visible_labels: List[str] = []
    visible_scores: List[float] = []
    palette_norm = normalize_palette_keys(palette)
    for n in nodes:
        d = G.nodes[n]
        node_lbl = normalize_combo(str(d.get(label_attr, "")).strip()) if label_attr else ""
        keep_focus = True
        if focus_label and label_attr:
            keep_focus = combo_contains_component(node_lbl, focus_label)
        if require_mag_pair:
            keep_focus = keep_focus and bool(d.get("has_mag_pair", False))

        raw_score = _safe_float(d.get(size_attr, 0.0), 0.0)
        has_signal = keep_focus and node_lbl and raw_score > 0
        c = palette_get(palette_norm, node_lbl, NOT_FOCUS_COLOR) if has_signal else NOT_FOCUS_COLOR
        node_colors.append(c)
        components = [
            part for part in re.split(r"\s*\+\s*", node_lbl)
            if part and palette_get(palette_norm, part, "") and has_signal
        ]
        node_components.append(components)

        s = isa_marker_area(raw_score, isa_scale)
        if not has_signal:
            s = 0.0
        if not np.isfinite(s) or s <= 0:
            s = 1.0
        node_sizes.append(s)
        if has_signal:
            # Multi-group indicators are drawn as pies, so their legend must
            # describe the component slices rather than the combined label.
            # Tracking only labels that reach the plot also prevents globally
            # configured but absent groups from appearing in the legend.
            visible_labels.extend(components if components else [node_lbl])
            visible_scores.append(raw_score)

    fig, ax = figure_ax((17, 14))

    # Edge widths are scaled by |weight| with a floor for visibility.
    e_w = edge_widths_from_weights(G, scale=edge_width_scale, min_width=0.25)
    nx.draw_networkx_edges(G, pos, edgelist=list(G.edges()), width=e_w, edge_color=NOT_FOCUS_COLOR, alpha=0.8, ax=ax)

    # Single-group indicators use their exact legend color. Multigroup
    # indicators are segmented into exact component colors instead of averaging
    # them into misleading earthy blends.
    single_nodes = [
        node for node, parts in zip(nodes, node_components)
        if len(parts) <= 1
    ]
    single_idx = [nodes.index(node) for node in single_nodes]
    nx.draw_networkx_nodes(
        G, pos,
        nodelist=single_nodes,
        node_color=[node_colors[idx] for idx in single_idx],
        node_size=[node_sizes[idx] for idx in single_idx],
        edgecolors="black", linewidths=0.25, alpha=1.0, ax=ax,
    )
    for node, size, parts in zip(nodes, node_sizes, node_components):
        if len(parts) <= 1:
            continue
        x, y = pos[node]
        angles = np.linspace(0.0, 360.0, len(parts) + 1)
        for index, part in enumerate(parts):
            radians = np.deg2rad(np.linspace(angles[index], angles[index + 1], 24))
            vertices = np.vstack([
                [0.0, 0.0],
                np.column_stack([np.cos(radians), np.sin(radians)]),
                [0.0, 0.0],
            ])
            codes = [MplPath.MOVETO] + [MplPath.LINETO] * len(radians) + [MplPath.CLOSEPOLY]
            ax.scatter(
                [x], [y], s=size, marker=MplPath(vertices, codes),
                facecolor=palette_get(palette_norm, part, NOT_FOCUS_COLOR),
                edgecolor="none", linewidth=0, zorder=3,
            )
        ax.scatter(
            [x], [y], s=size, marker="o", facecolor="none",
            edgecolor="black", linewidth=0.25, zorder=4,
        )

    if label:
        for n in nodes:
            if not is_not_focus_color(node_colors[nodes.index(n)]):
                x, y = pos[n]
                ax.text(x, y, G.nodes[n].get("Taxon", n),
                        fontsize=9, fontweight="bold",
                        ha="center", va="center")

    ordered_labels = ordered_present_labels(visible_labels, legend_order)
    if focus_label:
        ordered_labels = [
            lbl for lbl in ordered_labels
            if combo_contains_component(lbl, focus_label)
        ]
        if legend_order and "not_indicator" in normalize_order_list(legend_order) and "not_indicator" not in ordered_labels:
            ordered_labels.append("not_indicator")
    class_handles = [mpatches.Patch(color=palette_get(palette_norm, lbl, NOT_FOCUS_COLOR), label=lbl) for lbl in ordered_labels]
    if not class_handles:
        class_handles = [mpatches.Patch(color=NOT_FOCUS_COLOR, label="No significant ISA ASVs")]
    if require_mag_pair:
        class_handles.insert(0, mpatches.Patch(color=MAG_UNPAIRED_COLOR, label="No paired MAG"))

    size_legend_vals = build_isa_size_legend_values(visible_scores)
    size_handles = [ax.scatter([], [], s=isa_marker_area(v, isa_scale),
                                edgecolors="black", facecolors=NOT_FOCUS_COLOR, alpha=1,
                                label=f"ISA: {_format_legend_float(v, 2)}") for v in size_legend_vals]

    legend_handles = class_handles + size_handles
    legend = ax.legend(
        legend_handles,
        [h.get_label() for h in legend_handles],
        loc="upper left",
        bbox_to_anchor=(1.01, 1),
        frameon=False,
        title=legend_title,
        borderaxespad=0.0,
        scatterpoints=1,
        labelspacing=1.0,
    )

    finalize_network_axes(
        fig, ax, out_svg,
        title=title or "SPIEC-EASI Network",
        legends=[legend],
        right=0.72,
    )

def plot_type_isa(G, pos, out_svg, type_palette, isa_scale=500, label=False, title=None):
    plot_group_isa(
        G,
        pos,
        out_svg,
        type_palette,
        color_attr="type_color",
        size_attr="AxB_group1",
        isa_scale=isa_scale,
        label=label,
        title=title,
        legend_title="Group1 ISA",
    )


def plot_status_isa(G: nx.Graph, pos: Dict, out_svg: str, status_palette: Dict[str, str],
                    isa_scale: float, label: bool = False):
    plot_group_isa(
        G,
        pos,
        out_svg,
        status_palette,
        color_attr="status_color",
        size_attr="AxB_group2",
        isa_scale=isa_scale,
        label=label,
        title="SPIEC-EASI Network\nNode color: Group2 ISA | Node size: Indicator Species Strength",
        legend_title="Group2 ISA",
    )


def plot_phylum(
    G: nx.Graph,
    pos: Dict,
    out_svg: str,
    phylum_palette: Dict[str, str],
    size_attr: str,
    size_label: str,
    size_scale: float,
    edge_width_scale: float = 1.0,
    label: bool = False,
    abundance_min_area: float = 8.0,
    abundance_max_area: float = 420.0,
    abundance_scale_power: float = 1.6,
):
    fig, ax = figure_ax((17, 14))
    e_w = edge_widths_from_weights(G, scale=edge_width_scale, min_width=0.25)
    draw_edges_light(G, pos, alpha=1.0, edge_widths=e_w)

    def color_fn(n):
        if size_attr == "AxB_group1" or str(size_attr).startswith("isa_score__"):
            isa_val = _safe_float(G.nodes[n].get("AxB_group1", 0.0), 0.0)
            if str(size_attr).startswith("isa_score__"):
                isa_val = _safe_float(G.nodes[n].get(size_attr, 0.0), 0.0)
            if isa_val <= 0:
                return NOT_FOCUS_COLOR
        p = str(G.nodes[n].get('Phylum', '')).strip()
        return phylum_palette.get(p, NOT_FOCUS_COLOR)

    raw_vals_scaled = [
        _safe_float(G.nodes[n].get(size_attr, 0.0), 0.0) * size_scale
        for n in G.nodes()
    ]
    abund_mapper = None
    abund_legend_vals: List[float] = []
    if size_attr in {"mean", "median"}:
        abund_mapper, abund_legend_vals = build_log_size_scaler(
            raw_vals_scaled,
            min_area=abundance_min_area,
            max_area=abundance_max_area,
            scale_power=abundance_scale_power,
            n_legend=5,
        )

    def _size_to_area(raw_value: float) -> float:
        """Map raw values to marker area in a stable way for plotting."""
        if not np.isfinite(raw_value) or raw_value <= 0:
            return 1.0
        if size_attr in {"mean", "median"}:
            # Abundance ranges are wide; use log scaling for readability.
            return float(abund_mapper(raw_value)) if abund_mapper is not None else 1.0
        base_score = raw_value / max(float(size_scale), 1e-9)
        return isa_marker_area(base_score, size_scale)

    def size_fn(n):
        raw = _safe_float(G.nodes[n].get(size_attr, 0.0), 0.0) * size_scale
        return _size_to_area(raw)

    def alpha_fn(n): return 0.5 if is_not_focus_color(color_fn(n)) else 1.0

    draw_nodes_one_by_one(G, pos, color_fn, size_fn, alpha_fn=alpha_fn)

    if label:
        for n in G.nodes():
            if is_not_focus_color(color_fn(n)):
                continue
            x, y = pos[n]
            ax.text(x, y, str(G.nodes[n].get("Taxon", n)),
                    fontsize=8.5, fontweight="bold",
                    ha="center", va="center")

    # Build deterministic phylum legend from labels present in this graph.
    present_phyla = sorted({
        str(G.nodes[n].get('Phylum', '')).strip()
        for n in G.nodes()
        if str(G.nodes[n].get('Phylum', '')).strip()
        and not is_not_focus_color(color_fn(n))
        and phylum_palette.get(str(G.nodes[n].get('Phylum', '')).strip()) is not None
    })
    color_patches = [
        mpatches.Patch(color=phylum_palette[p], label=p)
        for p in present_phyla
    ]
    has_other = any(is_not_focus_color(color_fn(n)) for n in G.nodes())
    if has_other:
        color_patches.append(mpatches.Patch(color=NOT_FOCUS_COLOR, label="Unassigned/Other"))

    if size_attr in {"mean", "median"}:
        size_vals = abund_legend_vals if abund_legend_vals else [1.0]
        size_handles = [
            plt.scatter(
                [], [], s=_size_to_area(v),
                edgecolors="black", facecolors=NOT_FOCUS_COLOR, alpha=1,
                label=f"{size_label}: {_format_legend_num(v)}"
            )
            for v in size_vals
        ]
    else:
        size_vals = [0.1, 0.25, 0.5, 0.75, 1.0]
        size_handles = size_legend_handles(size_vals, size_label + ":", size_scale)

    phylum_ncol = 1 if len(color_patches) <= 14 else 2 if len(color_patches) <= 28 else 3
    leg_phylum = ax.legend(
        handles=color_patches,
        loc='upper left',
        bbox_to_anchor=(1.01, 1.0),
        title="Phylum",
        frameon=False,
        labelspacing=0.8,
        ncol=phylum_ncol
    )
    ax.add_artist(leg_phylum)

    size_legend = ax.legend(
        handles=size_handles,
        loc='upper left',
        bbox_to_anchor=(1.01, 0.45),
        title="Node size",
        frameon=False,
        scatterpoints=1,
        labelspacing=1.2
    )
    if size_attr == "mean":
        titletail = "Mean ASV Abundance (log-scaled)"
    elif size_attr == "median":
        titletail = "Median ASV Abundance (log-scaled)"
    elif size_attr == "AxB_group1":
        titletail = "ISA Strength (sig ISA only colored)"
    else:
        titletail = "Indicator Species Strength"
    finalize_network_axes(
        fig, ax, out_svg,
        title=f"SPIEC-EASI Network\nNode color: Phylum | Node size: {titletail}",
        legends=[leg_phylum, size_legend],
        right=0.72,
    )


def infer_group_name(path: str, fallback: str) -> str:
    stem = Path(path).name
    stem = re.sub(r"_indicator_species(_DULEG)?_summary\.tsv$", "", stem)
    stem = re.sub(r"[^0-9A-Za-z]+", "_", stem).strip("_")
    return stem if stem else fallback


def isa_summary_is_duleg(path: str) -> bool:
    return bool(re.search(r"_indicator_species_DULEG_summary\.tsv$", Path(path).name))


def slugify_group_name(name: str, fallback: str = "group") -> str:
    slug = re.sub(r"[^0-9A-Za-z._-]+", "_", str(name or "").strip()).strip("_")
    return slug or fallback


def collect_isa_summary_paths(
    isa_group_cols: List[str],
    group1_summary_path: str,
    group2_summary_path: str,
    group1_name: str,
    group2_name: str,
    summary_mode: Optional[str] = None,
) -> Dict[str, str]:
    variant_by_name: Dict[str, Dict[bool, str]] = {}
    for path in sorted(Path(".").glob("*_indicator_species*_summary.tsv")):
        name = infer_group_name(str(path), "")
        if name.startswith("stratified_"):
            try:
                preview = pd.read_csv(path, sep="\t", usecols=lambda col: col == "stratified_within_value")
            except Exception:
                preview = pd.DataFrame()
            # Pooled stratified tables combine several outer-state values and
            # cannot be represented honestly by one network color overlay.
            # Per-state files are retained as distinct overlays instead.
            if ("stratified_within_value" in preview and
                    preview["stratified_within_value"].dropna().astype(str).nunique() > 1):
                continue
        if name:
            variant_by_name.setdefault(name, {})[isa_summary_is_duleg(str(path))] = str(path.resolve())

    summary_by_name: Dict[str, str] = {}
    prefer_duleg = True if summary_mode == "duleg" else False if summary_mode == "default" else None
    for name, variants in variant_by_name.items():
        if prefer_duleg is True:
            chosen = variants.get(True) or variants.get(False)
        elif prefer_duleg is False:
            chosen = variants.get(False) or variants.get(True)
        else:
            chosen = variants.get(False) or variants.get(True)
        if chosen:
            summary_by_name[name] = chosen

    if group1_summary_path:
        summary_by_name[group1_name] = group1_summary_path
    if group2_summary_path:
        summary_by_name[group2_name] = group2_summary_path

    ordered = []
    for name in isa_group_cols:
        if name in summary_by_name and name not in ordered:
            ordered.append(name)
    for name in sorted(summary_by_name):
        if name not in ordered:
            ordered.append(name)
    return {name: summary_by_name[name] for name in ordered}


def isa_mode_variants(group_index: int) -> Dict[str, str]:
    prefix = f"group{group_index}"
    return {
        "isa": f"{prefix}_isa",
        "isa_labeled": f"{prefix}_isa_labeled",
        "isa_all": f"{prefix}_isa_all",
        "isa_all_labeled": f"{prefix}_isa_all_labeled",
        "isa_mag": f"{prefix}_isa_mag",
        "isa_mag_labeled": f"{prefix}_isa_mag_labeled",
        "isa_mag_all": f"{prefix}_isa_mag_all",
        "isa_mag_all_labeled": f"{prefix}_isa_mag_all_labeled",
        "isa_focus": f"{prefix}_isa_focus",
        "isa_focus_labeled": f"{prefix}_isa_focus_labeled",
        "isa_focus_all": f"{prefix}_isa_focus_all",
        "isa_focus_all_labeled": f"{prefix}_isa_focus_all_labeled",
    }


def auto_palette(labels: Iterable[str]) -> Dict[str, str]:
    vals = sorted({str(x).strip() for x in labels if pd.notna(x) and str(x).strip()})
    if not vals:
        return {}
    colors = sns.color_palette("tab20", n_colors=max(3, len(vals)))
    return {vals[i]: mcolors.to_hex(colors[i % len(colors)]) for i in range(len(vals))}


def load_modules_table(path: Optional[str], variant_name: str) -> pd.DataFrame:
    if not path:
        return pd.DataFrame(columns=["Taxon", "module_id", "module_label", "node_stability", "graph_variant"])
    if not os.path.exists(path):
        print(f"[WARN] Module table not found for {variant_name}: {path}")
        return pd.DataFrame(columns=["Taxon", "module_id", "module_label", "node_stability", "graph_variant"])
    df = load_table(path, sep="\t")
    if "Taxon" not in df.columns and "ASV_ID" in df.columns:
        df = df.rename(columns={"ASV_ID": "Taxon"})
    if "Taxon" not in df.columns:
        print(f"[WARN] Module table missing Taxon/ASV_ID for {variant_name}: {path}")
        return pd.DataFrame(columns=["Taxon", "module_id", "module_label", "node_stability", "graph_variant"])
    if "module_id" not in df.columns:
        print(f"[WARN] Module table missing module_id for {variant_name}: {path}")
        return pd.DataFrame(columns=["Taxon", "module_id", "module_label", "node_stability", "graph_variant"])
    if "module_label" not in df.columns:
        df["module_label"] = df["module_id"].apply(lambda x: f"M{int(x)}" if pd.notna(x) else "unassigned")
    if "node_stability" not in df.columns:
        df["node_stability"] = np.nan
    if "graph_variant" not in df.columns:
        df["graph_variant"] = variant_name
    keep = ["Taxon", "module_id", "module_label", "node_stability", "graph_variant"]
    return df[keep].drop_duplicates(subset=["Taxon"])


def load_asv_mag_pairing(path: Optional[str]) -> pd.DataFrame:
    if not path:
        return pd.DataFrame(columns=["ASV_ID", "has_mag_pair", "mag_pair_status", "best_genome_id", "mag_taxonomy_label"])
    if not os.path.exists(path):
        print(f"[WARN] ASV-MAG pairing table not found: {path}")
        return pd.DataFrame(columns=["ASV_ID", "has_mag_pair", "mag_pair_status", "best_genome_id", "mag_taxonomy_label"])
    df = load_table(path, sep="\t")
    if "ASV_ID" not in df.columns:
        print(f"[WARN] ASV-MAG pairing table missing ASV_ID: {path}")
        return pd.DataFrame(columns=["ASV_ID", "has_mag_pair", "mag_pair_status", "best_genome_id", "mag_taxonomy_label"])
    out = df.copy()
    if "best_genome_id" not in out.columns and "genome_id" in out.columns:
        out["best_genome_id"] = out["genome_id"]
    if "pairing_status" not in out.columns:
        if "best_genome_id" in out.columns:
            out["pairing_status"] = np.where(out["best_genome_id"].notna(), "paired", "unpaired")
        else:
            out["pairing_status"] = "unpaired"
    out["ASV_ID"] = out["ASV_ID"].astype(str).str.strip().str.split(";", n=1).str[0]
    out["pairing_status"] = out["pairing_status"].astype(str).str.strip()
    if "best_genome_id" not in out.columns:
        out["best_genome_id"] = pd.NA
    if "mag_phylum" not in out.columns:
        out["mag_phylum"] = pd.NA
    has_mimag_tier = "mag_mimag_tier" in out.columns
    if not has_mimag_tier:
        out["mag_mimag_tier"] = pd.NA
    out["mag_mimag_tier"] = (
        out["mag_mimag_tier"]
        .fillna("")
        .astype(str)
        .str.strip()
        .str.lower()
        .replace({"med": "medium"})
    )

    def summarize_asv(grp: pd.DataFrame) -> pd.Series:
        paired = grp.loc[grp["pairing_status"].ne("unpaired")].copy()
        if has_mimag_tier and paired["mag_mimag_tier"].ne("").any():
            paired = paired.loc[paired["mag_mimag_tier"].isin(ALLOWED_MAG_MIMAG_TIERS)].copy()
        genomes = sorted({str(x) for x in paired["best_genome_id"].dropna() if str(x).strip()})
        phyla = sorted({str(x).strip() for x in paired["mag_phylum"].dropna() if str(x).strip()})
        has_mag_pair = len(genomes) > 0
        if not has_mag_pair:
            status = "unpaired"
            best = pd.NA
            taxonomy = pd.NA
        elif len(genomes) == 1:
            status = "paired_unique"
            best = genomes[0]
            if len(phyla) == 1:
                taxonomy = phyla[0]
            elif len(phyla) > 1:
                taxonomy = MAG_MIXED_TAXONOMY
            else:
                taxonomy = MAG_UNKNOWN_TAXONOMY
        else:
            status = "paired_ambiguous"
            best = genomes[0]
            if len(phyla) == 1:
                taxonomy = phyla[0]
            elif len(phyla) > 1:
                taxonomy = MAG_MIXED_TAXONOMY
            else:
                taxonomy = MAG_UNKNOWN_TAXONOMY
        return pd.Series({
            "has_mag_pair": has_mag_pair,
            "mag_pair_status": status,
            "best_genome_id": best,
            "mag_taxonomy_label": taxonomy,
        })

    out = (
        out.groupby("ASV_ID", dropna=False, sort=False)
        .apply(summarize_asv)
        .reset_index()
    )
    return out


def build_module_palette(module_labels: Iterable[str]) -> Dict[str, str]:
    vals = sorted({str(x).strip() for x in module_labels if pd.notna(x) and str(x).strip()}, key=natural_sort_key)
    if not vals:
        return {}
    colors = sns.color_palette("tab20", n_colors=max(3, len(vals)))
    palette = {vals[i]: mcolors.to_hex(colors[i % len(colors)]) for i in range(len(vals))}
    for k in list(palette.keys()):
        if str(k).strip().lower() in {"unassigned", "mna", "na", "none"}:
            palette[k] = NOT_FOCUS_COLOR
    return palette


def select_best_modules(
    modules_df: pd.DataFrame,
    min_size: int,
    min_stability: float,
    graph: Optional[nx.Graph] = None,
    top_n: int = 8,
    ensure_one: bool = True
) -> Tuple[set, pd.DataFrame]:
    """Rank stable modules by topology and retain only the strongest few."""
    empty_stats = pd.DataFrame(columns=[
        "module_label", "n_nodes", "mean_node_stability", "internal_edges",
        "internal_weight", "edge_density", "conductance",
        "weighted_modularity_contribution", "module_quality_score",
        "quality_rank", "passes_hard_filters", "is_best",
    ])
    if modules_df is None or modules_df.empty:
        return set(), empty_stats

    x = modules_df.copy()
    x["module_label"] = x["module_label"].astype(str).str.strip()
    x["node_stability"] = pd.to_numeric(x["node_stability"], errors="coerce")

    # Never treat these as true biological modules.
    bad_labels = {"", "unassigned", "mna", "na", "none"}
    x = x[~x["module_label"].str.lower().isin(bad_labels)].copy()
    if x.empty:
        return set(), empty_stats

    stats = (
        x.groupby("module_label", dropna=False)
         .agg(
             n_nodes=("Taxon", "nunique"),
             mean_node_stability=("node_stability", "mean")
         )
         .reset_index()
    )
    for column in (
        "internal_edges", "internal_weight", "edge_density", "conductance",
        "weighted_modularity_contribution",
    ):
        stats[column] = np.nan

    if graph is not None and graph.number_of_nodes() > 0:
        taxon_to_node = {
            str(data.get("Taxon", data.get("name", node))): node
            for node, data in graph.nodes(data=True)
        }
        total_weight = sum(
            abs(_safe_float(data.get("weight", 1.0), 1.0))
            for _, _, data in graph.edges(data=True)
        )
        total_weight = max(float(total_weight), 1e-12)
        label_members = x.groupby("module_label")["Taxon"].apply(
            lambda values: {
                taxon_to_node[str(value)] for value in values
                if str(value) in taxon_to_node
            }
        )
        for idx, row in stats.iterrows():
            members = label_members.get(row["module_label"], set())
            subgraph = graph.subgraph(members)
            internal_weight = sum(
                abs(_safe_float(data.get("weight", 1.0), 1.0))
                for _, _, data in subgraph.edges(data=True)
            )
            cut_weight = sum(
                abs(_safe_float(data.get("weight", 1.0), 1.0))
                for node in members
                for neighbor, data in graph[node].items()
                if neighbor not in members
            )
            degree_weight = (2.0 * internal_weight) + cut_weight
            possible_edges = len(members) * (len(members) - 1) / 2.0
            stats.loc[idx, "internal_edges"] = subgraph.number_of_edges()
            stats.loc[idx, "internal_weight"] = internal_weight
            stats.loc[idx, "edge_density"] = (
                subgraph.number_of_edges() / possible_edges if possible_edges > 0 else 0.0
            )
            stats.loc[idx, "conductance"] = (
                cut_weight / degree_weight if degree_weight > 0 else 1.0
            )
            stats.loc[idx, "weighted_modularity_contribution"] = (
                (internal_weight / total_weight)
                - ((degree_weight / (2.0 * total_weight)) ** 2)
            )

    stats["passes_hard_filters"] = (
        (stats["n_nodes"] >= int(min_size))
        & (stats["mean_node_stability"].fillna(0.0) >= float(min_stability))
    )
    topology = stats["weighted_modularity_contribution"].fillna(0.0).clip(lower=0.0)
    separation = (1.0 - stats["conductance"].fillna(1.0)).clip(lower=0.0)
    stats["module_quality_score"] = (
        stats["mean_node_stability"].fillna(0.0)
        * np.log1p(stats["n_nodes"].fillna(0.0))
        * (1.0 + topology)
        * (0.5 + 0.5 * separation)
    )
    stats = stats.sort_values(
        ["passes_hard_filters", "module_quality_score", "mean_node_stability", "n_nodes", "module_label"],
        ascending=[False, False, False, False, True],
    ).reset_index(drop=True)
    stats["quality_rank"] = np.arange(1, len(stats) + 1)
    stats["is_best"] = False
    eligible = stats.index[stats["passes_hard_filters"]].tolist()[:max(1, int(top_n))]
    stats.loc[eligible, "is_best"] = True
    best = set(stats.loc[stats["is_best"], "module_label"].astype(str).tolist())

    if ensure_one and not best and not stats.empty:
        # If nothing passes thresholds, keep the single strongest module for visibility.
        pick = (
            stats.sort_values(
                by=["mean_node_stability", "n_nodes", "module_label"],
                ascending=[False, False, True]
            )
            .iloc[0]["module_label"]
        )
        best = {str(pick)}
        stats.loc[stats["module_label"] == pick, "is_best"] = True

    return best, stats


def _first_existing_col(df: pd.DataFrame, candidates: Iterable[str]) -> Optional[str]:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def _truthy_series(s: pd.Series) -> pd.Series:
    if s is None:
        return pd.Series([], dtype=bool)
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False)
    norm = s.astype(str).str.strip().str.lower()
    return norm.isin({"1", "true", "t", "yes", "y"})


def derive_significant_isa_hits(
    summary_df: pd.DataFrame,
    group_long_df: pd.DataFrame,
    *,
    min_stat: float = 0.25,
    max_q: float = 0.05,
    drop_all_combo: bool = True,
) -> pd.DataFrame:
    """
    Build per-ASV significant ISA calls with label/score, using index->label mapping.
    Returns columns: ASV_ID, isa_label, isa_score.
    """
    empty = pd.DataFrame(columns=["ASV_ID", "isa_label", "isa_score"])
    if summary_df is None or summary_df.empty:
        return empty

    asv_col = _first_existing_col(summary_df, ["ASV_ID", "ASV"])
    if not asv_col or "index" not in summary_df.columns:
        return empty

    stat_col = _first_existing_col(summary_df, ["stat", "Stat", "STAT"])
    q_col = _first_existing_col(summary_df, ["q.value", "q_value", "qvalue", "qval", "q"])
    sig_col = _first_existing_col(summary_df, ["significant", "is_significant", "sig"])

    x = summary_df[[asv_col, "index"]].copy()
    x = x.rename(columns={asv_col: "ASV_ID"})
    x["ASV_ID"] = x["ASV_ID"].astype(str).str.strip()
    x["index"] = pd.to_numeric(x["index"], errors="coerce").astype("Int64")
    sign_labels = labels_from_sign_columns(summary_df)
    if sign_labels is not None:
        x["isa_label"] = sign_labels
    else:
        index_map = infer_index_map_from_summary(summary_df)
        x["isa_label"] = x["index"].map(
            lambda v: normalize_combo(index_map.get(int(v))) if pd.notna(v) else None
        )

    mask = x["isa_label"].notna() & (x["isa_label"].astype(str).str.strip() != "")
    if stat_col:
        x["__stat"] = pd.to_numeric(summary_df[stat_col], errors="coerce")
        mask = mask & (x["__stat"] >= float(min_stat))
    else:
        x["__stat"] = np.nan

    if q_col:
        x["__q"] = pd.to_numeric(summary_df[q_col], errors="coerce")
        mask = mask & x["__q"].notna() & (x["__q"] <= float(max_q))
    if sig_col:
        mask = mask & _truthy_series(summary_df[sig_col])

    x = x.loc[mask, ["ASV_ID", "index", "isa_label", "__stat"]].copy()
    if drop_all_combo and not x.empty:
        s_cols = [str(c).strip() for c in summary_df.columns if str(c).strip().startswith("s.")]
        all_label = None
        if s_cols:
            groups = [c.split("s.", 1)[1].strip() for c in s_cols if c.split("s.", 1)[1].strip()]
            if groups:
                all_label = normalize_combo("+".join(groups))
        if all_label:
            x = x[x["isa_label"].map(normalize_combo) != all_label].copy()
    if x.empty:
        return empty

    if group_long_df is not None and not group_long_df.empty:
        g = group_long_df[["ASV_ID", "index", "Group", "AxB"]].copy()
        g["ASV_ID"] = g["ASV_ID"].astype(str).str.strip()
        g["index"] = pd.to_numeric(g["index"], errors="coerce").astype("Int64")
        g["__long_label"] = g["Group"].map(normalize_combo)
        g["__axb"] = pd.to_numeric(g["AxB"], errors="coerce")
        x = x.merge(g[["ASV_ID", "index", "__long_label", "__axb"]], on=["ASV_ID", "index"], how="left")
        x["isa_label"] = np.where(
            x["__long_label"].notna() & (x["__long_label"].astype(str).str.strip() != ""),
            x["__long_label"],
            x["isa_label"]
        )
        x["isa_score"] = x["__axb"].fillna(x["__stat"]).fillna(0.0)
    else:
        x["isa_score"] = x["__stat"].fillna(0.0)

    x["isa_score"] = pd.to_numeric(x["isa_score"], errors="coerce").fillna(0.0)
    x = (
        x[["ASV_ID", "isa_label", "isa_score"]]
        .sort_values(["ASV_ID", "isa_score"], ascending=[True, False])
        .drop_duplicates(subset=["ASV_ID"], keep="first")
    )
    return x.reset_index(drop=True)


def annotate_modules_with_isa(
    modules_df: pd.DataFrame,
    isa_hits_df: pd.DataFrame,
    isa_palette: Dict[str, str],
) -> pd.DataFrame:
    """
    Assign each module one ISA label/color based on strongest aggregate ISA score
    among ASVs in that module.
    """
    if modules_df is None or modules_df.empty:
        return modules_df

    out = modules_df.copy()
    out["module_label"] = out["module_label"].astype(str).str.strip()
    out["module_has_isa"] = False
    out["module_isa_label"] = ""
    out["module_isa_color"] = NOT_FOCUS_COLOR
    out["module_isa_legend"] = ""
    out["module_plot_keep"] = False

    if isa_hits_df is None or isa_hits_df.empty:
        return out

    hits = isa_hits_df.copy()
    hits["ASV_ID"] = hits["ASV_ID"].astype(str).str.strip()
    hits["isa_label"] = hits["isa_label"].astype(str).str.strip()
    hits["isa_score"] = pd.to_numeric(hits["isa_score"], errors="coerce").fillna(0.0)

    merged = out[["Taxon", "module_label"]].copy()
    merged["Taxon"] = merged["Taxon"].astype(str).str.strip()
    merged = merged.merge(hits.rename(columns={"ASV_ID": "Taxon"}), on="Taxon", how="left")
    merged = merged.dropna(subset=["isa_label"])
    merged = merged[merged["isa_label"].astype(str).str.strip() != ""]
    if merged.empty:
        return out

    bad_labels = {"", "unassigned", "mna", "na", "none"}
    merged = merged[~merged["module_label"].str.lower().isin(bad_labels)].copy()
    if merged.empty:
        return out

    module_scores = (
        merged.groupby(["module_label", "isa_label"], dropna=False)["isa_score"]
        .sum()
        .reset_index()
        .sort_values(["module_label", "isa_score", "isa_label"], ascending=[True, False, True])
    )
    best = module_scores.groupby("module_label", as_index=False).first()
    label_map = dict(zip(best["module_label"], best["isa_label"]))

    out["module_isa_label"] = out["module_label"].map(label_map).fillna("")
    out["module_has_isa"] = out["module_isa_label"].astype(str).str.strip() != ""
    out["module_isa_color"] = out["module_isa_label"].map(
        lambda lbl: palette_get(isa_palette, lbl, NOT_FOCUS_COLOR) if str(lbl).strip() else NOT_FOCUS_COLOR
    )
    out["module_isa_legend"] = np.where(
        out["module_has_isa"],
        out["module_label"].astype(str) + " (" + out["module_isa_label"].astype(str) + ")",
        ""
    )
    out["module_plot_keep"] = out["module_has_isa"]
    return out


def annotate_modules_with_renewal_tests(
    modules_df: pd.DataFrame,
    association_path: Optional[str],
    profiles_path: Optional[str],
    max_q: float,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Attach significant module-level renewal contrasts to every member ASV.

    These labels summarize ecological-module tests; they are deliberately not
    described as ASV-level renewal indicators. The source data retain the
    internal ``baseline`` name, whereas figure labels use ``stagnation``.
    """
    if modules_df is None or modules_df.empty:
        return modules_df, pd.DataFrame()
    out = modules_df.copy()
    out["renewal_module_label"] = "No supported renewal-phase contrast"
    out["renewal_module_color"] = NOT_FOCUS_COLOR
    out["renewal_module_supported_tests"] = ""
    out["renewal_module_min_q"] = np.nan
    out["active_renewal_module_label"] = "No supported renewal/post-renewal contrast"
    out["active_renewal_module_color"] = NOT_FOCUS_COLOR
    out["active_renewal_module_q"] = np.nan
    out["active_renewal_stagnation_mean"] = np.nan
    out["active_renewal_mean"] = np.nan
    out["active_renewal_minus_stagnation"] = np.nan
    out["active_renewal_to_stagnation_ratio"] = np.nan
    if not association_path or not os.path.isfile(association_path):
        return out, pd.DataFrame()
    tests = load_table(association_path)
    required = {"analysis", "module_label", "adjusted_pvalue"}
    missing = required - set(tests.columns)
    if missing:
        raise ValueError(
            "Module-renewal association table lacks required columns: "
            + ", ".join(sorted(missing))
        )
    tests = tests.copy()
    tests["module_label"] = tests["module_label"].astype(str)
    tests["adjusted_pvalue"] = pd.to_numeric(tests["adjusted_pvalue"], errors="coerce")
    significant = tests.loc[tests["adjusted_pvalue"].le(max_q)].copy()
    label_map: Dict[str, str] = {}
    tests_map: Dict[str, str] = {}
    q_map: Dict[str, float] = {}
    for module_label, frame in significant.groupby("module_label"):
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
        label_map[module_label] = label
        tests_map[module_label] = ";".join(sorted(analyses))
        q_map[module_label] = float(frame["adjusted_pvalue"].min())
    out["renewal_module_label"] = out["module_label"].astype(str).map(label_map).fillna(
        "No supported renewal-phase contrast"
    )
    colors = {
        "Renewal onset vs stagnation": "#0072B2",
        "Post-renewal vs stagnation": "#E69F00",
        "Renewal/post-renewal vs stagnation": "#56B4E9",
        "Renewal phase vs stagnation": "#CC79A7",
        "No supported renewal-phase contrast": NOT_FOCUS_COLOR,
    }
    out["renewal_module_color"] = out["renewal_module_label"].map(colors).fillna(NOT_FOCUS_COLOR)
    out["renewal_module_supported_tests"] = out["module_label"].astype(str).map(tests_map).fillna("")
    out["renewal_module_min_q"] = out["module_label"].astype(str).map(q_map)
    active = significant.loc[significant["analysis"].astype(str).eq("active_vs_baseline")].copy()
    active_q = active.groupby("module_label")["adjusted_pvalue"].min().to_dict()
    out["active_renewal_module_q"] = out["module_label"].astype(str).map(active_q)
    if profiles_path and os.path.isfile(profiles_path) and active_q:
        profiles = load_table(profiles_path)
        profile_required = {
            "module_label", "renewal_phase", "samples_n", "mean_module_relative_abundance"
        }
        profile_missing = profile_required - set(profiles.columns)
        if profile_missing:
            raise ValueError(
                "Module-renewal profile table lacks required columns: "
                + ", ".join(sorted(profile_missing))
            )
        profiles = profiles.copy()
        profiles["module_label"] = profiles["module_label"].astype(str)
        profiles["samples_n"] = pd.to_numeric(profiles["samples_n"], errors="coerce")
        profiles["mean_module_relative_abundance"] = pd.to_numeric(
            profiles["mean_module_relative_abundance"], errors="coerce"
        )
        direction_map: Dict[str, str] = {}
        stagnation_map: Dict[str, float] = {}
        active_mean_map: Dict[str, float] = {}
        delta_map: Dict[str, float] = {}
        ratio_map: Dict[str, float] = {}
        for module_label in active_q:
            frame = profiles.loc[profiles["module_label"].eq(module_label)]
            baseline = frame.loc[frame["renewal_phase"].astype(str).isin(["baseline", "stagnation"])]
            active_phases = frame.loc[
                frame["renewal_phase"].astype(str).isin(["renewal", "post-renewal"])
            ]
            if baseline.empty or active_phases.empty:
                continue
            stagnation_mean = float(np.average(
                baseline["mean_module_relative_abundance"], weights=baseline["samples_n"]
            ))
            active_mean = float(np.average(
                active_phases["mean_module_relative_abundance"], weights=active_phases["samples_n"]
            ))
            delta = active_mean - stagnation_mean
            direction_map[module_label] = (
                "Higher during renewal/post-renewal" if delta > 0
                else "Higher during stagnation" if delta < 0
                else "Equal pooled abundance"
            )
            stagnation_map[module_label] = stagnation_mean
            active_mean_map[module_label] = active_mean
            delta_map[module_label] = delta
            ratio_map[module_label] = active_mean / stagnation_mean if stagnation_mean > 0 else np.nan
        out["active_renewal_module_label"] = out["module_label"].astype(str).map(direction_map).fillna(
            "No supported renewal/post-renewal contrast"
        )
        active_colors = {
            "Higher during renewal/post-renewal": "#0072B2",
            "Higher during stagnation": "#E69F00",
            "Equal pooled abundance": "#7A5195",
        }
        out["active_renewal_module_color"] = out["active_renewal_module_label"].map(
            active_colors
        ).fillna(NOT_FOCUS_COLOR)
        out["active_renewal_stagnation_mean"] = out["module_label"].astype(str).map(stagnation_map)
        out["active_renewal_mean"] = out["module_label"].astype(str).map(active_mean_map)
        out["active_renewal_minus_stagnation"] = out["module_label"].astype(str).map(delta_map)
        out["active_renewal_to_stagnation_ratio"] = out["module_label"].astype(str).map(ratio_map)
    return out, significant


def plot_modules(
    G: nx.Graph,
    pos: Dict,
    out_svg: str,
    *,
    color_attr: str,
    label_attr: str,
    filter_attr: Optional[str] = None,
    degree_scale: float = 80.0,
    edge_width_scale: float = 1.0,
    label: bool = False,
    title: Optional[str] = None,
    legend_title: str = "Network Modules",
):
    if filter_attr:
        nodes_to_plot = [n for n in G.nodes() if bool(G.nodes[n].get(filter_attr, False))]
        H = G.subgraph(nodes_to_plot).copy()
        pos_h = {n: pos[n] for n in H.nodes() if n in pos}
    else:
        H = G
        pos_h = pos

    if H.number_of_nodes() == 0:
        print(f"[WARN] No nodes to plot for {out_svg} after applying filter.")
        return

    fig, ax = figure_ax((17, 14))
    e_w = edge_widths_from_weights(H, scale=edge_width_scale, min_width=0.25)
    draw_edges_light(H, pos_h, alpha=0.6, edge_widths=e_w)

    def color_fn(n):
        return H.nodes[n].get(color_attr, NOT_FOCUS_COLOR)

    def size_fn(n):
        return degree_marker_area(H.nodes[n].get("Degree", 0.0), degree_scale)

    def alpha_fn(n):
        stab = _safe_float(H.nodes[n].get("node_stability", np.nan), np.nan)
        if np.isfinite(stab):
            return max(0.35, min(1.0, stab))
        return 0.8 if not is_not_focus_color(color_fn(n)) else 0.35

    draw_nodes_one_by_one(H, pos_h, color_fn, size_fn, alpha_fn=alpha_fn)

    labels = []
    for n in H.nodes():
        lbl = H.nodes[n].get(label_attr, "")
        col = H.nodes[n].get(color_attr, NOT_FOCUS_COLOR)
        if lbl and not is_not_focus_color(col):
            labels.append((str(lbl), col))
    legend_items = {}
    for lbl, col in labels:
        legend_items.setdefault(lbl, col)
    sorted_items = sorted(legend_items.items(), key=lambda x: natural_sort_key(x[0]))
    patches = [mpatches.Patch(color=c, label=l) for l, c in sorted_items]
    legend = None
    if patches:
        legend_ncol = 1 if len(patches) <= 18 else 2 if len(patches) <= 36 else 3
        legend = ax.legend(handles=patches, loc='upper left', bbox_to_anchor=(1.01, 1),
                           title=legend_title, frameon=False, labelspacing=0.8,
                           borderaxespad=0.0, ncol=legend_ncol)

    if label:
        to_label = [n for n in H.nodes() if not is_not_focus_color(H.nodes[n].get(color_attr, NOT_FOCUS_COLOR))]
        label_selected(H, pos_h, to_label, text_attr='Taxon')

    finalize_network_axes(
        fig, ax, out_svg,
        title=title or "SPIEC-EASI Network\nNode color: Module assignment | Node size: Degree",
        legends=[legend] if legend else [],
        right=0.72,
    )


def _taxonomy_display(value: object, rank: str) -> str:
    """Return a stable display value for incomplete SILVA taxonomy."""
    if pd.isna(value):
        return f"Unclassified {rank.lower()}"
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none", "unassigned", "uncultured"}:
        return f"Unclassified {rank.lower()}"
    return text.replace("_", " ")


def _categorical_palette(labels: Iterable[object]) -> Dict[str, str]:
    """Build a deterministic, collision-free palette within one figure set."""
    clean = sorted({_taxonomy_display(value, "taxon") for value in labels}, key=natural_sort_key)
    classified = [label for label in clean if not label.lower().startswith("unclassified")]
    colors = sns.color_palette("husl", max(len(classified), 1))
    palette = {label: mcolors.to_hex(color) for label, color in zip(classified, colors)}
    for label in clean:
        if label.lower().startswith("unclassified"):
            palette[label] = "#BDBDBD"
    return palette


def render_module_taxonomy_subnetworks(
    G: nx.Graph,
    pos: Dict,
    modules: pd.DataFrame,
    outdir: str,
    *,
    anchor_top_n: int,
    prominence_metrics: List[str],
    prominence_threshold: float,
    prominence_label_top_n: int,
    prominence_min_area: float,
    prominence_max_area: float,
    edge_width_scale: float,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Render an induced, phylum-colored subnetwork for every ecological module.

    Complete-network SPIEC-EASI coordinates provide the structural reference
    for every induced module subnetwork. Each module is translated and scaled
    uniformly to fill its own plotting panel; no node-specific displacement or
    repulsion is applied. Numbered annotations are placed in ordered exterior
    lanes and connected back to their source nodes. A rank-based
    consortium-prominence index controls node area, and an eligible ASV--MAG
    link controls border width.
    """
    if modules is None or modules.empty:
        print("[WARN] No module assignments available for module subnetworks.")
        return pd.DataFrame(), pd.DataFrame()
    if anchor_top_n < 1:
        raise ValueError("anchor_top_n must be at least one")
    metric_attributes = {
        "max_relative_abundance": "max_relative_abundance",
        "eigenvector": "EigenCentral",
        "participation": "Participation",
        "degree": "Degree",
        "betweenness": "Betweenness",
    }
    invalid_metrics = [metric for metric in prominence_metrics if metric not in metric_attributes]
    if invalid_metrics:
        raise ValueError(
            "Unsupported module-subnetwork prominence metric(s): "
            + ", ".join(invalid_metrics)
        )
    if not prominence_metrics:
        raise ValueError("At least one module-subnetwork prominence metric is required")
    if prominence_label_top_n < 0:
        raise ValueError("prominence_label_top_n cannot be negative")
    if not 0.0 <= prominence_threshold <= 1.0:
        raise ValueError("prominence_threshold must be between zero and one")
    if prominence_min_area <= 0 or prominence_max_area < prominence_min_area:
        raise ValueError("Invalid module-subnetwork prominence marker-area range")

    root = Path(outdir) / "module_subnetworks"
    shutil.rmtree(root, ignore_errors=True)
    root.mkdir(parents=True, exist_ok=True)
    module_lookup = modules[["Taxon", "module_label"]].dropna().copy()
    module_lookup["Taxon"] = module_lookup["Taxon"].astype(str)
    module_lookup["module_label"] = module_lookup["module_label"].astype(str)
    module_by_taxon = module_lookup.drop_duplicates("Taxon").set_index("Taxon")["module_label"].to_dict()
    taxon_to_node = {
        str(G.nodes[node].get("Taxon", node)): node for node in G.nodes()
    }

    neighbor_modules: Dict[object, List[str]] = {node: [] for node in G.nodes()}
    for source, target in G.edges():
        source_taxon = str(G.nodes[source].get("Taxon", source))
        target_taxon = str(G.nodes[target].get("Taxon", target))
        if target_taxon in module_by_taxon:
            neighbor_modules[source].append(module_by_taxon[target_taxon])
        if source_taxon in module_by_taxon:
            neighbor_modules[target].append(module_by_taxon[source_taxon])

    def participation(node: object) -> float:
        labels = neighbor_modules.get(node, [])
        if not labels:
            return 0.0
        fractions = pd.Series(labels).value_counts().to_numpy(float) / len(labels)
        return float(1.0 - np.square(fractions).sum())

    def scaled_marker_area(score: float) -> float:
        bounded = min(max(_safe_float(score, 0.0), 0.0), 1.0)
        return float(prominence_min_area + bounded * (prominence_max_area - prominence_min_area))

    def pfg_taxonomy(attrs: Dict[str, object]) -> str:
        """Compact taxonomic label used in module-subnetwork figure keys."""
        return "; ".join([
            f"P: {_taxonomy_display(attrs.get('Phylum'), 'phylum')}",
            f"F: {_taxonomy_display(attrs.get('Family'), 'family')}",
            f"G: {_taxonomy_display(attrs.get('Genus'), 'genus')}",
        ])

    audit_rows: List[Dict[str, object]] = []
    manifest_rows: List[Dict[str, object]] = []
    for module_label in sorted(module_lookup["module_label"].unique(), key=natural_sort_key):
        taxa = set(module_lookup.loc[module_lookup["module_label"].eq(module_label), "Taxon"])
        nodes = [taxon_to_node[taxon] for taxon in taxa if taxon in taxon_to_node and taxon_to_node[taxon] in pos]
        if not nodes:
            continue
        H = G.subgraph(nodes).copy()
        pos_h = {node: pos[node] for node in H.nodes()}

        metrics = pd.DataFrame({
            "node": list(H.nodes()),
            "Degree": [_safe_float(G.nodes[node].get("Degree", 0.0)) for node in H.nodes()],
            "EigenCentral": [_safe_float(G.nodes[node].get("EigenCentral", 0.0)) for node in H.nodes()],
            "Betweenness": [get_betweenness_value(G.nodes[node]) for node in H.nodes()],
            "max_relative_abundance": [_safe_float(G.nodes[node].get("max_relative_abundance", 0.0)) for node in H.nodes()],
            "Participation": [participation(node) for node in H.nodes()],
        })
        anchor_flags = np.zeros(len(metrics), dtype=bool)
        for metric in ("Degree", "EigenCentral", "Betweenness"):
            metrics[f"anchor_rank_{metric.lower()}"] = metrics[metric].rank(
                method="min", ascending=False
            )
            anchor_flags |= metrics[f"anchor_rank_{metric.lower()}"].le(anchor_top_n).to_numpy()
        metrics["is_ecological_anchor"] = anchor_flags
        percentile_columns = []
        for metric in prominence_metrics:
            source_column = metric_attributes[metric]
            percentile_column = f"prominence_percentile_{metric}"
            values = pd.to_numeric(metrics[source_column], errors="coerce").fillna(0.0)
            if len(values) == 1:
                percentiles = pd.Series([1.0], index=values.index)
            else:
                ranks = values.rank(method="average", ascending=True)
                percentiles = (ranks - 1.0) / (len(values) - 1.0)
            metrics[percentile_column] = percentiles.clip(0.0, 1.0)
            percentile_columns.append(percentile_column)
        metrics["consortium_prominence_score"] = metrics[percentile_columns].mean(axis=1)
        metrics["consortium_prominence_rank"] = metrics["consortium_prominence_score"].rank(
            method="min", ascending=False
        ).astype(int)
        metrics["is_high_prominence_asv"] = metrics["consortium_prominence_score"].ge(
            prominence_threshold
        )
        score_by_node = metrics.set_index("node")["consortium_prominence_score"].to_dict()
        rank_by_node = metrics.set_index("node")["consortium_prominence_rank"].to_dict()
        high_prominence_by_node = metrics.set_index("node")["is_high_prominence_asv"].to_dict()
        marker_area_by_node = {
            node: scaled_marker_area(score_by_node[node]) for node in H.nodes()
        }
        ordered_nodes = sorted(
            H.nodes(), key=lambda node: natural_sort_key(str(G.nodes[node].get("Taxon", node)))
        )
        raw_coordinates = np.vstack([
            np.asarray(pos_h[node], dtype=float) for node in ordered_nodes
        ])
        raw_min = raw_coordinates.min(axis=0)
        raw_max = raw_coordinates.max(axis=0)
        raw_center = (raw_min + raw_max) / 2.0
        raw_span = float(np.max(raw_max - raw_min))
        if raw_span <= 1e-12:
            raw_span = 1.0
        # Uniform affine zoom retains orientation, relative distances, and
        # aspect ratio while making each module use the available panel.
        reference_coordinates = (raw_coordinates - raw_center) * (1.64 / raw_span)
        display_pos_h = {
            node: reference_coordinates[index].copy()
            for index, node in enumerate(ordered_nodes)
        }

        mag_nodes = [node for node in H.nodes() if bool(G.nodes[node].get("has_mag_pair", False))]
        high_prominence_nodes = [node for node in H.nodes() if bool(high_prominence_by_node[node])]
        top_nodes = [
            node for node in H.nodes()
            if rank_by_node[node] <= prominence_label_top_n
        ]
        label_nodes = sorted(
            set(mag_nodes + high_prominence_nodes + top_nodes),
            key=lambda node: natural_sort_key(str(G.nodes[node].get("Taxon", node))),
        )
        label_number = {node: index for index, node in enumerate(label_nodes, start=1)}

        for _, row in metrics.iterrows():
            node = row["node"]
            attrs = G.nodes[node]
            audit_rows.append({
                "ecological_module": module_label,
                "graph_node_id": str(node),
                "ASV_ID": str(attrs.get("Taxon", node)),
                "phylum": _taxonomy_display(attrs.get("Phylum"), "phylum"),
                "family": _taxonomy_display(attrs.get("Family"), "family"),
                "genus": _taxonomy_display(attrs.get("Genus"), "genus"),
                "degree": row["Degree"],
                "eigenvector_centrality": row["EigenCentral"],
                "betweenness": row["Betweenness"],
                "max_relative_abundance": row["max_relative_abundance"],
                "participation_coefficient": row["Participation"],
                "anchor_rank_degree": row["anchor_rank_degree"],
                "anchor_rank_eigenvector": row["anchor_rank_eigencentral"],
                "anchor_rank_betweenness": row["anchor_rank_betweenness"],
                "is_ecological_anchor": bool(row["is_ecological_anchor"]),
                **{column: row[column] for column in percentile_columns},
                "consortium_prominence_score": row["consortium_prominence_score"],
                "consortium_prominence_rank": row["consortium_prominence_rank"],
                "consortium_prominence_threshold": prominence_threshold,
                "is_high_prominence_asv": bool(row["is_high_prominence_asv"]),
                "consortium_prominence_metrics": ",".join(prominence_metrics),
                "has_eligible_mag_link": bool(attrs.get("has_mag_pair", False)),
                "plot_label_number": label_number.get(node, ""),
                "plot_label_text": (
                    f"{label_number[node]}: {str(attrs.get('Taxon', node))} - {pfg_taxonomy(attrs)}"
                    if node in label_number else ""
                ),
                "spieceasi_x": float(pos_h[node][0]),
                "spieceasi_y": float(pos_h[node][1]),
                "display_x": float(display_pos_h[node][0]),
                "display_y": float(display_pos_h[node][1]),
                "display_displacement": 0.0,
                "layout_scope": "uniform module zoom of SPIEC-EASI coordinates; nodes were not individually repositioned",
            })

        for rank, attr in (("phylum", "Phylum"),):
            rank_dir = root / rank
            rank_dir.mkdir(parents=True, exist_ok=True)
            labels = {
                node: _taxonomy_display(G.nodes[node].get(attr), rank)
                for node in H.nodes()
            }
            palette = _categorical_palette(labels.values())
            role_specs = [
                (
                    "Top-ranked prominence ASVs",
                    sorted(top_nodes, key=lambda node: natural_sort_key(str(G.nodes[node].get("Taxon", node)))),
                ),
                (
                    f"High-prominence ASVs (score >= {prominence_threshold:.2f})",
                    sorted(high_prominence_nodes, key=lambda node: natural_sort_key(str(G.nodes[node].get("Taxon", node)))),
                ),
                (
                    "Genome-linked ASVs (bold border)",
                    sorted(mag_nodes, key=lambda node: natural_sort_key(str(G.nodes[node].get("Taxon", node)))),
                ),
            ]
            # Long P/F/G labels are shown one ASV per row. Give every role its
            # own full-width legend band and place the genome-linked group
            # directly below the two prominence-based ASV groups.
            role_column_counts = [1 for _ in role_specs]
            role_row_counts = [
                max(1, int(math.ceil(max(1, len(role_nodes)) / column_count)))
                for (_, role_nodes), column_count in zip(role_specs, role_column_counts)
            ]
            role_band_heights = [1.20 + 0.46 * row_count for row_count in role_row_counts]
            network_band_height = 10.8
            fig_height = network_band_height + sum(role_band_heights) + 0.8
            fig = plt.figure(figsize=(20, fig_height))
            grid = fig.add_gridspec(
                4, 2,
                width_ratios=(4.6, 1.7),
                height_ratios=(network_band_height, *role_band_heights),
                wspace=0.04,
                hspace=0.25,
            )
            ax = fig.add_subplot(grid[0, 0])
            encoding_ax = fig.add_subplot(grid[0, 1])
            asv_key_axes = [fig.add_subplot(grid[index, :]) for index in range(1, 4)]
            encoding_ax.axis("off")
            for key_ax in asv_key_axes:
                key_ax.axis("off")
            edge_widths = edge_widths_from_weights(
                H, scale=edge_width_scale, min_width=0.30, max_width=1.65
            )
            network_edge_artist = nx.draw_networkx_edges(
                H, display_pos_h, ax=ax, width=edge_widths,
                edge_color="#4A4A4A", alpha=0.88,
            )
            set_edge_artist_zorder(network_edge_artist, 1)
            # Draw high-prominence and genome-linked ASVs last so both encodings remain visible.
            for is_high_prominence in (False, True):
                for has_mag in (False, True):
                    draw_nodes = [
                        node for node in H.nodes()
                        if bool(high_prominence_by_node[node]) == is_high_prominence
                        and bool(G.nodes[node].get("has_mag_pair", False)) == has_mag
                    ]
                    if not draw_nodes:
                        continue
                    draw_nodes = sorted(draw_nodes, key=lambda node: score_by_node[node])
                    node_artist = nx.draw_networkx_nodes(
                        H, display_pos_h, nodelist=draw_nodes, ax=ax,
                        node_color=[palette[labels[node]] for node in draw_nodes],
                        node_size=[marker_area_by_node[node] for node in draw_nodes],
                        node_shape="D" if is_high_prominence else "o",
                        edgecolors=["black" if has_mag else "#666666" for _ in draw_nodes],
                        linewidths=[2.8 if has_mag else 0.55 for _ in draw_nodes],
                        alpha=1.0,
                    )
                    node_artist.set_zorder(3)
            # Keep the network itself untouched. Numbered annotations occupy
            # deterministic exterior lanes, sorted in the same vertical order
            # as their source nodes so leader lines do not cross within a lane.
            external_label_positions: Dict[object, np.ndarray] = {}
            module_center = np.mean(
                np.vstack([display_pos_h[node] for node in H.nodes()]), axis=0
            )
            node_vectors = {
                node: np.asarray(display_pos_h[node], dtype=float) - module_center
                for node in H.nodes()
            }
            network_radius = max(
                (float(np.linalg.norm(vector)) for vector in node_vectors.values()),
                default=1.0,
            )
            if network_radius <= 1e-12:
                network_radius = 1.0
            left_labels = sorted(
                [node for node in label_nodes if display_pos_h[node][0] < module_center[0]],
                key=lambda node: display_pos_h[node][1],
            )
            right_labels = sorted(
                [node for node in label_nodes if display_pos_h[node][0] >= module_center[0]],
                key=lambda node: display_pos_h[node][1],
            )
            lane_x = 1.18 * network_radius
            lane_y_extent = 0.92 * network_radius
            for side, lane_nodes in ((-1.0, left_labels), (1.0, right_labels)):
                lane_y = (
                    [0.0]
                    if len(lane_nodes) == 1
                    else np.linspace(-lane_y_extent, lane_y_extent, len(lane_nodes))
                )
                for node, y_offset in zip(lane_nodes, lane_y):
                    external_label_positions[node] = module_center + np.asarray([
                        side * lane_x, float(y_offset)
                    ])

            plot_radius = network_radius * 1.34
            ax.set_xlim(module_center[0] - plot_radius, module_center[0] + plot_radius)
            ax.set_ylim(module_center[1] - plot_radius, module_center[1] + plot_radius)

            # Draw each leader all the way to its source coordinate beneath
            # the node layer. The node masks the final segment, producing an
            # unambiguous visual connection at the marker boundary without a
            # gap or a line across the marker face.
            for node in label_nodes:
                source = np.asarray(display_pos_h[node], dtype=float)
                target = np.asarray(external_label_positions[node], dtype=float)
                ax.plot(
                    [source[0], target[0]], [source[1], target[1]],
                    color="#2F2F2F", linewidth=1.0, alpha=0.95,
                    solid_capstyle="round", zorder=2,
                )

            for node in label_nodes:
                label_position = external_label_positions[node]
                ax.annotate(
                    str(label_number[node]),
                    xy=display_pos_h[node], xycoords="data",
                    xytext=label_position, textcoords="data",
                    ha="center", va="center", fontsize=10.5, fontweight="bold",
                    color="black", zorder=31, annotation_clip=False,
                    bbox=dict(
                        boxstyle="circle,pad=0.28", facecolor="white",
                        edgecolor="#2F2F2F", linewidth=0.9,
                    ),
                )

            tax_handles = [
                mpatches.Patch(facecolor=palette[label], edgecolor="none", label=label)
                for label in sorted(palette, key=natural_sort_key)
            ]
            enc_handles = [
                plt.Line2D([], [], marker="D", linestyle="None", markersize=9,
                           markerfacecolor="#BDBDBD", markeredgecolor="#666666",
                           markeredgewidth=0.55, label=f"High-prominence ASV (score >= {prominence_threshold:.2f})"),
                plt.Line2D([], [], marker="o", linestyle="None", markersize=9,
                           markerfacecolor="white", markeredgecolor="black",
                           markeredgewidth=2.8, label="Eligible ASV-genome link"),
            ]
            size_handles = [
                ax.scatter([], [], s=scaled_marker_area(value), facecolor="#BDBDBD", edgecolor="#666666",
                           linewidth=0.55, label=f"{value:.1f}")
                for value in (0.0, 0.5, 1.0)
            ]
            sidebar_handles = tax_handles + enc_handles + size_handles
            sidebar_legend = encoding_ax.legend(
                handles=sidebar_handles,
                title=f"{rank.title()} and node encoding",
                frameon=False,
                handletextpad=0.55,
                labelspacing=0.55,
                fontsize=8.0,
                loc="upper left",
                bbox_to_anchor=(0.0, 1.0),
            )
            role_legends = []
            for key_ax, (role_title, role_nodes), column_count in zip(
                asv_key_axes, role_specs, role_column_counts
            ):
                role_handles = [
                    plt.Line2D(
                        [], [], linestyle="None", marker=None,
                        label=(
                            f"{label_number[node]}: {str(G.nodes[node].get('Taxon', node))} - "
                            f"{pfg_taxonomy(G.nodes[node])}"
                        ),
                    )
                    for node in role_nodes
                ]
                if not role_handles:
                    role_handles = [
                        plt.Line2D([], [], linestyle="None", marker=None, label="None in this module")
                    ]
                role_legends.append(key_ax.legend(
                    handles=role_handles,
                    title=role_title,
                    frameon=False,
                    handlelength=0.0,
                    handletextpad=0.0,
                    columnspacing=1.6,
                    labelspacing=0.48,
                    fontsize=7.3,
                    loc="upper left",
                    bbox_to_anchor=(0.0, 0.0, 1.0, 1.0),
                    mode="expand" if column_count > 1 else None,
                    borderaxespad=0.0,
                    ncol=column_count,
                ))
            ax.set_aspect("equal", adjustable="box")
            ax.axis("off")
            ax.set_title(
                f"Ecological module {module_label}\n"
                f"Node color: {rank} | Node size: consortium prominence | Bold border: genome linked"
            )
            fig.subplots_adjust(left=0.035, right=0.985, top=0.93, bottom=0.035)
            stem = rank_dir / f"{module_label}_{rank}_prominence_mag_subnetwork"
            for extension in ("pdf", "png", "svg"):
                kwargs = {"dpi": 600} if extension == "png" else {}
                fig.savefig(
                    stem.with_suffix(f".{extension}"), bbox_inches="tight",
                    pad_inches=0.5,
                    **kwargs,
                )
            plt.close(fig)
            manifest_rows.append({
                "ecological_module": module_label,
                "taxonomy_rank": rank,
                "nodes_n": H.number_of_nodes(),
                "internal_edges_n": H.number_of_edges(),
                "labeled_top_prominence_n": len(top_nodes),
                "high_prominence_asvs_n": len(high_prominence_nodes),
                "labeled_total_n": len(label_nodes),
                "mag_linked_asvs_n": sum(bool(G.nodes[node].get("has_mag_pair", False)) for node in H.nodes()),
                "output_stem": str(stem.relative_to(Path(outdir))),
                "layout_scope": "uniform module zoom without node displacement; numbered labels use ordered exterior lanes",
            })

    audit = pd.DataFrame(audit_rows)
    manifest = pd.DataFrame(manifest_rows)
    if not audit.empty:
        audit.to_csv(root / "module_subnetwork_node_audit.tsv", sep="\t", index=False)
        high_prominence = audit.loc[audit["is_high_prominence_asv"]].sort_values(
            ["ecological_module", "consortium_prominence_rank", "ASV_ID"],
            kind="mergesort",
        )
        high_prominence.to_csv(root / "module_high_prominence_asvs.tsv", sep="\t", index=False)
    if not manifest.empty:
        manifest.to_csv(root / "module_subnetwork_plot_manifest.tsv", sep="\t", index=False)
    return audit, manifest


# ---------------------------- Main pipeline ----------------------------------
def main():
    p = argparse.ArgumentParser(
        description="Render SPIEC-EASI GraphML networks with ISA/metadata overlays."
    )
    # Base dir + defaults
    p.add_argument("--data-dir", default="/home/ryan/SeqData/SeqData/UBC/LMP_priority1",
                   help="Base data directory (used to resolve default paths).")
    p.add_argument("--outdir", required=True, help="Output directory for figures.")

    # Inputs (override if paths differ)
    p.add_argument("--graph-pos-all", default=None,
                   help="GraphML with all positive edges (default: <data-dir>/spark_combined_output/spieceasi/network_pos_all.graphml)")
    p.add_argument("--graph-pos-sub", default=None,
                   help="GraphML with thresholded positive edges (default: <data-dir>/spark_combined_output/spieceasi/network_pos_thr.graphml OR _pos_sub.graphml)")
    p.add_argument("--node-features", default=None,
                   help="Node features CSV (default: <data-dir>/spark_combined_output/spieceasi/node_features.csv)")
    p.add_argument("--asv-counts", default=None,
                   help="ASV count table (tsv; default: <data-dir>/spark_combined_output/ASVs/ASV_final.micro.tsv)")
    p.add_argument("--taxonomy", default=None,
                   help="taxonomy_updated.tsv (default: <data-dir>/spark_combined_output/metadata/taxonomy_updated.tsv)")
    p.add_argument("--group1-summary", "--type-summary", dest="group1_summary", default=None,
                   help="First ISA summary TSV (defaults to type_group_indicator_species_summary.tsv).")
    p.add_argument("--group2-summary", "--status-summary", dest="group2_summary", default=None,
                   help="Second ISA summary TSV (defaults to status_indicator_species_summary.tsv).")
    p.add_argument("--group1-name", default=None,
                   help="Display/slug name for first ISA summary; defaults inferred from filename.")
    p.add_argument("--group2-name", default=None,
                   help="Display/slug name for second ISA summary; defaults inferred from filename.")
    p.add_argument("--isa-group-cols", default="",
                   help="Comma-separated ISA grouping columns in configured order; used to resolve groupN overlay modes.")
    p.add_argument("--isa-palette-map-json", default="{}",
                   help="JSON object mapping ISA group name -> palette string.")
    p.add_argument("--isa-order-map-json", default="{}",
                   help="JSON object mapping ISA group name -> ordered label list.")
    p.add_argument("--isa-focus-map-json", default="{}",
                   help="JSON object mapping ISA group name -> focused ISA component label.")
    p.add_argument("--isa-summary-mode", choices=["auto", "default", "duleg"], default="auto",
                   help="Preferred ISA summary variant when both default and DULEG summaries are staged [default: auto].")
    p.add_argument("--metadata", default=None,
                   help="Optional metadata TSV/CSV used to derive ISA palettes from label/color columns.")
    p.add_argument("--sample-col", default="sampleID",
                   help="Sample ID column in metadata (for reference/validation).")
    p.add_argument("--color-col", default="Color",
                   help="Metadata color column for ISA palette derivation.")
    p.add_argument("--modules-sub", default=None,
                   help="Optional module assignment TSV for thresholded graph.")
    p.add_argument("--modules-all", default=None,
                   help="Optional module assignment TSV for all-edges graph.")
    p.add_argument("--asv-mag-pairing", default=None,
                   help="Optional ASV-to-MAG pairing TSV used for MAG-paired overlay plots.")
    p.add_argument("--module-best-min-size", type=int, default=5,
                   help="Minimum node count for a module to be considered best [default: 5].")
    p.add_argument("--module-best-min-stability", type=float, default=0.7,
                   help="Minimum mean node stability for a module to be considered best [default: 0.7].")
    p.add_argument("--module-best-top-n", type=int, default=8,
                   help="Maximum number of topologically strongest stable modules to highlight.")
    p.add_argument("--module-best-only", action="store_true",
                   help="Color/label only best modules in module plots; non-best modules are light gray.")
    p.add_argument("--module-isa-only", action="store_true",
                   help="Plot only modules that contain at least one significant ISA ASV.")
    p.add_argument("--module-color-by-isa", action="store_true",
                   help="Color each module by its dominant ISA group label instead of module ID.")
    p.add_argument("--module-isa-source", default="group1",
                   help="Which ISA group to use for ISA-associated module labeling (groupN alias or actual group name).")
    p.add_argument("--module-isa-min-stat", type=float, default=0.25,
                   help="Minimum ISA stat used to mark ASVs as significant for module ISA association.")
    p.add_argument("--module-isa-max-q", type=float, default=0.05,
                   help="Maximum ISA q-value used to mark ASVs as significant for module ISA association.")
    p.add_argument("--module-subnetworks", action="store_true",
                   help="Render one fixed-layout, phylum-colored subnetwork for every ecological module.")
    p.add_argument("--anchor-top-n", type=int, default=1,
                   help="Within-module top-N rank used for degree/eigenvector/betweenness anchor union [default: 1].")
    p.add_argument("--module-subnetwork-prominence-metrics", default="max_relative_abundance,eigenvector,participation",
                   help="Comma-separated properties combined as equal-weight within-module percentiles for node size.")
    p.add_argument("--module-subnetwork-prominence-threshold", type=float, default=0.75,
                   help="Composite-score threshold used to mark high-prominence ASVs with diamond nodes.")
    p.add_argument("--module-subnetwork-label-top-n", type=int, default=3,
                   help="Label this many highest-prominence ASVs per module in addition to every genome-linked ASV.")
    p.add_argument("--module-subnetwork-prominence-min-area", type=float, default=30.0)
    p.add_argument("--module-subnetwork-prominence-max-area", type=float, default=520.0)
    p.add_argument("--module-renewal-association", default=None,
                   help="Optional module-level renewal association table used for an indirect ASV-network overlay.")
    p.add_argument("--module-renewal-profiles", default=None,
                   help="Module abundance summaries used to assign direction to significant pooled renewal contrasts.")
    p.add_argument("--module-renewal-max-q", type=float, default=0.05,
                   help="Maximum adjusted p-value for module-level renewal overlay support [default: 0.05].")

    # Layout options
    p.add_argument("--layout-json-all", default=None, help="Cache/Load layout JSON for graph-pos-all.")
    p.add_argument("--layout-json-sub", default=None, help="Cache/Load layout JSON for graph-pos-sub.")
    p.add_argument("--layout-seed", type=int, default=42, help="Seed for spring_layout.")
    p.add_argument("--layout-scale", type=float, default=3.0, help="Scale factor applied to layout coordinates.")

    # Visual scales
    p.add_argument("--degree-scale", type=float, default=80.0, help="Base size multiplier for degree plots.")
    p.add_argument("--degree-size-mode", default="legacy", choices=["legacy", "linear"],
                   help="Degree size mapping mode.")
    p.add_argument("--degree-min-area", type=float, default=0.0,
                   help="Minimum marker area added in linear degree mode.")
    p.add_argument("--edge-width-scale", type=float, default=5.0, help="Edge width multiplier for |weight|.")
    p.add_argument("--isa-scale", type=float, default=700.0, help="Node size multiplier for ISA (AxB).")
    p.add_argument("--betweenness-threshold", type=float, default=0.05,
                   help="Threshold used to split dark/light betweenness classes [default: 0.05].")
    p.add_argument("--abundance-size-mode", default="legacy", choices=["legacy", "decade_doubling"],
                   help="Abundance size mapping mode.")
    p.add_argument("--abundance-reference", type=float, default=5000.0,
                   help="Reference abundance for decade_doubling mode.")
    p.add_argument("--abundance-reference-area", type=float, default=80.0,
                   help="Marker area assigned to abundance-reference in decade_doubling mode.")
    p.add_argument("--abundance-min-area", type=float, default=8.0,
                   help="Minimum marker area for abundance-scaled plots.")
    p.add_argument("--abundance-max-area", type=float, default=420.0,
                   help="Maximum marker area for abundance-scaled plots.")
    p.add_argument("--abundance-scale-power", type=float, default=1.6,
                   help="Power used to spread abundance sizes after log scaling.")
    p.add_argument("--max-labels", type=int, default=100,
                   help="Maximum collision-adjusted node labels per plot; 0 disables the limit [default: 100].")

    # Which plots to render
    p.add_argument("--modes", nargs="+", default=["all"],
                   help="Which figure(s) to render. Supports legacy fixed modes plus groupN ISA modes such as group3_isa_all.")

    args = p.parse_args()
    SIZE_CONFIG["degree_size_mode"] = str(args.degree_size_mode).strip().lower()
    SIZE_CONFIG["degree_min_area"] = _safe_float(args.degree_min_area, 0.0)
    SIZE_CONFIG["abundance_size_mode"] = str(args.abundance_size_mode).strip().lower()
    SIZE_CONFIG["abundance_reference"] = _safe_float(args.abundance_reference, 5000.0)
    SIZE_CONFIG["abundance_reference_area"] = _safe_float(args.abundance_reference_area, 80.0)
    LABEL_CONFIG["max_labels"] = max(0, int(args.max_labels))

    data_dir = args.data_dir
    os.makedirs(args.outdir, exist_ok=True)

    # Resolve defaults
    graph_all = args.graph_pos_all or os.path.join(data_dir, "spark_combined_output/spieceasi/network_pos_all.graphml")
    graph_sub = args.graph_pos_sub or os.path.join(data_dir, "spark_combined_output/spieceasi/network_pos_thr.graphml")
    node_features_path = args.node_features or os.path.join(data_dir, "spark_combined_output/spieceasi/node_features.csv")
    asv_counts_path = args.asv_counts or os.path.join(data_dir, "spark_combined_output/ASVs/ASV_final.micro.tsv")
    taxonomy_path = args.taxonomy or os.path.join(data_dir, "spark_combined_output/metadata/taxonomy_updated.tsv")
    group1_summary_path = args.group1_summary or os.path.join(data_dir, "spark_combined_output/indicspecies/type_group_indicator_species_summary.tsv")
    group2_summary_path = args.group2_summary or os.path.join(data_dir, "spark_combined_output/indicspecies/status_indicator_species_summary.tsv")
    group1_name = args.group1_name or infer_group_name(group1_summary_path, "group1")
    group2_name = args.group2_name or infer_group_name(group2_summary_path, "group2")
    group1_slug = re.sub(r"[^0-9A-Za-z._-]", "_", group1_name)
    group2_slug = re.sub(r"[^0-9A-Za-z._-]", "_", group2_name)
    isa_group_cols = [x for x in parse_csv_list(args.isa_group_cols) if x]
    if not isa_group_cols:
        isa_group_cols = [group1_name, group2_name]
    try:
        isa_palette_map = json.loads(args.isa_palette_map_json or "{}")
    except Exception:
        isa_palette_map = {}
    try:
        isa_order_map = json.loads(args.isa_order_map_json or "{}")
    except Exception:
        isa_order_map = {}
    try:
        isa_focus_map = json.loads(args.isa_focus_map_json or "{}")
    except Exception:
        isa_focus_map = {}
    summary_by_group = collect_isa_summary_paths(
        isa_group_cols,
        group1_summary_path if (args.group1_summary or not args.isa_group_cols) else "",
        group2_summary_path if (args.group2_summary or not args.isa_group_cols) else "",
        group1_name,
        group2_name,
        None if args.isa_summary_mode == "auto" else args.isa_summary_mode,
    )
    requested_groups = [name for name in isa_group_cols if name in summary_by_group]
    # Stratified ISA tables are discovered from the staged INDICSPECIES outputs.
    # Retain requested standalone overlays first, then append every stratified
    # comparison so downstream network figures do not silently discard them.
    isa_group_cols = requested_groups + [
        name for name in summary_by_group if name not in requested_groups
    ]
    modules_sub_path = args.modules_sub
    modules_all_path = args.modules_all
    asv_mag_pairing_path = args.asv_mag_pairing or os.path.join(data_dir, "asv_mag_link/tables/asv2mag_pairing.tsv")

    # Load inputs
    nf = load_table(node_features_path, sep=',', index_col=0)  # index = GraphML_ID
    ensure_cols(nf, ["Taxon", "Degree", "Betweenness", "Closeness", "EigenCentral"], "node_features")

    asv = load_table(asv_counts_path, sep='\t', index_col=0)
    asv_numeric = asv.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    sample_totals = asv_numeric.sum(axis=0).replace(0.0, np.nan)
    asv_relative = asv_numeric.divide(sample_totals, axis=1).fillna(0.0)
    asv_stack = asv_numeric.stack().reset_index()
    asv_stack.columns = ['ASV_ID', 'sample', 'count']
    abund_stats = asv_stack.groupby('ASV_ID')['count'].agg(['mean', 'median']).reset_index()
    abund_stats['mean'] = abund_stats['mean'].astype(float)
    abund_stats['median'] = abund_stats['median'].astype(float)
    abund_stats = abund_stats.set_index('ASV_ID')
    abund_stats['max_relative_abundance'] = asv_relative.max(axis=1)
    abund_stats['prevalence'] = asv_numeric.gt(0).mean(axis=1)
    abund_stats = abund_stats.reset_index()[['ASV_ID', 'mean', 'median', 'max_relative_abundance', 'prevalence']]

    tax = load_table(taxonomy_path, sep='\t')
    if "Feature ID" in tax.columns:
        tax['ASV_ID'] = tax['Feature ID'].astype(str).str.partition(';')[0]
    if 'ASV_ID' not in tax.columns or 'Taxon' not in tax.columns:
        # Handle case where ASV_ID is embedded; keep your original logic
        if 'ASV_ID' in tax.columns:
            pass
        else:
            die("taxonomy_updated.tsv must include 'ASV_ID' and 'Taxon' columns.")
    # normalize ASV_ID if it contains suffix like ';...' at end
    tax['ASV_ID'] = [x.rsplit(';', 1)[0] if isinstance(x, str) and ';' in x else x for x in tax['ASV_ID']]
    # expand taxonomy levels
    tax['CompleteTaxonomy'] = tax['Taxon'].fillna('Unassigned').astype(str)
    tdf = pd.DataFrame([split_taxa_string(x) for x in tax['Taxon']])
    tax = pd.concat([tax[['ASV_ID', 'CompleteTaxonomy']], tdf], axis=1).set_index('ASV_ID', drop=True)

    metadata_df = None
    if args.metadata and os.path.exists(args.metadata):
        metadata_df = load_table(args.metadata)
        if args.sample_col and args.sample_col not in metadata_df.columns:
            print(f"[WARN] metadata sample column not found: {args.sample_col}")
    isa_specs = []
    for idx, group_name in enumerate(isa_group_cols, start=1):
        summary_path = summary_by_group.get(group_name)
        if not summary_path or not os.path.exists(summary_path):
            print(f"[WARN] ISA summary not found for {group_name}; skipping.")
            continue
        group_sum = load_table(summary_path, sep='\t')
        if "ASV" not in group_sum.columns and "ASV_ID" not in group_sum.columns:
            print(f"[WARN] ISA summary missing ASV/ASV_ID column for {group_name}: {summary_path}")
            continue
        ensure_cols(group_sum, ['index'], f"{group_name}_summary")
        group_long = long_AB_for_group(group_sum.copy())
        group_all_combo = infer_all_combo_label(group_long["Group"].dropna().astype(str).tolist())
        group_sig_hits = derive_significant_isa_hits(
            group_sum, group_long,
            min_stat=args.module_isa_min_stat,
            max_q=args.module_isa_max_q,
            drop_all_combo=False
        )
        group_order = canonicalize_group_type_order([
            str(x).strip()
            for x in (isa_order_map.get(group_name, []) or [])
            if str(x).strip()
        ])
        focus_label = normalize_combo(str(isa_focus_map.get(group_name, "")).strip()) if isa_focus_map.get(group_name) else ""
        manual_palette = normalize_palette_keys(parse_mapping(str(isa_palette_map.get(group_name, ""))))
        group_palette: Dict[str, str] = {}
        if metadata_df is not None and group_name in metadata_df.columns and args.color_col in metadata_df.columns:
            group_palette.update(
                build_palette_from_metadata(metadata_df, group_name, args.color_col, group_order)
            )
        group_palette.update(manual_palette)
        if not group_palette:
            group_palette = auto_palette(group_long["Group"].dropna().unique().tolist())
        group_palette = augment_combo_palette(group_palette, group_long["Group"].dropna().astype(str).tolist())
        group_palette = ensure_palette_has_order_entries(group_palette, group_order)
        group_palette = normalize_palette_keys(group_palette)
        if idx == 2:
            group_palette = canonicalize_group2_palette_aliases(group_palette)
        focus_order: List[str] = []
        if focus_label:
            all_candidates = canonicalize_group_type_order(
                group_order +
                list(group_palette.keys()) +
                group_long["Group"].dropna().astype(str).tolist()
            )
            focus_order = [
                lbl for lbl in all_candidates
                if combo_contains_component(lbl, focus_label)
            ]
            if not focus_order:
                focus_order = [focus_label]

        slug = slugify_group_name(group_name, f"group{idx}")
        label_attr = f"isa_label__{slug}"
        score_attr = f"isa_score__{slug}"
        color_attr = f"isa_color__{slug}"

        nfeat_group = nf.reset_index().merge(
            group_long.set_index('ASV_ID'), left_on='Taxon', right_index=True, how='left'
        ).set_index('GraphML_ID')
        nfeat_group[label_attr] = nfeat_group["Group"].map(lambda x: normalize_combo(x) if pd.notna(x) else x)
        sig_label = group_sig_hits.set_index("ASV_ID")["isa_label"].to_dict()
        sig_score = group_sig_hits.set_index("ASV_ID")["isa_score"].to_dict()
        nfeat_group[label_attr] = nfeat_group["Taxon"].map(sig_label).fillna("not_indicator")
        nfeat_group[score_attr] = nfeat_group["Taxon"].map(sig_score).fillna(0.0)
        nfeat_group[color_attr] = nfeat_group[label_attr].map(lambda x: palette_get(group_palette, x, NOT_FOCUS_COLOR))
        nfeat_group = nfeat_group.reset_index().merge(
            tax.reset_index(), left_on="Taxon", right_on="ASV_ID", how="left"
        ).set_index("GraphML_ID")

        isa_specs.append({
            "index": idx,
            "name": group_name,
            "slug": slug,
            "summary_path": summary_path,
            "palette": group_palette,
            "order": group_order,
            "focus_label": focus_label,
            "focus_order": focus_order,
            "all_combo_label": group_all_combo or "",
            "sig_hits": group_sig_hits,
            "label_attr": label_attr,
            "score_attr": score_attr,
            "color_attr": color_attr,
            "nfeat": nfeat_group,
        })

    if len(isa_specs) < 2:
        die("Need at least two resolved ISA summaries for network overlay plotting.")

    render_manifest = {
        "isa_groups": [
            {
                "index": spec["index"],
                "name": spec["name"],
                "slug": spec["slug"],
                "summary_path": spec["summary_path"],
                "order": spec["order"],
                "focus_label": spec["focus_label"],
                "all_combo_label": spec["all_combo_label"],
                "palette": spec["palette"],
                "n_sig_hits": int(spec["sig_hits"].shape[0]) if spec["sig_hits"] is not None else 0,
            }
            for spec in isa_specs
        ],
        "module_isa_source": args.module_isa_source,
        "module_best_only": bool(args.module_best_only),
        "module_isa_only": bool(args.module_isa_only),
        "module_color_by_isa": bool(args.module_color_by_isa),
    }
    with open(os.path.join(args.outdir, "network_isa_render_manifest.json"), "w") as fh:
        json.dump(render_manifest, fh, indent=2)

    primary_isa_spec = isa_specs[0]
    secondary_isa_spec = isa_specs[1]
    group1_name = primary_isa_spec["name"]
    group2_name = secondary_isa_spec["name"]
    group1_slug = primary_isa_spec["slug"]
    group2_slug = secondary_isa_spec["slug"]
    group1_palette = primary_isa_spec["palette"]
    group2_palette = secondary_isa_spec["palette"]
    group1_order = primary_isa_spec["order"]
    group2_order = secondary_isa_spec["order"]
    focus_group1_label = primary_isa_spec["focus_label"]
    focus_group2_label = secondary_isa_spec["focus_label"]
    group1_focus_order = primary_isa_spec["focus_order"]
    group1_focus_all_combo = primary_isa_spec["all_combo_label"]
    group2_all_combo = secondary_isa_spec["all_combo_label"]
    group1_sig_hits = primary_isa_spec["sig_hits"]
    group2_sig_hits = secondary_isa_spec["sig_hits"]
    nfeat_group1 = primary_isa_spec["nfeat"]
    nfeat_group2 = secondary_isa_spec["nfeat"]
    nfeat_group1["group1_label"] = nfeat_group1[primary_isa_spec["label_attr"]]
    nfeat_group1["AxB_group1"] = pd.to_numeric(nfeat_group1[primary_isa_spec["score_attr"]], errors="coerce").fillna(0.0)
    nfeat_group1["group1_color"] = nfeat_group1[primary_isa_spec["color_attr"]]
    nfeat_group1["A_group1"] = 0.0
    nfeat_group1["B_group1"] = 0.0
    nfeat_group2["group2_label"] = nfeat_group2[secondary_isa_spec["label_attr"]]
    nfeat_group2["AxB_group2"] = pd.to_numeric(nfeat_group2[secondary_isa_spec["score_attr"]], errors="coerce").fillna(0.0)
    nfeat_group2["group2_color"] = nfeat_group2[secondary_isa_spec["color_attr"]]
    nfeat_group2["A_group2"] = 0.0
    nfeat_group2["B_group2"] = 0.0

    # abundance table with taxonomy
    nfeat_abund = nf.reset_index().merge(abund_stats, left_on='Taxon', right_on='ASV_ID', how='left').set_index('GraphML_ID')
    nfeat_abund = nfeat_abund.reset_index().merge(
        tax.reset_index(), left_on='Taxon', right_on='ASV_ID', how='left'
    ).set_index('GraphML_ID')

    # Module-quality statistics require topology, so load the graphs before
    # selecting and ranking modules.
    def load_graph(path: str) -> nx.Graph:
        if not os.path.exists(path):
            die(f"GraphML not found: {path}")
        return nx.read_graphml(path)

    G_all = load_graph(graph_all)
    G_sub = load_graph(graph_sub)
    for graph_name, graph_obj in (("all", G_all), ("sub", G_sub)):
        isolates = list(nx.isolates(graph_obj))
        if isolates:
            graph_obj.remove_nodes_from(isolates)
            print(f"[INFO] Removed {len(isolates)} unconnected nodes from {graph_name} network rendering.")

    # modules tables (optional)
    modules_sub = load_modules_table(modules_sub_path, "sub")
    modules_all = load_modules_table(modules_all_path, "all")
    asv_mag_pairing = load_asv_mag_pairing(asv_mag_pairing_path)
    module_palette = build_module_palette(
        pd.concat([modules_sub.get("module_label", pd.Series(dtype=str)),
                   modules_all.get("module_label", pd.Series(dtype=str))], ignore_index=True).dropna().tolist()
    )
    best_sub_labels, best_sub_stats = select_best_modules(
        modules_sub,
        min_size=args.module_best_min_size,
        min_stability=args.module_best_min_stability,
        graph=G_sub,
        top_n=args.module_best_top_n,
        ensure_one=True
    )
    best_all_labels, best_all_stats = select_best_modules(
        modules_all,
        min_size=args.module_best_min_size,
        min_stability=args.module_best_min_stability,
        graph=G_all,
        top_n=args.module_best_top_n,
        ensure_one=True
    )
    if not best_sub_stats.empty:
        best_sub_stats.to_csv(os.path.join(args.outdir, "network_modules_best_stats_sub.tsv"), sep="\t", index=False)
    if not best_all_stats.empty:
        best_all_stats.to_csv(os.path.join(args.outdir, "network_modules_best_stats_all.tsv"), sep="\t", index=False)

    def _apply_best_focus(df: pd.DataFrame, best_labels: set) -> pd.DataFrame:
        if df.empty:
            return df
        out = df.copy()
        if args.module_best_only:
            out["module_is_best"] = out["module_label"].astype(str).isin(best_labels)
            out["module_color_plot"] = np.where(out["module_is_best"], out["module_color"], NOT_FOCUS_COLOR)
            out["module_label_plot"] = np.where(out["module_is_best"], out["module_label"], "")
        else:
            out["module_is_best"] = True
            out["module_color_plot"] = out["module_color"]
            out["module_label_plot"] = out["module_label"]
        return out

    if not modules_sub.empty:
        modules_sub["module_color"] = modules_sub["module_label"].map(module_palette).fillna(NOT_FOCUS_COLOR)
        modules_sub = _apply_best_focus(modules_sub, best_sub_labels)
        modules_sub["module_plot_keep"] = True
    if not modules_all.empty:
        modules_all["module_color"] = modules_all["module_label"].map(module_palette).fillna(NOT_FOCUS_COLOR)
        modules_all = _apply_best_focus(modules_all, best_all_labels)
        modules_all["module_plot_keep"] = True

    isa_spec_by_name = {spec["name"]: spec for spec in isa_specs}
    isa_spec_by_alias = {f"group{spec['index']}": spec for spec in isa_specs}
    module_isa_spec = isa_spec_by_alias.get(args.module_isa_source) or isa_spec_by_name.get(args.module_isa_source) or primary_isa_spec
    isa_source_hits = module_isa_spec["sig_hits"]
    isa_source_palette = module_isa_spec["palette"]
    phylum_isa_size_attr = module_isa_spec["score_attr"]
    # Always annotate ISA-module membership, but only override plotting focus when ISA-module mode is requested.
    # This preserves module_best_only behavior unless module_isa_only/module_color_by_isa is enabled.
    if (not modules_sub.empty) or (not modules_all.empty):
        modules_sub = annotate_modules_with_isa(modules_sub, isa_source_hits, isa_source_palette)
        modules_all = annotate_modules_with_isa(modules_all, isa_source_hits, isa_source_palette)

        def _apply_isa_module_focus(df: pd.DataFrame) -> pd.DataFrame:
            if df.empty:
                return df
            out = df.copy()
            use_isa_focus = bool(args.module_isa_only or args.module_color_by_isa)
            if not use_isa_focus:
                return out

            isa_focus = pd.to_numeric(out["module_has_isa"], errors="coerce").fillna(False).astype(bool)
            if args.module_best_only and "module_is_best" in out.columns:
                best_focus = pd.to_numeric(out["module_is_best"], errors="coerce").fillna(False).astype(bool)
                focus = isa_focus & best_focus
            else:
                focus = isa_focus

            if args.module_color_by_isa:
                out["module_color_plot"] = np.where(focus, out["module_isa_color"], NOT_FOCUS_COLOR)
                out["module_label_plot"] = np.where(focus, out["module_isa_legend"], "")
            else:
                out["module_color_plot"] = np.where(focus, out["module_color"], NOT_FOCUS_COLOR)
                out["module_label_plot"] = np.where(focus, out["module_label"], "")

            # Keep all nodes unless explicit ISA-only filtering is requested.
            out["module_plot_keep"] = focus if args.module_isa_only else True
            return out

        modules_sub = _apply_isa_module_focus(modules_sub)
        modules_all = _apply_isa_module_focus(modules_all)

    modules_sub, renewal_tests_sub = annotate_modules_with_renewal_tests(
        modules_sub, args.module_renewal_association, args.module_renewal_profiles,
        args.module_renewal_max_q
    )
    modules_all, renewal_tests_all = annotate_modules_with_renewal_tests(
        modules_all, args.module_renewal_association, args.module_renewal_profiles,
        args.module_renewal_max_q
    )
    if not modules_all.empty:
        renewal_audit_columns = [
            "Taxon", "module_label", "renewal_module_label",
            "renewal_module_supported_tests", "renewal_module_min_q",
            "active_renewal_module_label", "active_renewal_module_q",
            "active_renewal_stagnation_mean", "active_renewal_mean",
            "active_renewal_minus_stagnation", "active_renewal_to_stagnation_ratio",
        ]
        modules_all[renewal_audit_columns].to_csv(
            os.path.join(args.outdir, "network_asv_module_renewal_overlay.tsv"),
            sep="\t", index=False,
        )
    if not renewal_tests_all.empty:
        renewal_tests_all.to_csv(
            os.path.join(args.outdir, "network_source_module_renewal_associations.tsv"),
            sep="\t", index=False,
        )

    if not modules_sub.empty and "module_plot_keep" in modules_sub.columns:
        keep_sub = int(pd.to_numeric(modules_sub["module_plot_keep"], errors="coerce").fillna(False).astype(bool).sum())
        kept_mod_sub = int(modules_sub.loc[modules_sub["module_plot_keep"].astype(bool), "module_label"].nunique()) if keep_sub > 0 else 0
        print(f"[INFO] module_sub nodes kept for plotting: {keep_sub}/{len(modules_sub)} (modules={kept_mod_sub})")
    if not modules_all.empty and "module_plot_keep" in modules_all.columns:
        keep_all = int(pd.to_numeric(modules_all["module_plot_keep"], errors="coerce").fillna(False).astype(bool).sum())
        kept_mod_all = int(modules_all.loc[modules_all["module_plot_keep"].astype(bool), "module_label"].nunique()) if keep_all > 0 else 0
        print(f"[INFO] module_all nodes kept for plotting: {keep_all}/{len(modules_all)} (modules={kept_mod_all})")
    nfeat_modules_sub = nf.reset_index().merge(
        modules_sub.set_index("Taxon"), left_on="Taxon", right_index=True, how="left"
    ).set_index("GraphML_ID")
    nfeat_modules_all = nf.reset_index().merge(
        modules_all.set_index("Taxon"), left_on="Taxon", right_index=True, how="left"
    ).set_index("GraphML_ID")
    nfeat_mag = nf.reset_index().merge(
        asv_mag_pairing.set_index("ASV_ID"), left_on="Taxon", right_index=True, how="left"
    ).set_index("GraphML_ID")
    if "has_mag_pair" in nfeat_mag.columns:
        nfeat_mag["has_mag_pair"] = nfeat_mag["has_mag_pair"].fillna(False).astype(bool)

    # Build deterministic Phylum palette. Avoid tab20 wraparound collisions for larger taxonomic sets.
    phyla = sorted({
        str(x).strip()
        for x in pd.concat([*(spec["nfeat"]["Phylum"] for spec in isa_specs), nfeat_abund['Phylum']]).dropna().tolist()
        if str(x).strip()
    })
    if len(phyla) <= 20:
        p_colors = sns.color_palette('tab20', len(phyla))
    else:
        p_colors = sns.color_palette('husl', len(phyla))
    phylum_palette = {p: mcolors.to_hex(c) for p, c in zip(phyla, p_colors)}

    keep_cols = [
        'Taxon', 'Degree', 'Betweenness', 'Betweenness_raw', 'Betweenness_norm', 'Closeness', 'EigenCentral',
        'A_group1', 'B_group1', 'AxB_group1', 'group1_label', 'group1_color',
        'A_group2', 'B_group2', 'AxB_group2', 'group2_label', 'group2_color',
        'module_id', 'module_label', 'module_color', 'module_is_best', 'module_color_plot', 'module_label_plot',
        'module_has_isa', 'module_isa_label', 'module_isa_color', 'module_isa_legend', 'module_plot_keep',
        'renewal_module_label', 'renewal_module_color', 'renewal_module_supported_tests', 'renewal_module_min_q',
        'active_renewal_module_label', 'active_renewal_module_color', 'active_renewal_module_q',
        'active_renewal_stagnation_mean', 'active_renewal_mean',
        'active_renewal_minus_stagnation', 'active_renewal_to_stagnation_ratio',
        'node_stability',
        'has_mag_pair', 'mag_pair_status', 'best_genome_id', 'mag_taxonomy_label',
        'CompleteTaxonomy', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species',
        'mean', 'median', 'max_relative_abundance', 'prevalence'
    ]
    for spec in isa_specs:
        keep_cols.extend([spec["label_attr"], spec["score_attr"], spec["color_attr"]])

    # Attach attributes (each plot function can use what it needs)
    for spec in isa_specs:
        add_node_attrs_from_df(G_all, spec["nfeat"], keep_cols)
    add_node_attrs_from_df(G_all, nfeat_abund, keep_cols)
    add_node_attrs_from_df(G_all, nfeat_modules_all, keep_cols)
    add_node_attrs_from_df(G_all, nfeat_mag, keep_cols)

    for spec in isa_specs:
        add_node_attrs_from_df(G_sub, spec["nfeat"], keep_cols)
    add_node_attrs_from_df(G_sub, nfeat_abund, keep_cols)
    add_node_attrs_from_df(G_sub, nfeat_modules_sub, keep_cols)
    add_node_attrs_from_df(G_sub, nfeat_mag, keep_cols)

    # Positions (cached)
    pos_all = spring_layout_cached(G_all, seed=args.layout_seed,
                                   scale_xy=args.layout_scale,
                                   layout_json=args.layout_json_all)
    layout_rows = []
    for node, coordinates in pos_all.items():
        attributes = G_all.nodes[node]
        asv_id = attributes.get("Taxon", attributes.get("name", node))
        asv_id = re.sub(r";.*$", "", str(asv_id).strip())
        layout_rows.append({
            "graph_node_id": str(node),
            "ASV_ID": asv_id,
            "spieceasi_x": float(coordinates[0]),
            "spieceasi_y": float(coordinates[1]),
            "layout_seed": int(args.layout_seed),
            "layout_scale": float(args.layout_scale),
        })
    layout_table = pd.DataFrame(layout_rows).sort_values("ASV_ID")
    if layout_table["ASV_ID"].duplicated().any():
        duplicated = sorted(layout_table.loc[
            layout_table["ASV_ID"].duplicated(keep=False), "ASV_ID"
        ].unique())
        raise ValueError(
            "SPIEC-EASI layout contains duplicate ASV identifiers: "
            + ", ".join(duplicated[:10])
        )
    layout_table.to_csv(
        os.path.join(args.outdir, "spieceasi_network_layout_all.tsv"),
        sep="\t", index=False,
    )
    if args.module_subnetworks:
        render_module_taxonomy_subnetworks(
            G_all,
            pos_all,
            modules_all,
            args.outdir,
            anchor_top_n=args.anchor_top_n,
            prominence_metrics=[value.strip() for value in args.module_subnetwork_prominence_metrics.split(",") if value.strip()],
            prominence_threshold=args.module_subnetwork_prominence_threshold,
            prominence_label_top_n=args.module_subnetwork_label_top_n,
            prominence_min_area=args.module_subnetwork_prominence_min_area,
            prominence_max_area=args.module_subnetwork_prominence_max_area,
            edge_width_scale=args.edge_width_scale,
        )
    same_graph_structure = (
        set(G_all.nodes()) == set(G_sub.nodes()) and
        {frozenset((u, v)) for u, v in G_all.edges()} == {frozenset((u, v)) for u, v in G_sub.edges()}
    )
    if same_graph_structure:
        pos_sub = dict(pos_all)
    else:
        pos_sub = spring_layout_cached(G_sub, seed=args.layout_seed,
                                       scale_xy=args.layout_scale,
                                       layout_json=args.layout_json_sub)

    # -------------------- Choose and render plots -----------------------------
    modes = set(args.modes)
    if "all" in modes:
        modes = {
            "degree_all", "degree_sub",
            "abundance_sub",
            "betweenness_all", "betweenness_all_labeled",
            "betweenness_mag_all", "betweenness_mag_all_labeled",
            "group1_isa", "group1_isa_labeled",
            "group1_isa_mag_all",
            "group1_isa_focus", "group1_isa_focus_labeled",
            "group2_isa", "group2_isa_labeled",
            "group2_isa_mag_all",
            "module_sub", "module_sub_labeled",
            "module_all",
            "mag_pair_all",
            "mag_pair_tax_all",
            "phylum_abund", "phylum_isa", "phylum_isa_labeled"
        }
        for spec in isa_specs[2:]:
            variants = isa_mode_variants(spec["index"])
            modes.update({
                variants["isa_all"],
                variants["isa_all_labeled"],
                variants["isa_mag_all"],
            })

    if asv_mag_pairing.empty:
        mag_modes = {mode for mode in modes if "mag" in mode.lower()}
        if mag_modes:
            print("[INFO] ASV-MAG linking is unavailable; omitting MAG-specific network modes.")
            modes.difference_update(mag_modes)

    # Backward-compatible aliases.
    if "type_isa" in modes:
        modes.add("group1_isa")
    if "type_isa_labeled" in modes:
        modes.add("group1_isa_labeled")
    if "status_isa" in modes:
        modes.add("group2_isa")
    if "status_isa_labeled" in modes:
        modes.add("group2_isa_labeled")
    if "type_venn" in modes or "type_venn_labeled" in modes:
        print("[WARN] Venn network modes are deprecated and ignored.")

    valid_static_modes = {
        "degree_all", "degree_sub",
        "abundance_sub", "abundance_all",
        "betweenness_sub", "betweenness_sub_labeled",
        "betweenness_all", "betweenness_all_labeled",
        "betweenness_mag_sub", "betweenness_mag_sub_labeled",
        "betweenness_mag_all", "betweenness_mag_all_labeled",
        "type_isa", "type_isa_labeled",
        "status_isa", "status_isa_labeled",
        "type_venn", "type_venn_labeled",
        "module_sub", "module_sub_labeled",
        "module_all", "module_all_labeled",
        "mag_pair_sub", "mag_pair_sub_labeled",
        "mag_pair_all", "mag_pair_all_labeled",
        "mag_pair_tax_sub", "mag_pair_tax_sub_labeled",
        "mag_pair_tax_all", "mag_pair_tax_all_labeled",
        "phylum_abund", "phylum_abund_mean", "phylum_abund_median",
        "phylum_isa", "phylum_isa_labeled",
        "phylum_abund_all", "phylum_abund_all_mean", "phylum_abund_all_median",
        "phylum_isa_all", "phylum_isa_all_labeled",
    }
    valid_dynamic_mode = re.compile(r"^group\d+_isa(?:_mag|_focus)?(?:_all)?(?:_labeled)?$")
    unknown_modes = sorted([m for m in modes if m not in valid_static_modes and not valid_dynamic_mode.match(m)])
    if unknown_modes:
        print(f"[WARN] Ignoring unsupported network modes: {', '.join(unknown_modes)}")
        modes = {m for m in modes if m not in unknown_modes}
    with open(os.path.join(args.outdir, "network_render_modes.txt"), "w") as fh:
        for mode in sorted(modes):
            fh.write(f"{mode}\n")

    # Degree (all edges; weighted widths)
    if "degree_all" in modes:
        out = os.path.join(args.outdir, "network_degree_POS_ALL.svg")
        plot_degree(G_all, pos_all, out, degree_scale=args.degree_scale,
                    edge_width_scale=args.edge_width_scale)

    # Degree (thresholded subgraph)
    if "degree_sub" in modes:
        out = os.path.join(args.outdir, "network_degree_POS_SUB.svg")
        plot_degree(G_sub, pos_sub, out, degree_scale=args.degree_scale,
                    edge_width_scale=args.edge_width_scale)

    # Abundance (thresholded subgraph)
    if "abundance_sub" in modes:
        out = os.path.join(args.outdir, "network_abundance.svg")
        plot_abundance(
            G_sub, pos_sub, out,
            edge_width_scale=args.edge_width_scale,
            abundance_min_area=args.abundance_min_area,
            abundance_max_area=args.abundance_max_area,
            abundance_scale_power=args.abundance_scale_power,
        )
    if "abundance_all" in modes:
        out = os.path.join(args.outdir, "network_abundance_POS_ALL.svg")
        plot_abundance(
            G_all, pos_all, out,
            edge_width_scale=args.edge_width_scale,
            abundance_min_area=args.abundance_min_area,
            abundance_max_area=args.abundance_max_area,
            abundance_scale_power=args.abundance_scale_power,
        )

    if "betweenness_sub" in modes:
        out = os.path.join(args.outdir, "network_betweenness_POS_SUB.svg")
        plot_betweenness(
            G_sub, pos_sub, out,
            degree_scale=args.degree_scale,
            edge_width_scale=args.edge_width_scale,
            betweenness_threshold=args.betweenness_threshold,
            require_mag_pair=False,
            label=False,
            title="SPIEC-EASI Network (POS_SUB)\nNode color: Betweenness class | Node size: Degree",
        )
    if "betweenness_sub_labeled" in modes:
        out = os.path.join(args.outdir, "network_betweenness_POS_SUB_LABELED.svg")
        plot_betweenness(
            G_sub, pos_sub, out,
            degree_scale=args.degree_scale,
            edge_width_scale=args.edge_width_scale,
            betweenness_threshold=args.betweenness_threshold,
            require_mag_pair=False,
            label=True,
            title="SPIEC-EASI Network (POS_SUB)\nNode color: Betweenness class | Node size: Degree (Labeled)",
        )
    if "betweenness_all" in modes:
        out = os.path.join(args.outdir, "network_betweenness_POS_ALL.svg")
        plot_betweenness(
            G_all, pos_all, out,
            degree_scale=args.degree_scale,
            edge_width_scale=args.edge_width_scale,
            betweenness_threshold=args.betweenness_threshold,
            require_mag_pair=False,
            label=False,
            title="SPIEC-EASI Network (POS_ALL)\nNode color: Betweenness class | Node size: Degree",
        )
    if "betweenness_all_labeled" in modes:
        out = os.path.join(args.outdir, "network_betweenness_POS_ALL_LABELED.svg")
        plot_betweenness(
            G_all, pos_all, out,
            degree_scale=args.degree_scale,
            edge_width_scale=args.edge_width_scale,
            betweenness_threshold=args.betweenness_threshold,
            require_mag_pair=False,
            label=True,
            title="SPIEC-EASI Network (POS_ALL)\nNode color: Betweenness class | Node size: Degree (Labeled)",
        )
    if "betweenness_mag_sub" in modes:
        out = os.path.join(args.outdir, "network_betweenness_MAG_POS_SUB.svg")
        plot_betweenness(
            G_sub, pos_sub, out,
            degree_scale=args.degree_scale,
            edge_width_scale=args.edge_width_scale,
            betweenness_threshold=args.betweenness_threshold,
            require_mag_pair=True,
            label=False,
            title="SPIEC-EASI Network (POS_SUB)\nNode color: Betweenness class (paired MAG ASVs only) | Node size: Degree",
        )
    if "betweenness_mag_sub_labeled" in modes:
        out = os.path.join(args.outdir, "network_betweenness_MAG_POS_SUB_LABELED.svg")
        plot_betweenness(
            G_sub, pos_sub, out,
            degree_scale=args.degree_scale,
            edge_width_scale=args.edge_width_scale,
            betweenness_threshold=args.betweenness_threshold,
            require_mag_pair=True,
            label=True,
            title="SPIEC-EASI Network (POS_SUB)\nNode color: Betweenness class (paired MAG ASVs only) | Node size: Degree (Labeled)",
        )
    if "betweenness_mag_all" in modes:
        out = os.path.join(args.outdir, "network_betweenness_MAG_POS_ALL.svg")
        plot_betweenness(
            G_all, pos_all, out,
            degree_scale=args.degree_scale,
            edge_width_scale=args.edge_width_scale,
            betweenness_threshold=args.betweenness_threshold,
            require_mag_pair=True,
            label=False,
            title="SPIEC-EASI Network (POS_ALL)\nNode color: Betweenness class (paired MAG ASVs only) | Node size: Degree",
        )
    if "betweenness_mag_all_labeled" in modes:
        out = os.path.join(args.outdir, "network_betweenness_MAG_POS_ALL_LABELED.svg")
        plot_betweenness(
            G_all, pos_all, out,
            degree_scale=args.degree_scale,
            edge_width_scale=args.edge_width_scale,
            betweenness_threshold=args.betweenness_threshold,
            require_mag_pair=True,
            label=True,
            title="SPIEC-EASI Network (POS_ALL)\nNode color: Betweenness class (paired MAG ASVs only) | Node size: Degree (Labeled)",
        )

    # Group1 ISA
    if "group1_isa" in modes:
        out = os.path.join(args.outdir, f"network_{group1_slug}_ISA.svg")
        plot_group_isa(
            G_sub, pos_sub, out, group1_palette,
            color_attr="group1_color", size_attr="AxB_group1", label_attr="group1_label",
            isa_scale=args.isa_scale, edge_width_scale=args.edge_width_scale, label=False,
            title=f"SPIEC-EASI Network\nNode color: {group1_name} ISA | Node size: Indicator Species Strength",
            legend_title=f"{group1_name} ISA",
            legend_order=group1_order
        )
    if "group1_isa_labeled" in modes:
        out = os.path.join(args.outdir, f"network_{group1_slug}_ISA_LABELED.svg")
        plot_group_isa(
            G_sub, pos_sub, out, group1_palette,
            color_attr="group1_color", size_attr="AxB_group1", label_attr="group1_label",
            isa_scale=args.isa_scale, edge_width_scale=args.edge_width_scale, label=True,
            title=f"SPIEC-EASI Network\nNode color: {group1_name} ISA | Node size: Indicator Species Strength (Labeled)",
            legend_title=f"{group1_name} ISA",
            legend_order=group1_order
        )
    if "group1_isa_all" in modes:
        out = os.path.join(args.outdir, f"network_{group1_slug}_ISA_POS_ALL.svg")
        plot_group_isa(
            G_all, pos_all, out, group1_palette,
            color_attr="group1_color", size_attr="AxB_group1", label_attr="group1_label",
            isa_scale=args.isa_scale, edge_width_scale=args.edge_width_scale, label=False,
            title=f"SPIEC-EASI Network (POS_ALL)\nNode color: {group1_name} ISA | Node size: Indicator Species Strength",
            legend_title=f"{group1_name} ISA",
            legend_order=group1_order
        )
    if "group1_isa_all_labeled" in modes:
        out = os.path.join(args.outdir, f"network_{group1_slug}_ISA_POS_ALL_LABELED.svg")
        plot_group_isa(
            G_all, pos_all, out, group1_palette,
            color_attr="group1_color", size_attr="AxB_group1", label_attr="group1_label",
            isa_scale=args.isa_scale, edge_width_scale=args.edge_width_scale, label=True,
            title=f"SPIEC-EASI Network (POS_ALL)\nNode color: {group1_name} ISA | Node size: Indicator Species Strength (Labeled)",
            legend_title=f"{group1_name} ISA",
            legend_order=group1_order
        )
    if "group1_isa_mag" in modes:
        out = os.path.join(args.outdir, f"network_{group1_slug}_ISA_MAG.svg")
        plot_group_isa(
            G_sub, pos_sub, out, group1_palette,
            color_attr="group1_color", size_attr="AxB_group1", label_attr="group1_label",
            isa_scale=args.isa_scale, edge_width_scale=args.edge_width_scale, label=False,
            title=f"SPIEC-EASI Network\nNode color: {group1_name} ISA (paired MAG ASVs only) | Node size: Indicator Species Strength",
            legend_title=f"{group1_name} ISA",
            legend_order=group1_order,
            require_mag_pair=True,
        )
    if "group1_isa_mag_labeled" in modes:
        out = os.path.join(args.outdir, f"network_{group1_slug}_ISA_MAG_LABELED.svg")
        plot_group_isa(
            G_sub, pos_sub, out, group1_palette,
            color_attr="group1_color", size_attr="AxB_group1", label_attr="group1_label",
            isa_scale=args.isa_scale, edge_width_scale=args.edge_width_scale, label=True,
            title=f"SPIEC-EASI Network\nNode color: {group1_name} ISA (paired MAG ASVs only) | Node size: Indicator Species Strength (Labeled)",
            legend_title=f"{group1_name} ISA",
            legend_order=group1_order,
            require_mag_pair=True,
        )
    if "group1_isa_mag_all" in modes:
        out = os.path.join(args.outdir, f"network_{group1_slug}_ISA_MAG_POS_ALL.svg")
        plot_group_isa(
            G_all, pos_all, out, group1_palette,
            color_attr="group1_color", size_attr="AxB_group1", label_attr="group1_label",
            isa_scale=args.isa_scale, edge_width_scale=args.edge_width_scale, label=False,
            title=f"SPIEC-EASI Network (POS_ALL)\nNode color: {group1_name} ISA (paired MAG ASVs only) | Node size: Indicator Species Strength",
            legend_title=f"{group1_name} ISA",
            legend_order=group1_order,
            require_mag_pair=True,
        )
    if "group1_isa_mag_all_labeled" in modes:
        out = os.path.join(args.outdir, f"network_{group1_slug}_ISA_MAG_POS_ALL_LABELED.svg")
        plot_group_isa(
            G_all, pos_all, out, group1_palette,
            color_attr="group1_color", size_attr="AxB_group1", label_attr="group1_label",
            isa_scale=args.isa_scale, edge_width_scale=args.edge_width_scale, label=True,
            title=f"SPIEC-EASI Network (POS_ALL)\nNode color: {group1_name} ISA (paired MAG ASVs only) | Node size: Indicator Species Strength (Labeled)",
            legend_title=f"{group1_name} ISA",
            legend_order=group1_order,
            require_mag_pair=True,
        )
    if "group1_isa_focus" in modes and focus_group1_label:
        out = os.path.join(args.outdir, f"network_{group1_slug}_ISA_FOCUS.svg")
        plot_group_isa(
            G_sub, pos_sub, out, group1_palette,
            color_attr="group1_color", size_attr="AxB_group1", label_attr="group1_label",
            isa_scale=args.isa_scale, edge_width_scale=args.edge_width_scale, label=False,
            title=f"SPIEC-EASI Network\nNode color: {group1_name} ISA ({focus_group1_label}) | Node size: Indicator Species Strength",
            legend_title=f"{group1_name} ISA ({focus_group1_label})",
            legend_order=group1_focus_order if group1_focus_order else group1_order,
            focus_label=focus_group1_label,
            all_combo_label=group1_focus_all_combo,
        )
    if "group1_isa_focus_labeled" in modes and focus_group1_label:
        out = os.path.join(args.outdir, f"network_{group1_slug}_ISA_FOCUS_LABELED.svg")
        plot_group_isa(
            G_sub, pos_sub, out, group1_palette,
            color_attr="group1_color", size_attr="AxB_group1", label_attr="group1_label",
            isa_scale=args.isa_scale, edge_width_scale=args.edge_width_scale, label=True,
            title=f"SPIEC-EASI Network\nNode color: {group1_name} ISA ({focus_group1_label}) | Node size: Indicator Species Strength (Labeled)",
            legend_title=f"{group1_name} ISA ({focus_group1_label})",
            legend_order=group1_focus_order if group1_focus_order else group1_order,
            focus_label=focus_group1_label,
            all_combo_label=group1_focus_all_combo,
        )
    if "group1_isa_focus_all" in modes and focus_group1_label:
        out = os.path.join(args.outdir, f"network_{group1_slug}_ISA_FOCUS_POS_ALL.svg")
        plot_group_isa(
            G_all, pos_all, out, group1_palette,
            color_attr="group1_color", size_attr="AxB_group1", label_attr="group1_label",
            isa_scale=args.isa_scale, edge_width_scale=args.edge_width_scale, label=False,
            title=f"SPIEC-EASI Network (POS_ALL)\nNode color: {group1_name} ISA ({focus_group1_label}) | Node size: Indicator Species Strength",
            legend_title=f"{group1_name} ISA ({focus_group1_label})",
            legend_order=group1_focus_order if group1_focus_order else group1_order,
            focus_label=focus_group1_label,
            all_combo_label=group1_focus_all_combo,
        )
    if "group1_isa_focus_all_labeled" in modes and focus_group1_label:
        out = os.path.join(args.outdir, f"network_{group1_slug}_ISA_FOCUS_POS_ALL_LABELED.svg")
        plot_group_isa(
            G_all, pos_all, out, group1_palette,
            color_attr="group1_color", size_attr="AxB_group1", label_attr="group1_label",
            isa_scale=args.isa_scale, edge_width_scale=args.edge_width_scale, label=True,
            title=f"SPIEC-EASI Network (POS_ALL)\nNode color: {group1_name} ISA ({focus_group1_label}) | Node size: Indicator Species Strength (Labeled)",
            legend_title=f"{group1_name} ISA ({focus_group1_label})",
            legend_order=group1_focus_order if group1_focus_order else group1_order,
            focus_label=focus_group1_label,
            all_combo_label=group1_focus_all_combo,
        )

    # Group2 ISA
    if "group2_isa" in modes:
        out = os.path.join(args.outdir, f"network_{group2_slug}_ISA.svg")
        plot_group_isa(
            G_sub, pos_sub, out, group2_palette,
            color_attr="group2_color", size_attr="AxB_group2", label_attr="group2_label",
            isa_scale=args.isa_scale, edge_width_scale=args.edge_width_scale, label=False,
            title=f"SPIEC-EASI Network\nNode color: {group2_name} ISA | Node size: Indicator Species Strength",
            legend_title=f"{group2_name} ISA",
            legend_order=group2_order,
            focus_label=focus_group2_label,
            all_combo_label=group2_all_combo or "",
        )
    if "group2_isa_labeled" in modes:
        out = os.path.join(args.outdir, f"network_{group2_slug}_ISA_LABELED.svg")
        plot_group_isa(
            G_sub, pos_sub, out, group2_palette,
            color_attr="group2_color", size_attr="AxB_group2", label_attr="group2_label",
            isa_scale=args.isa_scale, edge_width_scale=args.edge_width_scale, label=True,
            title=f"SPIEC-EASI Network\nNode color: {group2_name} ISA | Node size: Indicator Species Strength (Labeled)",
            legend_title=f"{group2_name} ISA",
            legend_order=group2_order,
            focus_label=focus_group2_label,
            all_combo_label=group2_all_combo or "",
        )
    if "group2_isa_all" in modes:
        out = os.path.join(args.outdir, f"network_{group2_slug}_ISA_POS_ALL.svg")
        plot_group_isa(
            G_all, pos_all, out, group2_palette,
            color_attr="group2_color", size_attr="AxB_group2", label_attr="group2_label",
            isa_scale=args.isa_scale, edge_width_scale=args.edge_width_scale, label=False,
            title=f"SPIEC-EASI Network (POS_ALL)\nNode color: {group2_name} ISA | Node size: Indicator Species Strength",
            legend_title=f"{group2_name} ISA",
            legend_order=group2_order,
            focus_label=focus_group2_label,
            all_combo_label=group2_all_combo or "",
        )
    if "group2_isa_all_labeled" in modes:
        out = os.path.join(args.outdir, f"network_{group2_slug}_ISA_POS_ALL_LABELED.svg")
        plot_group_isa(
            G_all, pos_all, out, group2_palette,
            color_attr="group2_color", size_attr="AxB_group2", label_attr="group2_label",
            isa_scale=args.isa_scale, edge_width_scale=args.edge_width_scale, label=True,
            title=f"SPIEC-EASI Network (POS_ALL)\nNode color: {group2_name} ISA | Node size: Indicator Species Strength (Labeled)",
            legend_title=f"{group2_name} ISA",
            legend_order=group2_order,
            focus_label=focus_group2_label,
            all_combo_label=group2_all_combo or "",
        )
    if "group2_isa_mag" in modes:
        out = os.path.join(args.outdir, f"network_{group2_slug}_ISA_MAG.svg")
        plot_group_isa(
            G_sub, pos_sub, out, group2_palette,
            color_attr="group2_color", size_attr="AxB_group2", label_attr="group2_label",
            isa_scale=args.isa_scale, edge_width_scale=args.edge_width_scale, label=False,
            title=f"SPIEC-EASI Network\nNode color: {group2_name} ISA (paired MAG ASVs only) | Node size: Indicator Species Strength",
            legend_title=f"{group2_name} ISA",
            legend_order=group2_order,
            focus_label=focus_group2_label,
            all_combo_label=group2_all_combo or "",
            require_mag_pair=True,
        )
    if "group2_isa_mag_labeled" in modes:
        out = os.path.join(args.outdir, f"network_{group2_slug}_ISA_MAG_LABELED.svg")
        plot_group_isa(
            G_sub, pos_sub, out, group2_palette,
            color_attr="group2_color", size_attr="AxB_group2", label_attr="group2_label",
            isa_scale=args.isa_scale, edge_width_scale=args.edge_width_scale, label=True,
            title=f"SPIEC-EASI Network\nNode color: {group2_name} ISA (paired MAG ASVs only) | Node size: Indicator Species Strength (Labeled)",
            legend_title=f"{group2_name} ISA",
            legend_order=group2_order,
            focus_label=focus_group2_label,
            all_combo_label=group2_all_combo or "",
            require_mag_pair=True,
        )
    if "group2_isa_mag_all" in modes:
        out = os.path.join(args.outdir, f"network_{group2_slug}_ISA_MAG_POS_ALL.svg")
        plot_group_isa(
            G_all, pos_all, out, group2_palette,
            color_attr="group2_color", size_attr="AxB_group2", label_attr="group2_label",
            isa_scale=args.isa_scale, edge_width_scale=args.edge_width_scale, label=False,
            title=f"SPIEC-EASI Network (POS_ALL)\nNode color: {group2_name} ISA (paired MAG ASVs only) | Node size: Indicator Species Strength",
            legend_title=f"{group2_name} ISA",
            legend_order=group2_order,
            focus_label=focus_group2_label,
            all_combo_label=group2_all_combo or "",
            require_mag_pair=True,
        )
    if "group2_isa_mag_all_labeled" in modes:
        out = os.path.join(args.outdir, f"network_{group2_slug}_ISA_MAG_POS_ALL_LABELED.svg")
        plot_group_isa(
            G_all, pos_all, out, group2_palette,
            color_attr="group2_color", size_attr="AxB_group2", label_attr="group2_label",
            isa_scale=args.isa_scale, edge_width_scale=args.edge_width_scale, label=True,
            title=f"SPIEC-EASI Network (POS_ALL)\nNode color: {group2_name} ISA (paired MAG ASVs only) | Node size: Indicator Species Strength (Labeled)",
            legend_title=f"{group2_name} ISA",
            legend_order=group2_order,
            focus_label=focus_group2_label,
            all_combo_label=group2_all_combo or "",
            require_mag_pair=True,
        )

    for spec in isa_specs[2:]:
        variants = isa_mode_variants(spec["index"])
        palette = spec["palette"]
        slug = spec["slug"]
        name = spec["name"]
        order = spec["order"]
        focus_label = spec["focus_label"]
        focus_order = spec["focus_order"]
        all_combo = spec["all_combo_label"]
        color_attr = spec["color_attr"]
        size_attr = spec["score_attr"]
        label_attr = spec["label_attr"]

        if variants["isa"] in modes:
            out = os.path.join(args.outdir, f"network_{slug}_ISA.svg")
            plot_group_isa(
                G_sub, pos_sub, out, palette,
                color_attr=color_attr, size_attr=size_attr, label_attr=label_attr,
                isa_scale=args.isa_scale, edge_width_scale=args.edge_width_scale, label=False,
                title=f"SPIEC-EASI Network\nNode color: {name} ISA | Node size: Indicator Species Strength",
                legend_title=f"{name} ISA",
                legend_order=order,
            )
        if variants["isa_labeled"] in modes:
            out = os.path.join(args.outdir, f"network_{slug}_ISA_LABELED.svg")
            plot_group_isa(
                G_sub, pos_sub, out, palette,
                color_attr=color_attr, size_attr=size_attr, label_attr=label_attr,
                isa_scale=args.isa_scale, edge_width_scale=args.edge_width_scale, label=True,
                title=f"SPIEC-EASI Network\nNode color: {name} ISA | Node size: Indicator Species Strength (Labeled)",
                legend_title=f"{name} ISA",
                legend_order=order,
            )
        if variants["isa_all"] in modes:
            out = os.path.join(args.outdir, f"network_{slug}_ISA_POS_ALL.svg")
            plot_group_isa(
                G_all, pos_all, out, palette,
                color_attr=color_attr, size_attr=size_attr, label_attr=label_attr,
                isa_scale=args.isa_scale, edge_width_scale=args.edge_width_scale, label=False,
                title=f"SPIEC-EASI Network (POS_ALL)\nNode color: {name} ISA | Node size: Indicator Species Strength",
                legend_title=f"{name} ISA",
                legend_order=order,
            )
        if variants["isa_all_labeled"] in modes:
            out = os.path.join(args.outdir, f"network_{slug}_ISA_POS_ALL_LABELED.svg")
            plot_group_isa(
                G_all, pos_all, out, palette,
                color_attr=color_attr, size_attr=size_attr, label_attr=label_attr,
                isa_scale=args.isa_scale, edge_width_scale=args.edge_width_scale, label=True,
                title=f"SPIEC-EASI Network (POS_ALL)\nNode color: {name} ISA | Node size: Indicator Species Strength (Labeled)",
                legend_title=f"{name} ISA",
                legend_order=order,
            )
        if variants["isa_mag"] in modes:
            out = os.path.join(args.outdir, f"network_{slug}_ISA_MAG.svg")
            plot_group_isa(
                G_sub, pos_sub, out, palette,
                color_attr=color_attr, size_attr=size_attr, label_attr=label_attr,
                isa_scale=args.isa_scale, edge_width_scale=args.edge_width_scale, label=False,
                title=f"SPIEC-EASI Network\nNode color: {name} ISA (paired MAG ASVs only) | Node size: Indicator Species Strength",
                legend_title=f"{name} ISA",
                legend_order=order,
                require_mag_pair=True,
            )
        if variants["isa_mag_labeled"] in modes:
            out = os.path.join(args.outdir, f"network_{slug}_ISA_MAG_LABELED.svg")
            plot_group_isa(
                G_sub, pos_sub, out, palette,
                color_attr=color_attr, size_attr=size_attr, label_attr=label_attr,
                isa_scale=args.isa_scale, edge_width_scale=args.edge_width_scale, label=True,
                title=f"SPIEC-EASI Network\nNode color: {name} ISA (paired MAG ASVs only) | Node size: Indicator Species Strength (Labeled)",
                legend_title=f"{name} ISA",
                legend_order=order,
                require_mag_pair=True,
            )
        if variants["isa_mag_all"] in modes:
            out = os.path.join(args.outdir, f"network_{slug}_ISA_MAG_POS_ALL.svg")
            plot_group_isa(
                G_all, pos_all, out, palette,
                color_attr=color_attr, size_attr=size_attr, label_attr=label_attr,
                isa_scale=args.isa_scale, edge_width_scale=args.edge_width_scale, label=False,
                title=f"SPIEC-EASI Network (POS_ALL)\nNode color: {name} ISA (paired MAG ASVs only) | Node size: Indicator Species Strength",
                legend_title=f"{name} ISA",
                legend_order=order,
                require_mag_pair=True,
            )
        if variants["isa_mag_all_labeled"] in modes:
            out = os.path.join(args.outdir, f"network_{slug}_ISA_MAG_POS_ALL_LABELED.svg")
            plot_group_isa(
                G_all, pos_all, out, palette,
                color_attr=color_attr, size_attr=size_attr, label_attr=label_attr,
                isa_scale=args.isa_scale, edge_width_scale=args.edge_width_scale, label=True,
                title=f"SPIEC-EASI Network (POS_ALL)\nNode color: {name} ISA (paired MAG ASVs only) | Node size: Indicator Species Strength (Labeled)",
                legend_title=f"{name} ISA",
                legend_order=order,
                require_mag_pair=True,
            )
        if variants["isa_focus"] in modes and focus_label:
            out = os.path.join(args.outdir, f"network_{slug}_ISA_FOCUS.svg")
            plot_group_isa(
                G_sub, pos_sub, out, palette,
                color_attr=color_attr, size_attr=size_attr, label_attr=label_attr,
                isa_scale=args.isa_scale, edge_width_scale=args.edge_width_scale, label=False,
                title=f"SPIEC-EASI Network\nNode color: {name} ISA ({focus_label}) | Node size: Indicator Species Strength",
                legend_title=f"{name} ISA ({focus_label})",
                legend_order=focus_order if focus_order else order,
                focus_label=focus_label,
                all_combo_label=all_combo,
            )
        if variants["isa_focus_labeled"] in modes and focus_label:
            out = os.path.join(args.outdir, f"network_{slug}_ISA_FOCUS_LABELED.svg")
            plot_group_isa(
                G_sub, pos_sub, out, palette,
                color_attr=color_attr, size_attr=size_attr, label_attr=label_attr,
                isa_scale=args.isa_scale, edge_width_scale=args.edge_width_scale, label=True,
                title=f"SPIEC-EASI Network\nNode color: {name} ISA ({focus_label}) | Node size: Indicator Species Strength (Labeled)",
                legend_title=f"{name} ISA ({focus_label})",
                legend_order=focus_order if focus_order else order,
                focus_label=focus_label,
                all_combo_label=all_combo,
            )
        if variants["isa_focus_all"] in modes and focus_label:
            out = os.path.join(args.outdir, f"network_{slug}_ISA_FOCUS_POS_ALL.svg")
            plot_group_isa(
                G_all, pos_all, out, palette,
                color_attr=color_attr, size_attr=size_attr, label_attr=label_attr,
                isa_scale=args.isa_scale, edge_width_scale=args.edge_width_scale, label=False,
                title=f"SPIEC-EASI Network (POS_ALL)\nNode color: {name} ISA ({focus_label}) | Node size: Indicator Species Strength",
                legend_title=f"{name} ISA ({focus_label})",
                legend_order=focus_order if focus_order else order,
                focus_label=focus_label,
                all_combo_label=all_combo,
            )
        if variants["isa_focus_all_labeled"] in modes and focus_label:
            out = os.path.join(args.outdir, f"network_{slug}_ISA_FOCUS_POS_ALL_LABELED.svg")
            plot_group_isa(
                G_all, pos_all, out, palette,
                color_attr=color_attr, size_attr=size_attr, label_attr=label_attr,
                isa_scale=args.isa_scale, edge_width_scale=args.edge_width_scale, label=True,
                title=f"SPIEC-EASI Network (POS_ALL)\nNode color: {name} ISA ({focus_label}) | Node size: Indicator Species Strength (Labeled)",
                legend_title=f"{name} ISA ({focus_label})",
                legend_order=focus_order if focus_order else order,
                focus_label=focus_label,
                all_combo_label=all_combo,
            )

    # Module overlays
    module_filter_attr = "module_plot_keep" if args.module_isa_only else None
    module_legend_title = (
        "ISA-associated ecological modules"
        if (args.module_isa_only or args.module_color_by_isa)
        else "Ecological modules"
    )
    if "module_sub" in modes:
        if modules_sub.empty:
            print("[WARN] module_sub requested, but no module assignments were loaded.")
        else:
            out = os.path.join(args.outdir, "network_modules_POS_SUB.svg")
            plot_modules(
                G_sub, pos_sub, out,
                color_attr="module_color_plot", label_attr="module_label_plot",
                filter_attr=module_filter_attr,
                degree_scale=args.degree_scale, edge_width_scale=args.edge_width_scale,
                label=False,
                title=f"SPIEC-EASI Network (POS_SUB)\nNode color: {module_legend_title} | Node size: Degree",
                legend_title=module_legend_title
            )
    if "module_sub_labeled" in modes:
        if modules_sub.empty:
            print("[WARN] module_sub_labeled requested, but no module assignments were loaded.")
        else:
            out = os.path.join(args.outdir, "network_modules_POS_SUB_LABELED.svg")
            plot_modules(
                G_sub, pos_sub, out,
                color_attr="module_color_plot", label_attr="module_label_plot",
                filter_attr=module_filter_attr,
                degree_scale=args.degree_scale, edge_width_scale=args.edge_width_scale,
                label=True,
                title=f"SPIEC-EASI Network (POS_SUB)\nNode color: {module_legend_title} | Node size: Degree (Labeled)",
                legend_title=module_legend_title
            )
    if "module_all" in modes:
        if modules_all.empty:
            print("[WARN] module_all requested, but no module assignments were loaded.")
        else:
            out = os.path.join(args.outdir, "network_modules_POS_ALL.svg")
            plot_modules(
                G_all, pos_all, out,
                color_attr="module_color_plot", label_attr="module_label_plot",
                filter_attr=module_filter_attr,
                degree_scale=args.degree_scale, edge_width_scale=args.edge_width_scale,
                label=False,
                title=f"SPIEC-EASI Network (POS_ALL)\nNode color: {module_legend_title} | Node size: Degree",
                legend_title=module_legend_title
            )
    if "module_all_labeled" in modes:
        if modules_all.empty:
            print("[WARN] module_all_labeled requested, but no module assignments were loaded.")
        else:
            out = os.path.join(args.outdir, "network_modules_POS_ALL_LABELED.svg")
            plot_modules(
                G_all, pos_all, out,
                color_attr="module_color_plot", label_attr="module_label_plot",
                filter_attr=module_filter_attr,
                degree_scale=args.degree_scale, edge_width_scale=args.edge_width_scale,
                label=True,
                title=f"SPIEC-EASI Network (POS_ALL)\nNode color: {module_legend_title} | Node size: Degree (Labeled)",
                legend_title=module_legend_title
            )

    # Renewal was tested at the ecological-module level, not as a direct ASV
    # indicator analysis.  The overlay therefore colors all ASVs belonging to
    # modules with supported renewal-phase contrasts and leaves other modules
    # gray.  This preserves the correct inferential unit in the figure itself.
    if args.module_renewal_association and not modules_all.empty:
        out = os.path.join(args.outdir, "network_ecological_module_renewal_phase_POS_ALL.svg")
        plot_modules(
            G_all, pos_all, out,
            color_attr="renewal_module_color", label_attr="renewal_module_label",
            filter_attr=None,
            degree_scale=args.degree_scale, edge_width_scale=args.edge_width_scale,
            label=False,
            title=("SPIEC-EASI Network (POS_ALL)\n"
                   "Node color: significant ecological-module renewal contrast "
                   "(baseline shown as stagnation) | Node size: Degree"),
            legend_title="Module-level renewal association",
        )
        active_out = os.path.join(
            args.outdir,
            "network_ecological_module_renewal_post_renewal_vs_stagnation_POS_ALL.svg",
        )
        plot_modules(
            G_all, pos_all, active_out,
            color_attr="active_renewal_module_color",
            label_attr="active_renewal_module_label",
            filter_attr=None,
            degree_scale=args.degree_scale,
            edge_width_scale=args.edge_width_scale,
            label=False,
            title=("SPIEC-EASI Network (POS_ALL)\n"
                   "Node color: direction of significant ecological-module difference between "
                   "renewal/post-renewal and stagnation | Node size: Degree"),
            legend_title="Direction of module-level difference",
        )

    if "mag_pair_sub" in modes:
        if asv_mag_pairing.empty:
            print("[WARN] mag_pair_sub requested, but no ASV-MAG pairing table was loaded.")
        else:
            out = os.path.join(args.outdir, "network_mag_pair_POS_SUB.svg")
            plot_mag_pairing(
                G_sub, pos_sub, out,
                degree_scale=args.degree_scale,
                edge_width_scale=args.edge_width_scale,
                label=False,
                title="SPIEC-EASI Network (POS_SUB)\nNode color: ASV with paired MAG | Node size: Degree",
            )
    if "mag_pair_sub_labeled" in modes:
        if asv_mag_pairing.empty:
            print("[WARN] mag_pair_sub_labeled requested, but no ASV-MAG pairing table was loaded.")
        else:
            out = os.path.join(args.outdir, "network_mag_pair_POS_SUB_LABELED.svg")
            plot_mag_pairing(
                G_sub, pos_sub, out,
                degree_scale=args.degree_scale,
                edge_width_scale=args.edge_width_scale,
                label=True,
                title="SPIEC-EASI Network (POS_SUB)\nNode color: ASV with paired MAG | Node size: Degree (Labeled)",
            )
    if "mag_pair_all" in modes:
        if asv_mag_pairing.empty:
            print("[WARN] mag_pair_all requested, but no ASV-MAG pairing table was loaded.")
        else:
            out = os.path.join(args.outdir, "network_mag_pair_POS_ALL.svg")
            plot_mag_pairing(
                G_all, pos_all, out,
                degree_scale=args.degree_scale,
                edge_width_scale=args.edge_width_scale,
                label=False,
                title="SPIEC-EASI Network (POS_ALL)\nNode color: ASV with paired MAG | Node size: Degree",
            )
    if "mag_pair_all_labeled" in modes:
        if asv_mag_pairing.empty:
            print("[WARN] mag_pair_all_labeled requested, but no ASV-MAG pairing table was loaded.")
        else:
            out = os.path.join(args.outdir, "network_mag_pair_POS_ALL_LABELED.svg")
            plot_mag_pairing(
                G_all, pos_all, out,
                degree_scale=args.degree_scale,
                edge_width_scale=args.edge_width_scale,
                label=True,
                title="SPIEC-EASI Network (POS_ALL)\nNode color: ASV with paired MAG | Node size: Degree (Labeled)",
            )
    if "mag_pair_tax_sub" in modes:
        if asv_mag_pairing.empty:
            print("[WARN] mag_pair_tax_sub requested, but no ASV-MAG pairing table was loaded.")
        else:
            out = os.path.join(args.outdir, "network_mag_pair_tax_POS_SUB.svg")
            plot_mag_pairing_taxonomy(
                G_sub, pos_sub, out,
                degree_scale=args.degree_scale,
                edge_width_scale=args.edge_width_scale,
                label=False,
                title="SPIEC-EASI Network (POS_SUB)\nNode color: paired MAG phylum | Node size: Degree",
            )
    if "mag_pair_tax_sub_labeled" in modes:
        if asv_mag_pairing.empty:
            print("[WARN] mag_pair_tax_sub_labeled requested, but no ASV-MAG pairing table was loaded.")
        else:
            out = os.path.join(args.outdir, "network_mag_pair_tax_POS_SUB_LABELED.svg")
            plot_mag_pairing_taxonomy(
                G_sub, pos_sub, out,
                degree_scale=args.degree_scale,
                edge_width_scale=args.edge_width_scale,
                label=True,
                title="SPIEC-EASI Network (POS_SUB)\nNode color: paired MAG phylum | Node size: Degree (Labeled)",
            )
    if "mag_pair_tax_all" in modes:
        if asv_mag_pairing.empty:
            print("[WARN] mag_pair_tax_all requested, but no ASV-MAG pairing table was loaded.")
        else:
            out = os.path.join(args.outdir, "network_mag_pair_tax_POS_ALL.svg")
            plot_mag_pairing_taxonomy(
                G_all, pos_all, out,
                degree_scale=args.degree_scale,
                edge_width_scale=args.edge_width_scale,
                label=False,
                title="SPIEC-EASI Network (POS_ALL)\nNode color: paired MAG phylum | Node size: Degree",
            )
    if "mag_pair_tax_all_labeled" in modes:
        if asv_mag_pairing.empty:
            print("[WARN] mag_pair_tax_all_labeled requested, but no ASV-MAG pairing table was loaded.")
        else:
            out = os.path.join(args.outdir, "network_mag_pair_tax_POS_ALL_LABELED.svg")
            plot_mag_pairing_taxonomy(
                G_all, pos_all, out,
                degree_scale=args.degree_scale,
                edge_width_scale=args.edge_width_scale,
                label=True,
                title="SPIEC-EASI Network (POS_ALL)\nNode color: paired MAG phylum | Node size: Degree (Labeled)",
            )

    # Phylum × {abundance, ISA}
    if "phylum_abund" in modes or "phylum_abund_mean" in modes:
        out = os.path.join(args.outdir, "network_phylum_ABUND.svg")
        plot_phylum(G_sub, pos_sub, out, phylum_palette, size_attr='mean', size_label='Mean abundance', size_scale=1.0, edge_width_scale=args.edge_width_scale, abundance_min_area=args.abundance_min_area, abundance_max_area=args.abundance_max_area, abundance_scale_power=args.abundance_scale_power)
        out_mean = os.path.join(args.outdir, "network_phylum_ABUND_MEAN.svg")
        plot_phylum(G_sub, pos_sub, out_mean, phylum_palette, size_attr='mean', size_label='Mean abundance', size_scale=1.0, edge_width_scale=args.edge_width_scale, abundance_min_area=args.abundance_min_area, abundance_max_area=args.abundance_max_area, abundance_scale_power=args.abundance_scale_power)
    if "phylum_abund_all" in modes or "phylum_abund_all_mean" in modes:
        out = os.path.join(args.outdir, "network_phylum_ABUND_POS_ALL.svg")
        plot_phylum(G_all, pos_all, out, phylum_palette, size_attr='mean', size_label='Mean abundance', size_scale=1.0, edge_width_scale=args.edge_width_scale, abundance_min_area=args.abundance_min_area, abundance_max_area=args.abundance_max_area, abundance_scale_power=args.abundance_scale_power)
        out_mean = os.path.join(args.outdir, "network_phylum_ABUND_POS_ALL_MEAN.svg")
        plot_phylum(G_all, pos_all, out_mean, phylum_palette, size_attr='mean', size_label='Mean abundance', size_scale=1.0, edge_width_scale=args.edge_width_scale, abundance_min_area=args.abundance_min_area, abundance_max_area=args.abundance_max_area, abundance_scale_power=args.abundance_scale_power)

    if "phylum_abund" in modes or "phylum_abund_median" in modes:
        out_median = os.path.join(args.outdir, "network_phylum_ABUND_MEDIAN.svg")
        plot_phylum(G_sub, pos_sub, out_median, phylum_palette, size_attr='median', size_label='Median abundance', size_scale=1.0, edge_width_scale=args.edge_width_scale, abundance_min_area=args.abundance_min_area, abundance_max_area=args.abundance_max_area, abundance_scale_power=args.abundance_scale_power)
    if "phylum_abund_all" in modes or "phylum_abund_all_median" in modes:
        out_median = os.path.join(args.outdir, "network_phylum_ABUND_POS_ALL_MEDIAN.svg")
        plot_phylum(G_all, pos_all, out_median, phylum_palette, size_attr='median', size_label='Median abundance', size_scale=1.0, edge_width_scale=args.edge_width_scale, abundance_min_area=args.abundance_min_area, abundance_max_area=args.abundance_max_area, abundance_scale_power=args.abundance_scale_power)

    if "phylum_isa" in modes:
        out = os.path.join(args.outdir, "network_phylum_ISA.svg")
        plot_phylum(
            G_sub, pos_sub, out, phylum_palette,
            size_attr=phylum_isa_size_attr, size_label='ISA', size_scale=args.isa_scale, edge_width_scale=args.edge_width_scale, label=False
        )
    if "phylum_isa_labeled" in modes:
        out = os.path.join(args.outdir, "network_phylum_ISA_LABELED.svg")
        plot_phylum(
            G_sub, pos_sub, out, phylum_palette,
            size_attr=phylum_isa_size_attr, size_label='ISA', size_scale=args.isa_scale, edge_width_scale=args.edge_width_scale, label=True
        )
    if "phylum_isa_all" in modes:
        out = os.path.join(args.outdir, "network_phylum_ISA_POS_ALL.svg")
        plot_phylum(
            G_all, pos_all, out, phylum_palette,
            size_attr=phylum_isa_size_attr, size_label='ISA', size_scale=args.isa_scale, edge_width_scale=args.edge_width_scale, label=False
        )
    if "phylum_isa_all_labeled" in modes:
        out = os.path.join(args.outdir, "network_phylum_ISA_POS_ALL_LABELED.svg")
        plot_phylum(
            G_all, pos_all, out, phylum_palette,
            size_attr=phylum_isa_size_attr, size_label='ISA', size_scale=args.isa_scale, edge_width_scale=args.edge_width_scale, label=True
        )

    ok("All done.")


if __name__ == "__main__":
    main()
