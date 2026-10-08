#!/usr/bin/env python3
"""
data_loss_sankey.py
Build Sankey diagrams for read/ASV flow with flexible I/O, grouping, and colors.

Two modes:
A) COMPUTE from pipeline TSVs (default)
B) MANUAL via --steps/--lmp-in/--lmp-out

Examples
--------
# A) Compute from files (defaults mirror your script paths/columns)
python data_loss_sankey.py \
  --data-dir /path/to/project \
  --sub-dir spark_combined_output \
  --metadata /path/to/project/ref_db/spark_metadata.tsv \
  --group1-col type_group \
  --samp-col lmp_id \
  --keep-types "Oral Rinse,Lung Brush,BAL,Skin Brush,Scope Flush" \
  --fastq-stats stats/fastq_stats.tsv \
  --filtered-stats stats/filtered_fastqs.tsv \
  --asv-raw ASVs/ASV_counts.tsv \
  --asv-decon ASVs/ASV_target.decon.tsv \
  --asv-micro ASVs/ASV_target.micro.tsv \
  --palette "Scope Flush:#E69F00,Skin Brush:#CC79A7,Lung Brush:#009E73,BAL:#0072B2,Oral Rinse:#6A3D9A,Failed-QC:lightgray" \
  --title "Data Loss Flow" \
  --output-prefix metadata/data_loss_sankey --make-labeled --make-unlabeled

# B) Manual counts
python data_loss_sankey.py \
  --steps "Quality Control:123456,Error Correction:110000,Decontamination:98000,Off-Target Filtering:82000,Finished Data:76000" \
  --lmp-in "Oral Rinse:40000,Lung Brush:35000,BAL:28000,Skin Brush:12000,Scope Flush:8400" \
  --lmp-out "Oral Rinse:18000,Lung Brush:22000,BAL:24000,Skin Brush:9000,Scope Flush:5100" \
  --palette "Oral Rinse:#6A3D9A,Lung Brush:#009E73,BAL:#0072B2,Skin Brush:#CC79A7,Scope Flush:#E69F00" \
  --output-prefix out/sankey --make-labeled --make-unlabeled
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Sequence, Optional
from xml.sax.saxutils import escape

import pandas as pd

# =========================
# Utility parsers / helpers
# =========================
def parse_kv_csv(s: str, val_cast=int) -> Dict[str, object]:
    """
    Parse 'A:1,B:2' into dict. Whitespace tolerated. Empty string -> {}.
    """
    out: Dict[str, object] = {}
    if not s:
        return out
    for item in s.split(','):
        item = item.strip()
        if not item:
            continue
        if ':' not in item:
            raise ValueError(f"Expected key:value pair, got '{item}'")
        k, v = item.split(':', 1)
        k = k.strip()
        v = v.strip()
        out[k] = val_cast(v) if val_cast is not None else v
    return out


def parse_steps_csv(s: str) -> Tuple[List[str], List[int]]:
    """
    Parse 'StepA:100,StepB:90,...' -> (['StepA','StepB',...],[100,90,...])
    """
    d = parse_kv_csv(s, val_cast=int)
    return list(d.keys()), list(d.values())


def extract_sample_id_from_path(path_str: str) -> str:
    """
    Extract a sample id from a file path.
    Returns the basename without common sequencing extensions.
    """
    base = os.path.basename(path_str)
    # Fallback: strip extensions
    stem = base
    for ext in ('.fastq.gz', '.fq.gz', '.fastq', '.fq',
                '.fasta.gz', '.fasta', '.fa.gz', '.fa',
                '.gz', '.tsv', '.csv', '.txt'):
        if stem.endswith(ext):
            stem = stem[: -len(ext)]
    if (('.filtered' in path_str) or
        ('.merged' in path_str) or
        ('.trimmed' in path_str)
        ):
        stem = re.sub(r'(\.filtered|\.merged|\.trimmed)$', '', stem)
    else:
        stem = re.sub(r'(_R[12]|_[12])?(_001)?$', '', stem)
    stem = re.sub(r'(-)$', '_', stem)

    return stem


def safe_int(x) -> int:
    try:
        return int(x)
    except Exception:
        return 0


# =========================
# I/O readers (compute mode)
# =========================
def read_metadata(path: Path, samp_col: str, group_col: str,
                  keep_groups: Optional[Sequence[str]]) -> pd.DataFrame:
    df = pd.read_csv(path, sep='\t', header=0)
    if keep_groups:
        df = df[df[group_col].isin(keep_groups)].copy()
    # Make sure sample ids are strings
    df[samp_col] = df[samp_col].astype(str)
    return df


def load_sample_manifest(path: Path) -> Dict[str, str]:
    """
    Build a lookup from FASTQ file path (or basename) to sample ID.
    Manifest columns: sample_id, fastq_r1, fastq_r2 (no header).
    """
    df = pd.read_csv(path, sep='\t', header=None, names=['sample_id', 'r1', 'r2'])
    mapping: Dict[str, str] = {}
    for _, row in df.iterrows():
        sample_id = str(row['sample_id']).strip()
        if not sample_id:
            continue
        for col in ('r1', 'r2'):
            fastq_path = str(row[col]).strip()
            if not fastq_path or fastq_path.lower() == 'nan':
                continue
            candidates = {
                fastq_path,
                os.path.basename(fastq_path),
            }
            try:
                candidates.add(str(Path(fastq_path).resolve()))
            except Exception:
                pass
            for cand in candidates:
                if cand in mapping and mapping[cand] != sample_id:
                    raise ValueError(
                        f"FASTQ '{cand}' maps to multiple sample IDs ({mapping[cand]} vs {sample_id})"
                    )
                mapping[cand] = sample_id
    if not mapping:
        raise ValueError(f"No FASTQ entries were parsed from manifest: {path}")
    return mapping


def read_fastq_stats(path: Path, samp_col: str,
                     manifest_map: Optional[Dict[str, str]] = None) -> pd.DataFrame:
    """
    Expects columns: file, num_seqs
    Collapses replicates by sample ID via groupby+sum.
    """
    df = pd.read_csv(path, sep='\t', header=0)
    if 'file' not in df or 'num_seqs' not in df:
        raise ValueError(f"{path} must contain columns: file, num_seqs")
    stats_dir = path.parent

    def lookup_sample(file_path: str) -> str:
        if manifest_map:
            candidates = [
                file_path,
                os.path.basename(file_path),
                extract_sample_id_from_path(file_path),
            ]
            rel_path = (stats_dir / file_path)
            candidates.append(str(rel_path))
            candidates.append(os.path.basename(rel_path))
            try:
                candidates.append(str(rel_path.resolve()))
            except Exception:
                pass
            for cand in candidates:
                if cand in manifest_map:
                    return manifest_map[cand]
            print(candidates)
            print(manifest_map)
            raise ValueError(f"File '{file_path}' not found in manifest")
        return extract_sample_id_from_path(file_path)

    df[samp_col] = df['file'].apply(lookup_sample)
    out = df.groupby(samp_col, as_index=False)['num_seqs'].sum()
    return out


def read_asv_matrix(path: Path, samp_col: str) -> pd.DataFrame:
    """
    Input: wide matrix (rows=ASVs, columns=samples), counts.
    Returns long: [ASV_ID, samp_col, count] with count>0
    """
    df = pd.read_csv(path, sep='\t', header=0, index_col=0)
    long_df = df.stack().reset_index()
    long_df.columns = ['ASV_ID', 'sample_raw', 'count']
    long_df = long_df[long_df['count'] > 0].copy()
    long_df[samp_col] = long_df['sample_raw'].astype(str)
    long_df.drop(columns=['sample_raw'], inplace=True)
    return long_df


def group_counts_by_group(long_counts: pd.DataFrame, metadata: pd.DataFrame,
                          samp_col: str, group_col: str) -> pd.DataFrame:
    """
    Merge counts with metadata and sum by the chosen grouping column.
    Replicates with the same sample ID and group are naturally summed.
    """
    merged = long_counts.merge(metadata[[samp_col, group_col]], on=samp_col, how='inner')
    merged[group_col] = merged[group_col].astype(str)
    grp = merged.groupby(group_col, as_index=False)['count'].sum()
    grp.rename(columns={'count': 'num_reads'}, inplace=True)
    return grp


# =========================
# Sankey construction
# =========================
def build_sankey(steps: List[str], counts: List[int],
                 lmp_in: Dict[str, int], lmp_out: Dict[str, int],
                 palette: Dict[str, str], title: str,
                 output_html: Path, labeled: bool,
                 arrangement: str = "snap", loss_groups=None) -> None:
    """Export the same count-conserving, ordered layout to SVG and interactive HTML."""
    import json
    import textwrap
    from html import escape as esc
    if len(steps) != len(counts) or not steps:
        raise ValueError('Sankey stages and totals must be nonempty and aligned')
    if any(v < 0 for v in counts) or any(b > a for a, b in zip(counts, counts[1:])):
        raise ValueError('Sankey stage counts must be nonnegative and nonincreasing')
    lmp_in = {k: int(v) for k, v in lmp_in.items() if v > 0}
    lmp_out = {k: int(v) for k, v in lmp_out.items() if v > 0}
    if sum(lmp_in.values()) != counts[0] or sum(lmp_out.values()) != counts[-1]:
        raise ValueError('Sankey endpoint groups must sum to their stage totals')
    total = max(counts[0], 1)
    group_gaps = 36 * max(len(lmp_in)-1, len(lmp_out)-1, 0)
    retained_height = max(360, group_gaps + 100)
    scale = min(260, retained_height - group_gaps) / total
    divider = 135 + retained_height + 25
    loss_top = divider + 40
    width, height = (len(steps)+2)*230+160, 880
    nodes, links = [], []

    def add_node(label, value, column, y, color, lane, stage=False):
        idx = len(nodes)
        nodes.append(dict(label=label, value=int(value), x=45+column*230,
                          y=float(y), initial_y=float(y), h=max(2., value*scale),
                          color=color, lane=lane, column=column, stage=stage))
        return idx

    def add_link(source, target, value, color='#98A2AD'):
        if value > 0:
            links.append(dict(source=source, target=target, value=int(value), color=color))

    source_nodes, output_nodes = [], []
    for groups, column, dest in [(lmp_in, 0, source_nodes), (lmp_out, len(steps)+1, output_nodes)]:
        cursor = 135.
        for label, value in groups.items():
            i = add_node(label, value, column, cursor, palette.get(label, '#6B7280'), 'retained')
            dest.append(i); cursor += nodes[i]['h'] + 36
    stage_nodes = [add_node(label, value, i+1, 135, '#253D52', 'retained', True)
                   for i, (label, value) in enumerate(zip(steps, counts))]
    for i in source_nodes: add_link(i, stage_nodes[0], nodes[i]['value'], nodes[i]['color'])
    for j, (source, target) in enumerate(zip(stage_nodes, stage_nodes[1:])):
        add_link(source, target, counts[j+1])
        difference = counts[j] - counts[j+1]
        groups = (loss_groups or {}).get(j, [(f'{steps[j+1]} removed', difference)])
        if sum(int(value) for _, value in groups) != difference:
            raise ValueError(f'Loss accounting does not balance at {steps[j+1]}')
        cursor = float(loss_top)
        for label, value in groups:
            if not value: continue
            i = add_node(label, value, j+2, cursor, '#CBD2D9', 'removed')
            add_link(source, i, value, '#CBD2D9'); cursor += nodes[i]['h'] + 50
    for i in output_nodes: add_link(stage_nodes[-1], i, nodes[i]['value'], nodes[i]['color'])

    height = max(height, int(max(n["y"] + n["h"] for n in nodes)) + 90)

    # Sort ports by the actual opposite endpoint, including the final group split.
    # Geometry is identical in both exports and does not depend on browser size.
    for index, node in enumerate(nodes):
        for key, port in [('source', 'sy'), ('target', 'ty')]:
            edge_ids = [i for i, edge in enumerate(links) if edge[key] == index]
            other = 'target' if key == 'source' else 'source'
            edge_ids.sort(key=lambda i: (nodes[links[i][other]]['y'], links[i][other]))
            cursor = node['y']
            for i in edge_ids:
                links[i][port] = cursor - node['y']
                cursor += links[i]['value'] * scale

    def ribbon(edge):
        a,b = nodes[edge['source']],nodes[edge['target']]
        x1,x2=a['x']+18,b['x']; mid=(x1+x2)/2
        y1,y2=a['y']+edge['sy'],b['y']+edge['ty']; h=edge['value']*scale
        return f'M{x1},{y1} C{mid},{y1} {mid},{y2} {x2},{y2} L{x2},{y2+h} C{mid},{y2+h} {mid},{y1+h} {x1},{y1+h} Z'

    parts=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" width="{width}" height="{height}" role="img" aria-label="{esc(title)}">',
           '<style>text{font-family:Arial,sans-serif;fill:#243447;font-size:12px}.node{cursor:grab}.ribbon{fill-opacity:.55}.ribbon:hover{fill-opacity:.85}</style>',
           '<rect width="100%" height="100%" fill="white"/>',
           f'<text x="30" y="30" style="font-size:23px;font-weight:bold">{esc(title)}</text>',
           '<text x="30" y="54">Read pairs / merged-read equivalents; band widths represent counts.</text>',
           f'<path d="M25,{divider} H{width-25}" stroke="#CBD2D9" stroke-dasharray="5 5"/>',
           f'<text x="30" y="{divider+22}" style="font-weight:bold">Removed from downstream biological analysis</text>']
    for i, edge in enumerate(links):
        a,b=nodes[edge['source']],nodes[edge['target']]
        parts.append(f'<path id="edge-{i}" class="ribbon" d="{ribbon(edge)}" fill="{esc(edge["color"])}"><title>{esc(a["label"])} → {esc(b["label"])}: {edge["value"]:,}</title></path>')
    for i, node in enumerate(nodes):
        parts.append(f'<g id="node-{i}" class="node" data-index="{i}" transform="translate({node["x"]},{node["y"]})"><title>{esc(node["label"])}: {node["value"]:,}</title><rect width="18" height="{node["h"]}" fill="{esc(node["color"])}" stroke="#566573"/>')
        if labeled:
            lines=textwrap.wrap(node['label'], 26)
            if node['stage']:
                y=-16-len(lines)*15
                for line in lines:
                    parts.append(f'<text x="0" y="{y}" style="font-weight:bold">{esc(line)}</text>');y+=15
                parts.append(f'<text x="0" y="-7">{node["value"]:,}</text>')
            else:
                for j,line in enumerate(lines):parts.append(f'<text x="26" y="{12+j*14}">{esc(line)}</text>')
                parts.append(f'<text x="26" y="{12+len(lines)*14}">{node["value"]:,}</text>')
        parts.append('</g>')
    parts.append('</svg>'); svg=''.join(parts)
    payload=dict(nodes=nodes,links=links,scale=scale,steps=steps,counts=counts,arrangement=arrangement,loss_top=loss_top,retained_bottom=divider-20,removed_bottom=height-55)
    encoded=json.dumps(payload).replace('<','\\u003c')
    javascript=r'''
const flow=JSON.parse(document.getElementById('flow-data').textContent),svg=document.querySelector('svg');
let active=null,startY=0,originY=0;
function point(e){return new DOMPoint(e.clientX,e.clientY).matrixTransform(svg.getScreenCTM().inverse()).y;}
function redraw(){
 flow.nodes.forEach((n,i)=>document.getElementById('node-'+i).setAttribute('transform',`translate(${n.x},${n.y})`));
 flow.links.forEach((e,i)=>{let a=flow.nodes[e.source],b=flow.nodes[e.target],x=a.x+18,z=b.x,m=(x+z)/2,y=a.y+e.sy,t=b.y+e.ty,h=e.value*flow.scale;
 document.getElementById('edge-'+i).setAttribute('d',`M${x},${y} C${m},${y} ${m},${t} ${z},${t} L${z},${t+h} C${m},${t+h} ${m},${y+h} ${x},${y+h} Z`);});
}
svg.addEventListener('pointerdown',e=>{let g=e.target.closest('.node');if(!g||flow.arrangement==='fixed')return;
 active=+g.dataset.index;startY=point(e);originY=flow.nodes[active].y;svg.setPointerCapture(e.pointerId);e.preventDefault();});
svg.addEventListener('pointermove',e=>{if(active===null)return;let n=flow.nodes[active],lo=n.lane==='removed'?flow.loss_top:135,hi=(n.lane==='removed'?flow.removed_bottom:flow.retained_bottom)-n.h;
 flow.nodes.forEach((o,i)=>{if(i===active||o.column!==n.column||o.lane!==n.lane)return;if(o.initial_y<n.initial_y)lo=Math.max(lo,o.y+o.h+10);else hi=Math.min(hi,o.y-n.h-10);});
 n.y=Math.max(lo,Math.min(hi,originY+point(e)-startY));redraw();});
function release(){if(active!==null&&flow.arrangement==='snap')flow.nodes[active].y=flow.nodes[active].initial_y;active=null;redraw();}
svg.addEventListener('pointerup',release);svg.addEventListener('pointercancel',release);
document.getElementById('reset').onclick=()=>{flow.nodes.forEach(n=>n.y=n.initial_y);redraw();};
'''
    html=f'<!doctype html><html><head><meta charset="utf-8"><title>{esc(title)}</title><style>body{{margin:0;font:14px Arial;color:#243447}}.toolbar{{padding:10px;background:#f3f6f8;position:sticky;top:0}}.canvas{{overflow:auto}}svg{{display:block;touch-action:none}}button{{margin-right:12px}}</style></head><body><div class="toolbar"><button id="reset">Reset layout</button>Drag nodes vertically within their retained/removed area. Scroll horizontally to view all stages.</div><div class="canvas">{svg}</div><script id="flow-data" type="application/json">{encoded}</script><script>{javascript}</script></body></html>'
    output_html.parent.mkdir(parents=True,exist_ok=True)
    output_html.write_text(html)
    output_html.with_suffix('.svg').write_text(svg+'\n')
    output_html.with_suffix('.flow.json').write_text(json.dumps(payload,indent=2)+'\n')
    print(f'✔ Sankey HTML, SVG and flow audit saved: {output_html}')


# =========================
# CLI
# =========================
def get_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Generate Sankey diagrams for data loss/flow.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # --- Compute mode inputs
    io = p.add_argument_group("Compute Mode Inputs")
    io.add_argument("--data-dir", type=Path, help="Project root (used to resolve defaults)")
    io.add_argument("--sub-dir", default="spark_combined_output", help="Subdir under data-dir for outputs/stats")
    io.add_argument("--metadata", type=Path, help="TSV with sample metadata")
    io.add_argument("--sample-manifest", type=Path,
                    help="TSV with columns: sample_id, fastq_r1, fastq_r2")
    io.add_argument("--samp-col", default="lmp_id", help="Sample column name in metadata")
    io.add_argument("--group1-col", default="group1", help="Grouping column in metadata")
    io.add_argument("--color-col", default="Color", help="Color column in metadata")
    io.add_argument("--keep-types", default="",
                    help="Comma-separated list; if empty, keep all types")
    io.add_argument("--all-samples", action="store_true",
                    help="Use all metadata samples as the sample universe instead of only samples present in the microbial ASV table.")

    io.add_argument("--fastq-stats", default="stats/fastq_stats.tsv",
                    help="Path (relative to sub-dir or absolute) to raw fastq stats TSV")

    io.add_argument("--filtered-stats", default="stats/filtered_fastqs.tsv",
                    help="Path to filtered fastq stats TSV")

    io.add_argument("--sample-qc", type=Path, help="CONTROL_DECONTAM sample_qc.tsv for explicit class/depth dropouts")
    io.add_argument("--asv-cleaned", type=Path, help="Biological counts after TECH/BIO ASV removal")
    io.add_argument("--asv-final", type=Path, help="Final filtered ASV_target.tsv before batch correction")
    io.add_argument("--asv-raw", default="ASVs/ASV_counts.tsv", help="Wide ASV counts matrix")

    io.add_argument("--asv-decon", default="ASVs/ASV_target.decon.tsv", help="Wide ASV after decontamination")
    io.add_argument("--asv-micro", default="ASVs/ASV_target.micro.tsv", help="Wide ASV microbial (finished)")

    # --- Appearance / output
    out = p.add_argument_group("Output")
    out.add_argument("--title", default="Data Loss Flow", help="Plot title")
    out.add_argument("--output-prefix", default="metadata/data_loss_sankey",
                     help="Output prefix ('.html' appended automatically)")
    out.add_argument("--make-labeled", action="store_true", help="Create labeled-node HTML")
    out.add_argument("--make-unlabeled", action="store_true", help="Create unlabeled-node HTML")
    out.add_argument(
        "--arrangement",
        default="snap",
        choices=["snap", "perpendicular", "freeform", "fixed"],
        help="Layout interaction: freeform/perpendicular allow vertical movement within ordered lanes; snap resets on release; fixed disables dragging.",
    )
    out.add_argument(
        "--vertical-order",
        default="",
        help="Comma-separated group order from top to bottom for input/output Sankey nodes. Unlisted groups are appended.",
    )

    # --- Misc
    p.add_argument("--verbose", action="store_true", help="Verbose logs")

    return p


def main():
    args = get_parser().parse_args()

    # ---- Compute mode ----
    if not args.data_dir:
        raise SystemExit("--data-dir is required in compute mode")
    data_dir: Path = args.data_dir

    # Resolve default paths if relative
    def resolve(rel_or_abs: str) -> Path:
        p = Path(rel_or_abs)
        if p.is_absolute():
            return p
        return data_dir / args.sub_dir / rel_or_abs

    metadata_path = args.metadata or (data_dir / "ref_db" / "spark_metadata.tsv")
    manifest_path = args.sample_manifest or (data_dir / "ref_db" / "sample_manifest.tsv")
    fastq_stats_path = resolve(args.fastq_stats)
    filtered_stats_path = resolve(args.filtered_stats)
    asv_raw_path = resolve(args.asv_raw)
    asv_decon_path = resolve(args.asv_decon)
    asv_micro_path = resolve(args.asv_micro)
    asv_final_path = resolve(args.asv_final) if args.asv_final else asv_micro_path
    final_long = read_asv_matrix(asv_final_path, args.samp_col)

    keep_types = [t.strip() for t in args.keep_types.split(',')] if args.keep_types.strip() else None

    if args.verbose:
        print(f"[i] Metadata: {metadata_path}")
        print(f"[i] Sample manifest: {manifest_path}")
        print(f"[i] Raw fastq stats: {fastq_stats_path}")
        print(f"[i] Filtered stats: {filtered_stats_path}")
        print(f"[i] ASV raw: {asv_raw_path}")
        print(f"[i] ASV decon: {asv_decon_path}")
        print(f"[i] ASV micro: {asv_micro_path}")

    meta = read_metadata(metadata_path, args.samp_col, args.group1_col, keep_types)
    meta[args.group1_col] = meta[args.group1_col].astype(str)

    # ASV matrices -> long -> merge -> sum
    # Use consistent sample ID parsing across all ASV matrices
    asv_micro_long = read_asv_matrix(
        asv_micro_path,
        args.samp_col,
    )

    if args.all_samples:
        sample_list = meta[args.samp_col].astype(str).unique().tolist()
    else:
        sample_list = final_long[args.samp_col].unique().tolist()

    asv_raw_long = read_asv_matrix(
        asv_raw_path,
        args.samp_col,
    )
    asv_raw_long = asv_raw_long[asv_raw_long[args.samp_col].isin(sample_list)].copy()

    asv_decon_long = read_asv_matrix(
        asv_decon_path,
        args.samp_col,
    )
    asv_decon_long = asv_decon_long[asv_decon_long[args.samp_col].isin(sample_list)].copy()

    # Restrict metadata to the active sample universe
    meta = meta[meta[args.samp_col].isin(sample_list)].copy()

    manifest_map = load_sample_manifest(manifest_path)
    filter_map = {v: v for k, v in manifest_map.items()}

    # Raw reads (pairs): sum num_seqs across files, then /2, with replicates collapsed
    raw_df = read_fastq_stats(
        fastq_stats_path,
        args.samp_col,
        manifest_map,
    )
    raw_df = raw_df[raw_df[args.samp_col].isin(sample_list)].copy()

    # Filtered reads (already single-end counts in your script), replicates collapsed
    filt_df = read_fastq_stats(
        filtered_stats_path,
        args.samp_col,
        filter_map,
    )
    filt_df = filt_df[filt_df[args.samp_col].isin(sample_list)].copy()

    # Build palette from metadata: grouping column -> color, with string keys
    palette = {str(t): str(c) for t, c in zip(meta[args.group1_col], meta[args.color_col])}

    # Sort palette deterministically (numeric if possible, else lexical)
    try:
        palette = dict(sorted(palette.items(), key=lambda x: float(x[0])))
    except (ValueError, TypeError):
        palette = dict(sorted(palette.items()))

    # Sum by group — this implicitly respects keep_types and drops samples without metadata
    raw_by_type = raw_df.merge(meta[[args.samp_col, args.group1_col]],
                               on=args.samp_col, how='inner') \
                        .groupby(args.group1_col, as_index=False)['num_seqs'].sum()
    raw_by_type['num_reads'] = (raw_by_type['num_seqs'] // 2).astype(int)
    
    filt_by_type = filt_df.merge(meta[[args.samp_col, args.group1_col]],
                                 on=args.samp_col, how='inner') \
                          .groupby(args.group1_col, as_index=False)['num_seqs'].sum()
    filt_by_type['num_reads'] = filt_by_type['num_seqs'].astype(int)
    
    # Override step totals so node labels use exactly the subset represented in the ribbons
    raw_reads_total = int(raw_by_type['num_reads'].sum())
    filt_reads_total = int(filt_by_type['num_reads'].sum())

    asv_raw_by_type = group_counts_by_group(asv_raw_long, meta, args.samp_col, args.group1_col)
    asv_decon_by_type = group_counts_by_group(asv_decon_long, meta, args.samp_col, args.group1_col)
    asv_micro_by_type = group_counts_by_group(asv_micro_long, meta, args.samp_col, args.group1_col)

    # Totals for remaining steps
    asv_raw_reads = int(asv_raw_by_type['num_reads'].sum())
    asv_decon_reads = int(asv_decon_by_type['num_reads'].sum())
    asv_micro_reads = int(asv_micro_by_type['num_reads'].sum())

    final_by_type = group_counts_by_group(final_long, meta, args.samp_col, args.group1_col)
    final_reads = int(final_by_type['num_reads'].sum())
    steps = ['Input pairs', 'Read QC', 'ASV inference']
    counts = [raw_reads_total, filt_reads_total, asv_raw_reads]
    losses = {}
    if args.sample_qc:
        if not args.asv_cleaned:
            raise ValueError('--sample-qc requires --asv-cleaned')
        qc = pd.read_csv(resolve(args.sample_qc), sep='\t', dtype={args.samp_col:str})
        qc = qc[qc[args.samp_col].isin(sample_list)]
        raw_per_sample = asv_raw_long.groupby(args.samp_col)['count'].sum()
        qc = qc.set_index(args.samp_col)
        if not set(raw_per_sample.index) <= set(qc.index):
            raise ValueError('Sample QC audit does not cover all ASV samples')
        qc['reads'] = raw_per_sample.reindex(qc.index).fillna(0).astype(int)
        eligible = qc['biological_pass'].astype(str).str.lower().eq('true')
        cohort_reads = int(qc.loc[eligible,'reads'].sum())
        classes = [('technical','TECH controls (used for decontam)'),
                   ('bio_control','BIO controls (used for decontam)'),
                   ('positive','Positive controls (QC only)')]
        cohort_losses = [(label,int(qc.loc[qc.decontam_class.eq(cls),'reads'].sum())) for cls,label in classes]
        cohort_losses.append(('Biological samples below depth cutoff',int(qc.loc[qc.decontam_class.eq('biological') & ~eligible,'reads'].sum())))
        accounted = sum(value for _,value in cohort_losses)
        remaining = asv_raw_reads - cohort_reads - accounted
        if remaining < 0: raise ValueError('Overlapping sample classes in cohort audit')
        if remaining: cohort_losses.append(('Other samples excluded',remaining))
        losses[2] = cohort_losses
        steps.append('Eligible biological cohort'); counts.append(cohort_reads)
        cleaned = read_asv_matrix(resolve(args.asv_cleaned), args.samp_col)
        cleaned_reads = int(cleaned.loc[cleaned[args.samp_col].isin(sample_list),'count'].sum())
        losses[3] = [('TECH/BIO contaminant ASVs removed',cohort_reads-cleaned_reads)]
        steps.append('After TECH/BIO ASV filtering'); counts.append(cleaned_reads)
    steps.extend(['After non-target screening','Finished biological data'])
    counts.extend([asv_micro_reads,final_reads])
    losses[len(counts)-3] = [('Host / mitochondrial reads removed', counts[-3]-counts[-2])]
    losses[len(counts)-2] = [('Abundance / prevalence / taxonomy removed', counts[-2]-counts[-1])]

    vertical_order = [t.strip() for t in args.vertical_order.split(',') if t.strip()]
    if keep_types:
        types = keep_types
    else:
        unique_types = [str(t) for t in meta[args.group1_col].dropna().unique()]
        # Sort numerically if possible, otherwise alphabetically
        try:
            types = sorted(unique_types, key=lambda x: float(x))
        except (ValueError, TypeError):
            types = sorted(unique_types)
    if vertical_order:
        present = {str(t) for t in types}
        ordered = [t for t in vertical_order if t in present]
        ordered.extend([str(t) for t in types if str(t) not in set(ordered)])
        types = ordered

    # Input and output dicts for sankey ends (string keys to match palette)
    lmp_in = {
        str(t): int(raw_by_type.loc[raw_by_type[args.group1_col] == t, 'num_reads'].sum())
        for t in types
    }
    lmp_out = {
        str(t): int(final_by_type.loc[final_by_type[args.group1_col] == t, 'num_reads'].sum())
        for t in types
    }

    # Preserve explicit top-to-bottom order when requested; otherwise sort deterministically.
    if not vertical_order:
        try:
            lmp_in = dict(sorted(lmp_in.items(), key=lambda x: float(x[0])))
            lmp_out = dict(sorted(lmp_out.items(), key=lambda x: float(x[0])))
        except (ValueError, TypeError):
            lmp_in = dict(sorted(lmp_in.items()))
            lmp_out = dict(sorted(lmp_out.items()))

    if args.verbose:
        print("[i] Steps:")
        for s, c in zip(steps, counts):
            print(f"  - {s}: {c}")
        print("[i] Inputs by type:", lmp_in)
        print("[i] Outputs by type:", lmp_out)

    # Outputs
    out_pref = (args.data_dir / args.sub_dir / args.output_prefix) if args.data_dir else Path(args.output_prefix)
    # default: generate both if none chosen
    if not args.make_labeled and not args.make_unlabeled:
        args.make_labeled = True
        args.make_unlabeled = True

    if args.make_labeled:
        build_sankey(
            steps, counts, lmp_in, lmp_out, palette,
            args.title, out_pref.with_suffix(".label.html"), True,
            arrangement=args.arrangement, loss_groups=losses
        )
    if args.make_unlabeled:
        build_sankey(
            steps, counts, lmp_in, lmp_out, palette,
            args.title, out_pref.with_suffix(".html"), False,
            arrangement=args.arrangement, loss_groups=losses
        )


if __name__ == "__main__":
    main()
