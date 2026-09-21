#!/usr/bin/env python3
"""Export ASV association networks annotated with sequence-derived MAG links."""

from __future__ import annotations

import argparse
import json
import math
import re
from collections import Counter, defaultdict
from itertools import permutations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.offsetbox import AnnotationBbox, DrawingArea, HPacker, TextArea, VPacker
import sys as _sys
_sys.path.insert(0, str(Path(__file__).resolve().parent))
_sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shared_plot_style import install_publication_style
install_publication_style()
from mag_grouping_diagnostics import run_mag_grouping_diagnostics
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns


RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
RANK_PREFIXES = {
    "d": "domain",
    "k": "domain",
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
}
TITAN_DIRECTION_COLORS = {"z-": "#BDBDBD", "z+": "#000000"}
TITAN_Z_LEGEND_VALUES = (5.0, 10.0, 20.0, 30.0)


def titan_marker_area(value: object) -> float:
    """Map absolute TITAN z scores to the shared bounded marker-area scale."""
    score = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if not np.isfinite(score):
        return 20.0
    return float(45.0 + 20.0 * min(abs(float(score)), TITAN_Z_LEGEND_VALUES[-1]))


def species_marker_map(taxa: list[str]) -> dict[str, str]:
    """Assign stable species shapes; the first two are circle and square."""
    markers = ("o", "s", "D", "^", "v", "P", "X")
    ordered = sorted(set(str(taxon) for taxon in taxa))
    return {taxon: markers[index % len(markers)] for index, taxon in enumerate(ordered)}


def node_outline_color(facecolor: str) -> str:
    """Give black nodes a thin light rim while retaining dark rims elsewhere."""
    rgb = np.asarray(matplotlib.colors.to_rgb(facecolor), dtype=float)
    return "#E6E6E6" if float(rgb.max()) <= 0.15 else "#000000"


def info(msg: str) -> None:
    print(f"[i] {msg}")


def die(msg: str) -> None:
    raise SystemExit(msg)


def read_table(path: Path | None) -> pd.DataFrame:
    if path is None or not path.is_file() or path.stat().st_size == 0:
        return pd.DataFrame()
    sep = "\t" if path.suffix.lower() in {".tsv", ".tab", ".txt"} else ","
    return pd.read_csv(path, sep=sep, low_memory=False)


def write_table(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def clean_taxon(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none", "unclassified", "uncultured", "unknown"}:
        return ""
    text = re.sub(r"^[dkpcofgs]__?", "", text, flags=re.IGNORECASE)
    text = re.sub(r"^[a-z]__", "", text, flags=re.IGNORECASE)
    text = text.replace("_", " ")
    text = re.sub(r"\s+", " ", text).strip()
    if text.lower() in {"", "na", "n/a", "uncultured", "unclassified", "unknown"}:
        return ""
    # Preserve the source taxonomy capitalization for display. Comparisons
    # canonicalize case explicitly where needed.
    return text


def module_sort_key(value: object) -> tuple[int, str]:
    text = str(value).strip()
    match = re.fullmatch(r"[Mm](\d+)", text)
    return (int(match.group(1)), text) if match else (10**9, text)


def strip_sequence_suffix(value: str) -> str:
    text = str(value).strip()
    changed = True
    while changed:
        changed = False
        for suffix in (".gz", ".gff3", ".gff", ".fasta", ".fna", ".fa", ".fas", ".ffn", ".faa"):
            if text.endswith(suffix):
                text = text[: -len(suffix)]
                changed = True
                break
    text = re.sub(r"\.full$", "", text, flags=re.IGNORECASE)
    return text


def standardize_asv_id(value: object) -> str:
    """Remove only the trailing VSEARCH abundance annotation from an ASV ID."""
    if value is None or pd.isna(value):
        return ""
    return re.sub(r";size=[0-9]+.*$", "", str(value).strip(), flags=re.IGNORECASE)


def standardize_mag_id(value: object, mode: str) -> str:
    if value is None or pd.isna(value):
        return ""
    text = strip_sequence_suffix(str(value).strip())
    if "::" in text:
        text = text.split("::", 1)[1]
    if mode == "suffix_after_double_underscore" and "__" in text:
        return text.rsplit("__", 1)[-1]
    return text


def pairing_mag_match_id(row: pd.Series, mode: str) -> str:
    for col in ("mag_native_genome_id_raw", "Genome_Id", "native_genome_id", "genome_id"):
        if col in row.index:
            value = standardize_mag_id(row.get(col), mode)
            if value:
                return value
    return ""


def parse_taxonomy_string(value: object) -> dict[str, str]:
    if value is None or pd.isna(value):
        return {}
    parts = [part.strip() for part in re.split(r";|\|", str(value)) if part.strip()]
    out: dict[str, str] = {}
    positional = 0
    for part in parts:
        match = re.match(r"^([dkpcofgs])__?(.*)$", part, flags=re.IGNORECASE)
        if match:
            rank = RANK_PREFIXES.get(match.group(1).lower())
            value_clean = clean_taxon(match.group(2))
        else:
            rank = RANKS[positional] if positional < len(RANKS) else ""
            value_clean = clean_taxon(part)
        if rank and value_clean:
            out[rank] = value_clean
        positional += 1
    return out


def load_asv_taxonomy(
    path: Path | None,
    min_confidence: float | None = None,
) -> pd.DataFrame:
    df = read_table(path)
    if df.empty:
        return pd.DataFrame(columns=[
            "ASV_ID",
            *[f"asv_{rank}" for rank in RANKS],
            "asv_taxonomy",
            "asv_taxonomy_confidence",
            "passes_asv_taxonomy_confidence",
        ])
    id_col = next((c for c in df.columns if str(c).strip().lower() in {"feature id", "feature_id", "asv_id", "asv", "id"}), df.columns[0])
    tax_col = next((c for c in df.columns if str(c).strip().lower() in {"taxon", "taxonomy", "consensus"}), None)
    confidence_col = next(
        (
            c for c in df.columns
            if str(c).strip().lower() in {
                "confidence", "consensus", "taxonomy_confidence",
                "taxonomic_confidence", "confidence_score",
            }
            and c != tax_col
        ),
        None,
    )
    if min_confidence is not None and confidence_col is None:
        die(
            "ASV taxonomy confidence filtering was requested, but the taxonomy "
            "table has no recognized confidence column."
        )
    rows = []
    for _, row in df.iterrows():
        asv_id = standardize_asv_id(row.get(id_col, ""))
        tax = parse_taxonomy_string(row.get(tax_col, "")) if tax_col else {}
        confidence = (
            pd.to_numeric(pd.Series([row.get(confidence_col)]), errors="coerce").iloc[0]
            if confidence_col else np.nan
        )
        passes_confidence = (
            True
            if min_confidence is None
            else bool(pd.notna(confidence) and float(confidence) >= min_confidence)
        )
        record = {
            "ASV_ID": asv_id,
            "asv_taxonomy": row.get(tax_col, "") if tax_col else "",
            "asv_taxonomy_confidence": confidence,
            "passes_asv_taxonomy_confidence": passes_confidence,
        }
        for rank in RANKS:
            record[f"asv_{rank}"] = tax.get(rank, "")
        rows.append(record)
    return pd.DataFrame(rows).drop_duplicates("ASV_ID")


def mag_taxonomy_from_row(row: pd.Series) -> dict[str, str]:
    out: dict[str, str] = {}
    for rank in RANKS:
        candidates = [
            f"mag_{rank}",
            rank,
            rank.capitalize(),
            f"gtdb_{rank}",
            f"mag_gtdb_{rank}",
        ]
        for col in candidates:
            if col in row.index:
                value = clean_taxon(row.get(col))
                if value:
                    out[rank] = value
                    break
    joined_cols = [c for c in row.index if str(c).lower() in {"taxonomy", "mag_taxonomy", "gtdb_taxonomy", "classification"}]
    for col in joined_cols:
        parsed = parse_taxonomy_string(row.get(col))
        for rank, value in parsed.items():
            out.setdefault(rank, value)
    return out


def normalize_taxonomy_source(value: str) -> str:
    text = (value or "auto").strip().lower()
    if text in {"", "auto", "unknown", "none"}:
        return "unknown"
    if text in {"ncbi", "silva", "rdp", "pr2"}:
        return "ncbi_like"
    if text in {"gtdb", "gtdbtk", "gtdb-tk"}:
        return "gtdb"
    return text


def validate_pair(row: pd.Series, asv_source: str = "unknown", mag_source: str = "unknown") -> dict[str, object]:
    comparable: list[str] = []
    conflicts: list[str] = []
    agreements: list[str] = []
    lowest = ""
    crossdb = asv_source != "unknown" and mag_source != "unknown" and asv_source != mag_source
    for rank in RANKS:
        asv_value = clean_taxon(row.get(f"asv_{rank}", ""))
        mag_value = clean_taxon(row.get(f"mag_tax_{rank}", ""))
        if not asv_value or not mag_value:
            continue
        comparable.append(rank)
        lowest = rank
        if asv_value.casefold() == mag_value.casefold():
            agreements.append(rank)
        else:
            conflicts.append(rank)
    if conflicts and crossdb:
        status = "taxonomy_unvalidated_crossdb"
        agreement = pd.NA
    elif conflicts:
        status = "taxonomy_rejected"
        agreement = False
    elif comparable:
        status = "taxonomy_accepted"
        agreement = True
    else:
        status = "taxonomy_unvalidated"
        agreement = pd.NA
    return {
        "taxonomy_validation_status": status,
        "taxonomy_agreement": agreement,
        "taxonomy_cross_database": crossdb,
        "asv_taxonomy_source": asv_source,
        "mag_taxonomy_source": mag_source,
        "lowest_comparable_rank": lowest,
        "comparable_ranks": ",".join(comparable),
        "conflicting_ranks": ",".join(conflicts),
        "agreeing_ranks": ",".join(agreements),
    }


def normalize_pairing(
    pairing: pd.DataFrame,
    asv_tax: pd.DataFrame,
    asv_taxonomy_source: str,
    mag_taxonomy_source: str,
    mag_id_mode: str,
) -> pd.DataFrame:
    if pairing.empty:
        return pd.DataFrame()
    out = pairing.copy()
    if "ASV_ID" not in out.columns:
        id_col = next((c for c in out.columns if str(c).lower() in {"asv", "asv_id", "taxon"}), None)
        if id_col:
            out = out.rename(columns={id_col: "ASV_ID"})
    if "genome_id" not in out.columns:
        genome_col = next((c for c in out.columns if str(c).lower() in {"mag", "mag_id", "genome", "genome_id"}), None)
        if genome_col:
            out = out.rename(columns={genome_col: "genome_id"})
    if "ASV_ID" not in out.columns or "genome_id" not in out.columns:
        die("ASV-MAG pairing table must contain ASV_ID and genome_id columns.")
    out["ASV_ID_raw"] = out["ASV_ID"].astype(str).str.strip()
    out["ASV_ID"] = out["ASV_ID"].map(standardize_asv_id)
    out["mag_match_id"] = out.apply(lambda row: pairing_mag_match_id(row, mag_id_mode), axis=1)

    for rank in RANKS:
        values = []
        for _, row in out.iterrows():
            values.append(mag_taxonomy_from_row(row).get(rank, ""))
        out[f"mag_tax_{rank}"] = values

    out = out.merge(asv_tax, on="ASV_ID", how="left")
    asv_source = normalize_taxonomy_source(asv_taxonomy_source)
    mag_source = normalize_taxonomy_source(mag_taxonomy_source)
    validations = out.apply(lambda row: validate_pair(row, asv_source, mag_source), axis=1, result_type="expand")
    out = pd.concat([out, validations], axis=1)
    passes_taxonomy_confidence = (
        out["passes_asv_taxonomy_confidence"].fillna(False).astype(bool)
        if "passes_asv_taxonomy_confidence" in out.columns
        else pd.Series(True, index=out.index)
    )
    out.loc[
        ~passes_taxonomy_confidence,
        "taxonomy_validation_status",
    ] = "asv_taxonomy_below_confidence"

    status = out.get("pairing_status", pd.Series("", index=out.index)).astype(str)
    out["mapping_class"] = np.select(
        [
            status.eq("unpaired") | out["genome_id"].isna() | out["genome_id"].astype(str).str.strip().eq(""),
            ~passes_taxonomy_confidence,
            status.str.contains("ambiguous", case=False, na=False),
            out["taxonomy_validation_status"].eq("taxonomy_rejected"),
            status.str.contains("unique", case=False, na=False),
        ],
        [
            "no_match",
            "asv_taxonomy_below_confidence",
            "ambiguous_match",
            "taxonomy_rejected",
            "unique_match",
        ],
        default=status.replace({"": "unknown"}),
    )
    out["accepted_paper_pair"] = out["mapping_class"].eq("unique_match") & out["taxonomy_validation_status"].isin(
        ["taxonomy_accepted", "taxonomy_unvalidated", "taxonomy_unvalidated_crossdb"]
    )
    return out.sort_values(["ASV_ID", "genome_id"], na_position="last")


def graph_edges_to_table(graph: nx.Graph) -> pd.DataFrame:
    rows = []
    for source, target, attrs in sorted(graph.edges(data=True), key=lambda item: (str(item[0]), str(item[1]))):
        weight = attrs.get("weight", attrs.get("Weight", attrs.get("partial_correlation", "")))
        try:
            weight_float = float(weight)
        except (TypeError, ValueError):
            weight_float = np.nan
        rows.append(
            {
                "source": source,
                "target": target,
                "edge_type": attrs.get("edge_type", "asv_association"),
                "sign": "positive" if not np.isfinite(weight_float) or weight_float >= 0 else "negative",
                "weight": weight_float if np.isfinite(weight_float) else weight,
                "association_interpretation": "SPIEC-EASI statistical association; not a direct biological interaction",
            }
        )
    return pd.DataFrame(rows)


def read_node_features(path: Path | None) -> pd.DataFrame:
    df = read_table(path)
    if df.empty:
        return pd.DataFrame(columns=["ASV_ID"])
    if "Taxon" in df.columns:
        df = df.rename(columns={"Taxon": "ASV_ID"})
    elif "ASV_ID" not in df.columns:
        df = df.rename(columns={df.columns[0]: "ASV_ID"})
    return df.drop_duplicates("ASV_ID")


def annotate_ecological_modules(
    node_features: pd.DataFrame,
    modules: pd.DataFrame,
    anchor_top_n: int,
    module_selection: pd.DataFrame | None = None,
    display_all_modules: bool = False,
) -> pd.DataFrame:
    """Attach ecological-module membership and within-module anchor ranks."""
    nodes = node_features.copy()
    if modules.empty:
        return nodes
    asv_col = "Taxon" if "Taxon" in modules else "ASV_ID" if "ASV_ID" in modules else None
    if asv_col is None or "module_label" not in modules:
        return nodes
    annotation_cols = [
        col for col in
        [asv_col, "module_label", "module_id", "node_stability", "is_best"]
        if col in modules
    ]
    annotation = modules[annotation_cols].copy()
    annotation["ASV_ID"] = annotation[asv_col].astype(str).str.strip()
    annotation = annotation.drop(columns=[asv_col], errors="ignore").drop_duplicates("ASV_ID")
    nodes["ASV_ID"] = nodes["ASV_ID"].astype(str)
    nodes = nodes.merge(annotation, on="ASV_ID", how="left")
    rank_specs = {
        "Degree": "anchor_rank_degree",
        "EigenCentral": "anchor_rank_eigencentral",
        "Betweenness": "anchor_rank_betweenness",
    }
    anchor_columns = []
    for metric, rank_column in rank_specs.items():
        values = pd.to_numeric(nodes.get(metric), errors="coerce")
        nodes[metric] = values
        nodes[rank_column] = (
            nodes.assign(_metric=values)
            .groupby("module_label")["_metric"]
            .rank(method="dense", ascending=False, na_option="bottom")
            .astype("Int64")
        )
        anchor_flag = rank_column.replace("rank", "is")
        nodes[anchor_flag] = nodes[rank_column].le(anchor_top_n).fillna(False)
        anchor_columns.append(anchor_flag)
    nodes["is_ecological_anchor"] = nodes[anchor_columns].any(axis=1)
    selected_modules: set[str] = set()
    if not display_all_modules and module_selection is not None and not module_selection.empty:
        selection = module_selection.copy()
        if "is_best" in selection and "module_label" in selection:
            is_best = selection["is_best"].astype(str).str.lower().isin(
                {"true", "1", "yes"}
            )
            selected_modules = set(
                selection.loc[is_best, "module_label"].astype(str)
            )
    if not selected_modules:
        selected_modules = set(nodes["module_label"].dropna().astype(str))
    nodes["is_display_ecological_module"] = (
        nodes["module_label"].astype(str).isin(selected_modules)
    )
    nodes["display_module_label"] = nodes["module_label"].where(
        nodes["is_display_ecological_module"], "Other modules"
    )
    nodes["is_ecological_anchor"] = (
        nodes["is_ecological_anchor"] & nodes["is_display_ecological_module"]
    )
    return nodes


def summarize_functional_annotations(paths: list[Path], module_min_fraction: float = 0.5, mag_id_mode: str = "exact") -> pd.DataFrame:
    summaries: dict[str, Counter] = defaultdict(Counter)
    module_rows: list[dict[str, object]] = []
    if not paths:
        return pd.DataFrame(columns=["genome_id"])
    for path in paths:
        df = read_table(path)
        if df.empty:
            continue
        lower_cols = {str(c).lower(): c for c in df.columns}
        genome_col = next((c for c in df.columns if str(c).lower() in {"genome_id", "genome", "mag", "mag_id", "bin_id", "bin id"}), None)
        if genome_col is None:
            continue
        if {"module", "module_name", "logic_score", "complete"}.issubset(lower_cols):
            module_id_col = lower_cols["module"]
            module_name_col = lower_cols["module_name"]
            logic_score_col = lower_cols["logic_score"]
            complete_col = lower_cols["complete"]
            observed_col = lower_cols.get("observed_module_kos")
            total_col = lower_cols.get("total_module_kos")
            for _, row in df.iterrows():
                genome = str(row.get(genome_col, "")).strip()
                mag_match_id = standardize_mag_id(genome, mag_id_mode)
                module_id = str(row.get(module_id_col, "")).strip()
                if not genome or not mag_match_id or not module_id:
                    continue
                logic_score = pd.to_numeric(
                    pd.Series([row.get(logic_score_col)]), errors="coerce"
                ).iloc[0]
                complete_value = str(row.get(complete_col, "")).strip().lower()
                module_complete = complete_value in {"true", "1", "yes", "y"}
                if module_complete and not (
                    pd.notna(logic_score) and np.isclose(float(logic_score), 1.0)
                ):
                    raise ValueError(
                        f"KEGG module {module_id} for {genome} is marked complete "
                        f"but has logic_score={logic_score!r}"
                    )
                module_rows.append(
                    {
                        "genome_id": genome,
                        "mag_match_id": mag_match_id,
                        "module_id": module_id,
                        "module_name": row.get(module_name_col, ""),
                        "observed_ko_count": row.get(observed_col, pd.NA) if observed_col else pd.NA,
                        "total_feature_kos": row.get(total_col, pd.NA) if total_col else pd.NA,
                        "fraction_covered": logic_score,
                        "logic_score": logic_score,
                        "module_complete": module_complete,
                        "present": module_complete,
                        "presence_definition": "KEGG logical definition complete (logic_score = 1.0)",
                    }
                )
                if module_complete:
                    summaries[mag_match_id][
                        f"{module_id}:{row.get(module_name_col, '')}"
                    ] += 1
            continue
        if {"feature_id", "feature_name", "fraction_covered"}.issubset(lower_cols):
            feature_id_col = lower_cols["feature_id"]
            feature_name_col = lower_cols["feature_name"]
            fraction_col = lower_cols["fraction_covered"]
            observed_col = lower_cols.get("observed_ko_count")
            total_col = lower_cols.get("total_feature_kos")
            for _, row in df.iterrows():
                genome = str(row.get(genome_col, "")).strip()
                mag_match_id = standardize_mag_id(genome, mag_id_mode)
                module_id = str(row.get(feature_id_col, "")).strip()
                if not genome or not mag_match_id or not module_id:
                    continue
                fraction = pd.to_numeric(pd.Series([row.get(fraction_col)]), errors="coerce").iloc[0]
                present = bool(pd.notna(fraction) and fraction >= module_min_fraction)
                module_rows.append(
                    {
                        "genome_id": genome,
                        "mag_match_id": mag_match_id,
                        "module_id": module_id,
                        "module_name": row.get(feature_name_col, ""),
                        "observed_ko_count": row.get(observed_col, pd.NA) if observed_col else pd.NA,
                        "total_feature_kos": row.get(total_col, pd.NA) if total_col else pd.NA,
                        "fraction_covered": fraction,
                        "logic_score": pd.NA,
                        "module_complete": bool(np.isclose(float(fraction), 1.0)) if pd.notna(fraction) else False,
                        "present": present,
                        "presence_definition": f"fraction_covered >= {module_min_fraction}",
                    }
                )
                if present:
                    summaries[genome][f"{module_id}:{row.get(feature_name_col, '')}"] += 1
            continue
        ann_cols = [c for c in df.columns if str(c).lower() in {"ko", "kegg", "ec", "pfam", "cog", "pathway", "module"}]
        if not ann_cols:
            ann_cols = [c for c in df.columns if c != genome_col][:3]
        for _, row in df.iterrows():
            genome = str(row.get(genome_col, "")).strip()
            mag_match_id = standardize_mag_id(genome, mag_id_mode)
            if not genome or not mag_match_id:
                continue
            for col in ann_cols:
                value = row.get(col)
                if pd.isna(value) or not str(value).strip():
                    continue
                summaries[mag_match_id][f"{col}:{value}"] += 1
    rows = []
    for mag_match_id, counts in summaries.items():
        top = [item for item, _count in counts.most_common(25)]
        module_records = [row for row in module_rows if row["mag_match_id"] == mag_match_id]
        present_modules = [row for row in module_records if row.get("present")]
        rows.append(
            {
                "mag_match_id": mag_match_id,
                "genome_id": next((row["genome_id"] for row in module_records), mag_match_id),
                "functional_feature_count": sum(counts.values()),
                "functional_summary": "|".join(top),
                "module_count": len(module_records),
                "module_present_count": len(present_modules),
                "module_present_ids": "|".join(str(row["module_id"]) for row in present_modules),
            }
        )
    out = pd.DataFrame(rows)
    if module_rows:
        module_df = pd.DataFrame(module_rows)
        out.attrs["module_table"] = module_df
    return out


def select_well_covered_cycle_modules(
    module_table: pd.DataFrame,
    minimum_fraction: float,
) -> tuple[list[str], pd.DataFrame]:
    """Select environmentally relevant N/S KEGG modules supported by ≥1 linked MAG."""
    if module_table.empty:
        return [], pd.DataFrame()
    work = module_table.copy()
    work["module_name"] = work["module_name"].fillna("").astype(str)
    work["fraction_covered"] = pd.to_numeric(
        work["fraction_covered"], errors="coerce"
    ).fillna(0.0)
    nitrogen_pattern = re.compile(
        r"nitrogen fixation|nitrification|denitrification|"
        r"(?:assimilatory|dissimilatory) nitrate reduction|"
        r"nitrate assimilation|anammox",
        flags=re.IGNORECASE,
    )
    sulfur_pattern = re.compile(
        r"(?:assimilatory|dissimilatory) sulfate reduction|"
        r"sulfur oxidation|sulfide oxidation|sulfur reduction|"
        r"dimethylsulfoniopropionate .*degradation|DMSP degradation|"
        r"sulfoquinovose degradation",
        flags=re.IGNORECASE,
    )
    work["cycle"] = np.select(
        [
            work["module_name"].str.contains(nitrogen_pattern, regex=True),
            work["module_name"].str.contains(sulfur_pattern, regex=True),
        ],
        ["N", "S"],
        default="",
    )
    candidates = work.loc[work["cycle"].ne("")].copy()
    if candidates.empty:
        return [], pd.DataFrame()
    per_mag = (
        candidates.groupby(
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


def load_mag_abundance(
    path: Path,
    fmt: str,
    genome_col: str,
    sample_col: str,
    value_col: str,
    mag_id_mode: str,
) -> pd.DataFrame:
    mag = read_table(path)
    if mag.empty:
        return pd.DataFrame()
    fmt = fmt.lower()
    if fmt == "auto":
        lower = {c.lower(): c for c in mag.columns}
        fmt = "long" if genome_col.lower() in lower and sample_col.lower() in lower and value_col.lower() in lower else "wide"
    if fmt == "long":
        lower = {c.lower(): c for c in mag.columns}
        gcol = lower.get(genome_col.lower(), genome_col)
        scol = lower.get(sample_col.lower(), sample_col)
        vcol = lower.get(value_col.lower(), value_col)
        missing = [c for c in [gcol, scol, vcol] if c not in mag.columns]
        if missing:
            die(f"MAG abundance long table missing columns: {', '.join(missing)}")
        mag = mag[[gcol, scol, vcol]].copy()
        converted = pd.to_numeric(mag[vcol], errors="coerce")
        invalid = converted.isna() & mag[vcol].notna()
        if invalid.any():
            die(
                "MAG abundance long table contains nonnumeric values in "
                f"{vcol}: rows {', '.join(map(str, mag.index[invalid][:5]))}"
            )
        mag[vcol] = converted.fillna(0.0)
        mag["_mag_match_id"] = mag[gcol].map(lambda value: standardize_mag_id(value, mag_id_mode))
        mag = mag.loc[mag["_mag_match_id"].astype(bool)].copy()
        return mag.pivot_table(index="_mag_match_id", columns=scol, values=vcol, aggfunc="sum").fillna(0.0)
    mag = mag.set_index(mag.columns[0])
    mag.index = mag.index.map(lambda value: standardize_mag_id(value, mag_id_mode))
    mag = mag.loc[mag.index.astype(bool)]
    converted = mag.apply(pd.to_numeric, errors="coerce")
    invalid = converted.isna() & mag.notna()
    if invalid.any().any():
        row, column = invalid.stack().loc[lambda values: values].index[0]
        die(
            "MAG abundance wide table contains a nonnumeric value at "
            f"genome={row}, sample={column}"
        )
    return converted.fillna(0.0)


def relative_by_sample(table: pd.DataFrame) -> pd.DataFrame:
    numeric = table.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    denom = numeric.sum(axis=0).replace(0, np.nan)
    return numeric.div(denom, axis=1).fillna(0.0)


def median_ratio_normalize(table: pd.DataFrame) -> tuple[pd.DataFrame, pd.Series]:
    """Normalize nonexclusive counts without treating their column sum as a composition."""
    numeric = table.apply(pd.to_numeric, errors="coerce").fillna(0.0).clip(lower=0.0)
    positive = numeric.where(numeric > 0)
    geometric_means = np.exp(np.log(positive).mean(axis=1, skipna=True))
    geometric_means = geometric_means.replace([np.inf, -np.inf, 0], np.nan)
    ratios = numeric.div(geometric_means, axis=0).where(numeric > 0)
    size_factors = ratios.median(axis=0, skipna=True)
    valid = size_factors.replace([np.inf, -np.inf], np.nan).dropna()
    valid = valid.loc[valid > 0]
    if valid.empty:
        die("MAG median-ratio normalization could not estimate any positive sample size factors.")
    center = float(np.exp(np.log(valid).mean()))
    size_factors = size_factors / center
    missing = ~np.isfinite(size_factors) | size_factors.le(0)
    if missing.any():
        fallback = numeric.sum(axis=0)
        fallback_positive = fallback.loc[fallback > 0]
        fallback_center = float(np.exp(np.log(fallback_positive).mean()))
        size_factors.loc[missing] = fallback.loc[missing] / fallback_center
    if (~np.isfinite(size_factors) | size_factors.le(0)).any():
        bad = size_factors.index[(~np.isfinite(size_factors) | size_factors.le(0))].tolist()
        die(f"MAG normalization failed for samples: {', '.join(map(str, bad[:5]))}")
    return numeric.div(size_factors, axis=1), size_factors


def parse_seqkit_paired_libraries(
    seqkit: pd.DataFrame,
    *,
    file_col: str = "file",
    count_col: str = "num_seqs",
) -> tuple[pd.Series, pd.DataFrame]:
    """Return input read-pair fragments after strict SeqKit mate validation."""
    required = {file_col, count_col}
    if not required.issubset(seqkit):
        raise ValueError(
            "SeqKit table lacks required columns: "
            + ", ".join(sorted(required - set(seqkit)))
        )
    rows = []
    for source_row, row in seqkit[[file_col, count_col]].iterrows():
        filename = Path(str(row[file_col]).strip()).name
        stem = re.sub(
            r"(?i)\.(?:fastq|fq)(?:\.(?:gz|bz2|xz|zst))?$", "", filename
        )
        match = re.fullmatch(
            r"(?i)(?P<sample>.+?)(?:_pe[._-]?|[._-]R?)(?P<mate>[12])"
            r"(?:_001)?",
            stem,
        )
        if match is None:
            raise ValueError(
                "Could not parse paired-end sample and mate from SeqKit file "
                f"'{row[file_col]}'. Expected suffixes such as _pe.1/_pe.2, "
                "_R1/_R2, or _1/_2."
            )
        count = pd.to_numeric(pd.Series([row[count_col]]), errors="coerce").iloc[0]
        if not np.isfinite(count) or count <= 0 or float(count) != int(count):
            raise ValueError(
                f"SeqKit {count_col} must be a positive integer for '{filename}'"
            )
        rows.append({
            "seqkit_source_row": int(source_row),
            "seqkit_file": str(row[file_col]),
            "seqkit_filename": filename,
            "sample": match.group("sample"),
            "mate": int(match.group("mate")),
            "num_seqs": int(count),
        })
    if not rows:
        raise ValueError("SeqKit table contains no paired-end FASTQ rows")
    audit = pd.DataFrame(rows)
    audit["duplicate_source_rows"] = audit.groupby(
        ["sample", "mate"]
    )["seqkit_source_row"].transform("size")
    duplicates = audit["duplicate_source_rows"].gt(1)
    conflicting = (
        audit.loc[duplicates]
        .groupby(["sample", "mate"])
        .agg(
            filenames=("seqkit_filename", "nunique"),
            counts=("num_seqs", "nunique"),
        )
    )
    conflicting = conflicting.loc[
        conflicting["filenames"].gt(1) | conflicting["counts"].gt(1)
    ]
    if not conflicting.empty:
        raise ValueError(
            "SeqKit table contains conflicting duplicate sample/mate rows: "
            + ", ".join(
                f"{sample}/R{mate}" for sample, mate in conflicting.index
            )
        )
    audit = audit.drop_duplicates(["sample", "mate"], keep="first").copy()
    pairs = audit.pivot(index="sample", columns="mate", values="num_seqs")
    missing = pairs.index[pairs.reindex(columns=[1, 2]).isna().any(axis=1)]
    if len(missing):
        raise ValueError(
            "SeqKit paired libraries lack R1 or R2 for: "
            + ", ".join(map(str, missing[:10]))
        )
    mismatch = pairs.index[pairs[1].ne(pairs[2])]
    if len(mismatch):
        details = ", ".join(
            f"{sample} (R1={int(pairs.loc[sample, 1])}, "
            f"R2={int(pairs.loc[sample, 2])})"
            for sample in mismatch[:10]
        )
        raise ValueError(f"SeqKit R1/R2 sequence counts disagree: {details}")
    input_fragments = pairs[1].astype(float).rename("input_fragments")
    audit["input_fragments"] = audit["sample"].map(input_fragments)
    audit["pair_validation"] = np.where(
        audit["duplicate_source_rows"].gt(1),
        "matched_duplicate_reference_collapsed",
        "matched",
    )
    return input_fragments, audit.sort_values(
        ["sample", "mate"]
    ).reset_index(drop=True)


def normalize_mag_recruitment(
    table: pd.DataFrame,
    *,
    method: str,
    seqkit_path: Path | None = None,
    seqkit_file_col: str = "file",
    seqkit_count_col: str = "num_seqs",
) -> tuple[pd.DataFrame, pd.Series, str, pd.DataFrame]:
    """Normalize MAG recruitment, preferring input-fragment FPM when available."""
    converted = table.apply(pd.to_numeric, errors="coerce")
    selected = method.lower()
    if selected in {"provided_fpkm", "provided_tpm"}:
        invalid = converted.isna() & table.notna()
        if invalid.any().any():
            row, column = invalid.stack().loc[lambda x: x].index[0]
            raise ValueError(
                f"{selected} input contains a nonnumeric value at "
                f"genome={row}, sample={column}"
            )
        numeric = converted.fillna(0.0)
        values = numeric.to_numpy(dtype=float)
        if not np.isfinite(values).all():
            raise ValueError(f"{selected} input contains nonfinite values")
        if (values < 0).any():
            raise ValueError(f"{selected} input contains negative values")
        # These values were normalized by the upstream quantification workflow.
        # Preserve them exactly; absent genome/sample combinations are zeros.
        return numeric, pd.Series(dtype=float), selected, pd.DataFrame()

    numeric = converted.fillna(0.0).clip(lower=0.0)
    input_fragments = None
    seqkit_audit = pd.DataFrame()
    if seqkit_path is not None:
        input_fragments, seqkit_audit = parse_seqkit_paired_libraries(
            read_table(seqkit_path),
            file_col=seqkit_file_col,
            count_col=seqkit_count_col,
        )
    if selected == "auto":
        if input_fragments is None:
            raise ValueError(
                "auto MAG normalization requires a paired-end SeqKit library "
                "table; select median_ratio explicitly to use the legacy "
                "normalization"
            )
        selected = "input_fragment_fpm"
    if selected == "input_fragment_fpm":
        if input_fragments is None:
            raise ValueError(
                "input_fragment_fpm normalization requires a SeqKit library table"
            )
        denominators = pd.to_numeric(
            input_fragments.reindex(numeric.columns), errors="coerce"
        )
        missing = denominators.index[
            denominators.isna() | ~np.isfinite(denominators) | denominators.le(0)
        ]
        if len(missing):
            raise ValueError(
                "SeqKit input-fragment counts are missing or invalid for MAG "
                "recruitment samples: " + ", ".join(map(str, missing[:10]))
            )
        normalized = numeric.div(denominators, axis=1).mul(1_000_000.0)
        return normalized, denominators, selected, seqkit_audit
    if selected == "median_ratio":
        normalized, factors = median_ratio_normalize(numeric)
        return normalized, factors, "median_ratio_positive_counts", seqkit_audit
    if selected == "raw":
        return numeric, pd.Series(dtype=float), "raw_counts", seqkit_audit
    if selected == "relative":
        return (
            relative_by_sample(numeric),
            numeric.sum(axis=0),
            "selected_catalog_relative_abundance",
            seqkit_audit,
        )
    raise ValueError(
        "MAG normalization must be auto, input_fragment_fpm, provided_fpkm, "
        "provided_tpm, median_ratio, raw, or relative"
    )


def standardize_profiles(table: pd.DataFrame) -> pd.DataFrame:
    """Center and scale each feature across samples for profile-shape comparison."""
    means = table.mean(axis=1)
    scales = table.std(axis=1, ddof=0).replace(0, np.nan)
    return table.sub(means, axis=0).div(scales, axis=0).fillna(0.0)


def harmonize_abundance_samples(
    asv: pd.DataFrame,
    mag: pd.DataFrame,
    sample_metadata: pd.DataFrame | None = None,
    target_sample_col: str | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    """Harmonize assay IDs to biological sample IDs using metadata when available."""
    exact = sorted(set(asv.columns).intersection(mag.columns))
    regex_renamed = {
        str(sample): re.sub(r"_plate[0-9]+$", "", str(sample), flags=re.IGNORECASE)
        for sample in asv.columns
    }
    normalized_overlap = sorted(set(regex_renamed.values()).intersection(mag.columns))
    best_renamed = regex_renamed
    best_overlap = len(normalized_overlap)
    best_mode = "strip_terminal_plate"
    metadata_source_col = ""
    if (
        isinstance(sample_metadata, pd.DataFrame)
        and target_sample_col
        and target_sample_col in sample_metadata.columns
    ):
        target = sample_metadata[target_sample_col].astype(str)
        asv_ids = set(map(str, asv.columns))
        for source_col in sample_metadata.columns:
            source = sample_metadata[source_col].astype(str)
            matched = source.isin(asv_ids) & target.notna()
            if not matched.any():
                continue
            mapping = dict(zip(source.loc[matched], target.loc[matched]))
            candidate = {
                str(sample): mapping.get(str(sample), str(sample))
                for sample in asv.columns
            }
            overlap = len(set(candidate.values()).intersection(mag.columns))
            if overlap > best_overlap:
                best_renamed = candidate
                best_overlap = overlap
                best_mode = "metadata_mapping"
                metadata_source_col = str(source_col)
    mode = "exact"
    collisions = 0
    if best_overlap > len(exact):
        mode = best_mode
        before = len(asv.columns)
        asv = asv.rename(columns=best_renamed).T.groupby(level=0, sort=False).sum().T
        collisions = before - len(asv.columns)
    shared = sorted(set(asv.columns).intersection(mag.columns))
    audit = {
        "sample_harmonization_mode": mode,
        "asv_sample_count": int(len(asv.columns)),
        "mag_sample_count": int(len(mag.columns)),
        "exact_shared_sample_count": int(len(exact)),
        "normalized_shared_sample_count": int(best_overlap),
        "shared_sample_count_used": int(len(shared)),
        "normalization_collisions_collapsed": int(collisions),
        "metadata_mapping_source_col": metadata_source_col,
        "metadata_mapping_target_col": target_sample_col if mode == "metadata_mapping" else "",
    }
    return asv, mag, audit


def abundance_agreement(
    pairing: pd.DataFrame,
    mag_abundance_path: Path | None,
    asv_counts_path: Path | None,
    mag_abundance_format: str,
    mag_genome_col: str,
    mag_sample_col: str,
    mag_value_col: str,
    mag_id_mode: str,
    min_shared_samples: int,
    transform: str,
    mag_normalization: str,
    seqkit_path: Path | None = None,
    seqkit_file_col: str = "file",
    seqkit_count_col: str = "num_seqs",
    sample_metadata: pd.DataFrame | None = None,
    target_sample_col: str | None = None,
    include_ambiguous: bool = False,
) -> tuple[pd.DataFrame, dict[str, object]]:
    if mag_abundance_path is None or not mag_abundance_path.is_file() or asv_counts_path is None or not asv_counts_path.is_file():
        return pd.DataFrame(), {}
    mag = load_mag_abundance(mag_abundance_path, mag_abundance_format, mag_genome_col, mag_sample_col, mag_value_col, mag_id_mode)
    asv = pd.read_csv(asv_counts_path, sep="\t", index_col=0)
    if mag.empty or asv.empty:
        return pd.DataFrame(), {}
    asv.index = asv.index.astype(str)
    asv.columns = asv.columns.astype(str)
    mag.index = mag.index.astype(str)
    mag.columns = mag.columns.astype(str)
    mag_all = mag.copy()
    (
        mag_all_normalized,
        mag_all_size_factors,
        mag_all_normalization_used,
        seqkit_audit,
    ) = normalize_mag_recruitment(
        mag_all,
        method=mag_normalization,
        seqkit_path=seqkit_path,
        seqkit_file_col=seqkit_file_col,
        seqkit_count_col=seqkit_count_col,
    )
    mag_all_transformed = (
        np.log1p(mag_all_normalized)
        if transform == "log1p"
        else mag_all_normalized.copy()
    )
    # ASPIRE count tables are usually features x samples; transpose only if needed.
    pair_samples = set(asv.columns).intersection(mag.columns)
    if len(pair_samples) < 2 and set(asv.index).intersection(mag.columns):
        asv = asv.T
    asv, mag, sample_audit = harmonize_abundance_samples(
        asv,
        mag,
        sample_metadata=sample_metadata,
        target_sample_col=target_sample_col,
    )
    pair_samples = set(asv.columns).intersection(mag.columns)
    samples = sorted(pair_samples)
    asv_all_relative = relative_by_sample(asv)
    asv = asv_all_relative.loc[:, samples]
    mag = mag.loc[:, samples]
    mag, mag_size_factors, normalization_used, _ = normalize_mag_recruitment(
        mag,
        method=mag_normalization,
        seqkit_path=seqkit_path,
        seqkit_file_col=seqkit_file_col,
        seqkit_count_col=seqkit_count_col,
    )
    sample_audit["mag_normalization"] = normalization_used
    sample_audit["mag_normalization_factor_min"] = (
        float(mag_size_factors.min()) if not mag_size_factors.empty else np.nan
    )
    sample_audit["mag_normalization_factor_max"] = (
        float(mag_size_factors.max()) if not mag_size_factors.empty else np.nan
    )
    unit_names = {
        "input_fragment_fpm": "fragments_per_million_input_fragments",
        "provided_fpkm": "fragments_per_kilobase_per_million",
        "provided_tpm": "transcripts_per_million",
    }
    sample_audit["mag_abundance_units"] = unit_names.get(
        normalization_used, normalization_used
    )
    mag_normalized = mag.copy()
    if transform == "log1p":
        asv = np.log1p(asv)
        mag = np.log1p(mag)
    mag_transformed = mag.copy()
    asv = standardize_profiles(asv)
    mag = standardize_profiles(mag)
    sample_audit["profile_standardization"] = "per_feature_zscore_across_samples"
    rows = []
    if include_ambiguous:
        eligible = pairing["accepted_paper_pair"] | (
            pairing["mapping_class"].eq("ambiguous_match")
            & pairing["passes_paper_thresholds"]
            & ~pairing["taxonomy_validation_status"].eq("taxonomy_rejected")
        )
    else:
        eligible = pairing["accepted_paper_pair"]
    accepted = pairing.loc[eligible].dropna(subset=["genome_id"])
    for _, row in accepted.iterrows():
        asv_id = str(row["ASV_ID"])
        genome_id = str(row["genome_id"])
        mag_match_id = str(row.get("mag_match_id", genome_id))
        if asv_id not in asv.index or mag_match_id not in mag.index or len(samples) < min_shared_samples:
            continue
        x = pd.to_numeric(asv.loc[asv_id, samples], errors="coerce")
        y = pd.to_numeric(mag.loc[mag_match_id, samples], errors="coerce")
        valid = x.notna() & y.notna()
        if valid.sum() < min_shared_samples:
            continue
        rows.append(
            {
                "ASV_ID": asv_id,
                "genome_id": genome_id,
                "mag_match_id": mag_match_id,
                "shared_samples": int(valid.sum()),
                "abundance_metric": f"standardized_ASV_relative_abundance_vs_MAG_{sample_audit['mag_normalization']}",
                "transform": f"{transform}+per_feature_zscore",
                "spearman_rho": float(x[valid].corr(y[valid], method="spearman")),
            }
        )
    matrices = {
        "asv": asv.loc[:, samples],
        "asv_all_relative": asv_all_relative,
        "mag": mag.loc[:, samples],
        "mag_normalized": mag_normalized.loc[:, samples],
        "mag_transformed": mag_transformed.loc[:, samples],
        "mag_all_normalized": mag_all_normalized,
        "mag_all_size_factors": mag_all_size_factors,
        "mag_all_normalization": mag_all_normalization_used,
        "mag_all_transformed": mag_all_transformed,
        "mag_all_raw": mag_all,
        "mag_all_samples": list(mag_all_transformed.columns),
        "seqkit_audit": seqkit_audit,
        "samples": samples,
        "sample_audit": sample_audit,
    }
    return pd.DataFrame(rows), matrices


def ambiguity_outputs(
    pairing: pd.DataFrame,
    agreement: pd.DataFrame,
    matrices: dict[str, object],
    module_table: pd.DataFrame,
    outdir: Path,
    prefix: str,
    target_modules: list[str],
    formats: list[str],
    cruise_metadata: pd.DataFrame | None = None,
    profile_only: bool = False,
) -> None:
    """Plot every eligible linked ASV while diagnosing one-to-many mappings separately."""
    def flag_mask(values: pd.Series) -> pd.Series:
        return values.fillna(False).map(
            lambda value: bool(value)
            if isinstance(value, (bool, np.bool_))
            else str(value).strip().lower() in {"1", "true", "yes"}
        )

    if "analysis_eligible_pair" in pairing.columns:
        eligible_mask = flag_mask(pairing["analysis_eligible_pair"])
    else:
        eligible_mask = flag_mask(pairing["accepted_paper_pair"]) | (
            pairing["mapping_class"].eq("ambiguous_match")
            & flag_mask(pairing["passes_paper_thresholds"])
            & ~pairing["taxonomy_validation_status"].eq("taxonomy_rejected")
        )
    profile_linkages = pairing.loc[eligible_mask].copy()
    profile_linkages = profile_linkages.dropna(
        subset=["ASV_ID", "genome_id", "mag_match_id"]
    )
    profile_linkages = profile_linkages.drop_duplicates(["ASV_ID", "mag_match_id"])

    credible = pairing.loc[
        pairing["mapping_class"].eq("ambiguous_match")
        & pairing["passes_paper_thresholds"]
        & ~pairing["taxonomy_validation_status"].eq("taxonomy_rejected")
    ].copy()
    credible = credible.dropna(subset=["ASV_ID", "genome_id", "mag_match_id"])
    credible = credible.drop_duplicates(["ASV_ID", "mag_match_id"])
    counts = credible.groupby("ASV_ID")["mag_match_id"].nunique()
    credible = credible.loc[credible["ASV_ID"].isin(counts[counts >= 2].index)].copy()
    if not profile_only:
        write_table(credible, outdir / "mapping" / f"{prefix}_credible_ambiguous_mappings.tsv")
    if not profile_only and not agreement.empty:
        write_table(
            agreement.loc[agreement["ASV_ID"].isin(credible["ASV_ID"])],
            outdir / "abundance" / f"{prefix}_ambiguous_asv_mag_abundance_agreement.tsv",
        )

    target = module_table.loc[module_table.get("module_id", pd.Series(dtype=str)).astype(str).isin(target_modules)].copy() if not module_table.empty else pd.DataFrame()
    if not profile_only and not target.empty:
        target = target.loc[target["mag_match_id"].astype(str).isin(credible["mag_match_id"].astype(str))]
        write_table(target, outdir / "functional" / f"{prefix}_ambiguous_mag_target_modules.tsv")

    mag_matrix = matrices.get("mag")
    asv_matrix = matrices.get("asv")
    samples = matrices.get("samples", [])
    pairwise_rows = []
    ranking_rows = []
    for asv_id, grp in credible.groupby("ASV_ID", sort=False):
        mags = sorted(grp["mag_match_id"].astype(str).unique())
        correlations = []
        if isinstance(mag_matrix, pd.DataFrame):
            for idx, mag_a in enumerate(mags):
                for mag_b in mags[idx + 1:]:
                    if mag_a not in mag_matrix.index or mag_b not in mag_matrix.index:
                        continue
                    rho = mag_matrix.loc[mag_a, samples].corr(mag_matrix.loc[mag_b, samples], method="spearman")
                    correlations.append(rho)
                    pairwise_rows.append({"ASV_ID": asv_id, "mag_a": mag_a, "mag_b": mag_b, "spearman_rho": rho})
        function_vectors = {}
        for mag in mags:
            vals = []
            for module in target_modules:
                hit = target.loc[(target["mag_match_id"].astype(str) == mag) & (target["module_id"].astype(str) == module)] if not target.empty else pd.DataFrame()
                vals.append(float(pd.to_numeric(hit.get("fraction_covered"), errors="coerce").max()) if not hit.empty else 0.0)
            function_vectors[mag] = vals
        functional_divergence = max(
            (float(np.max(np.abs(np.asarray(function_vectors[a]) - np.asarray(function_vectors[b])))) for i, a in enumerate(mags) for b in mags[i + 1:]),
            default=0.0,
        )
        finite_corr = [float(x) for x in correlations if pd.notna(x)]
        abundance_divergence = 1.0 - min(finite_corr) if finite_corr else 0.0
        ranking_rows.append({
            "ASV_ID": asv_id,
            "candidate_mag_count": len(mags),
            "minimum_pairwise_mag_spearman_rho": min(finite_corr) if finite_corr else np.nan,
            "abundance_divergence_score": abundance_divergence,
            "maximum_target_module_difference": functional_divergence,
            "informative_score": abundance_divergence + functional_divergence,
        })
    ranking_columns = [
        "ASV_ID",
        "candidate_mag_count",
        "minimum_pairwise_mag_spearman_rho",
        "abundance_divergence_score",
        "maximum_target_module_difference",
        "informative_score",
    ]
    ranking = pd.DataFrame(ranking_rows, columns=ranking_columns)
    if not ranking.empty:
        ranking = ranking.sort_values(
            ["informative_score", "candidate_mag_count"], ascending=False
        )
    if not profile_only:
        write_table(pd.DataFrame(pairwise_rows), outdir / "abundance" / f"{prefix}_ambiguous_mag_pairwise_abundance.tsv")
        write_table(ranking, outdir / "qc" / f"{prefix}_ambiguous_asv_ranking.tsv")

    if not isinstance(mag_matrix, pd.DataFrame) or not isinstance(asv_matrix, pd.DataFrame):
        return
    mag_all_raw = matrices.get("mag_all_raw")
    mag_all_normalized = matrices.get("mag_all_normalized")
    mag_all_size_factors = matrices.get("mag_all_size_factors")
    mag_all_normalization = str(
        matrices.get("mag_all_normalization", "unknown")
    )
    if isinstance(mag_all_size_factors, pd.Series) and not mag_all_size_factors.empty:
        write_table(
            mag_all_size_factors.rename("normalization_factor")
            .rename_axis("sample")
            .reset_index()
            .assign(normalization=mag_all_normalization),
            outdir / "abundance" / f"{prefix}_mag_normalization_factors.tsv",
        )
    cruise_dates: dict[int, pd.Timestamp] = {}
    if isinstance(cruise_metadata, pd.DataFrame) and not cruise_metadata.empty:
        lower_columns = {str(column).lower(): column for column in cruise_metadata.columns}
        cruise_col = lower_columns.get("cruise")
        date_col = next(
            (lower_columns[name] for name in ["date", "profile_date"] if name in lower_columns),
            None,
        )
        year_col = lower_columns.get("year")
        month_col = lower_columns.get("month")
        day_col = lower_columns.get("day")
        if cruise_col is not None:
            metadata = cruise_metadata.copy()
            metadata["_cruise_number"] = pd.to_numeric(metadata[cruise_col], errors="coerce")
            if date_col is not None:
                metadata["_cruise_date"] = pd.to_datetime(metadata[date_col], errors="coerce")
            elif year_col is not None and month_col is not None:
                metadata["_cruise_date"] = pd.to_datetime({
                    "year": pd.to_numeric(metadata[year_col], errors="coerce"),
                    "month": pd.to_numeric(metadata[month_col], errors="coerce"),
                    "day": (
                        pd.to_numeric(metadata[day_col], errors="coerce")
                        if day_col is not None else 1
                    ),
                }, errors="coerce")
            else:
                metadata["_cruise_date"] = pd.NaT
            metadata = metadata.dropna(subset=["_cruise_number", "_cruise_date"])
            cruise_dates = {
                int(cruise): pd.Timestamp(rows["_cruise_date"].min())
                for cruise, rows in metadata.groupby("_cruise_number")
            }
    if isinstance(mag_all_raw, pd.DataFrame) and isinstance(mag_all_normalized, pd.DataFrame):
        all_mag_samples = list(mag_all_raw.columns)
        for asv_id, grp in profile_linkages.groupby("ASV_ID", sort=False):
            species_col = "mag_tax_species"
            if species_col not in grp.columns:
                continue
            species_groups = {
                species: sorted(set(rows["mag_match_id"].astype(str)).intersection(mag_all_raw.index))
                for species, rows in grp.dropna(subset=[species_col]).groupby(species_col)
                if clean_taxon(species)
            }
            species_groups = {
                clean_taxon(species): mags
                for species, mags in species_groups.items()
                if mags
            }
            if not species_groups:
                continue
            profile_rows = []
            for species, mags in sorted(species_groups.items()):
                for sample in all_mag_samples:
                    cruise_match = re.match(r"^(SI\d+)", str(sample), flags=re.IGNORECASE)
                    depth_match = re.search(r"_(\d+(?:\.\d+)?)m(?:_|$)", str(sample), flags=re.IGNORECASE)
                    if cruise_match is None:
                        continue
                    cruise_number = int(re.search(r"\d+", cruise_match.group(1)).group())
                    cruise_date = cruise_dates.get(cruise_number)
                    if cruise_date is None:
                        continue
                    values = pd.to_numeric(
                        mag_all_raw.loc[mags, sample], errors="coerce"
                    ).dropna()
                    for mag_id, raw_count in values.items():
                        normalized_count = pd.to_numeric(
                            pd.Series([mag_all_normalized.at[mag_id, sample]]),
                            errors="coerce",
                        ).iloc[0]
                        profile_rows.append({
                            "ASV_ID": str(asv_id),
                            "species": species,
                            "mag_match_id": str(mag_id),
                            "cruise": cruise_match.group(1).upper(),
                            "cruise_number": cruise_number,
                            "date": cruise_date.date().isoformat(),
                            "month": int(cruise_date.month),
                            "sample": str(sample),
                            "depth_m": float(depth_match.group(1)) if depth_match else np.nan,
                            "raw_recruitment_count": float(raw_count),
                            "normalized_recruitment": float(normalized_count),
                            "normalization": mag_all_normalization,
                        })
            profile = pd.DataFrame(profile_rows)
            if profile.empty:
                continue
            safe_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(asv_id))[:80]
            write_table(
                profile,
                outdir / "abundance" / f"{prefix}_{safe_id}_species_recruitment_by_sample.tsv",
            )
            mag_cruise_profile = (
                profile.groupby(
                    [
                        "ASV_ID", "species", "mag_match_id",
                        "cruise", "cruise_number", "date", "month",
                    ],
                    as_index=False,
                )
                .agg(
                    sampled_depth_count=("sample", "nunique"),
                    total_raw_recruitment_count=("raw_recruitment_count", "sum"),
                    total_normalized_recruitment=(
                        "normalized_recruitment", "sum"
                    ),
                )
            )
            write_table(
                mag_cruise_profile,
                outdir / "abundance" / f"{prefix}_{safe_id}_mag_recruitment_by_cruise.tsv",
            )
            cruise_profile = (
                mag_cruise_profile.groupby(
                    [
                        "ASV_ID", "species",
                        "cruise", "cruise_number", "date", "month",
                    ],
                    as_index=False,
                )
                .agg(
                    linked_mag_count=("mag_match_id", "nunique"),
                    mean_mag_total_raw_recruitment_count=(
                        "total_raw_recruitment_count", "mean"
                    ),
                    mean_mag_total_normalized_recruitment=(
                        "total_normalized_recruitment", "mean"
                    ),
                )
            )
            write_table(
                cruise_profile,
                outdir / "abundance" / f"{prefix}_{safe_id}_lineage_recruitment_by_cruise.tsv",
            )
            species_order = sorted(profile["species"].unique())
            lineage_grays = (
                [0.0]
                if len(species_order) == 1
                else np.linspace(0.0, 1.0, len(species_order))
            )
            species_styles = {
                species: {
                    "face": matplotlib.colors.to_hex((gray, gray, gray)),
                    "edge": "white" if np.isclose(gray, 0.0) else "black",
                }
                for species, gray in zip(species_order, lineage_grays)
            }
            fig, ax = plt.subplots(figsize=(18, 9))
            month_positions = np.arange(1, 13, dtype=float)
            month_labels = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                            "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
            smooth_grid = np.linspace(1.0, 13.0, 300, endpoint=False)
            wrapped_grid = ((smooth_grid - 1.0) % 12.0) + 1.0
            bandwidth_months = 1.15
            asv_legend_handles = []

            def cyclic_kernel_smooth(months: np.ndarray, values: np.ndarray) -> np.ndarray:
                distance = np.abs(wrapped_grid[:, None] - months[None, :])
                distance = np.minimum(distance, 12.0 - distance)
                weights = np.exp(-0.5 * (distance / bandwidth_months) ** 2)
                return (weights @ values) / weights.sum(axis=1)

            asv_all_relative = matrices.get("asv_all_relative")
            if (
                isinstance(asv_all_relative, pd.DataFrame)
                and str(asv_id) in asv_all_relative.index
            ):
                asv_rows = []
                for sample, relative_abundance in pd.to_numeric(
                    asv_all_relative.loc[str(asv_id)], errors="coerce"
                ).items():
                    cruise_match = re.match(
                        r"^(SI\d+)", str(sample), flags=re.IGNORECASE
                    )
                    if cruise_match is None or pd.isna(relative_abundance):
                        continue
                    cruise_number = int(
                        re.search(r"\d+", cruise_match.group(1)).group()
                    )
                    cruise_date = cruise_dates.get(cruise_number)
                    if cruise_date is None:
                        continue
                    asv_rows.append({
                        "ASV_ID": str(asv_id),
                        "sample": str(sample),
                        "cruise": cruise_match.group(1).upper(),
                        "cruise_number": cruise_number,
                        "date": cruise_date.date().isoformat(),
                        "month": int(cruise_date.month),
                        "asv_relative_abundance": float(relative_abundance),
                    })
                asv_sample_profile = pd.DataFrame(asv_rows)
                if not asv_sample_profile.empty:
                    asv_cruise_profile = (
                        asv_sample_profile.groupby(
                            [
                                "ASV_ID", "cruise", "cruise_number",
                                "date", "month",
                            ],
                            as_index=False,
                        )
                        .agg(
                            sampled_depth_count=("sample", "nunique"),
                            mean_asv_relative_abundance=(
                                "asv_relative_abundance", "mean"
                            ),
                        )
                    )
                    write_table(
                        asv_cruise_profile,
                        outdir / "abundance"
                        / f"{prefix}_{safe_id}_asv_relative_abundance_by_cruise.tsv",
                    )
                    asv_months = asv_cruise_profile["month"].to_numpy(dtype=float)
                    asv_values = (
                        100.0
                        * asv_cruise_profile[
                            "mean_asv_relative_abundance"
                        ].to_numpy(dtype=float)
                    )
                    asv_smooth = cyclic_kernel_smooth(asv_months, asv_values)
                    asv_order = np.argsort(wrapped_grid)
                    asv_x = wrapped_grid[asv_order]
                    asv_axis = ax.twinx()
                    asv_axis.set_zorder(ax.get_zorder() - 1)
                    ax.patch.set_alpha(0.0)
                    asv_axis.plot(
                        asv_x,
                        asv_smooth[asv_order],
                        color="#666666",
                        linestyle="--",
                        linewidth=2.6,
                        alpha=0.8,
                        zorder=1,
                    )
                    asv_axis.scatter(
                        asv_months,
                        asv_values,
                        marker="D",
                        s=30,
                        facecolor="#BDBDBD",
                        edgecolor="#666666",
                        linewidth=0.7,
                        alpha=0.65,
                        zorder=2,
                    )
                    asv_axis.set_ylabel(
                        "ASV relative abundance (%)", color="#666666"
                    )
                    asv_axis.tick_params(axis="y", colors="#666666")
                    asv_axis.spines["top"].set_visible(False)
                    asv_axis.spines["right"].set_color("#666666")
                    asv_axis.set_ylim(
                        0.0,
                        max(
                            0.01,
                            float(np.nanmax(asv_values)) * 1.08,
                        ),
                    )
                    asv_legend_handles.append(
                        plt.Line2D(
                            [], [], color="#666666", linestyle="--",
                            marker="D", markerfacecolor="#BDBDBD",
                            markeredgecolor="#666666", linewidth=2.6,
                            markersize=7, label=f"{asv_id} relative abundance",
                        )
                    )

            for species_index, species in enumerate(species_order):
                species_data = cruise_profile.loc[
                    cruise_profile["species"].eq(species)
                ].copy()
                jitter = -0.06 if species_index == 0 else 0.06
                style = species_styles[species]
                ax.scatter(
                    species_data["month"] + jitter,
                    species_data[
                        "mean_mag_total_normalized_recruitment"
                    ],
                    s=38,
                    facecolor=style["face"],
                    edgecolor=style["edge"],
                    linewidth=0.65,
                    alpha=0.48,
                    zorder=3,
                )
                observed_months = species_data["month"].to_numpy(dtype=float)
                observed_values = species_data[
                    "mean_mag_total_normalized_recruitment"
                ].to_numpy(dtype=float)
                smooth = cyclic_kernel_smooth(observed_months, observed_values)
                order = np.argsort(wrapped_grid)
                x_smooth = wrapped_grid[order]
                smooth = smooth[order]
                if style["face"].lower() == "#ffffff":
                    ax.plot(x_smooth, smooth, color="black", linewidth=5.0, zorder=4)
                    ax.plot(x_smooth, smooth, color="white", linewidth=3.0, zorder=5)
                elif style["face"].lower() == "#000000":
                    ax.plot(x_smooth, smooth, color="black", linewidth=3.0, zorder=5)
                else:
                    ax.plot(x_smooth, smooth, color="black", linewidth=5.0, zorder=4)
                    ax.plot(
                        x_smooth,
                        smooth,
                        color=style["face"],
                        linewidth=3.0,
                        zorder=5,
                    )
            ax.set_xticks(month_positions, month_labels)
            ax.set_xlim(0.6, 12.4)
            ax.set_xlabel("Calendar month")
            ax.set_ylabel("Mean normalized lineage recruitment")
            ax.grid(axis="y", color="#D9D9D9", linestyle="--", linewidth=0.8, alpha=0.7)
            ax.spines[["top", "right"]].set_visible(False)
            legend_handles = []
            for species in species_order:
                style = species_styles[species]
                legend_handles.append(
                    plt.Line2D(
                        [], [], marker="o", linestyle="-",
                        color="black",
                        markerfacecolor=style["face"],
                        markeredgecolor=style["edge"],
                        linewidth=3.0,
                        markersize=10,
                        label=species,
                    )
                )
            legend_handles.extend(asv_legend_handles)
            ax.legend(
                handles=legend_handles,
                frameon=False,
                loc="upper left",
                bbox_to_anchor=(1.22, 1.0),
            )
            # Reserve separate right-side lanes for the secondary ASV axis and
            # the legend; placing both immediately at x=1 caused their text to
            # overlap in multi-lineage profiles.
            fig.tight_layout(rect=[0.0, 0.0, 0.80, 1.0])
            for fmt in formats:
                fig.savefig(
                    outdir / "abundance" / f"{prefix}_{safe_id}_species_recruitment_monthly_profile.{fmt}",
                    dpi=300,
                    bbox_inches="tight",
                )
            plt.close(fig)

    if profile_only:
        return

    display_names = {"M00529": "Denitrification", "M00973": "Anammox", "M00595": "SOX", "M00984": "S4I"}
    for asv_id in ranking["ASV_ID"]:
        if asv_id not in asv_matrix.index:
            continue
        grp = credible.loc[credible["ASV_ID"] == asv_id]
        mags = [m for m in sorted(grp["mag_match_id"].astype(str).unique()) if m in mag_matrix.index]
        if len(mags) < 2:
            continue
        # Paginate dense mappings so every eligible MAG remains legible; pagination
        # is a layout operation and never filters or ranks candidates.
        page_size = 10
        page_count = math.ceil(len(mags) / page_size)
        palette = plt.get_cmap("tab10")
        colors = {mag: palette(i % 10) for i, mag in enumerate(mags)}
        order = pd.to_numeric(asv_matrix.loc[asv_id, samples], errors="coerce").sort_values().index.tolist()
        x = np.arange(len(order))
        asv_values = pd.to_numeric(asv_matrix.loc[asv_id, order], errors="coerce").to_numpy(dtype=float)
        safe_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(asv_id))[:80]
        for page_index in range(page_count):
            page_mags = mags[page_index * page_size:(page_index + 1) * page_size]
            fig = plt.figure(figsize=(24, max(10, 0.85 * len(page_mags) + 5)), constrained_layout=True)
            grid = fig.add_gridspec(1, 3, width_ratios=[1.2, 2.8, 1.8])
            ax_link = fig.add_subplot(grid[0, 0])
            ax_abund = fig.add_subplot(grid[0, 1])
            ax_func = fig.add_subplot(grid[0, 2])

            ax_link.scatter([0], [0.5], s=900, marker="o", facecolor="white", edgecolor="black", linewidth=2.0, zorder=3)
            ax_link.text(0, 0.5, "ASV", ha="center", va="center", fontweight="bold")
            y_positions = np.linspace(0.88, 0.12, len(page_mags))
            for y, mag in zip(y_positions, page_mags):
                ax_link.plot([0.12, 0.88], [0.5, y], color="0.55", lw=1.5, zorder=1)
                ax_link.scatter([1], [y], s=650, marker="s", facecolor=colors[mag], edgecolor="black", linewidth=1.2, zorder=3)
                ax_link.text(1.13, y, mag, ha="left", va="center", fontsize=22)
            ax_link.set_xlim(-0.35, 3.3)
            ax_link.set_ylim(0, 1)
            ax_link.axis("off")
            ax_link.set_title("Candidate genome links", pad=14)

            ax_abund.plot(x, asv_values, color="black", lw=2.8, ls="--", label="ASV")
            for mag in page_mags:
                vals = pd.to_numeric(mag_matrix.loc[mag, order], errors="coerce").to_numpy(dtype=float)
                ax_abund.plot(x, vals, color=colors[mag], lw=2.2, alpha=0.9, label=mag)
            ax_abund.set_xlabel("Samples ordered by ASV relative abundance")
            ax_abund.set_ylabel("Standardized abundance profile (z-score)")
            ax_abund.set_title("Distinct abundance profiles", pad=14)
            ax_abund.tick_params(axis="x", labelbottom=False)
            ax_abund.legend(frameon=False, fontsize=22, loc="upper left")
            ax_abund.spines[["top", "right"]].set_visible(False)

            matrix = np.zeros((len(page_mags), len(target_modules)))
            for i, mag in enumerate(page_mags):
                for j, module in enumerate(target_modules):
                    hit = target.loc[(target["mag_match_id"].astype(str) == mag) & (target["module_id"].astype(str) == module)] if not target.empty else pd.DataFrame()
                    if not hit.empty:
                        matrix[i, j] = float(pd.to_numeric(hit["fraction_covered"], errors="coerce").max())
            image = ax_func.imshow(matrix, vmin=0, vmax=1, cmap="Greys", aspect="auto")
            ax_func.set_xticks(range(len(target_modules)), [display_names.get(m, m) for m in target_modules], rotation=45, ha="right")
            ax_func.set_yticks(range(len(page_mags)), page_mags)
            ax_func.set_title("Metabolic-module completeness", pad=14)
            for i in range(len(page_mags)):
                for j in range(len(target_modules)):
                    ax_func.text(j, i, f"{matrix[i, j] * 100:.0f}%", ha="center", va="center", color="white" if matrix[i, j] > 0.55 else "black", fontsize=18)
            colorbar = fig.colorbar(image, ax=ax_func, fraction=0.05, pad=0.04)
            colorbar.set_label("Completeness")
            for letter, axis in zip("ABC", [ax_link, ax_abund, ax_func]):
                axis.text(-0.08, 1.04, letter, transform=axis.transAxes, fontsize=28, fontweight="bold", va="bottom", ha="right")
            page_suffix = f"_page_{page_index + 1:02d}_of_{page_count:02d}" if page_count > 1 else ""
            stem = f"{prefix}_ambiguous_{safe_id}{page_suffix}"
            for fmt in formats:
                fig.savefig(outdir / "mapping" / f"{stem}.{fmt}", dpi=300, bbox_inches="tight")
            plt.close(fig)


def cytoscape_json(nodes: pd.DataFrame, edges: pd.DataFrame) -> dict:
    elements = {"nodes": [], "edges": []}
    for _, row in nodes.iterrows():
        data = {k: (None if pd.isna(v) else v) for k, v in row.to_dict().items()}
        data["id"] = str(data.get("id", data.get("node_id", data.get("ASV_ID", data.get("genome_id", "")))))
        elements["nodes"].append({"data": data})
    for idx, row in edges.iterrows():
        data = {k: (None if pd.isna(v) else v) for k, v in row.to_dict().items()}
        data["id"] = str(data.get("id", f"e{idx}"))
        elements["edges"].append({"data": data})
    return {"elements": elements}


def add_table_attrs(graph: nx.Graph, node_table: pd.DataFrame, id_col: str) -> None:
    for _, row in node_table.iterrows():
        node_id = str(row.get(id_col, ""))
        if node_id not in graph:
            continue
        for key, value in row.items():
            if key == id_col or pd.isna(value):
                continue
            graph.nodes[node_id][str(key)] = str(value) if isinstance(value, (list, dict)) else value


def make_paper_outputs(
    graph: nx.Graph,
    node_features: pd.DataFrame,
    validation: pd.DataFrame,
    functional: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, nx.Graph]:
    accepted = validation.loc[validation["analysis_eligible_pair"]].copy()
    if not functional.empty:
        functional_for_merge = functional.drop(columns=["genome_id"], errors="ignore")
        accepted = accepted.merge(functional_for_merge, on="mag_match_id", how="left")
    mag_cols = [
        c
        for c in accepted.columns
        if c.startswith("mag_")
        or c
        in {
            "genome_id",
            "link_pident",
            "link_qcov",
            "link_bitscore",
            "taxonomy_validation_status",
            "taxonomy_agreement",
            "lowest_comparable_rank",
            "functional_feature_count",
            "functional_summary",
        }
    ]
    multiplicity = accepted.groupby("ASV_ID").agg(
        eligible_mag_count=("mag_match_id", "nunique"),
        eligible_mag_ids=("mag_match_id", lambda values: "|".join(sorted(set(values.dropna().astype(str))))),
    ).reset_index()
    accepted_one = accepted.sort_values(["ASV_ID", "genome_id"]).drop_duplicates("ASV_ID", keep="first")
    # Validation tables already carry eligible_mag_count. Recompute both
    # multiplicity fields here from the exact accepted subset used for this
    # output, and remove stale copies before merging to avoid pandas suffixing
    # them to eligible_mag_count_x/eligible_mag_count_y.
    accepted_one = accepted_one.drop(columns=["eligible_mag_count", "eligible_mag_ids"], errors="ignore")
    accepted_one = accepted_one.merge(multiplicity, on="ASV_ID", how="left")
    mag_cols.extend(["eligible_mag_count", "eligible_mag_ids"])
    nodes = pd.DataFrame({"ASV_ID": sorted(str(n) for n in graph.nodes())})
    nodes = nodes.merge(node_features, on="ASV_ID", how="left")
    # Preserve every parsed ASV rank used by the heterogeneous-network label
    # formatter. Previously only genus/species and the unsplit taxonomy string
    # were copied, so valid phylum and family assignments appeared as
    # "unclassified" in plot callouts.
    asv_tax_cols = [
        col
        for col in [
            *[f"asv_{rank}" for rank in RANKS],
            "asv_taxonomy",
            "asv_taxonomy_confidence",
            "passes_asv_taxonomy_confidence",
        ]
        if col in validation.columns
    ]
    if asv_tax_cols:
        asv_taxonomy = (
            validation[["ASV_ID", *asv_tax_cols]]
            .drop_duplicates("ASV_ID")
        )
        nodes = nodes.merge(asv_taxonomy, on="ASV_ID", how="left")
    if mag_cols:
        nodes = nodes.merge(accepted_one[["ASV_ID", *mag_cols]], on="ASV_ID", how="left")
    nodes["node_type"] = "ASV"
    nodes["has_accepted_mag"] = nodes["genome_id"].notna() if "genome_id" in nodes.columns else False
    edges = graph_edges_to_table(graph)
    paper_graph = graph.copy()
    for node in paper_graph.nodes:
        paper_graph.nodes[node]["node_type"] = "ASV"
    add_table_attrs(paper_graph, nodes, "ASV_ID")
    for u, v, attrs in paper_graph.edges(data=True):
        attrs["edge_type"] = "asv_association"
        weight = attrs.get("weight", attrs.get("Weight", 0))
        try:
            attrs["sign"] = "positive" if float(weight) >= 0 else "negative"
        except (TypeError, ValueError):
            attrs["sign"] = "unknown"
        attrs["association_interpretation"] = "SPIEC-EASI statistical association; not a direct biological interaction"
    return nodes, edges, paper_graph


def make_heterogeneous_outputs(
    graph: nx.Graph,
    paper_nodes: pd.DataFrame,
    paper_edges: pd.DataFrame,
    validation: pd.DataFrame,
    agreement: pd.DataFrame,
    abundance_matrices: dict[str, object],
    mag_knn: int = 3,
    association_modality: str = "metagenome",
) -> tuple[pd.DataFrame, pd.DataFrame, nx.Graph, pd.DataFrame, pd.DataFrame]:
    accepted_all = validation.loc[validation["analysis_eligible_pair"]].dropna(subset=["genome_id"]).copy()
    graph_nodes = sorted(str(n) for n in graph.nodes())
    accepted = accepted_all.loc[accepted_all["ASV_ID"].astype(str).isin(graph_nodes)].copy()
    excluded = accepted_all.loc[~accepted_all["ASV_ID"].astype(str).isin(graph_nodes)].copy()
    if not excluded.empty:
        excluded["heterogeneous_exclusion_reason"] = "ASV_not_present_in_SPIEC_EASI_network"
    paired_asvs = set(accepted["ASV_ID"].dropna().astype(str))
    asv_nodes = pd.DataFrame({"id": graph_nodes, "node_type": "ASV"})
    asv_nodes = asv_nodes.merge(paper_nodes.rename(columns={"ASV_ID": "id"}), on=["id", "node_type"], how="left")
    asv_nodes["has_accepted_mag"] = asv_nodes["id"].isin(paired_asvs)
    mag_records = []
    for genome_id, grp in accepted.groupby("genome_id", dropna=True):
        row = {"id": str(genome_id), "node_type": "MAG", "genome_id": str(genome_id), "linked_asv_count": grp["ASV_ID"].nunique()}
        for rank in RANKS:
            values = sorted({clean_taxon(v) for v in grp.get(f"mag_tax_{rank}", pd.Series(dtype=object)).dropna() if clean_taxon(v)})
            if values:
                row[f"mag_tax_{rank}"] = values[0]
        for col in [c for c in grp.columns if c.startswith("mag_") and c not in row]:
            vals = [v for v in grp[col].dropna().astype(str).unique() if v and v.lower() != "nan"]
            if vals:
                row[col] = vals[0]
        mag_match_ids = sorted(set(grp["mag_match_id"].dropna().astype(str)))
        row["mag_match_id"] = "|".join(mag_match_ids)
        mag_matrix = abundance_matrices.get("mag")
        if isinstance(mag_matrix, pd.DataFrame):
            present_ids = [value for value in mag_match_ids if value in mag_matrix.index]
            if present_ids:
                values = mag_matrix.loc[present_ids].sum(axis=0)
                row["abundance_samples"] = int((values > 0).sum())
                row["mean_transformed_normalized_abundance"] = float(values.mean())
                row["max_transformed_normalized_abundance"] = float(values.max())
        mag_records.append(row)
    mag_nodes = pd.DataFrame(mag_records)
    nodes = pd.concat([asv_nodes, mag_nodes], ignore_index=True, sort=False)

    seq_edges = accepted.loc[accepted["ASV_ID"].astype(str).isin(graph_nodes)].copy()
    if not agreement.empty:
        # Inputs reused from an earlier published run may already carry these
        # columns. Replace them with the current abundance calculation rather
        # than allowing pandas to create ambiguous _x/_y copies.
        seq_edges = seq_edges.drop(
            columns=["shared_samples", "abundance_metric", "transform", "spearman_rho"],
            errors="ignore",
        )
        agreement_cols = [
            "ASV_ID", "genome_id", "mag_match_id", "shared_samples",
            "abundance_metric", "transform", "spearman_rho",
        ]
        seq_edges = seq_edges.merge(
            agreement[[c for c in agreement_cols if c in agreement]].drop_duplicates(
                ["ASV_ID", "genome_id", "mag_match_id"]
            ),
            on=["ASV_ID", "genome_id", "mag_match_id"],
            how="left",
        )
    seq_edges = seq_edges.assign(
        source=seq_edges["ASV_ID"].astype(str),
        target=seq_edges["genome_id"].astype(str),
        edge_type="sequence_match",
        evidence_type="near_exact_16S_alignment",
        identity=seq_edges.get("link_pident", pd.NA),
        coverage=seq_edges.get("link_qcov", pd.NA),
    )
    seq_edge_cols = [
        "source",
        "target",
        "edge_type",
        "evidence_type",
        "identity",
        "coverage",
        "taxonomy_agreement",
        "taxonomy_validation_status",
        "lowest_comparable_rank",
        "shared_samples",
        "abundance_metric",
        "transform",
        "spearman_rho",
    ]
    mag_matrix = abundance_matrices.get("mag")
    samples = abundance_matrices.get("samples", [])
    pairwise_rows = []
    if isinstance(mag_matrix, pd.DataFrame):
        for asv_id, group in seq_edges.groupby("ASV_ID", sort=False):
            candidates = (
                group[["genome_id", "mag_match_id"]]
                .dropna()
                .drop_duplicates()
            )
            candidates = candidates.loc[candidates["mag_match_id"].astype(str).isin(mag_matrix.index)]
            records = candidates.to_dict("records")
            for left_index, left in enumerate(records):
                for right in records[left_index + 1:]:
                    left_match = str(left["mag_match_id"])
                    right_match = str(right["mag_match_id"])
                    rho = mag_matrix.loc[left_match, samples].corr(
                        mag_matrix.loc[right_match, samples], method="spearman",
                    )
                    pairwise_rows.append({
                        "ASV_ID": str(asv_id),
                        "source": str(left["genome_id"]),
                        "target": str(right["genome_id"]),
                        "source_mag_match_id": left_match,
                        "target_mag_match_id": right_match,
                        "spearman_rho": rho,
                        "shared_samples": len(samples),
                        "data_modality": association_modality,
                    })
    mag_pairwise = pd.DataFrame(pairwise_rows)
    mag_layout_edges = pd.DataFrame()
    if not mag_pairwise.empty:
        mag_pairwise["spearman_rho"] = pd.to_numeric(
            mag_pairwise["spearman_rho"], errors="coerce"
        )
        finite_pairwise = mag_pairwise.loc[np.isfinite(mag_pairwise["spearman_rho"])].copy()
        directed = pd.concat([
            finite_pairwise.assign(anchor=finite_pairwise["source"], neighbor=finite_pairwise["target"]),
            finite_pairwise.assign(anchor=finite_pairwise["target"], neighbor=finite_pairwise["source"]),
        ], ignore_index=True)
        positive = directed.loc[directed["spearman_rho"] > 0].sort_values(
            ["ASV_ID", "anchor", "spearman_rho", "neighbor"],
            ascending=[True, True, False, True],
        )
        negative = directed.loc[directed["spearman_rho"] < 0].sort_values(
            ["ASV_ID", "anchor", "spearman_rho", "neighbor"],
            ascending=[True, True, True, True],
        )
        per_sign_n = max(1, int(mag_knn))
        selected = pd.concat([
            positive.groupby(["ASV_ID", "anchor"], sort=False).head(per_sign_n),
            negative.groupby(["ASV_ID", "anchor"], sort=False).head(per_sign_n),
        ], ignore_index=True)
        selected["_edge_key"] = selected.apply(
            lambda row: "||".join(sorted([str(row["anchor"]), str(row["neighbor"])])),
            axis=1,
        )
        mag_layout_edges = (
            selected.sort_values(["spearman_rho", "ASV_ID"], ascending=[False, True])
            .drop_duplicates(["ASV_ID", "_edge_key"])
            .assign(layout_source=lambda frame: frame["anchor"],
                    layout_target=lambda frame: frame["neighbor"])
        )
        mag_layout_edges = mag_layout_edges.assign(
            edge_type="mag_abundance_association",
            evidence_type=(
                "MAG_normalized_expression_spearman"
                if association_modality == "metatranscriptome"
                else "MAG_normalized_abundance_spearman"
            ),
            data_modality=association_modality,
        )
        mag_layout_edges = mag_layout_edges[[
            "layout_source", "layout_target", "edge_type", "evidence_type", "ASV_ID",
            "spearman_rho", "shared_samples", "data_modality",
        ]].rename(columns={"layout_source": "source", "layout_target": "target"})
        expected_multi = set(
            seq_edges.groupby("ASV_ID")["genome_id"].nunique().loc[lambda values: values >= 2].index.astype(str)
        )
        represented_multi = set(mag_layout_edges["ASV_ID"].dropna().astype(str))
        missing_multi = sorted(expected_multi - represented_multi)
        if missing_multi:
            die(
                "MAG abundance-subnetwork integrity failure: no finite MAG-MAG "
                f"abundance edge was available for {len(missing_multi)} multi-MAG ASVs. "
                f"Examples: {', '.join(missing_multi[:5])}"
            )
    edges = pd.concat(
        [
            paper_edges,
            seq_edges[[c for c in seq_edge_cols if c in seq_edges.columns]],
            mag_layout_edges,
        ],
        ignore_index=True,
        sort=False,
    )
    hetero = nx.Graph()
    for _, row in nodes.iterrows():
        node_id = str(row["id"])
        attrs = {k: v for k, v in row.to_dict().items() if k != "id" and not pd.isna(v)}
        hetero.add_node(node_id, **attrs)
    for _, row in edges.iterrows():
        source = str(row["source"])
        target = str(row["target"])
        attrs = {k: v for k, v in row.to_dict().items() if k not in {"source", "target"} and not pd.isna(v)}
        if (
            hetero.has_edge(source, target)
            and attrs.get("edge_type") == "mag_abundance_association"
            and hetero.edges[source, target].get("edge_type") == "mag_abundance_association"
        ):
            prior = hetero.edges[source, target]
            asv_ids = set(str(prior.get("ASV_IDs", prior.get("ASV_ID", ""))).split("|"))
            asv_ids.add(str(attrs.get("ASV_ID", "")))
            prior.update(attrs)
            prior["ASV_IDs"] = "|".join(sorted(value for value in asv_ids if value))
            continue
        if attrs.get("edge_type") == "mag_abundance_association":
            attrs["ASV_IDs"] = str(attrs.get("ASV_ID", ""))
        hetero.add_edge(source, target, **attrs)
    sequence_pairs = {
        tuple(sorted((str(row["source"]), str(row["target"]))))
        for _, row in seq_edges.iterrows()
    }
    orphan_mags = [
        node for node, data in hetero.nodes(data=True)
        if data.get("node_type") == "MAG"
        and not any(tuple(sorted((str(node), str(neighbor)))) in sequence_pairs for neighbor in hetero.neighbors(node))
    ]
    if orphan_mags:
        die(
            f"Heterogeneous-network integrity failure: {len(orphan_mags)} MAG nodes "
            f"lack an ASV-MAG sequence edge. Examples: {', '.join(orphan_mags[:5])}"
        )
    sequence_groups: dict[str, set[str]] = defaultdict(set)
    for left, right, attrs in hetero.edges(data=True):
        if attrs.get("edge_type") != "sequence_match":
            continue
        if hetero.nodes[left].get("node_type") == "ASV":
            sequence_groups[str(left)].add(str(right))
        else:
            sequence_groups[str(right)].add(str(left))
    disconnected_local_mags = []
    for asv_id, mags in sequence_groups.items():
        if len(mags) < 2:
            continue
        for mag_id in mags:
            has_local_association = any(
                hetero.edges[mag_id, neighbor].get("edge_type")
                == "mag_abundance_association"
                and neighbor in mags
                and asv_id in set(
                    str(
                        hetero.edges[mag_id, neighbor].get(
                            "ASV_IDs",
                            hetero.edges[mag_id, neighbor].get("ASV_ID", ""),
                        )
                    ).split("|")
                )
                for neighbor in hetero.neighbors(mag_id)
            )
            if not has_local_association:
                disconnected_local_mags.append((asv_id, mag_id))
    if disconnected_local_mags:
        preview = ", ".join(
            f"{asv_id}/{mag_id}"
            for asv_id, mag_id in disconnected_local_mags[:5]
        )
        die(
            "Heterogeneous-network integrity failure: "
            f"{len(disconnected_local_mags)} MAG nodes in multi-MAG subnetworks "
            "lack a retained within-subnetwork abundance association. "
            f"Examples: {preview}"
        )
    return nodes, edges, hetero, mag_pairwise, excluded


def plot_heterogeneous_network(
    graph: nx.Graph,
    outdir: Path,
    prefix: str,
    formats: list[str],
    overlay_mode: str = "base",
    functional_module_id: str | None = None,
    functional_module_name: str = "",
    minimum_fraction: float = 0.5,
    biochemical_grouping: str = "",
    biochemical_associations: dict[str, str] | None = None,
    biochemical_palette: dict[str, str] | None = None,
    biochemical_order: list[str] | None = None,
    titan_responses: dict[str, dict[str, object]] | None = None,
    titan_variable: str = "",
) -> None:
    """Render all evidence layers on one identical heterogeneous-network layout."""
    if graph.number_of_nodes() == 0:
        return
    asv_nodes_all = [
        node for node, data in graph.nodes(data=True)
        if data.get("node_type") == "ASV"
    ]
    mag_nodes = [
        node for node, data in graph.nodes(data=True)
        if data.get("node_type") == "MAG"
    ]
    module_levels = sorted({
        str(graph.nodes[node].get("display_module_label"))
        for node in asv_nodes_all
        if str(graph.nodes[node].get("display_module_label", "")).lower()
        not in {"", "nan", "none", "other modules"}
    }, key=module_sort_key)
    module_palette = {
        module: matplotlib.colors.to_hex(plt.get_cmap("tab20")(index % 20))
        for index, module in enumerate(module_levels)
    }
    biochemical_associations = biochemical_associations or {}
    biochemical_palette = biochemical_palette or {}
    biochemical_order = biochemical_order or []
    titan_responses = titan_responses or {}

    def pale_color(color: str, white_fraction: float = 0.68) -> str:
        rgb = np.asarray(matplotlib.colors.to_rgb(color), dtype=float)
        return matplotlib.colors.to_hex(
            rgb * (1.0 - white_fraction) + white_fraction
        )
    qualifying_mags = {
        str(node) for node in mag_nodes
        if functional_module_id
        and node_cycle_coverage(graph.nodes[node], functional_module_id) >= minimum_fraction
    }
    asv_graph = graph.subgraph(asv_nodes_all).copy()
    # This layout is computed solely from ASV-ASV edges. MAGs therefore cannot
    # alter the ASV network geometry.
    asv_pos = nx.spring_layout(asv_graph, seed=42)
    asv_pos = {node: np.asarray(coords, dtype=float) * 0.72 for node, coords in asv_pos.items()}

    mag_groups: dict[str, list[str]] = defaultdict(list)
    for left, right, attrs in graph.edges(data=True):
        if attrs.get("edge_type") != "sequence_match":
            continue
        if graph.nodes[left].get("node_type") == "ASV":
            asv_id, mag_id = str(left), str(right)
        else:
            asv_id, mag_id = str(right), str(left)
        mag_groups[asv_id].append(mag_id)
    mag_groups = {
        asv_id: sorted(set(mags))
        for asv_id, mags in mag_groups.items()
        if asv_id in asv_pos and mags
    }
    qualifying_asvs = {
        asv_id for asv_id, mags in mag_groups.items()
        if any(str(mag) in qualifying_mags for mag in mags)
    }
    paired_asv_nodes = sorted(mag_groups)
    unpaired_asv_nodes = [node for node in asv_nodes_all if node not in mag_groups]

    taxon_by_mag = {}
    for node in mag_nodes:
        species = clean_taxon(graph.nodes[node].get("mag_tax_species"))
        taxon_by_mag[node] = species or "unclassified species"

    def compact_taxonomy(node: str, node_type: str) -> str:
        data = graph.nodes[node]
        prefix = "asv_" if node_type == "ASV" else "mag_tax_"
        values = []
        ranks = [("P", "phylum"), ("F", "family"), ("G", "genus")]
        if node_type == "MAG":
            ranks.append(("S", "species"))
        for abbreviation, rank in ranks:
            value = clean_taxon(data.get(f"{prefix}{rank}"))
            values.append(f"{abbreviation}: {value or 'unclassified'}")
        return "; ".join(values)

    # Allocate callouts separately across the upper and lower network border.
    # Within each half, preserve the ASVs' left-to-right order. Thus upper
    # ASVs always connect to upper callouts, lower ASVs to lower callouts, and
    # the straight anchors do not cross because their order is unchanged.
    center = np.mean(np.asarray(list(asv_pos.values())), axis=0)
    central_network_radius = max(
        float(np.linalg.norm(coords - center))
        for coords in asv_pos.values()
    )
    angular_order = sorted(
        paired_asv_nodes,
        key=lambda node: math.atan2(
            asv_pos[node][1] - center[1],
            asv_pos[node][0] - center[0],
        ),
    )
    upper_nodes = sorted(
        [node for node in paired_asv_nodes if asv_pos[node][1] >= center[1]],
        key=lambda node: asv_pos[node][0],
    )
    lower_nodes = sorted(
        [node for node in paired_asv_nodes if asv_pos[node][1] < center[1]],
        key=lambda node: asv_pos[node][0],
    )
    border_margin = math.radians(18.0)
    assigned_angles: dict[str, float] = {}

    def segments_cross(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray) -> bool:
        def orientation(p: np.ndarray, q: np.ndarray, r: np.ndarray) -> float:
            qp = q - p
            rp = r - p
            return float(qp[0] * rp[1] - qp[1] * rp[0])

        return (
            orientation(a, b, c) * orientation(a, b, d) < 0
            and orientation(c, d, a) * orientation(c, d, b) < 0
        )

    def optimize_border_assignment(nodes: list[str], slot_angles: np.ndarray) -> None:
        if not nodes:
            return
        if len(nodes) == 1:
            assigned_angles[nodes[0]] = float(slot_angles[0])
            return
        slot_radius = central_network_radius + 0.72
        slot_centers = [
            center + slot_radius * np.asarray([math.cos(angle), math.sin(angle)])
            for angle in slot_angles
        ]
        length_cost = np.asarray([
            [
                float(np.sum((slot_center - asv_pos[node]) ** 2))
                for slot_center in slot_centers
            ]
            for node in nodes
        ])
        crossing_lookup: dict[tuple[int, int, int, int], bool] = {}
        for left_index in range(len(nodes)):
            for right_index in range(left_index + 1, len(nodes)):
                for left_slot in range(len(slot_centers)):
                    for right_slot in range(len(slot_centers)):
                        crossing_lookup[
                            (left_index, left_slot, right_index, right_slot)
                        ] = segments_cross(
                            asv_pos[nodes[left_index]], slot_centers[left_slot],
                            asv_pos[nodes[right_index]], slot_centers[right_slot],
                        )
        # Exhaustively optimize the small callout sets used here. For unusually
        # large sets, retain the monotone assignment to avoid factorial cost.
        candidate_orders = permutations(range(len(nodes))) if len(nodes) <= 8 else [range(len(nodes))]
        best_score = math.inf
        best_order = tuple(range(len(nodes)))
        for order in candidate_orders:
            crossings = sum(
                crossing_lookup[
                    (left_index, order[left_index], right_index, order[right_index])
                ]
                for left_index in range(len(nodes))
                for right_index in range(left_index + 1, len(nodes))
            )
            squared_length = sum(
                length_cost[node_index, slot_index]
                for node_index, slot_index in enumerate(order)
            )
            score = crossings * 1_000_000.0 + squared_length
            if score < best_score:
                best_score = score
                best_order = tuple(order)
        for node, slot_index in zip(nodes, best_order):
            assigned_angles[node] = float(slot_angles[slot_index])

    if upper_nodes:
        if len(upper_nodes) == 1:
            natural = math.atan2(
                asv_pos[upper_nodes[0]][1] - center[1],
                asv_pos[upper_nodes[0]][0] - center[0],
            )
            upper_angles = np.asarray([
                np.clip(natural, border_margin, math.pi - border_margin)
            ])
        else:
            upper_angles = np.linspace(
                math.pi - border_margin, border_margin, len(upper_nodes)
            )
        optimize_border_assignment(upper_nodes, upper_angles)
    if lower_nodes:
        if len(lower_nodes) == 1:
            natural = math.atan2(
                asv_pos[lower_nodes[0]][1] - center[1],
                asv_pos[lower_nodes[0]][0] - center[0],
            )
            lower_angles = np.asarray([
                np.clip(natural, -math.pi + border_margin, -border_margin)
            ])
        else:
            lower_angles = np.linspace(
                -math.pi + border_margin, -border_margin, len(lower_nodes)
            )
        optimize_border_assignment(lower_nodes, lower_angles)
    callouts = {}
    occupied: list[tuple[np.ndarray, float, float]] = []
    placement_order = sorted(
        angular_order,
        key=lambda asv_id: len(mag_groups[asv_id]),
        reverse=True,
    )
    for asv_id in placement_order:
        desired_angle = assigned_angles[asv_id]
        n_mags = len(mag_groups[asv_id])
        is_singleton = n_mags == 1
        width = 0.10 if is_singleton else min(0.96, 0.25 + 0.085 * math.sqrt(n_mags))
        height = 0.10 if is_singleton else min(0.74, 0.21 + 0.065 * math.sqrt(n_mags))
        radial = np.asarray([math.cos(desired_angle), math.sin(desired_angle)])
        ellipse_inward_extent = math.sqrt(
            ((width / 2) * radial[0]) ** 2 +
            ((height / 2) * radial[1]) ** 2
        )
        minimum_radius = central_network_radius + ellipse_inward_extent + 0.16
        callout_center = None
        for radius in np.arange(minimum_radius, 4.81, 0.08):
            candidate = center + float(radius) * radial
            overlaps = any(
                abs(candidate[0] - other_center[0]) < (width + other_width) / 2 + 0.06
                and abs(candidate[1] - other_center[1]) < (height + other_height) / 2 + 0.06
                for other_center, other_width, other_height in occupied
            )
            if not overlaps:
                callout_center = candidate
                break
        if callout_center is None:
            callout_center = center + 4.88 * radial
        occupied.append((callout_center, width, height))
        callouts[asv_id] = {
            "center": callout_center,
            "radial": radial,
            "width": width,
            "height": height,
            "singleton": is_singleton,
        }

    # These are a network plus its key, not independent analytical panels.
    # Use manually positioned axes so the shared publication-style hook does
    # not add subplot letters, and devote the canvas only to actual content.
    if overlay_mode in {"modules", "ecological", "biochemical"}:
        fig = plt.figure(figsize=(24, 18))
        ax = fig.add_axes([0.02, 0.04, 0.67, 0.92])
        legend_ax = fig.add_axes([0.71, 0.04, 0.28, 0.92])
    else:
        fig = plt.figure(figsize=(34, 18))
        ax = fig.add_axes([0.02, 0.04, 0.57, 0.92])
        legend_ax = fig.add_axes([0.61, 0.04, 0.38, 0.92])
    legend_ax.set_axis_off()
    taxonomy_legend_groups = []
    association_edges = [
        (u, v) for u, v, data in asv_graph.edges(data=True)
        if data.get("edge_type") == "asv_association"
    ]
    background_edges = nx.draw_networkx_edges(
        asv_graph, asv_pos, edgelist=association_edges,
        edge_color="#BDBDBD", width=0.8, alpha=0.55, ax=ax,
    )
    if background_edges is not None:
        if isinstance(background_edges, list):
            for artist in background_edges:
                artist.set_zorder(1.0)
        else:
            background_edges.set_zorder(1.0)
    def asv_face(node: str, paired: bool) -> str:
        if overlay_mode in {"modules", "ecological", "biochemical"}:
            return module_palette.get(
                str(graph.nodes[node].get("display_module_label")), "#D9D9D9"
            )
        if overlay_mode == "functional" and node in qualifying_asvs:
            return module_palette.get(
                str(graph.nodes[node].get("display_module_label")), "#808080"
            )
        if overlay_mode == "functional":
            module = str(graph.nodes[node].get("display_module_label"))
            return pale_color(module_palette[module]) if module in module_palette else "#E6E6E6"
        if overlay_mode == "titan":
            direction = str(titan_responses.get(str(node), {}).get("response_direction", ""))
            return TITAN_DIRECTION_COLORS.get(direction, "#F2F2F2")
        return "#BDBDBD" if paired else "#F2F2F2"

    def asv_edge(node: str, paired: bool) -> str:
        if overlay_mode == "biochemical":
            association = biochemical_associations.get(str(node), "")
            if not association:
                return "#BDBDBD"
            if "+" in association:
                return "#000000"
            return biochemical_palette.get(association, "#4D4D4D")
        if overlay_mode in {"modules", "ecological"}:
            if (
                overlay_mode == "ecological"
                and paired
                and bool(graph.nodes[node].get("is_ecological_anchor"))
            ):
                return "black"
            return "#A6A6A6"
        if overlay_mode == "functional" and node in qualifying_asvs:
            return "black"
        if overlay_mode == "titan":
            face = asv_face(node, paired)
            if matplotlib.colors.to_hex(face).lower() == "#000000":
                return node_outline_color(face)
            return "black" if str(node) in titan_responses or paired else "#D9D9D9"
        return "black" if paired else "#D9D9D9"

    def asv_size(node: str, paired: bool) -> float:
        if overlay_mode == "titan":
            return titan_marker_area(titan_responses.get(str(node), {}).get("z_score"))
        if overlay_mode in {"modules", "ecological", "biochemical"}:
            return 300.0 if (
                overlay_mode == "ecological"
                and bool(graph.nodes[node].get("is_ecological_anchor"))
            ) else 70.0
        return 520.0 if paired else 65.0

    def asv_linewidth(node: str, paired: bool) -> float:
        if overlay_mode == "titan":
            return 0.70 if matplotlib.colors.to_hex(asv_face(node, paired)).lower() == "#000000" else (1.25 if paired else 0.55)
        if overlay_mode in {"modules", "ecological"}:
            return 2.4 if (
                overlay_mode == "ecological" and paired
                and bool(graph.nodes[node].get("is_ecological_anchor"))
            ) else 0.35
        if overlay_mode == "biochemical":
            return 2.2 if biochemical_associations.get(str(node), "") else 0.35
        if overlay_mode == "functional":
            return 2.8 if node in qualifying_asvs else (0.8 if paired else 0.3)
        return 0.8 if paired else 0.3

    unpaired_collection = nx.draw_networkx_nodes(
        asv_graph, asv_pos, nodelist=unpaired_asv_nodes,
        node_shape="o",
        node_size=[asv_size(node, False) for node in unpaired_asv_nodes],
        node_color=[asv_face(node, False) for node in unpaired_asv_nodes],
        edgecolors=[asv_edge(node, False) for node in unpaired_asv_nodes],
        linewidths=[asv_linewidth(node, False) for node in unpaired_asv_nodes], ax=ax,
    )
    if unpaired_collection is not None:
        unpaired_collection.set_zorder(2.0)
    paired_collection = nx.draw_networkx_nodes(
        asv_graph, asv_pos, nodelist=paired_asv_nodes,
        node_shape="o",
        node_size=[asv_size(node, True) for node in paired_asv_nodes],
        node_color=[asv_face(node, True) for node in paired_asv_nodes],
        edgecolors=[asv_edge(node, True) for node in paired_asv_nodes],
        linewidths=[asv_linewidth(node, True) for node in paired_asv_nodes], ax=ax,
    )
    if paired_collection is not None:
        # Linked ASVs remain visible even where the fixed network layout places
        # another node or association immediately beneath them.
        paired_collection.set_zorder(7.0)

    label_positions: list[np.ndarray] = []
    for group_index, asv_id in enumerate(angular_order):
        if overlay_mode in {"modules", "ecological", "biochemical"}:
            continue
        details = callouts[asv_id]
        callout_center = details["center"]
        radial = details["radial"]
        width = details["width"]
        height = details["height"]
        is_singleton = bool(details["singleton"])
        if not is_singleton:
            highlights_function = (
                overlay_mode == "functional" and asv_id in qualifying_asvs
            )
            ellipse = matplotlib.patches.Ellipse(
                callout_center, width=width, height=height,
                facecolor="white",
                edgecolor="black" if highlights_function else "#4D4D4D",
                linewidth=2.4 if highlights_function else 0.9,
                zorder=4.0,
            )
            ax.add_patch(ellipse)

        # Straight radial anchor from the ASV to its own callout boundary.
        anchor = asv_pos[asv_id]
        toward_asv = anchor - callout_center
        norm = np.linalg.norm(toward_asv)
        unit = toward_asv / norm if norm > 0 else -radial
        if is_singleton:
            boundary = callout_center + unit * 0.025
        else:
            boundary_scale = 1.0 / math.sqrt(
                (unit[0] / (width / 2)) ** 2 +
                (unit[1] / (height / 2)) ** 2
            )
            boundary = callout_center + unit * boundary_scale
        ax.plot(
            [anchor[0], boundary[0]],
            [anchor[1], boundary[1]],
            color=(
                "black"
                if overlay_mode == "functional" and asv_id in qualifying_asvs
                else "#4D4D4D"
            ),
            linewidth=(
                2.4
                if overlay_mode == "functional" and asv_id in qualifying_asvs
                else 1.0
            ),
            # Keep paired-ASV connectors legible above the background
            # SPIEC-EASI network. MAG callout ellipses and their contents are
            # drawn on still higher layers, so a connector cannot obscure a
            # callout that happens to lie along its route.
            zorder=3.0,
            solid_capstyle="round",
        )

        mags = mag_groups[asv_id]
        local_taxa = sorted(set(taxon_by_mag[node] for node in mags))
        linked_asv_color = asv_face(asv_id, True)
        local_taxon_colors = {taxon: linked_asv_color for taxon in local_taxa}
        local_taxon_shapes = species_marker_map(local_taxa)
        local_graph = nx.Graph()
        local_graph.add_nodes_from(mags)
        for left, right, attrs in graph.edges(data=True):
            if attrs.get("edge_type") != "mag_abundance_association":
                continue
            edge_asv_ids = set(
                str(attrs.get("ASV_IDs", attrs.get("ASV_ID", ""))).split("|")
            )
            if asv_id not in edge_asv_ids:
                continue
            if left in local_graph and right in local_graph:
                rho = pd.to_numeric(attrs.get("spearman_rho"), errors="coerce")
                local_graph.add_edge(
                    left, right,
                    spearman_rho=rho,
                    layout_weight=0.05 + 2.5 * max(float(rho), 0.0) if pd.notna(rho) else 0.05,
                )
        if len(mags) == 1:
            local_raw = {mags[0]: np.zeros(2)}
        elif local_graph.number_of_edges() > 0:
            local_raw = nx.spring_layout(
                local_graph, seed=42 + group_index,
                weight="layout_weight",
                k=max(0.35, 1.8 / math.sqrt(len(mags))),
                iterations=400,
            )
        else:
            local_raw = nx.circular_layout(local_graph)
        local_array = np.asarray(list(local_raw.values()), dtype=float)
        # Rescale each local solution to use the available ellipse interior.
        # Dense, highly correlated subnetworks otherwise collapse into a tiny
        # central knot whose edges disappear underneath the MAG markers.
        x_scale = max(float(np.ptp(local_array[:, 0])), 0.05)
        y_scale = max(float(np.ptp(local_array[:, 1])), 0.05)
        local_pos = {
            node: callout_center + np.asarray([
                0.56 * width * coords[0] / x_scale,
                0.56 * height * coords[1] / y_scale,
            ])
            for node, coords in local_raw.items()
        }
        # Resolve near-identical projected coordinates without changing the
        # abundance-derived ordering more than necessary. This is especially
        # important when many genomes have nearly identical abundance profiles.
        if len(mags) > 1:
            local_rng = np.random.default_rng(4200 + group_index)
            normalized = {
                node: np.asarray([
                    (coords[0] - callout_center[0]) / (0.38 * width),
                    (coords[1] - callout_center[1]) / (0.38 * height),
                ])
                for node, coords in local_pos.items()
            }
            min_separation = 0.12 if len(mags) > 40 else 0.16
            for _ in range(180):
                movement = {node: np.zeros(2) for node in mags}
                changed = False
                for left_index, left in enumerate(mags):
                    for right in mags[left_index + 1:]:
                        delta = normalized[right] - normalized[left]
                        distance = float(np.linalg.norm(delta))
                        if distance >= min_separation:
                            continue
                        if distance < 1e-9:
                            delta = local_rng.normal(size=2)
                            distance = float(np.linalg.norm(delta))
                        push = 0.5 * (min_separation - distance) * delta / distance
                        movement[left] -= push
                        movement[right] += push
                        changed = True
                if not changed:
                    break
                for node in mags:
                    normalized[node] += 0.35 * movement[node]
                    radius = float(np.linalg.norm(normalized[node]))
                    if radius > 0.90:
                        normalized[node] *= 0.90 / radius
            local_pos = {
                node: callout_center + np.asarray([
                    coords[0] * 0.38 * width,
                    coords[1] * 0.38 * height,
                ])
                for node, coords in normalized.items()
            }
        positive = [
            (u, v) for u, v, data in local_graph.edges(data=True)
            if pd.notna(data.get("spearman_rho")) and float(data["spearman_rho"]) >= 0
        ]
        negative = [
            (u, v) for u, v, data in local_graph.edges(data=True)
            if pd.notna(data.get("spearman_rho")) and float(data["spearman_rho"]) < 0
        ]
        association_styles = (
            (positive, "solid", "#D55E00"),
            (negative, "dashed", "#0072B2"),
        )
        for edges, style, color in association_styles:
            if edges:
                # White underlay separates association edges from grayscale
                # species nodes and dense neighboring edges.
                white_edges = nx.draw_networkx_edges(
                    local_graph, local_pos, edgelist=edges,
                    edge_color="white", width=2.8, style=style,
                    alpha=0.75, ax=ax,
                )
                color_edges = nx.draw_networkx_edges(
                    local_graph, local_pos, edgelist=edges,
                    edge_color=color, width=1.4, style=style,
                    alpha=0.65, ax=ax,
                )
                if white_edges is not None:
                    white_edges.set_zorder(5.0)
                if color_edges is not None:
                    color_edges.set_zorder(5.1)
        for taxon in local_taxa:
            taxon_nodes = [node for node in mags if taxon_by_mag[node] == taxon]
            local_nodes = nx.draw_networkx_nodes(
                local_graph, local_pos, nodelist=taxon_nodes,
                node_shape=local_taxon_shapes[taxon],
                node_size=72,
                node_color=[local_taxon_colors[taxon]] * len(taxon_nodes),
                edgecolors=[node_outline_color(local_taxon_colors[taxon])] * len(taxon_nodes),
                linewidths=[0.65] * len(taxon_nodes), ax=ax,
            )
            if local_nodes is not None:
                local_nodes.set_zorder(6.0)

        # Keep the plot callout itself compact; the detailed ASV/MAG taxonomy
        # key is assembled in a dedicated panel to the right.
        if is_singleton:
            outward_extent = 0.04
        else:
            outward_extent = 1.0 / math.sqrt(
                (radial[0] / (width / 2)) ** 2
                + (radial[1] / (height / 2)) ** 2
            )
        # Start immediately outside the owning ellipse; collision handling may
        # move the label farther outward, but never laterally toward another
        # callout where its ownership would become ambiguous.
        label_position = callout_center + radial * (outward_extent + 0.045)
        # A label belongs outside its own callout and must not cover any other
        # MAG ellipse. Move it farther along the same straight radial direction
        # until its padded center clears every callout and prior label.
        for _ in range(80):
            overlaps_callout = any(
                abs(label_position[0] - other["center"][0])
                < other["width"] / 2 + 0.24
                and abs(label_position[1] - other["center"][1])
                < other["height"] / 2 + 0.12
                for other_asv, other in callouts.items()
                if other_asv != asv_id
            )
            overlaps_label = any(
                abs(label_position[0] - prior[0]) < 0.34
                and abs(label_position[1] - prior[1]) < 0.13
                for prior in label_positions
            )
            if not overlaps_callout and not overlaps_label:
                break
            label_position = label_position + radial * 0.06
        label_positions.append(label_position.copy())
        asv_taxon = compact_taxonomy(asv_id, "ASV")
        ecological_module = str(
            graph.nodes[asv_id].get("module_label", "unassigned")
        )
        linked_asv_label = f"{asv_id} [{ecological_module}]"
        legend_rows = [
            TextArea(
                f"{asv_id} — {asv_taxon}",
                textprops={"fontsize": 16, "fontweight": "bold"},
            ),
        ]
        for taxon in local_taxa:
            representative_mag = next(
                node for node in mags if taxon_by_mag[node] == taxon
            )
            mag_taxon = compact_taxonomy(representative_mag, "MAG")
            marker = DrawingArea(14, 14, 0, 0)
            if local_taxon_shapes[taxon] == "s":
                marker_patch = matplotlib.patches.Rectangle(
                    (2.5, 2.5), 9, 9,
                    facecolor=local_taxon_colors[taxon],
                    edgecolor=node_outline_color(local_taxon_colors[taxon]),
                    linewidth=0.8,
                )
            else:
                marker_patch = matplotlib.patches.Circle(
                    (7, 7), radius=4.5,
                    facecolor=local_taxon_colors[taxon],
                    edgecolor=node_outline_color(local_taxon_colors[taxon]),
                    linewidth=0.8,
                )
            marker.add_artist(marker_patch)
            legend_rows.append(HPacker(
                children=[
                    marker,
                    TextArea(f" MAG: {mag_taxon}", textprops={"fontsize": 12}),
                ],
                align="center", pad=0, sep=1,
            ))
        taxonomy_legend_groups.append((
            ecological_module,
            asv_id,
            VPacker(children=legend_rows, align="left", pad=0, sep=1),
            len(legend_rows),
        ))
        horizontal = "left" if radial[0] > 0.2 else "right" if radial[0] < -0.2 else "center"
        vertical = "bottom" if radial[1] > 0.2 else "top" if radial[1] < -0.2 else "center"
        ax.text(
            label_position[0], label_position[1], linked_asv_label,
            ha=horizontal, va=vertical, fontsize=16, fontweight="bold",
            bbox={"boxstyle": "round,pad=0.12", "facecolor": "white", "edgecolor": "none"},
            zorder=10,
        )
    network_handles = [
        plt.Line2D([], [], marker="o", linestyle="", markerfacecolor="#BDBDBD",
                   markeredgecolor="black", markersize=14, label="MAG-paired ASV"),
        plt.Line2D([], [], marker="o", linestyle="", markerfacecolor="#F2F2F2",
                   markeredgecolor="#D9D9D9", markersize=7, label="Unpaired ASV"),
        plt.Line2D([], [], color="#BDBDBD", linewidth=2, label="SPIEC-EASI ASV association"),
        plt.Line2D([], [], color="#666666", linewidth=1.2, label="ASV-to-MAG-subnetwork anchor"),
        plt.Line2D([], [], color="#D55E00", linestyle="solid", linewidth=2.5,
                   label="Positive MAG abundance association"),
        plt.Line2D([], [], color="#0072B2", linestyle="dashed", linewidth=2.5,
                   label="Negative MAG abundance association"),
    ]
    if overlay_mode in {"modules", "ecological"}:
        network_handles = [
            *[
                plt.Line2D(
                    [], [], marker="o", linestyle="",
                    markerfacecolor=color, markeredgecolor="black",
                    markersize=10, label=module,
                )
                for module, color in module_palette.items()
            ],
            plt.Line2D(
                [], [], marker="o", linestyle="",
                markerfacecolor="#D9D9D9", markeredgecolor="black",
                markersize=8, label="Other modules",
            ),
            *(
                [
                    plt.Line2D(
                        [], [], marker="o", linestyle="",
                        markerfacecolor="#808080", markeredgecolor="#A6A6A6",
                        markeredgewidth=0.5, markersize=17,
                        label="Anchor",
                    ),
                    plt.Line2D(
                        [], [], marker="o", linestyle="",
                        markerfacecolor="#808080", markeredgecolor="black",
                        markeredgewidth=2.4, markersize=17,
                        label="Anchor with paired MAG",
                    ),
                ]
                if overlay_mode == "ecological"
                else []
            ),
            plt.Line2D(
                [], [], color="#BDBDBD", linewidth=2,
                label="SPIEC-EASI edge",
            ),
        ]
    elif overlay_mode == "biochemical":
        observed_levels = {
            value
            for value in biochemical_associations.values()
            if value and "+" not in value
        }
        ordered_levels = [
            level for level in biochemical_order if level in observed_levels
        ] + sorted(observed_levels - set(biochemical_order))
        network_handles = [
            *[
                plt.Line2D(
                    [], [], marker="o", linestyle="",
                    markerfacecolor=color, markeredgecolor="#BDBDBD",
                    markeredgewidth=0.5, markersize=10,
                    label=f"Fill: {module}",
                )
                for module, color in module_palette.items()
            ],
            plt.Line2D(
                [], [], marker="o", linestyle="",
                markerfacecolor="#D9D9D9", markeredgecolor="#BDBDBD",
                markeredgewidth=0.5, markersize=10,
                label="Fill: other modules",
            ),
            *[
                plt.Line2D(
                    [], [], marker="o", linestyle="",
                    markerfacecolor="#D9D9D9",
                    markeredgecolor=biochemical_palette.get(level, "#4D4D4D"),
                    markeredgewidth=2.4, markersize=10, label=f"Rim: {level}",
                )
                for level in ordered_levels
            ],
            *(
                [
                    plt.Line2D(
                        [], [], marker="o", linestyle="",
                        markerfacecolor="#D9D9D9", markeredgecolor="black",
                        markeredgewidth=2.4, markersize=10,
                        label="Rim: multigroup indicator",
                    )
                ]
                if any("+" in value for value in biochemical_associations.values())
                else []
            ),
            plt.Line2D(
                [], [], marker="o", linestyle="",
                markerfacecolor="#D9D9D9", markeredgecolor="#BDBDBD",
                markeredgewidth=0.5, markersize=10,
                label="Rim: none",
            ),
        ]
    elif overlay_mode == "functional":
        network_handles = [
            *[
                plt.Line2D(
                    [], [], marker="o", linestyle="",
                    markerfacecolor=module_palette[module],
                    markeredgecolor="#A6A6A6", markeredgewidth=0.5,
                    markersize=10, label=module,
                )
                for module in module_levels
            ],
            plt.Line2D(
                [], [], marker="o", linestyle="",
                markerfacecolor="white",
                markeredgecolor="black", markeredgewidth=2.4,
                markersize=11, label="KEGG-contributing ASV",
            ),
            plt.Line2D(
                [], [], color="#D55E00", linestyle="solid", linewidth=2.5,
                label="Positive MAG association",
            ),
            plt.Line2D(
                [], [], color="#0072B2", linestyle="dashed", linewidth=2.5,
                label="Negative MAG association",
            ),
        ]
    elif overlay_mode == "titan":
        titan_missing_present = any(str(node) not in titan_responses for node in asv_nodes_all)
        network_handles = [
            plt.Line2D(
                [], [], marker="o", linestyle="",
                markerfacecolor=TITAN_DIRECTION_COLORS["z-"],
                markeredgecolor="black", markersize=10,
                label="z−: declining response",
            ),
            plt.Line2D(
                [], [], marker="o", linestyle="",
                markerfacecolor=TITAN_DIRECTION_COLORS["z+"],
                markeredgecolor="#E6E6E6", markeredgewidth=0.7, markersize=10,
                label="z+: increasing response",
            ),
            *(
                [
                    plt.Line2D(
                        [], [], marker="o", linestyle="",
                        markerfacecolor="#F2F2F2", markeredgecolor="#D9D9D9",
                        markersize=6, label="No TITAN response",
                    )
                ]
                if titan_missing_present else []
            ),
            *[
                plt.Line2D(
                    [], [], marker="o", linestyle="",
                    markerfacecolor="white", markeredgecolor="black",
                    markersize=np.sqrt(titan_marker_area(value)),
                    label=f"|z| = {value:g}",
                )
                for value in TITAN_Z_LEGEND_VALUES
            ],
            plt.Line2D([], [], color="#BDBDBD", linewidth=2,
                       label="SPIEC-EASI ASV association"),
            plt.Line2D([], [], color="#4D4D4D", linewidth=1.2,
                       label="ASV-to-MAG-subnetwork anchor"),
            plt.Line2D([], [], color="#D55E00", linestyle="solid", linewidth=2.5,
                       label="Positive MAG abundance association"),
            plt.Line2D([], [], color="#0072B2", linestyle="dashed", linewidth=2.5,
                       label="Negative MAG abundance association"),
        ]
    legend_title = (
        "Ecological modules and anchors"
        if overlay_mode == "ecological"
        else "Ecological modules"
        if overlay_mode == "modules"
        else biochemical_grouping
        if overlay_mode == "biochemical"
        else f"{functional_module_id}: {functional_module_name}"
        if overlay_mode == "functional"
        else f"TITAN response to {titan_variable}"
        if overlay_mode == "titan"
        else "Network features"
    )
    network_legend = legend_ax.legend(
        handles=network_handles, loc="upper left",
        bbox_to_anchor=(0.0, 1.0), frameon=False, title=legend_title,
        ncol=2 if overlay_mode in {"modules", "ecological", "functional", "biochemical", "titan"} else 1,
    )
    legend_ax.add_artist(network_legend)
    def asv_numeric_key(value: str) -> tuple[int, str]:
        match = re.search(r"(\d+)", value)
        return (int(match.group(1)) if match else 10**12, value)

    ordered_taxonomy_groups = sorted(
        taxonomy_legend_groups,
        key=lambda item: (module_sort_key(item[0]), asv_numeric_key(item[1])),
    )
    if overlay_mode not in {"modules", "ecological", "biochemical"}:
        taxonomy_children = []
        current_module = None
        for module, _asv_id, group_box, _row_count in ordered_taxonomy_groups:
            if module != current_module:
                taxonomy_children.append(TextArea(
                    module if module else "Unassigned module",
                    textprops={"fontsize": 17, "fontweight": "bold"},
                ))
                current_module = module
            taxonomy_children.append(group_box)
        taxonomy_column = VPacker(
            children=taxonomy_children,
            align="left",
            pad=0,
            sep=8,
        )
        taxonomy_key = VPacker(
            children=[
                TextArea("ASV and local MAG taxonomy", textprops={"fontsize": 18, "fontweight": "bold"}),
                taxonomy_column,
            ],
            align="left",
            pad=0,
            sep=5,
        )
        legend_ax.add_artist(AnnotationBbox(
            taxonomy_key,
            (
                0.0,
                max(0.35, 0.78 - 0.025 * len(network_handles))
                if overlay_mode in {"functional", "titan"}
                else 0.79,
            ),
            xycoords=legend_ax.transAxes,
            box_alignment=(0.0, 1.0),
            frameon=False,
            pad=0,
            zorder=10,
        ))
    if overlay_mode in {"modules", "ecological", "biochemical"}:
        x_min = min(coords[0] for coords in asv_pos.values())
        x_max = max(coords[0] for coords in asv_pos.values())
        y_min = min(coords[1] for coords in asv_pos.values())
        y_max = max(coords[1] for coords in asv_pos.values())
    else:
        x_min = min(
            min(coords[0] for coords in asv_pos.values()),
            min(details["center"][0] - details["width"] / 2 for details in callouts.values()),
        )
        x_max = max(
            max(coords[0] for coords in asv_pos.values()),
            max(details["center"][0] + details["width"] / 2 for details in callouts.values()),
        )
        y_min = min(
            min(coords[1] for coords in asv_pos.values()),
            min(details["center"][1] - details["height"] / 2 for details in callouts.values()),
        )
        y_max = max(
            max(coords[1] for coords in asv_pos.values()),
            max(details["center"][1] + details["height"] / 2 for details in callouts.values()),
        )
    x_padding = max(0.22, 0.08 * (x_max - x_min))
    y_padding = max(0.22, 0.08 * (y_max - y_min))
    ax.set_xlim(x_min - x_padding, x_max + x_padding)
    ax.set_ylim(y_min - y_padding, y_max + y_padding)
    ax.set_aspect("equal", adjustable="box")
    ax.set_axis_off()
    if overlay_mode == "modules":
        output_stem = f"{prefix}_ecological_modules"
    elif overlay_mode == "ecological":
        output_stem = f"{prefix}_ecological_modules_and_anchors"
    elif overlay_mode == "functional":
        output_stem = f"{prefix}_kegg_{functional_module_id}_contributors"
    elif overlay_mode == "biochemical":
        safe_grouping = re.sub(r"[^A-Za-z0-9_.-]+", "_", biochemical_grouping)
        output_stem = f"{prefix}_ecological_modules_by_{safe_grouping}"
    elif overlay_mode == "titan":
        safe_variable = re.sub(r"[^A-Za-z0-9_.-]+", "_", titan_variable).strip("_").lower()
        output_stem = f"{prefix}_titan_{safe_variable}_heterogeneous_network"
    else:
        output_stem = f"{prefix}_heterogeneous_network"
    for fmt in formats:
        fig.savefig(
            outdir / "network" / f"{output_stem}.{fmt}",
            dpi=300,
            bbox_inches="tight",
        )
    plt.close(fig)


def annotate_heterogeneous_functions(
    graph: nx.Graph,
    module_table: pd.DataFrame,
    selected_modules: list[str],
) -> None:
    """Attach selected functional-module coverage to MAG nodes in place."""
    if module_table.empty:
        return
    work = module_table.loc[
        module_table["module_id"].astype(str).isin(selected_modules)
    ].copy()
    work["fraction_covered"] = pd.to_numeric(
        work["fraction_covered"], errors="coerce"
    ).fillna(0.0)
    by_match = {}
    for row in work.itertuples():
        match_id = str(row.mag_match_id)
        module_id = str(row.module_id)
        by_match.setdefault(match_id, {})
        by_match[match_id][module_id] = max(
            by_match[match_id].get(module_id, 0.0),
            float(row.fraction_covered),
        )
    for node, attrs in graph.nodes(data=True):
        if attrs.get("node_type") != "MAG":
            continue
        match_ids = [
            value for value in str(attrs.get("mag_match_id", "")).split("|")
            if value
        ]
        coverage = {
            module_id: max(
                (by_match.get(match_id, {}).get(module_id, 0.0) for match_id in match_ids),
                default=0.0,
            )
            for module_id in selected_modules
        }
        attrs["selected_cycle_module_coverage"] = "|".join(
            f"{module_id}={coverage[module_id]:.6g}" for module_id in selected_modules
        )
        attrs["selected_cycle_modules_present"] = "|".join(
            module_id for module_id, value in coverage.items() if value > 0
        )


def load_biochemical_indicator_associations(
    isa_dir: Path,
    grouping: str,
    q_threshold: float,
) -> dict[str, str]:
    """Return significant ISA group memberships for each ASV."""
    path = isa_dir / f"{grouping}_indicator_species_summary.tsv"
    table = read_table(path)
    if table.empty:
        return {}
    asv_col = (
        "ASV" if "ASV" in table
        else "ASV_ID" if "ASV_ID" in table
        else table.columns[0]
    )
    q_col = (
        "q.value" if "q.value" in table
        else "p.value" if "p.value" in table
        else None
    )
    if q_col is None:
        return {}
    sign_cols = [column for column in table if str(column).startswith("s.")]
    associations: dict[str, set[str]] = defaultdict(set)
    significant = table.loc[
        pd.to_numeric(table[q_col], errors="coerce").le(q_threshold)
    ]
    for _, row in significant.iterrows():
        levels = [
            str(column)[2:]
            for column in sign_cols
            if pd.to_numeric(row.get(column), errors="coerce") == 1
        ]
        if not levels and "index" in table and pd.notna(row.get("index")):
            levels = [str(row["index"])]
        associations[str(row[asv_col]).strip()].update(levels)
    return {
        asv_id: "+".join(sorted(levels))
        for asv_id, levels in associations.items()
        if levels
    }


def fixed_heterogeneous_layout(graph: nx.Graph) -> dict[str, np.ndarray]:
    """One deterministic layout reused by every evidence overlay for a graph."""
    return {
        str(node): np.asarray(coords, dtype=float)
        for node, coords in nx.spring_layout(graph, seed=42, weight="weight").items()
    }


def node_cycle_coverage(attrs: dict, module_id: str) -> float:
    values = {
        item.split("=", 1)[0]: float(item.split("=", 1)[1])
        for item in str(attrs.get("selected_cycle_module_coverage", "")).split("|")
        if "=" in item
    }
    return float(values.get(module_id, 0.0))


def functional_contributor_table(
    graph: nx.Graph,
    module_id: str,
    minimum_fraction: float,
    modality: str,
) -> pd.DataFrame:
    rows = []
    for left, right, attrs in graph.edges(data=True):
        if attrs.get("edge_type") != "sequence_match":
            continue
        asv = left if graph.nodes[left].get("node_type") == "ASV" else right
        mag = right if asv == left else left
        coverage = node_cycle_coverage(graph.nodes[mag], module_id)
        if coverage < minimum_fraction:
            continue
        rows.append({
            "data_modality": modality,
            "functional_module_id": module_id,
            "ASV_ID": str(asv),
            "ecological_module": graph.nodes[asv].get("module_label", ""),
            "is_ecological_anchor": bool(
                graph.nodes[asv].get("is_ecological_anchor", False)
            ),
            "genome_id": str(mag),
            "mag_match_id": graph.nodes[mag].get("mag_match_id", ""),
            "fraction_covered": coverage,
            "link_pident": attrs.get("identity", ""),
            "link_qcov": attrs.get("coverage", ""),
            "mapping_class": attrs.get("mapping_class", ""),
        })
    return pd.DataFrame(rows)


def plot_heterogeneous_evidence_overlay(
    graph: nx.Graph,
    positions: dict[str, np.ndarray],
    outdir: Path,
    prefix: str,
    formats: list[str],
    mode: str,
    functional_module_id: str | None = None,
    functional_module_name: str = "",
    minimum_fraction: float = 0.5,
) -> None:
    """Plot module/anchor or KEGG evidence on the unchanged heterogeneous graph."""
    if graph.number_of_nodes() == 0:
        return
    asvs = [n for n, d in graph.nodes(data=True) if d.get("node_type") == "ASV"]
    mags = [n for n, d in graph.nodes(data=True) if d.get("node_type") == "MAG"]
    module_levels = sorted({
        str(graph.nodes[node].get("module_label"))
        for node in asvs
        if str(graph.nodes[node].get("module_label", "")).lower() not in {"", "nan", "none"}
    }, key=module_sort_key)
    palette = {
        module: matplotlib.colors.to_hex(plt.get_cmap("tab20")(index % 20))
        for index, module in enumerate(module_levels)
    }
    sequence_edges = [
        (u, v) for u, v, data in graph.edges(data=True)
        if data.get("edge_type") == "sequence_match"
    ]
    asv_edges = [
        (u, v) for u, v, data in graph.edges(data=True)
        if data.get("edge_type") == "asv_association"
    ]
    mag_edges = [
        (u, v, data) for u, v, data in graph.edges(data=True)
        if data.get("edge_type") == "mag_abundance_association"
    ]
    positive_mags: set[str] = set()
    if functional_module_id:
        for node in mags:
            if node_cycle_coverage(graph.nodes[node], functional_module_id) >= minimum_fraction:
                positive_mags.add(str(node))
    positive_asvs = {
        str(u if graph.nodes[u].get("node_type") == "ASV" else v)
        for u, v in sequence_edges
        if str(u) in positive_mags or str(v) in positive_mags
    }

    fig, ax = plt.subplots(figsize=(18, 16))
    nx.draw_networkx_edges(
        graph, positions, edgelist=asv_edges, edge_color="#D9D9D9",
        width=0.7, alpha=0.45, ax=ax,
    )
    for sign, style, color in (
        ("positive", "solid", "#D55E00"),
        ("negative", "dashed", "#0072B2"),
    ):
        edges = [
            (u, v) for u, v, data in mag_edges
            if (float(data.get("spearman_rho", 0)) >= 0) == (sign == "positive")
        ]
        if edges:
            nx.draw_networkx_edges(
                graph, positions, edgelist=edges, edge_color=color,
                style=style, width=1.2, alpha=0.45, ax=ax,
            )
    if mode == "ecological":
        sequence_colors = [
            palette.get(
                str(graph.nodes[u if graph.nodes[u].get("node_type") == "ASV" else v].get("module_label")),
                "#808080",
            )
            for u, v in sequence_edges
        ]
    else:
        sequence_colors = [
            "#D55E00" if str(u) in positive_mags or str(v) in positive_mags else "#BDBDBD"
            for u, v in sequence_edges
        ]
    nx.draw_networkx_edges(
        graph, positions, edgelist=sequence_edges, edge_color=sequence_colors,
        width=1.4, alpha=0.8, ax=ax,
    )

    if mode == "ecological":
        asv_colors = [
            palette.get(str(graph.nodes[node].get("module_label")), "#E6E6E6")
            for node in asvs
        ]
        asv_edges_colors = [
            "#D55E00" if bool(graph.nodes[node].get("is_ecological_anchor")) else "black"
            for node in asvs
        ]
        asv_widths = [
            2.4 if bool(graph.nodes[node].get("is_ecological_anchor")) else 0.5
            for node in asvs
        ]
        mag_colors = []
        for node in mags:
            linked = [
                neighbor for neighbor in graph.neighbors(node)
                if graph.nodes[neighbor].get("node_type") == "ASV"
            ]
            module = str(graph.nodes[linked[0]].get("module_label")) if linked else ""
            mag_colors.append(palette.get(module, "#E6E6E6"))
    else:
        asv_colors = ["#4D4D4D" if str(node) in positive_asvs else "#F2F2F2" for node in asvs]
        asv_edges_colors = ["black"] * len(asvs)
        asv_widths = [0.7] * len(asvs)
        mag_colors = ["#D55E00" if str(node) in positive_mags else "white" for node in mags]
    if mode == "ecological":
        asv_sizes = [
            220 if bool(graph.nodes[n].get("is_ecological_anchor")) else 85
            for n in asvs
        ]
    else:
        linked_positive_counts = Counter()
        for left, right in sequence_edges:
            if str(left) in positive_mags or str(right) in positive_mags:
                asv = left if graph.nodes[left].get("node_type") == "ASV" else right
                linked_positive_counts[str(asv)] += 1
        asv_sizes = [
            100 + 35 * linked_positive_counts.get(str(node), 0)
            if str(node) in positive_asvs else 55
            for node in asvs
        ]
    nx.draw_networkx_nodes(
        graph, positions, nodelist=asvs, node_shape="o", node_color=asv_colors,
        node_size=asv_sizes,
        edgecolors=asv_edges_colors, linewidths=asv_widths, ax=ax,
    )
    if mode == "ecological":
        mag_sizes = [65] * len(mags)
    else:
        mag_sizes = [
            55 + 150 * node_cycle_coverage(graph.nodes[node], functional_module_id)
            if str(node) in positive_mags else 35
            for node in mags
        ]
    nx.draw_networkx_nodes(
        graph, positions, nodelist=mags, node_shape="o", node_color=mag_colors,
        node_size=mag_sizes, edgecolors="black", linewidths=0.7, ax=ax,
    )
    if mode == "ecological":
        label_nodes = [
            node for node in asvs
            if bool(graph.nodes[node].get("is_ecological_anchor"))
            and any(graph.nodes[neighbor].get("node_type") == "MAG" for neighbor in graph.neighbors(node))
        ]
    else:
        qualifying_links = Counter()
        for left, right in sequence_edges:
            if str(left) in positive_mags or str(right) in positive_mags:
                asv = left if graph.nodes[left].get("node_type") == "ASV" else right
                qualifying_links[str(asv)] += 1
        label_nodes = [
            node for node, _count in
            sorted(qualifying_links.items(), key=lambda item: (-item[1], item[0]))[:8]
        ]
    if label_nodes:
        all_xy = np.asarray([positions[node] for node in graph.nodes])
        x_min, y_min = all_xy.min(axis=0)
        x_max, y_max = all_xy.max(axis=0)
        x_span = max(float(x_max - x_min), 1.0)
        y_span = max(float(y_max - y_min), 1.0)
        center_x = float((x_min + x_max) / 2)
        sides = {
            "left": sorted(
                [node for node in label_nodes if positions[node][0] < center_x],
                key=lambda node: positions[node][1],
            ),
            "right": sorted(
                [node for node in label_nodes if positions[node][0] >= center_x],
                key=lambda node: positions[node][1],
            ),
        }
        for side, nodes_on_side in sides.items():
            if not nodes_on_side:
                continue
            label_ys = np.linspace(
                y_min + 0.08 * y_span,
                y_max - 0.08 * y_span,
                len(nodes_on_side),
            )
            label_x = x_min - 0.04 * x_span if side == "left" else x_max + 0.04 * x_span
            for node, label_y in zip(nodes_on_side, label_ys):
                ax.annotate(
                    str(node),
                    xy=positions[node],
                    xytext=(label_x, float(label_y)),
                    textcoords="data",
                    ha="right" if side == "left" else "left",
                    va="center",
                    fontsize=22,
                    bbox={"boxstyle": "round,pad=0.15", "facecolor": "white", "edgecolor": "none", "alpha": 0.9},
                    arrowprops={"arrowstyle": "-", "color": "#4D4D4D", "lw": 0.7},
                    zorder=20,
                )
    handles = []
    if mode == "ecological":
        handles.extend(
            Line2D([], [], marker="o", linestyle="", markerfacecolor=color,
                   markeredgecolor="black", markersize=9, label=module)
            for module, color in palette.items()
        )
        handles.append(
            Line2D([], [], marker="o", linestyle="", markerfacecolor="white",
                   markeredgecolor="#D55E00", markeredgewidth=2.4,
                   markersize=11, label="Ecological-module anchor")
        )
        output_stem = f"{prefix}_ecological_modules_and_anchors"
    else:
        handles.extend([
            Line2D([], [], marker="o", linestyle="", markerfacecolor="#D55E00",
                   markeredgecolor="black", markersize=9, label="MAG meeting coverage threshold"),
            Line2D([], [], marker="o", linestyle="", markerfacecolor="#4D4D4D",
                   markeredgecolor="black", markersize=9, label="ASV linked to qualifying MAG"),
        ])
        output_stem = f"{prefix}_kegg_{functional_module_id}_contributors"
    ax.legend(handles=handles, frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    ax.set_axis_off()
    ax.set_title(
        "Ecological modules and anchors"
        if mode == "ecological"
        else f"{functional_module_id}: {functional_module_name}",
        loc="left",
    )
    fig.tight_layout()
    for fmt in formats:
        fig.savefig(
            outdir / "network" / f"{output_stem}.{fmt}",
            dpi=300,
            bbox_inches="tight",
        )
    plt.close(fig)


def plot_summary(summary: pd.DataFrame, outdir: Path) -> None:
    plots_dir = outdir / "qc"
    plots_dir.mkdir(parents=True, exist_ok=True)
    if "mapping_class" in summary.columns:
        fig, ax = plt.subplots(figsize=(7, 4))
        sns.barplot(data=summary, x="mapping_class", y="n_asvs", color="#4C78A8", ax=ax)
        ax.set_xlabel("Mapping class")
        ax.set_ylabel("ASVs")
        ax.set_title("ASV-MAG network mapping classes")
        ax.tick_params(axis="x", rotation=25)
        fig.tight_layout()
        fig.savefig(plots_dir / "asv_mag_network_mapping_classes.png", dpi=300)
        fig.savefig(plots_dir / "asv_mag_network_mapping_classes.svg")
        plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph", required=True)
    parser.add_argument("--node-features", required=True)
    parser.add_argument(
        "--ecological-modules",
        help="Network-module membership table used for exact-graph module and anchor overlays.",
    )
    parser.add_argument(
        "--ecological-module-selection",
        help="Module-quality table; rows with is_best=true define displayed modules.",
    )
    parser.add_argument(
        "--display-all-ecological-modules",
        action="store_true",
        help=(
            "Display every inferred ecological module with its own label and color, "
            "irrespective of the visualization-focused is_best ranking."
        ),
    )
    parser.add_argument(
        "--anchor-top-n", type=int, default=1,
        help="Within-module top rank retained for degree, eigenvector centrality, or betweenness.",
    )
    parser.add_argument("--asv-mag-pairing", required=True)
    parser.add_argument("--taxonomy", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--prefix", default="asv_mag_network")
    parser.add_argument("--genome-summary")
    parser.add_argument("--reference-catalog")
    parser.add_argument("--asv-counts")
    parser.add_argument("--mag-abundance")
    parser.add_argument(
        "--cruise-metadata",
        help="Optional table containing Cruise and Date or Year/Month/Day columns for monthly MAG profiles.",
    )
    parser.add_argument("--mag-id-mode", default="exact", choices=["exact", "suffix_after_double_underscore"])
    parser.add_argument("--mag-abundance-format", default="auto", choices=["auto", "long", "wide"])
    parser.add_argument("--mag-abundance-genome-col", default="genome_id")
    parser.add_argument("--mag-abundance-sample-col", default="sample_id")
    parser.add_argument("--mag-abundance-value-col", default="read_count")
    parser.add_argument("--mag-abundance-seqkit")
    parser.add_argument("--mag-abundance-seqkit-file-col", default="file")
    parser.add_argument("--mag-abundance-seqkit-count-col", default="num_seqs")
    parser.add_argument(
        "--mag-transcript-abundance",
        help="Optional MAG metatranscriptome recruitment table used for a parallel set of lineage profiles.",
    )
    parser.add_argument(
        "--mag-transcript-abundance-format",
        default="auto",
        choices=["auto", "long", "wide"],
    )
    parser.add_argument("--mag-transcript-abundance-genome-col", default="genome_id")
    parser.add_argument("--mag-transcript-abundance-sample-col", default="sample_id")
    parser.add_argument("--mag-transcript-abundance-value-col", default="read_count")
    parser.add_argument("--mag-transcript-abundance-seqkit")
    parser.add_argument(
        "--mag-transcript-abundance-seqkit-file-col", default="file"
    )
    parser.add_argument(
        "--mag-transcript-abundance-seqkit-count-col", default="num_seqs"
    )
    parser.add_argument(
        "--mag-transcript-abundance-normalization",
        default="auto",
        choices=[
            "auto", "input_fragment_fpm", "provided_fpkm", "provided_tpm",
            "median_ratio", "raw", "relative",
        ],
    )
    parser.add_argument("--min-shared-samples", type=int, default=5)
    parser.add_argument(
        "--mag-knn",
        type=int,
        default=3,
        help="Strongest positive and strongest negative abundance-correlation neighbors retained per MAG and sign.",
    )
    parser.add_argument("--abundance-transform", default="log1p", choices=["none", "log1p"])
    parser.add_argument(
        "--mag-abundance-normalization",
        default="auto",
        choices=[
            "auto", "input_fragment_fpm", "provided_fpkm", "provided_tpm",
            "median_ratio", "raw", "relative",
        ],
    )
    parser.add_argument("--functional-annotation", action="append", default=[])
    parser.add_argument("--functional-module-min-fraction", type=float, default=0.5)
    parser.add_argument("--ambiguity-target-modules", default="M00529,M00973,M00595,M00984")
    parser.add_argument("--ambiguity-formats", default="pdf,png,svg")
    parser.add_argument("--isa-dir")
    parser.add_argument(
        "--biochemical-groupings",
        default="",
        help="Comma-separated ISA groupings rendered as rims over ecological-module nodes.",
    )
    parser.add_argument("--group-palette-map-json", default="{}")
    parser.add_argument("--group-order-map-json", default="{}")
    parser.add_argument("--isa-q-threshold", type=float, default=0.05)
    parser.add_argument("--sample-metadata")
    parser.add_argument("--grouping-diagnostic-sample-col", default="sample")
    parser.add_argument("--grouping-diagnostic-cruise-col", default="Cruise")
    parser.add_argument(
        "--grouping-diagnostic-groups",
        default="o2_subcompartment_final,cruise_group",
    )
    parser.add_argument(
        "--grouping-diagnostic-cruise-level-groups",
        default="Season,cruise_group",
        help="Comma-separated grouping columns evaluated after aggregation to one profile per cruise.",
    )
    parser.add_argument("--grouping-diagnostic-permutations", type=int, default=999)
    parser.add_argument("--font-family", default="Times New Roman")
    parser.add_argument("--asv-taxonomy-source", default="ncbi")
    parser.add_argument(
        "--asv-taxonomy-min-confidence",
        type=float,
        help="Minimum inclusive ASV taxonomy confidence on a 0-1 scale.",
    )
    parser.add_argument("--mag-taxonomy-source", default="gtdb")
    parser.add_argument("--min-pident", type=float, default=99.5)
    parser.add_argument("--min-qcov", type=float, default=100.0)
    args = parser.parse_args()
    if args.mag_knn < 1:
        parser.error("--mag-knn must be at least 1")
    if (
        args.asv_taxonomy_min_confidence is not None
        and not 0 <= args.asv_taxonomy_min_confidence <= 1
    ):
        parser.error("--asv-taxonomy-min-confidence must be between 0 and 1")
    try:
        biochemical_palettes = json.loads(args.group_palette_map_json or "{}")
        biochemical_orders = json.loads(args.group_order_map_json or "{}")
    except json.JSONDecodeError as exc:
        parser.error(f"Invalid biochemical group palette/order JSON: {exc}")

    plt.rcParams.update({
        "font.family": args.font_family,
        "font.size": 22,
        "axes.titlesize": 24,
        "axes.labelsize": 22,
        "xtick.labelsize": 18,
        "ytick.labelsize": 18,
        "legend.fontsize": 18,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    outdir = Path(args.outdir)
    for subdir in ["validation", "filtering", "spieceasi", "rrna", "mapping", "network", "abundance", "functional", "qc", "group_diagnostics"]:
        (outdir / subdir).mkdir(parents=True, exist_ok=True)
    for obsolete in (outdir / "group_diagnostics").glob(f"{args.prefix}_*"):
        if obsolete.is_file() or obsolete.is_symlink():
            obsolete.unlink()
    # Remove only obsolete, modality-unlabeled outputs from earlier versions.
    # All evidence is regenerated below under explicit metagenome/metatranscriptome names.
    legacy_patterns = {
        "network": [
            f"{args.prefix}_heterogeneous*",
            f"{args.prefix}_metagenome_ecological_modules*",
            f"{args.prefix}_metatranscriptome_ecological_modules*",
        ],
        "abundance": [
            f"{args.prefix}_ASV*_species_recruitment*",
            f"{args.prefix}_ASV*_mag_recruitment*",
            f"{args.prefix}_ASV*_lineage_recruitment*",
            f"{args.prefix}_mag_mag_abundance_associations.tsv",
            f"{args.prefix}_asv_mag_abundance_agreement.tsv",
            f"{args.prefix}_deseq2_size_factors.tsv",
            f"{args.prefix}_shared_metag_metat_lineage_y_axis_ranges.tsv",
        ],
        "qc": [
            f"{args.prefix}_abundance_sample_harmonization.tsv",
            f"{args.prefix}_heterogeneous_excluded_links.tsv",
        ],
    }
    for subdir, patterns in legacy_patterns.items():
        for pattern in patterns:
            for legacy_path in (outdir / subdir).glob(pattern):
                if legacy_path.is_file() or legacy_path.is_symlink():
                    legacy_path.unlink()

    graph_path = Path(args.graph)
    graph = nx.read_graphml(graph_path)
    graph = nx.relabel_nodes(graph, {node: str(data.get("name", node)) for node, data in graph.nodes(data=True)})
    graph.remove_nodes_from(list(nx.isolates(graph)))
    for u, v, attrs in graph.edges(data=True):
        attrs["edge_type"] = "asv_association"

    node_features = read_node_features(Path(args.node_features))
    ecological_modules = (
        read_table(Path(args.ecological_modules))
        if args.ecological_modules
        else pd.DataFrame()
    )
    ecological_module_selection = (
        read_table(Path(args.ecological_module_selection))
        if args.ecological_module_selection
        else pd.DataFrame()
    )
    node_features = annotate_ecological_modules(
        node_features,
        ecological_modules,
        args.anchor_top_n,
        ecological_module_selection,
        args.display_all_ecological_modules,
    )
    asv_tax = load_asv_taxonomy(
        Path(args.taxonomy),
        min_confidence=args.asv_taxonomy_min_confidence,
    )
    pairing = normalize_pairing(
        read_table(Path(args.asv_mag_pairing)),
        asv_tax,
        args.asv_taxonomy_source,
        args.mag_taxonomy_source,
        args.mag_id_mode,
    )
    if pairing.empty:
        die("ASV-MAG pairing table is empty; cannot build ASV_MAG_NETWORK outputs.")
    pairing["passes_paper_thresholds"] = (
        pd.to_numeric(pairing.get("link_pident"), errors="coerce").ge(args.min_pident)
        & pd.to_numeric(pairing.get("link_qcov"), errors="coerce").ge(args.min_qcov)
    )
    pairing["accepted_paper_pair"] = pairing["accepted_paper_pair"] & pairing["passes_paper_thresholds"]
    pairing["analysis_eligible_pair"] = pairing["passes_paper_thresholds"] & ~pairing["taxonomy_validation_status"].eq("taxonomy_rejected") & pairing["mapping_class"].isin(["unique_match", "ambiguous_match"])
    eligible_counts = pairing.loc[pairing["analysis_eligible_pair"]].groupby("ASV_ID")["mag_match_id"].transform("nunique")
    pairing["eligible_mag_count"] = eligible_counts.reindex(pairing.index).fillna(0).astype(int)
    pairing["asv_link_weight"] = np.where(
        pairing["analysis_eligible_pair"] & pairing["eligible_mag_count"].gt(0),
        1.0 / pairing["eligible_mag_count"],
        0.0,
    )

    functional = summarize_functional_annotations(
        [Path(p) for p in args.functional_annotation],
        args.functional_module_min_fraction,
        args.mag_id_mode,
    )
    module_table = pd.DataFrame()
    if not functional.empty:
        write_table(functional, outdir / "functional" / f"{args.prefix}_functional_summary.tsv")
        module_table = functional.attrs.get("module_table")
        if isinstance(module_table, pd.DataFrame) and not module_table.empty:
            accepted_genomes = set(pairing.loc[pairing["analysis_eligible_pair"], "mag_match_id"].dropna().astype(str))
            if accepted_genomes:
                accepted_module_table = module_table.loc[module_table["mag_match_id"].astype(str).isin(accepted_genomes)].copy()
            else:
                accepted_module_table = module_table.iloc[0:0].copy()
            write_table(accepted_module_table, outdir / "functional" / f"{args.prefix}_functional_modules.tsv")
    requested_target_modules = [
        value.strip()
        for value in args.ambiguity_target_modules.split(",")
        if value.strip()
    ]
    cycle_module_audit = pd.DataFrame()
    if any(value.lower() == "auto_n_s" for value in requested_target_modules):
        target_modules, cycle_module_audit = select_well_covered_cycle_modules(
            module_table,
            args.functional_module_min_fraction,
        )
        write_table(
            cycle_module_audit,
            outdir / "functional" / f"{args.prefix}_nitrogen_sulfur_module_selection.tsv",
        )
        if not target_modules:
            die(
                "Automatic N/S KEGG-module selection found no module with at least "
                f"one MAG at fraction_covered >= {args.functional_module_min_fraction}."
            )
    else:
        target_modules = requested_target_modules

    sample_metadata = (
        read_table(Path(args.sample_metadata))
        if args.sample_metadata
        else None
    )
    agreement, abundance_matrices = abundance_agreement(
        pairing,
        Path(args.mag_abundance) if args.mag_abundance else None,
        Path(args.asv_counts) if args.asv_counts else None,
        args.mag_abundance_format,
        args.mag_abundance_genome_col,
        args.mag_abundance_sample_col,
        args.mag_abundance_value_col,
        args.mag_id_mode,
        args.min_shared_samples,
        args.abundance_transform,
        args.mag_abundance_normalization,
        seqkit_path=(
            Path(args.mag_abundance_seqkit)
            if args.mag_abundance_seqkit else None
        ),
        seqkit_file_col=args.mag_abundance_seqkit_file_col,
        seqkit_count_col=args.mag_abundance_seqkit_count_col,
        sample_metadata=sample_metadata,
        target_sample_col=args.grouping_diagnostic_sample_col,
        include_ambiguous=True,
    )
    if not agreement.empty:
        if "analysis_eligible_pair" in pairing.columns:
            eligible_pair_mask = pairing["analysis_eligible_pair"].fillna(False).map(
                lambda value: bool(value)
                if isinstance(value, (bool, np.bool_))
                else str(value).strip().lower() in {"1", "true", "yes"}
            )
        else:
            eligible_pair_mask = pairing["accepted_paper_pair"].fillna(False).map(
                lambda value: bool(value)
                if isinstance(value, (bool, np.bool_))
                else str(value).strip().lower() in {"1", "true", "yes"}
            )
        eligible_pairs = pairing.loc[
            eligible_pair_mask,
            ["ASV_ID", "genome_id", "mag_match_id"],
        ].dropna(subset=["ASV_ID", "genome_id", "mag_match_id"]).drop_duplicates()
        eligible_genome_counts = eligible_pairs.groupby("ASV_ID")["mag_match_id"].nunique()
        single_eligible_pairs = eligible_pairs.loc[
            eligible_pairs["ASV_ID"].isin(
                eligible_genome_counts.loc[eligible_genome_counts.eq(1)].index
            )
        ]
        single_eligible_agreement = agreement.merge(
            single_eligible_pairs,
            on=["ASV_ID", "genome_id", "mag_match_id"],
            how="inner",
        )
        write_table(
            single_eligible_agreement,
            outdir / "abundance" / f"{args.prefix}_metagenome_asv_mag_abundance_agreement.tsv",
        )
        pairing = pairing.merge(agreement, on=["ASV_ID", "genome_id", "mag_match_id"], how="left")
    sample_audit = abundance_matrices.get("sample_audit", {})
    if sample_audit:
        write_table(
            pd.DataFrame([sample_audit]),
            outdir / "qc" / f"{args.prefix}_metagenome_abundance_sample_harmonization.tsv",
        )
    metagenome_seqkit_audit = abundance_matrices.get("seqkit_audit")
    if (
        isinstance(metagenome_seqkit_audit, pd.DataFrame)
        and not metagenome_seqkit_audit.empty
    ):
        write_table(
            metagenome_seqkit_audit,
            outdir / "qc" / f"{args.prefix}_metagenome_seqkit_library_audit.tsv",
        )
    if args.mag_abundance:
        expected_links = {
            (str(row.ASV_ID), str(row.genome_id), str(row.mag_match_id))
            for row in pairing.loc[
                pairing["analysis_eligible_pair"]
                & pairing["ASV_ID"].astype(str).isin(graph.nodes)
            ].itertuples()
        }
        observed_links = {
            (str(row.ASV_ID), str(row.genome_id), str(row.mag_match_id))
            for row in agreement.loc[agreement["spearman_rho"].notna()].itertuples()
        } if not agreement.empty else set()
        missing_links = expected_links - observed_links
        if missing_links:
            preview = ", ".join("/".join(values) for values in sorted(missing_links)[:5])
            die(
                f"MAG abundance was configured, but {len(missing_links)} eligible ASV-MAG links "
                f"lack a finite abundance correlation. Examples: {preview}"
            )

    cruise_metadata = (
        read_table(Path(args.cruise_metadata))
        if args.cruise_metadata
        else None
    )
    transcript_abundance_matrices: dict[str, object] = {}
    if args.mag_transcript_abundance:
        transcript_raw = load_mag_abundance(
            Path(args.mag_transcript_abundance),
            args.mag_transcript_abundance_format,
            args.mag_transcript_abundance_genome_col,
            args.mag_transcript_abundance_sample_col,
            args.mag_transcript_abundance_value_col,
            args.mag_id_mode,
        )
        (
            transcript_normalized,
            transcript_size_factors,
            transcript_normalization_used,
            transcript_seqkit_audit,
        ) = normalize_mag_recruitment(
            transcript_raw,
            method=args.mag_transcript_abundance_normalization,
            seqkit_path=(
                Path(args.mag_transcript_abundance_seqkit)
                if args.mag_transcript_abundance_seqkit else None
            ),
            seqkit_file_col=args.mag_transcript_abundance_seqkit_file_col,
            seqkit_count_col=args.mag_transcript_abundance_seqkit_count_col,
        )
        if not transcript_seqkit_audit.empty:
            write_table(
                transcript_seqkit_audit,
                outdir / "qc"
                / f"{args.prefix}_metatranscriptome_seqkit_library_audit.tsv",
            )
        transcript_transformed = (
            np.log1p(transcript_normalized)
            if args.abundance_transform == "log1p"
            else transcript_normalized.copy()
        )
        transcript_samples = list(transcript_transformed.columns)
        transcript_abundance_matrices = {
            "mag": standardize_profiles(transcript_transformed),
            "asv": abundance_matrices.get("asv", pd.DataFrame()),
            "asv_all_relative": abundance_matrices.get(
                "asv_all_relative", pd.DataFrame()
            ),
            "samples": transcript_samples,
            "mag_all_raw": transcript_raw,
            "mag_all_normalized": transcript_normalized,
            "mag_all_size_factors": transcript_size_factors,
            "mag_all_normalization": transcript_normalization_used,
            "mag_all_transformed": transcript_transformed,
        }
    if sample_metadata is not None:
        diagnostic_modalities: dict[str, pd.DataFrame] = {}
        asv_features = abundance_matrices.get("asv_all_relative")
        metag_features = abundance_matrices.get("mag_all_normalized")
        metat_features = transcript_abundance_matrices.get("mag_all_normalized")
        base_modalities = {
            "asv": asv_features,
            "metagenome": metag_features,
            "metatranscriptome": metat_features,
        }
        for modality, features in base_modalities.items():
            if isinstance(features, pd.DataFrame) and not features.empty:
                current = features.copy()
                current.index = [
                    f"{modality}::{feature}" for feature in current.index
                ]
                diagnostic_modalities[modality] = current
        joint_specs = [
            ("joint_asv_metagenome", ["asv", "metagenome"]),
            ("joint_asv_metatranscriptome", ["asv", "metatranscriptome"]),
            (
                "joint_asv_metagenome_metatranscriptome",
                ["asv", "metagenome", "metatranscriptome"],
            ),
        ]
        for joint_name, components in joint_specs:
            if not all(component in diagnostic_modalities for component in components):
                continue
            shared_samples = sorted(set.intersection(*[
                set(diagnostic_modalities[component].columns)
                for component in components
            ]))
            if len(shared_samples) < 4:
                continue
            diagnostic_modalities[joint_name] = pd.concat([
                diagnostic_modalities[component].loc[:, shared_samples]
                for component in components
            ])
        run_mag_grouping_diagnostics(
            diagnostic_modalities,
            sample_metadata,
            outdir,
            args.prefix,
            args.grouping_diagnostic_sample_col,
            args.grouping_diagnostic_cruise_col,
            [
                value.strip()
                for value in args.grouping_diagnostic_groups.split(",")
                if value.strip()
            ],
            [
                value.strip()
                for value in args.grouping_diagnostic_cruise_level_groups.split(",")
                if value.strip()
            ],
            args.grouping_diagnostic_permutations,
            42,
            [
                value.strip().lower()
                for value in args.ambiguity_formats.split(",")
                if value.strip()
            ],
        )
    ambiguity_formats = [
        value.strip().lower()
        for value in args.ambiguity_formats.split(",")
        if value.strip()
    ]
    ambiguity_outputs(
        pairing,
        agreement,
        abundance_matrices,
        module_table,
        outdir,
        f"{args.prefix}_metagenome",
        target_modules,
        ambiguity_formats,
        cruise_metadata,
    )
    if transcript_abundance_matrices:
        ambiguity_outputs(
            pairing,
            pd.DataFrame(),
            transcript_abundance_matrices,
            pd.DataFrame(),
            outdir,
            f"{args.prefix}_metatranscriptome",
            target_modules,
            ambiguity_formats,
            cruise_metadata,
            profile_only=True,
        )

    write_table(pairing, outdir / "validation" / f"{args.prefix}_taxonomy_validation.tsv")
    write_table(pairing.loc[pairing["mapping_class"].eq("ambiguous_match")], outdir / "mapping" / f"{args.prefix}_ambiguous_mappings.tsv")
    write_table(pairing.loc[pairing["accepted_paper_pair"]], outdir / "mapping" / f"{args.prefix}_accepted_mappings.tsv")
    write_table(pairing.loc[pairing["analysis_eligible_pair"]], outdir / "mapping" / f"{args.prefix}_analysis_eligible_mappings.tsv")

    if args.reference_catalog:
        ref = read_table(Path(args.reference_catalog))
        if not ref.empty:
            write_table(ref, outdir / "rrna" / f"{args.prefix}_mag_16s_reference_catalog.tsv")
            genome_ids = set(pairing["genome_id"].dropna().astype(str))
            lacking = ref.loc[~ref["genome_id"].astype(str).isin(genome_ids)].copy() if "genome_id" in ref.columns else pd.DataFrame()
            write_table(lacking, outdir / "rrna" / f"{args.prefix}_mags_without_accepted_asv.tsv")

    paper_nodes, paper_edges, paper_graph = make_paper_outputs(graph, node_features, pairing, functional)
    hetero_nodes, hetero_edges, hetero_graph, mag_pairwise, excluded_network_links = make_heterogeneous_outputs(
        graph, paper_nodes, paper_edges, pairing, agreement, abundance_matrices,
        mag_knn=args.mag_knn,
        association_modality="metagenome",
    )
    annotate_heterogeneous_functions(hetero_graph, module_table, target_modules)
    for column in [
        "selected_cycle_module_coverage",
        "selected_cycle_modules_present",
    ]:
        hetero_nodes[column] = hetero_nodes["id"].map(
            lambda node: hetero_graph.nodes[str(node)].get(column, "")
        )

    write_table(paper_nodes, outdir / "network" / f"{args.prefix}_paper_nodes.tsv")
    write_table(paper_edges, outdir / "network" / f"{args.prefix}_paper_edges.tsv")
    nx.write_graphml(paper_graph, outdir / "network" / f"{args.prefix}_paper.graphml")
    (outdir / "network" / f"{args.prefix}_paper.cyjs").write_text(json.dumps(cytoscape_json(paper_nodes.rename(columns={"ASV_ID": "id"}), paper_edges), indent=2))

    metagenome_prefix = f"{args.prefix}_metagenome"
    write_table(hetero_nodes, outdir / "network" / f"{metagenome_prefix}_heterogeneous_nodes.tsv")
    write_table(hetero_edges, outdir / "network" / f"{metagenome_prefix}_heterogeneous_edges.tsv")
    write_table(
        mag_pairwise,
        outdir / "abundance" / f"{metagenome_prefix}_mag_mag_abundance_associations.tsv",
    )
    write_table(
        excluded_network_links,
        outdir / "qc" / f"{metagenome_prefix}_heterogeneous_excluded_links.tsv",
    )
    nx.write_graphml(hetero_graph, outdir / "network" / f"{metagenome_prefix}_heterogeneous.graphml")
    (outdir / "network" / f"{metagenome_prefix}_heterogeneous.cyjs").write_text(json.dumps(cytoscape_json(hetero_nodes, hetero_edges), indent=2))
    plot_heterogeneous_network(
        hetero_graph,
        outdir,
        metagenome_prefix,
        [x.strip().lower() for x in args.ambiguity_formats.split(",") if x.strip()],
    )
    overlay_formats = [
        x.strip().lower()
        for x in args.ambiguity_formats.split(",")
        if x.strip()
    ]
    plot_heterogeneous_network(
        hetero_graph,
        outdir,
        args.prefix,
        overlay_formats,
        overlay_mode="modules",
        minimum_fraction=args.functional_module_min_fraction,
    )
    plot_heterogeneous_network(
        hetero_graph,
        outdir,
        args.prefix,
        overlay_formats,
        overlay_mode="ecological",
        minimum_fraction=args.functional_module_min_fraction,
    )
    if args.isa_dir:
        for grouping in [
            value.strip()
            for value in args.biochemical_groupings.split(",")
            if value.strip()
        ]:
            associations = load_biochemical_indicator_associations(
                Path(args.isa_dir), grouping, args.isa_q_threshold
            )
            palette_value = biochemical_palettes.get(grouping, {})
            if isinstance(palette_value, str):
                palette = {
                    item.split("=", 1)[0].strip(): item.split("=", 1)[1].strip()
                    for item in palette_value.split(",")
                    if "=" in item
                }
            else:
                palette = dict(palette_value)
            order_value = biochemical_orders.get(grouping, [])
            order = (
                [value.strip() for value in order_value.split(",") if value.strip()]
                if isinstance(order_value, str)
                else [str(value) for value in order_value]
            )
            association_table = pd.DataFrame([
                {
                    "ASV_ID": asv_id,
                    "ecological_module": hetero_graph.nodes[asv_id].get(
                        "display_module_label", "Other modules"
                    ) if asv_id in hetero_graph else np.nan,
                    "biochemical_grouping": grouping,
                    "indicator_group_membership": membership,
                    "association_scope": (
                        "multigroup" if "+" in membership else "single_group"
                    ),
                    "isa_q_threshold": args.isa_q_threshold,
                }
                for asv_id, membership in associations.items()
            ])
            write_table(
                association_table,
                outdir / "network"
                / f"{args.prefix}_ecological_modules_by_{grouping}_node_evidence.tsv",
            )
            plot_heterogeneous_network(
                hetero_graph,
                outdir,
                args.prefix,
                overlay_formats,
                overlay_mode="biochemical",
                biochemical_grouping=grouping,
                biochemical_associations=associations,
                biochemical_palette=palette,
                biochemical_order=order,
            )
    module_names = (
        module_table[["module_id", "module_name"]]
        .drop_duplicates("module_id")
        .set_index("module_id")["module_name"]
        .astype(str)
        .to_dict()
        if not module_table.empty
        else {}
    )
    for module_id in target_modules:
        write_table(
            functional_contributor_table(
                hetero_graph,
                module_id,
                args.functional_module_min_fraction,
                "metagenome",
            ),
            outdir / "functional" / f"{metagenome_prefix}_kegg_{module_id}_contributors.tsv",
        )
        plot_heterogeneous_network(
            hetero_graph,
            outdir,
            metagenome_prefix,
            overlay_formats,
            overlay_mode="functional",
            functional_module_id=module_id,
            functional_module_name=module_names.get(module_id, ""),
            minimum_fraction=args.functional_module_min_fraction,
        )
    if transcript_abundance_matrices:
        (
            transcript_nodes,
            transcript_edges,
            transcript_graph,
            transcript_pairwise,
            transcript_excluded_links,
        ) = make_heterogeneous_outputs(
            graph,
            paper_nodes,
            paper_edges,
            pairing,
            pd.DataFrame(),
            transcript_abundance_matrices,
            mag_knn=args.mag_knn,
            association_modality="metatranscriptome",
        )
        annotate_heterogeneous_functions(
            transcript_graph,
            module_table,
            target_modules,
        )
        for column in [
            "selected_cycle_module_coverage",
            "selected_cycle_modules_present",
        ]:
            transcript_nodes[column] = transcript_nodes["id"].map(
                lambda node: transcript_graph.nodes[str(node)].get(column, "")
            )
        metatranscriptome_prefix = f"{args.prefix}_metatranscriptome"
        write_table(
            transcript_nodes,
            outdir / "network" / f"{metatranscriptome_prefix}_heterogeneous_nodes.tsv",
        )
        write_table(
            transcript_edges,
            outdir / "network" / f"{metatranscriptome_prefix}_heterogeneous_edges.tsv",
        )
        write_table(
            transcript_pairwise,
            outdir / "abundance" / f"{metatranscriptome_prefix}_mag_mag_expression_associations.tsv",
        )
        write_table(
            transcript_excluded_links,
            outdir / "qc" / f"{metatranscriptome_prefix}_heterogeneous_excluded_links.tsv",
        )
        nx.write_graphml(
            transcript_graph,
            outdir / "network" / f"{metatranscriptome_prefix}_heterogeneous.graphml",
        )
        (
            outdir / "network" / f"{metatranscriptome_prefix}_heterogeneous.cyjs"
        ).write_text(
            json.dumps(cytoscape_json(transcript_nodes, transcript_edges), indent=2)
        )
        plot_heterogeneous_network(
            transcript_graph,
            outdir,
            metatranscriptome_prefix,
            [
                x.strip().lower()
                for x in args.ambiguity_formats.split(",")
                if x.strip()
            ],
        )
        for module_id in target_modules:
            write_table(
                functional_contributor_table(
                    transcript_graph,
                    module_id,
                    args.functional_module_min_fraction,
                    "metatranscriptome",
                ),
                outdir / "functional" / f"{metatranscriptome_prefix}_kegg_{module_id}_contributors.tsv",
            )
            plot_heterogeneous_network(
                transcript_graph,
                outdir,
                metatranscriptome_prefix,
                overlay_formats,
                overlay_mode="functional",
                functional_module_id=module_id,
                functional_module_name=module_names.get(module_id, ""),
                minimum_fraction=args.functional_module_min_fraction,
            )

    panel_modalities = ["metagenome"]
    if transcript_abundance_matrices:
        panel_modalities.append("metatranscriptome")
    cycle_lookup = (
        cycle_module_audit.set_index("module_id")["cycle"].astype(str).to_dict()
        if not cycle_module_audit.empty
        else {}
    )
    panel_rows = []
    panel_rows.extend([
        {
            "modality": "shared",
            "panel_order": 1,
            "panel_type": "ecological_modules",
            "functional_module_id": "",
            "functional_module_name": "",
            "cycle": "",
            "file_stem": f"network/{args.prefix}_ecological_modules",
        },
        {
            "modality": "shared",
            "panel_order": 2,
            "panel_type": "ecological_modules_and_anchors",
            "functional_module_id": "",
            "functional_module_name": "",
            "cycle": "",
            "file_stem": f"network/{args.prefix}_ecological_modules_and_anchors",
        },
    ])
    for panel_offset, grouping in enumerate(
        [
            value.strip()
            for value in args.biochemical_groupings.split(",")
            if value.strip()
        ],
        start=3,
    ):
        panel_rows.append({
            "modality": "shared",
            "panel_order": panel_offset,
            "panel_type": "ecological_module_biochemical_overlay",
            "functional_module_id": "",
            "functional_module_name": grouping,
            "cycle": "",
            "file_stem": (
                f"network/{args.prefix}_ecological_modules_by_{grouping}"
            ),
        })
    biochemical_panel_count = len([
        value for value in args.biochemical_groupings.split(",") if value.strip()
    ])
    for modality in panel_modalities:
        modality_prefix = f"{args.prefix}_{modality}"
        panel_rows.extend([
            {
                "modality": modality,
                "panel_order": 3 + biochemical_panel_count,
                "panel_type": "base_heterogeneous_network",
                "functional_module_id": "",
                "functional_module_name": "",
                "cycle": "",
                "file_stem": f"network/{modality_prefix}_heterogeneous_network",
            },
        ])
        panel_rows.extend({
            "modality": modality,
            "panel_order": 4 + biochemical_panel_count,
            "panel_type": "kegg_module_contributors",
            "functional_module_id": module_id,
            "functional_module_name": module_names.get(module_id, ""),
            "cycle": cycle_lookup.get(module_id, ""),
            "file_stem": f"network/{modality_prefix}_kegg_{module_id}_contributors",
        } for module_id in target_modules)
    write_table(
        pd.DataFrame(panel_rows),
        outdir / "network" / f"{args.prefix}_network_panel_manifest.tsv",
    )

    summary = pairing.groupby("mapping_class", dropna=False)["ASV_ID"].nunique().reset_index(name="n_asvs")
    summary["fraction_asvs"] = summary["n_asvs"] / max(1, pairing["ASV_ID"].nunique())
    write_table(summary, outdir / "qc" / f"{args.prefix}_summary.tsv")
    plot_summary(summary, outdir)

    params = {
        "graph": str(graph_path),
        "node_features": args.node_features,
        "asv_mag_pairing": args.asv_mag_pairing,
        "taxonomy": args.taxonomy,
        "min_pident": args.min_pident,
        "min_qcov": args.min_qcov,
        "asv_taxonomy_source": args.asv_taxonomy_source,
        "asv_taxonomy_min_confidence": args.asv_taxonomy_min_confidence,
        "mag_taxonomy_source": args.mag_taxonomy_source,
        "mag_abundance": args.mag_abundance,
        "mag_abundance_seqkit": args.mag_abundance_seqkit,
        "mag_abundance_normalization": args.mag_abundance_normalization,
        "mag_transcript_abundance": args.mag_transcript_abundance,
        "mag_transcript_abundance_seqkit": args.mag_transcript_abundance_seqkit,
        "mag_transcript_abundance_normalization": args.mag_transcript_abundance_normalization,
        "mag_id_mode": args.mag_id_mode,
        "mag_abundance_format": args.mag_abundance_format,
        "mag_abundance_value_col": args.mag_abundance_value_col,
        "min_shared_samples": args.min_shared_samples,
        "abundance_transform": args.abundance_transform,
        "abundance_sample_harmonization": sample_audit,
        "functional_annotation": args.functional_annotation,
        "functional_module_min_fraction": args.functional_module_min_fraction,
        "ambiguity_target_modules": args.ambiguity_target_modules,
        "note": "SPIEC-EASI edges are ASV-only statistical associations. ASV-MAG edges are sequence-derived links; abundance correlation is supporting evidence attached to those links and does not establish MAG identity.",
    }
    (outdir / "qc" / f"{args.prefix}_parameters.json").write_text(json.dumps(params, indent=2))
    info(f"ASV_MAG_NETWORK wrote outputs to {outdir}")


if __name__ == "__main__":
    main()
