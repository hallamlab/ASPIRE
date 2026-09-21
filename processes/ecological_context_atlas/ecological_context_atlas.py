#!/usr/bin/env python3
"""Build a normalized, post hoc ecological evidence atlas for retained ASVs.

All inference products are consumed read-only.  This module never changes ASV,
TITAN, compartment, network, or ASV--MAG results.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/aspire_matplotlib")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import networkx as nx
import numpy as np
import pandas as pd
import scipy

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "asv_mag_network"))
from shared_plot_style import install_publication_style
from asv_mag_network import plot_heterogeneous_network

install_publication_style()


TITAN_DIRECTION_COLORS = {"z-": "#BDBDBD", "z+": "#000000"}
TITAN_Z_LEGEND_VALUES = (5.0, 10.0, 20.0, 30.0)


def titan_marker_area(value: object) -> float:
    """Map absolute TITAN z scores to a bounded, visually distinct marker area."""
    score = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if not np.isfinite(score):
        return 20.0
    return float(45.0 + 20.0 * min(abs(float(score)), TITAN_Z_LEGEND_VALUES[-1]))


def read_table(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix.lower() in {".tsv", ".tab", ".txt"} else ","
    return pd.read_csv(path, sep=sep, low_memory=False)


def write_table(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, na_rep="")


def clean_asv(value: object) -> str:
    return re.sub(r";size=.*$", "", str(value).strip())


def split_values(value: str) -> list[str]:
    return [x.strip() for x in re.split(r"[,|]", value or "") if x.strip()]


def first_existing(frame: pd.DataFrame, names: list[str]) -> str | None:
    lower = {str(c).lower(): c for c in frame.columns}
    for name in names:
        if name in frame.columns:
            return name
        if name.lower() in lower:
            return lower[name.lower()]
    return None


def taxonomy_table(path: Path) -> pd.DataFrame:
    raw = read_table(path)
    id_col = first_existing(raw, ["ASV_ID", "Feature ID", "ASV", "Taxon"])
    tax_col = first_existing(raw, ["taxonomy", "Taxon", "lineage"])
    if id_col is None or tax_col is None:
        raise ValueError("Taxonomy table requires an ASV identifier and taxonomy column")
    out = pd.DataFrame({"ASV_ID": raw[id_col].map(clean_asv), "taxonomy": raw[tax_col].fillna("").astype(str)})
    prefixes = {"d": "domain", "p": "phylum", "c": "class", "o": "order", "f": "family", "g": "genus", "s": "species"}
    for prefix, label in prefixes.items():
        pattern = rf"(?:^|;\s*){prefix}__([^;]*)"
        out[label] = out.taxonomy.str.extract(pattern, expand=False).fillna("").str.strip()
    return out.drop_duplicates("ASV_ID")


def load_network(modules_path: Path, features_path: Path, layout_path: Path) -> pd.DataFrame:
    modules = read_table(modules_path)
    mid = first_existing(modules, ["ASV_ID", "Taxon"])
    modules = modules.rename(columns={mid: "ASV_ID"})
    modules["ASV_ID"] = modules.ASV_ID.map(clean_asv)
    features = read_table(features_path)
    fid = first_existing(features, ["ASV_ID", "Taxon"])
    features = features.rename(columns={fid: "ASV_ID"})
    features["ASV_ID"] = features.ASV_ID.map(clean_asv)
    layout = read_table(layout_path)
    lid = first_existing(layout, ["ASV_ID", "Taxon"])
    layout = layout.rename(columns={lid: "ASV_ID"})
    layout["ASV_ID"] = layout.ASV_ID.map(clean_asv)
    return modules.merge(features, on="ASV_ID", how="outer", suffixes=("", "_network")).merge(layout, on="ASV_ID", how="left")


def weighted_quantile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    valid = np.isfinite(values) & np.isfinite(weights) & (weights >= 0)
    values, weights = values[valid], weights[valid]
    if not len(values) or weights.sum() <= 0:
        return np.nan
    order = np.argsort(values)
    values, weights = values[order], weights[order]
    return float(np.interp(q * weights.sum(), np.cumsum(weights), values))


def make_sample_context(ra: pd.DataFrame, metadata: pd.DataFrame, microbial: pd.DataFrame,
                        taxonomy: pd.DataFrame, network: pd.DataFrame, sample_col: str) -> pd.DataFrame:
    if sample_col not in ra.columns:
        sample_col = ra.columns[0]
    long = ra.rename(columns={sample_col: "sample_ID"}).melt("sample_ID", var_name="ASV_ID", value_name="relative_abundance")
    long["ASV_ID"] = long.ASV_ID.map(clean_asv)
    long["relative_abundance"] = pd.to_numeric(long.relative_abundance, errors="coerce").fillna(0.0)
    long["detected"] = long.relative_abundance > 0
    meta_sample = first_existing(metadata, [sample_col, "sampleID", "Sample", "sample_code"])
    meta = metadata.rename(columns={meta_sample: "sample_ID"}).drop_duplicates("sample_ID")
    mc_sample = first_existing(microbial, [sample_col, "sampleID", "Sample", "sample_ID"])
    mc = microbial.rename(columns={mc_sample: "sample_ID"}).drop_duplicates("sample_ID")
    long = long.merge(meta, on="sample_ID", how="left", validate="many_to_one")
    long = long.merge(mc, on="sample_ID", how="left", validate="many_to_one", suffixes=("", "_mc"))
    long = long.merge(taxonomy, on="ASV_ID", how="left", validate="many_to_one")
    keep_net = [c for c in network.columns if c not in {"taxonomy", "domain", "phylum", "class", "order", "family", "genus", "species"}]
    return long.merge(network[keep_net].drop_duplicates("ASV_ID"), on="ASV_ID", how="left", validate="many_to_one")


def context_profiles(sample: pd.DataFrame, context_cols: list[str], cruise_col: str,
                     depth_col: str, date_col: str) -> pd.DataFrame:
    rows = []
    for col in context_cols:
        if col not in sample.columns:
            continue
        frame = sample.loc[sample[col].notna() & (sample[col].astype(str).str.strip() != "")].copy()
        for (asv, level), group in frame.groupby(["ASV_ID", col], observed=True):
            ra = group.relative_abundance.to_numpy(float)
            detected = ra > 0
            row = {
                "ASV_ID": asv, "context_type": col, "context_level": level,
                "number_of_available_samples": len(group), "number_of_detected_samples": int(detected.sum()),
                "prevalence": float(detected.mean()), "mean_relative_abundance": float(np.mean(ra)),
                "median_relative_abundance": float(np.median(ra)), "maximum_relative_abundance": float(np.max(ra)),
                "mean_relative_abundance_when_detected": float(np.mean(ra[detected])) if detected.any() else 0.0,
            }
            if cruise_col in group:
                row["number_of_cruises"] = group[cruise_col].nunique(dropna=True)
                row["number_of_detected_cruises"] = group.loc[detected, cruise_col].nunique(dropna=True)
            if depth_col in group:
                depths = pd.to_numeric(group.loc[detected, depth_col], errors="coerce")
                row.update({"detected_depth_min": depths.min(), "detected_depth_median": depths.median(), "detected_depth_max": depths.max()})
            if date_col in group:
                dates = pd.to_datetime(group.loc[detected, date_col], errors="coerce")
                row.update({"first_detected_date": dates.min(), "last_detected_date": dates.max()})
            rows.append(row)
    result = pd.DataFrame(rows)
    if len(result):
        totals = sample.groupby("ASV_ID").relative_abundance.sum().rename("asv_total_relative_abundance")
        subtotals = sample.groupby(["ASV_ID"] + [], observed=True).relative_abundance.sum()
        result = result.merge(totals, on="ASV_ID", how="left")
        level_sum = sample.melt(id_vars=["ASV_ID", "relative_abundance"], value_vars=[c for c in context_cols if c in sample], var_name="context_type", value_name="context_level")
        level_sum = level_sum.dropna(subset=["context_level"]).groupby(["ASV_ID", "context_type", "context_level"], observed=True).relative_abundance.sum().rename("context_relative_abundance_sum").reset_index()
        result = result.merge(level_sum, on=["ASV_ID", "context_type", "context_level"], how="left")
        result["share_of_asv_total_abundance"] = result.context_relative_abundance_sum / result.asv_total_relative_abundance.replace(0, np.nan)
    return result


def environmental_response(sample: pd.DataFrame, titan: pd.DataFrame, correlations: pd.DataFrame,
                           cruise_col: str, depth_col: str, date_col: str) -> pd.DataFrame:
    titan = titan.copy()
    titan["ASV_ID"] = titan.ASV_ID.map(clean_asv)
    corr = correlations.copy()
    corr["ASV_ID"] = corr.ASV_ID.map(clean_asv)
    corr = corr.rename(columns={"measurement": "environmental_variable", "rho": "spearman_rho", "p_value": "spearman_p_value", "q_value": "spearman_q_value", "n": "spearman_n"})
    result = titan.merge(corr, on=["ASV_ID", "environmental_variable"], how="outer", suffixes=("", "_spearman"))
    support = []
    for row in result[["ASV_ID", "environmental_variable", "change_point"]].itertuples(index=False):
        variable = row.environmental_variable
        if variable not in sample.columns:
            support.append({"ASV_ID": row.ASV_ID, "environmental_variable": variable})
            continue
        support_cols = [c for c in [cruise_col, depth_col, date_col] if c in sample]
        group = sample.loc[sample.ASV_ID.eq(row.ASV_ID), ["relative_abundance", variable] + support_cols].copy()
        group[variable] = pd.to_numeric(group[variable], errors="coerce")
        group = group.dropna(subset=[variable])
        values = group[variable].to_numpy(float); weights = group.relative_abundance.to_numpy(float)
        cp = pd.to_numeric(pd.Series([row.change_point]), errors="coerce").iloc[0]
        below = group[variable] < cp if np.isfinite(cp) else pd.Series(False, index=group.index)
        above = group[variable] >= cp if np.isfinite(cp) else pd.Series(False, index=group.index)
        item = {
            "ASV_ID": row.ASV_ID, "environmental_variable": variable,
            "observed_measurement_min": np.min(values) if len(values) else np.nan,
            "observed_measurement_q05": np.quantile(values, .05) if len(values) else np.nan,
            "observed_measurement_median": np.median(values) if len(values) else np.nan,
            "observed_measurement_q95": np.quantile(values, .95) if len(values) else np.nan,
            "observed_measurement_max": np.max(values) if len(values) else np.nan,
            "abundance_weighted_measurement_q25": weighted_quantile(values, weights, .25),
            "abundance_weighted_measurement_median": weighted_quantile(values, weights, .5),
            "abundance_weighted_measurement_q75": weighted_quantile(values, weights, .75),
            "samples_below_change_point": int(below.sum()), "samples_at_or_above_change_point": int(above.sum()),
            "detected_samples_below_change_point": int(((group.relative_abundance > 0) & below).sum()),
            "detected_samples_at_or_above_change_point": int(((group.relative_abundance > 0) & above).sum()),
        }
        if cruise_col in group:
            item["cruises_below_change_point"] = group.loc[below, cruise_col].nunique(dropna=True)
            item["cruises_at_or_above_change_point"] = group.loc[above, cruise_col].nunique(dropna=True)
        detected = group.relative_abundance > 0
        for label, mask in [("below", below & detected), ("at_or_above", above & detected)]:
            if depth_col in group:
                depths = pd.to_numeric(group.loc[mask, depth_col], errors="coerce")
                item[f"detected_depth_min_{label}_change_point"] = depths.min()
                item[f"detected_depth_median_{label}_change_point"] = depths.median()
                item[f"detected_depth_max_{label}_change_point"] = depths.max()
            if date_col in group:
                dates = pd.to_datetime(group.loc[mask, date_col], errors="coerce")
                item[f"first_detected_date_{label}_change_point"] = dates.min()
                item[f"last_detected_date_{label}_change_point"] = dates.max()
        support.append(item)
    result = result.merge(pd.DataFrame(support), on=["ASV_ID", "environmental_variable"], how="left")
    titan_sign = result.response_direction.map({"z+": 1, "z-": -1})
    rho_sign = np.sign(pd.to_numeric(result.spearman_rho, errors="coerce"))
    result["titan_spearman_direction_concordant"] = np.where(titan_sign.notna() & rho_sign.ne(0), titan_sign.eq(rho_sign), pd.NA)
    return result


def _unique_text(values: pd.Series) -> str:
    """Return sorted, pipe-delimited nonempty values for compact bridge columns."""
    cleaned = {
        str(value).strip() for value in values.dropna()
        if str(value).strip() and str(value).strip().lower() not in {"nan", "none"}
    }
    return "|".join(sorted(cleaned))


def _boolean_series(values: pd.Series) -> pd.Series:
    """Parse native or text booleans without treating the string 'False' as true."""
    if pd.api.types.is_bool_dtype(values):
        return values.fillna(False).astype(bool)
    return values.fillna("").astype(str).str.strip().str.lower().isin(
        {"true", "t", "1", "yes", "y"}
    )


def aggregate_mag_links(mappings: pd.DataFrame) -> pd.DataFrame:
    """Collapse analysis-eligible ASV--MAG links without discarding one-to-many links."""
    if mappings.empty:
        return pd.DataFrame(columns=[
            "ASV_ID", "has_analysis_eligible_mag_link",
            "analysis_eligible_mag_link_count", "analysis_eligible_mag_species_count",
            "analysis_eligible_mag_match_ids", "analysis_eligible_mag_genome_ids",
            "analysis_eligible_mag_species", "analysis_eligible_mag_genera",
            "analysis_eligible_mag_phyla", "asv_mag_pairing_statuses",
        ])
    frame = mappings.copy()
    if "ASV_ID" not in frame:
        raise ValueError("Analysis-eligible ASV--MAG mappings require ASV_ID")
    frame["ASV_ID"] = frame.ASV_ID.map(clean_asv)
    match_col = first_existing(frame, ["mag_match_id", "mag_genome_id", "genome_id"])
    genome_col = first_existing(frame, ["genome_id", "mag_genome_id", "mag_match_id"])
    species_col = first_existing(frame, ["mag_tax_species", "mag_species"])
    genus_col = first_existing(frame, ["mag_tax_genus", "mag_genus"])
    phylum_col = first_existing(frame, ["mag_tax_phylum", "mag_phylum"])
    status_col = first_existing(frame, ["pairing_status", "mapping_class"])
    rows = []
    for asv, group in frame.groupby("ASV_ID", observed=True, sort=True):
        rows.append({
            "ASV_ID": asv,
            "has_analysis_eligible_mag_link": True,
            "analysis_eligible_mag_link_count": (
                int(group[match_col].dropna().astype(str).nunique()) if match_col else len(group)
            ),
            "analysis_eligible_mag_species_count": (
                int(group[species_col].dropna().astype(str).nunique()) if species_col else 0
            ),
            "analysis_eligible_mag_match_ids": _unique_text(group[match_col]) if match_col else "",
            "analysis_eligible_mag_genome_ids": _unique_text(group[genome_col]) if genome_col else "",
            "analysis_eligible_mag_species": _unique_text(group[species_col]) if species_col else "",
            "analysis_eligible_mag_genera": _unique_text(group[genus_col]) if genus_col else "",
            "analysis_eligible_mag_phyla": _unique_text(group[phylum_col]) if phylum_col else "",
            "asv_mag_pairing_statuses": _unique_text(group[status_col]) if status_col else "",
        })
    return pd.DataFrame(rows)


def ecological_linkage_crosswalk(
    network: pd.DataFrame,
    taxonomy: pd.DataFrame,
    response: pd.DataFrame,
    module_compartment: pd.DataFrame,
    module_measurement: pd.DataFrame,
    member_support: pd.DataFrame,
    mappings: pd.DataFrame,
    q_threshold: float,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Join ASV, module, compartment, measurement, and genome evidence.

    Module membership is a consensus-network assignment rather than a null-
    hypothesis test.  Its evidential quantity is therefore node/module
    stability, not a fabricated p-value.  The remaining links retain their
    native p- and q-values.
    """
    required_compartment = {
        "grouping", "ecological_module", "effect_eta_squared", "p_value",
        "q_value", "adjusted_r_squared", "highest_group_level",
    }
    required_measurement = {
        "ecological_module", "measurement", "rho", "p_value", "q_value",
    }
    required_member = {
        "ecological_module", "measurement", "enrichment_p_value",
        "enrichment_q_value",
    }
    for label, frame, required in [
        ("module-compartment association", module_compartment, required_compartment),
        ("module-measurement association", module_measurement, required_measurement),
        ("module-member support", member_support, required_member),
    ]:
        missing = sorted(required.difference(frame.columns))
        if missing:
            raise ValueError(f"{label} table is missing required columns: {', '.join(missing)}")

    module_col = first_existing(network, ["module_label", "ecological_module"])
    if module_col is None:
        raise ValueError("Network node table requires module_label or ecological_module")
    network_fields = [
        "ASV_ID", module_col, "module_id", "node_stability", "graph_variant",
        "method", "Degree", "Betweenness", "Closeness", "EigenCentral",
        "participation_coefficient", "is_module_anchor", "is_anchor_any",
    ]
    members = network[[c for c in network_fields if c in network.columns]].copy()
    members = members.rename(columns={module_col: "ecological_module"})
    members["ASV_ID"] = members.ASV_ID.map(clean_asv)
    members = members.loc[
        members.ecological_module.notna()
        & members.ecological_module.astype(str).str.strip().ne("")
    ].drop_duplicates("ASV_ID")
    members = members.merge(taxonomy, on="ASV_ID", how="left", validate="one_to_one")

    comp = module_compartment.rename(columns={
        "grouping": "compartment_strategy",
        "effect_eta_squared": "module_compartment_eta_squared",
        "adjusted_r_squared": "module_compartment_adjusted_r_squared",
        "p_value": "module_compartment_p_value",
        "q_value": "module_compartment_q_value",
        "highest_group_level": "module_highest_compartment",
        "highest_group_mean": "module_highest_compartment_mean",
        "n_units": "module_compartment_n_units",
        "n_levels": "module_compartment_n_levels",
    }).copy()
    comp["module_compartment_significant_q"] = (
        pd.to_numeric(comp.module_compartment_q_value, errors="coerce") <= q_threshold
    )

    measurement = module_measurement.rename(columns={
        "measurement": "environmental_variable",
        "rho": "module_measurement_spearman_rho",
        "p_value": "module_measurement_p_value",
        "q_value": "module_measurement_q_value",
        "samples_measured": "module_measurement_samples",
        "cruises_measured": "module_measurement_cruises",
        "depths_measured": "module_measurement_depths",
    }).copy()
    measurement["module_measurement_significant_q"] = (
        pd.to_numeric(measurement.module_measurement_q_value, errors="coerce") <= q_threshold
    )

    member = member_support.rename(columns={
        "measurement": "environmental_variable",
        "enrichment_p_value": "module_member_enrichment_p_value",
        "enrichment_q_value": "module_member_enrichment_q_value",
        "enrichment_odds_ratio": "module_member_enrichment_odds_ratio",
        "significant_asvs": "module_significant_member_asvs",
        "significant_fraction": "module_significant_member_fraction",
        "significant_positive_asvs": "module_significant_positive_asvs",
        "significant_negative_asvs": "module_significant_negative_asvs",
        "sign_concordance": "module_member_sign_concordance",
        "median_rho_all_asvs": "module_member_median_rho_all",
        "median_rho_significant_asvs": "module_member_median_rho_significant",
    }).copy()
    member["module_member_enrichment_significant_q"] = (
        pd.to_numeric(member.module_member_enrichment_q_value, errors="coerce") <= q_threshold
    )

    direct_fields = [
        "ASV_ID", "environmental_variable", "spearman_rho", "spearman_p_value",
        "spearman_q_value", "spearman_n", "change_point", "response_direction",
        "indicator_score", "z_score", "purity", "reliability",
        "bootstrap_cp_q05", "bootstrap_cp_q50", "bootstrap_cp_q95",
        "passes_purity", "passes_reliability", "passes_purity_and_reliability",
        "titan_spearman_direction_concordant",
    ]
    direct = response[[c for c in direct_fields if c in response.columns]].copy()
    direct["ASV_ID"] = direct.ASV_ID.map(clean_asv)
    direct["asv_measurement_spearman_significant_q"] = (
        pd.to_numeric(direct.spearman_q_value, errors="coerce") <= q_threshold
    )
    if "passes_purity_and_reliability" in direct:
        titan_supported = _boolean_series(direct.passes_purity_and_reliability)
    else:
        titan_supported = pd.Series(False, index=direct.index)
    direct["asv_measurement_titan_supported"] = titan_supported
    direct["asv_measurement_supported_by_either_method"] = (
        direct.asv_measurement_spearman_significant_q | titan_supported
    )

    full = members.merge(comp, on="ecological_module", how="inner", validate="many_to_many")
    full = full.merge(measurement, on="ecological_module", how="inner", validate="many_to_many")
    full = full.merge(
        member,
        on=["ecological_module", "environmental_variable"],
        how="left", validate="many_to_one", suffixes=("", "_member"),
    )
    full = full.merge(
        direct,
        on=["ASV_ID", "environmental_variable"],
        how="left", validate="many_to_one",
    )
    mag = aggregate_mag_links(mappings)
    full = full.merge(mag, on="ASV_ID", how="left", validate="many_to_one")
    full["has_analysis_eligible_mag_link"] = _boolean_series(full.has_analysis_eligible_mag_link)
    for col in ["analysis_eligible_mag_link_count", "analysis_eligible_mag_species_count"]:
        full[col] = pd.to_numeric(full[col], errors="coerce").fillna(0).astype(int)
    for col in [
        "analysis_eligible_mag_match_ids", "analysis_eligible_mag_genome_ids",
        "analysis_eligible_mag_species", "analysis_eligible_mag_genera",
        "analysis_eligible_mag_phyla", "asv_mag_pairing_statuses",
    ]:
        full[col] = full[col].fillna("")

    module_rho = pd.to_numeric(full.module_measurement_spearman_rho, errors="coerce")
    asv_rho = pd.to_numeric(full.spearman_rho, errors="coerce")
    full["asv_module_measurement_direction_concordant"] = np.where(
        module_rho.notna() & asv_rho.notna() & module_rho.ne(0) & asv_rho.ne(0),
        np.sign(module_rho).eq(np.sign(asv_rho)), pd.NA,
    )
    full["complete_significant_chain"] = (
        full.module_compartment_significant_q
        & full.module_measurement_significant_q
        & full.asv_measurement_spearman_significant_q
    )
    full["complete_significant_enriched_chain"] = (
        full.complete_significant_chain
        & full.module_member_enrichment_significant_q.fillna(False)
    )
    full["linkage_q_threshold"] = q_threshold
    sort_cols = ["ASV_ID", "ecological_module", "compartment_strategy", "environmental_variable"]
    full = full.sort_values(sort_cols, kind="stable").reset_index(drop=True)
    supported = full.loc[full.complete_significant_chain].reset_index(drop=True)
    genome_linked = supported.loc[supported.has_analysis_eligible_mag_link].reset_index(drop=True)

    definitions = pd.DataFrame([
        {
            "linkage": "ASV to ecological module",
            "evidence_fields": "ecological_module; node_stability",
            "support_rule": "Consensus module assignment retained; node_stability reported",
            "inferential_statistic": "None",
            "interpretation": "Network membership, not a null-hypothesis test; no p-value is applicable",
        },
        {
            "linkage": "Ecological module to compartment strategy",
            "evidence_fields": "module_compartment_eta_squared; module_compartment_adjusted_r_squared; module_compartment_p_value; module_compartment_q_value",
            "support_rule": f"module_compartment_q_value <= {q_threshold:g}",
            "inferential_statistic": "Cruise-restricted permutation test with Benjamini-Hochberg correction",
            "interpretation": "Omnibus difference in module abundance among levels; module_highest_compartment identifies the largest mean but is not a separate pairwise test",
        },
        {
            "linkage": "Ecological module to environmental measurement",
            "evidence_fields": "module_measurement_spearman_rho; module_measurement_p_value; module_measurement_q_value",
            "support_rule": f"module_measurement_q_value <= {q_threshold:g}",
            "inferential_statistic": "Spearman rank correlation with Benjamini-Hochberg correction",
            "interpretation": "Bivariate association between module abundance and observed measurement values",
        },
        {
            "linkage": "Module-member enrichment for environmental measurement",
            "evidence_fields": "module_member_enrichment_odds_ratio; module_member_enrichment_p_value; module_member_enrichment_q_value",
            "support_rule": f"module_member_enrichment_q_value <= {q_threshold:g}",
            "inferential_statistic": "Member enrichment test with Benjamini-Hochberg correction",
            "interpretation": "Tests whether significant ASV associations are overrepresented within the module relative to the network cohort",
        },
        {
            "linkage": "ASV to environmental measurement",
            "evidence_fields": "spearman_rho; spearman_p_value; spearman_q_value; TITAN change_point; response_direction; purity; reliability",
            "support_rule": f"Strict chain: spearman_q_value <= {q_threshold:g}; complementary threshold support: TITAN passes purity and reliability",
            "inferential_statistic": "Spearman rank correlation with Benjamini-Hochberg correction; independent TITAN bootstrap support",
            "interpretation": "Direct ASV response to the same measurement represented by the module-level association",
        },
        {
            "linkage": "ASV to analysis-eligible MAGs",
            "evidence_fields": "analysis_eligible_mag_match_ids; analysis_eligible_mag_species; analysis_eligible_mag_link_count",
            "support_rule": "Exact downstream sequence-link criteria defined by the ASV-MAG process",
            "inferential_statistic": "None",
            "interpretation": "Sequence identity/coverage mapping, not a null-hypothesis test; one-to-many links are retained",
        },
    ])
    return full, supported, genome_linked, definitions


def asv_ecological_module_mag_crosswalk(
    network: pd.DataFrame,
    taxonomy: pd.DataFrame,
    mappings: pd.DataFrame,
) -> pd.DataFrame:
    """Return one auditable row per connected ASV and ecological module.

    The statistical linkage table deliberately repeats an ASV across compartment
    strategies and measurements.  This companion inventory instead provides the
    direct module-membership and genome-link lookup needed to enumerate and
    interpret the ASVs assigned to each ecological module.
    """
    network_cols = [c for c in [
        "ASV_ID", "module_label", "module_id", "node_stability",
        "graph_variant", "method", "Degree", "Degree_thresholded",
        "Betweenness", "Closeness", "EigenCentral",
    ] if c in network]
    members = network[network_cols].drop_duplicates("ASV_ID").copy()
    members = members.rename(columns={"module_label": "ecological_module"})
    members["module_asv_count"] = members.groupby(
        "ecological_module", observed=True
    )["ASV_ID"].transform("size")

    tax_cols = [c for c in [
        "ASV_ID", "taxonomy", "domain", "phylum", "class", "order",
        "family", "genus", "species",
    ] if c in taxonomy]
    if tax_cols:
        members = members.merge(
            taxonomy[tax_cols].drop_duplicates("ASV_ID"),
            on="ASV_ID", how="left", validate="one_to_one",
        )

    mag = aggregate_mag_links(mappings)
    members = members.merge(mag, on="ASV_ID", how="left", validate="one_to_one")
    members["has_analysis_eligible_mag_link"] = _boolean_series(
        members.get("has_analysis_eligible_mag_link", False)
    )
    for col in ["analysis_eligible_mag_link_count", "analysis_eligible_mag_species_count"]:
        members[col] = pd.to_numeric(members.get(col), errors="coerce").fillna(0).astype(int)
    for col in [
        "analysis_eligible_mag_match_ids", "analysis_eligible_mag_genome_ids",
        "analysis_eligible_mag_species", "analysis_eligible_mag_genera",
        "analysis_eligible_mag_phyla", "asv_mag_pairing_statuses",
    ]:
        members[col] = members.get(col, pd.Series(index=members.index, dtype=object)).fillna("")
    members["module_mag_linked_asv_count"] = members.groupby(
        "ecological_module", observed=True
    )["has_analysis_eligible_mag_link"].transform("sum").astype(int)
    members["module_mag_linked_asv_fraction"] = (
        members.module_mag_linked_asv_count / members.module_asv_count
    )
    return members.sort_values(
        ["module_id", "ASV_ID"], kind="stable"
    ).reset_index(drop=True)


def indicator_evidence(directory: Path) -> pd.DataFrame:
    rows = []
    for path in sorted(directory.glob("*_indicator_species_summary.tsv")):
        frame = read_table(path)
        aid = first_existing(frame, ["ASV_ID", "ASV"])
        if aid is None:
            continue
        grouping = path.name.removesuffix("_indicator_species_summary.tsv")
        for rec in frame.to_dict("records"):
            rows.append({
                "ASV_ID": clean_asv(rec.get(aid, "")), "grouping_analysis": grouping,
                "indicator_statistic": rec.get("stat"), "p_value": rec.get("p.value"),
                "q_value": rec.get("q.value"), "significant": rec.get("significant"),
                "association_group_count": rec.get("association_group_count"),
                "association_groups": rec.get("association_groups"),
                "association_scope": rec.get("association_scope"), "indicator_class": rec.get("indicator_class"),
            })
    return pd.DataFrame(rows)


def mag_function_links(mappings: pd.DataFrame, functional: pd.DataFrame, network: pd.DataFrame,
                       taxonomy: pd.DataFrame) -> pd.DataFrame:
    mappings = mappings.copy(); mappings["ASV_ID"] = mappings.ASV_ID.map(clean_asv)
    # mag_match_id is the normalized cross-process identifier.  genome_id retains
    # source-specific spelling and is intentionally not used as a compound key.
    keys = ["mag_match_id"] if "mag_match_id" in mappings and "mag_match_id" in functional else ["genome_id"]
    mapping_cols = [c for c in [
        "ASV_ID", "genome_id", "mag_match_id", "link_pident", "link_qcov",
        "link_bitscore", "mapping_class", "asv_link_weight", "shared_samples",
        "spearman_rho", "mag_completeness", "mag_contamination", "mag_qscore",
        "mag_mimag_tier", "mag_tax_domain", "mag_tax_phylum", "mag_tax_class",
        "mag_tax_order", "mag_tax_family", "mag_tax_genus", "mag_tax_species",
    ] if c in mappings]
    merged = mappings[mapping_cols].drop_duplicates().merge(functional, on=keys, how="left", suffixes=("", "_function"))
    net_cols = [c for c in ["ASV_ID", "module_id", "module_label", "Degree", "Betweenness", "Closeness", "EigenCentral"] if c in network]
    return merged.merge(network[net_cols].drop_duplicates("ASV_ID"), on="ASV_ID", how="left").merge(taxonomy, on="ASV_ID", how="left", suffixes=("", "_asv"))


def multiomics_context(abundance_dir: Path, sample: pd.DataFrame) -> pd.DataFrame:
    frames = []
    pattern = re.compile(r".*_(metagenome|metatranscriptome)_(ASV[^_]+)_species_recruitment_by_sample\.tsv$")
    for path in sorted(abundance_dir.glob("*_species_recruitment_by_sample.tsv")):
        match = pattern.match(path.name)
        if not match:
            continue
        frame = read_table(path); frame["data_modality"] = match.group(1); frame["ASV_ID"] = match.group(2)
        sample_key = first_existing(frame, ["sample_ID", "sample", "sampleID"])
        if sample_key is not None:
            frame = frame.rename(columns={sample_key: "sample_join_key"})
        frames.append(frame)
    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, ignore_index=True, sort=False)
    context = sample.copy()
    # Use the explicit environmental sample field shared with recruitment tables;
    # do not infer correspondence from row order or edit identifiers heuristically.
    context_key = "sample" if "sample" in context.columns else "sample_ID"
    context["sample_join_key"] = context[context_key].astype(str)
    categorical = [c for c in ["microbial_compartment", "o2_compartment", "gmm_component", "o2_subcompartment_final", "Season", "renewal_phase", "cruise_group"] if c in context]
    keys = ["ASV_ID", "sample_join_key"]
    if context.duplicated(keys).any():
        grouped = context.groupby(keys, dropna=False, observed=True)
        collapsed = grouped.relative_abundance.mean().rename("relative_abundance").reset_index()
        collapsed = collapsed.merge(grouped.size().rename("number_of_asv_sample_records").reset_index(), on=keys)
        for col in categorical:
            first = grouped[col].first().rename(col).reset_index()
            unique = grouped[col].nunique(dropna=True).rename("_nunique").reset_index()
            first = first.merge(unique, on=keys); first.loc[first._nunique > 1, col] = pd.NA
            collapsed = collapsed.merge(first[keys + [col]], on=keys)
        context = collapsed
    else:
        context = context[keys + ["relative_abundance"] + categorical].copy()
        context["number_of_asv_sample_records"] = 1
    cols = [c for c in ["ASV_ID", "sample_join_key", "number_of_asv_sample_records", "relative_abundance", "microbial_compartment", "o2_compartment", "gmm_component", "o2_subcompartment_final", "Season", "renewal_phase", "cruise_group"] if c in context]
    out = out.merge(context[cols], on=["ASV_ID", "sample_join_key"], how="left")
    out["sample_context_matched"] = out["relative_abundance"].notna()
    return out


def normalize_mag_join_id(value: object, mode: str) -> str:
    """Normalize an upstream quantification ID to the ASV-MAG join key."""
    if value is None or pd.isna(value):
        return ""
    text = re.sub(r"\.(?:fa|fna|fasta)(?:\.gz)?$", "", str(value).strip())
    if mode == "suffix_after_double_underscore" and "__" in text:
        text = text.rsplit("__", 1)[-1]
    return text


def load_wide_genome_metric(
    path: Path | None, value_name: str, mag_id_mode: str
) -> pd.DataFrame:
    """Convert a genome-by-sample metric matrix to stable long-form join keys."""
    if path is None:
        return pd.DataFrame(columns=["mag_match_id", "sample_join_key", value_name])
    table = read_table(path)
    if table.empty or len(table.columns) < 2:
        return pd.DataFrame(columns=["mag_match_id", "sample_join_key", value_name])
    genome_col = "genome_id" if "genome_id" in table else table.columns[0]
    long = table.melt(
        id_vars=[genome_col], var_name="sample_join_key", value_name=value_name
    )
    long["mag_match_id"] = long[genome_col].map(
        lambda value: normalize_mag_join_id(value, mag_id_mode)
    )
    long[value_name] = pd.to_numeric(long[value_name], errors="coerce")
    return long[["mag_match_id", "sample_join_key", value_name]].drop_duplicates(
        ["mag_match_id", "sample_join_key"]
    )


def paired_multiomics(
    multiomics: pd.DataFrame,
    *,
    log2_ratio_path: Path | None = None,
    log2_se_path: Path | None = None,
    mag_id_mode: str = "exact",
    max_log2_se: float = 0.3,
) -> pd.DataFrame:
    """Pair DNA/RNA recruitment only at identical ASV--MAG--sample keys."""
    if multiomics.empty:
        return pd.DataFrame()
    keys = [c for c in ["ASV_ID", "mag_match_id", "species", "sample_join_key", "cruise", "date", "depth_m"] if c in multiomics]
    values = multiomics.pivot_table(index=keys, columns="data_modality", values="normalized_recruitment", aggfunc="first").reset_index()
    values.columns.name = None
    context_cols = [c for c in ["number_of_asv_sample_records", "relative_abundance", "microbial_compartment", "o2_compartment", "gmm_component", "o2_subcompartment_final", "Season", "renewal_phase", "cruise_group", "sample_context_matched"] if c in multiomics]
    if context_cols:
        context = multiomics[["ASV_ID", "sample_join_key"] + context_cols].drop_duplicates(["ASV_ID", "sample_join_key"])
        values = values.merge(context, on=["ASV_ID", "sample_join_key"], how="left", validate="many_to_one")
    modes = {}
    if "normalization" in multiomics:
        for modality, group in multiomics.groupby("data_modality", observed=True):
            observed = sorted(set(group["normalization"].dropna().astype(str)))
            modes[str(modality)] = observed[0] if len(observed) == 1 else "mixed"
    values["metagenome_units"] = modes.get("metagenome", pd.NA)
    values["metatranscriptome_units"] = modes.get("metatranscriptome", pd.NA)

    if log2_ratio_path is not None:
        ratio = load_wide_genome_metric(
            log2_ratio_path, "log2_rna_to_dna_recruitment_ratio", mag_id_mode
        )
        values = values.merge(
            ratio, on=["mag_match_id", "sample_join_key"], how="left",
            validate="many_to_one",
        )
        if log2_se_path is not None:
            uncertainty = load_wide_genome_metric(
                log2_se_path, "log2_rna_to_dna_standard_error", mag_id_mode
            )
            values = values.merge(
                uncertainty, on=["mag_match_id", "sample_join_key"], how="left",
                validate="many_to_one",
            )
            passes = (
                values["log2_rna_to_dna_standard_error"].notna()
                & values["log2_rna_to_dna_standard_error"].le(max_log2_se)
            )
            values.loc[~passes, "log2_rna_to_dna_recruitment_ratio"] = np.nan
            values["paired_tpm_ratio_passes_se"] = passes
        values["rna_to_dna_recruitment_ratio"] = np.exp2(
            values["log2_rna_to_dna_recruitment_ratio"]
        )
        values["rna_dna_ratio_basis"] = "paired_tpm_over_tpm"
    elif (
        "metagenome" in values
        and "metatranscriptome" in values
        and modes.get("metagenome") == modes.get("metatranscriptome")
    ):
        values["rna_to_dna_recruitment_ratio"] = (
            values.metatranscriptome / values.metagenome.replace(0, np.nan)
        )
        positive = values.rna_to_dna_recruitment_ratio.where(
            values.rna_to_dna_recruitment_ratio > 0
        )
        values["log2_rna_to_dna_recruitment_ratio"] = np.log2(positive)
        values["rna_dna_ratio_basis"] = "identically_normalized_values"
    if "metagenome" in values and "metatranscriptome" in values:
        values["paired_dna_rna_available"] = (
            values.metagenome.notna() & values.metatranscriptome.notna()
        )
    return values


def evidence_index(sample: pd.DataFrame, response: pd.DataFrame, indicators: pd.DataFrame,
                   network: pd.DataFrame, mappings: pd.DataFrame, functional: pd.DataFrame) -> pd.DataFrame:
    base = sample.groupby("ASV_ID").agg(number_of_samples=("sample_ID", "size"), detected_samples=("detected", "sum"), mean_relative_abundance=("relative_abundance", "mean"), maximum_relative_abundance=("relative_abundance", "max")).reset_index()
    base["prevalence"] = base.detected_samples / base.number_of_samples
    summaries = {
        "number_of_titan_variables": response.groupby("ASV_ID").environmental_variable.nunique(),
        "number_of_reliable_titan_variables": response.loc[response.get("passes_purity_and_reliability", False).astype(str).str.lower().eq("true")].groupby("ASV_ID").environmental_variable.nunique() if "passes_purity_and_reliability" in response else pd.Series(dtype=float),
        "number_of_significant_spearman_variables": response.loc[pd.to_numeric(response.get("spearman_q_value"), errors="coerce") <= .05].groupby("ASV_ID").environmental_variable.nunique(),
        "number_of_indicator_analyses": indicators.groupby("ASV_ID").grouping_analysis.nunique() if len(indicators) else pd.Series(dtype=float),
        "number_of_significant_indicator_analyses": indicators.loc[indicators.significant.astype(str).str.lower().eq("true")].groupby("ASV_ID").grouping_analysis.nunique() if len(indicators) else pd.Series(dtype=float),
        "number_of_linked_mags": mappings.groupby("ASV_ID").genome_id.nunique() if len(mappings) else pd.Series(dtype=float),
    }
    for name, series in summaries.items():
        base = base.merge(series.rename(name), on="ASV_ID", how="left")
    net_cols = [c for c in ["ASV_ID", "module_label", "Degree", "Betweenness", "Closeness", "EigenCentral"] if c in network]
    base = base.merge(network[net_cols].drop_duplicates("ASV_ID"), on="ASV_ID", how="left")
    base["in_spieceasi_network"] = base.module_label.notna() if "module_label" in base else False
    base["has_mag_link"] = base.number_of_linked_mags.fillna(0) > 0
    if len(functional):
        key = "mag_match_id" if "mag_match_id" in functional and "mag_match_id" in mappings else "genome_id"
        function_counts = functional.loc[functional.present.astype(str).str.lower().eq("true")].groupby(key).module_id.nunique()
        map_fun = mappings[["ASV_ID", key]].drop_duplicates().merge(function_counts.rename("n"), on=key, how="left").groupby("ASV_ID").n.sum()
        base = base.merge(map_fun.rename("number_of_present_mag_function_modules"), on="ASV_ID", how="left")
    return base.fillna({c: 0 for c in base.columns if c.startswith("number_of_")})


def savefig(fig: plt.Figure, path_base: Path, formats: list[str]) -> None:
    path_base.parent.mkdir(parents=True, exist_ok=True)
    for fmt in formats:
        fig.savefig(path_base.with_suffix(f".{fmt}"), bbox_inches="tight", dpi=300)
    plt.close(fig)


def plot_network_overlays(graph_path: Path, layout: pd.DataFrame, response: pd.DataFrame,
                          variables: list[str], outdir: Path, formats: list[str]) -> None:
    graph = nx.read_graphml(graph_path)
    id_to_asv = {str(row.get("graph_node_id", row.get("GraphML_ID", ""))): row.ASV_ID for _, row in layout.iterrows()}
    asv_pos = {row.ASV_ID: (row.spieceasi_x, row.spieceasi_y) for _, row in layout.dropna(subset=["spieceasi_x", "spieceasi_y"]).iterrows()}
    edges = [(id_to_asv.get(str(a), str(a)), id_to_asv.get(str(b), str(b))) for a, b in graph.edges()]
    for variable in variables:
        data = response.loc[response.environmental_variable.astype(str).str.casefold().eq(variable.casefold())].set_index("ASV_ID")
        if data.empty:
            continue
        fig, ax = plt.subplots(figsize=(13, 9))
        for a, b in edges:
            if a in asv_pos and b in asv_pos:
                ax.plot([asv_pos[a][0], asv_pos[b][0]], [asv_pos[a][1], asv_pos[b][1]], color="#D0D0D0", lw=.55, zorder=1)
        absent = [a for a in asv_pos if a not in data.index or str(data.loc[a].get("response_direction", "")) not in TITAN_DIRECTION_COLORS]
        ax.scatter([asv_pos[a][0] for a in absent], [asv_pos[a][1] for a in absent], s=20, c="#F2F2F2", edgecolors="#D9D9D9", linewidths=.25, zorder=2)
        for direction in ["z-", "z+"]:
            ids = [a for a in asv_pos if a in data.index and data.loc[a].get("response_direction") == direction]
            z = np.array([abs(float(data.loc[a].get("z_score", 0) or 0)) for a in ids])
            ax.scatter(
                [asv_pos[a][0] for a in ids], [asv_pos[a][1] for a in ids],
                s=[titan_marker_area(value) for value in z],
                c=TITAN_DIRECTION_COLORS[direction], marker="o",
                edgecolors="black", linewidths=.55, zorder=3,
            )
        direction_handles = [
            Line2D([], [], marker="o", linestyle="", markerfacecolor=TITAN_DIRECTION_COLORS[direction],
                   markeredgecolor="black", markersize=10,
                   label=f"{direction}: {'declining' if direction == 'z-' else 'increasing'} response")
            for direction in ["z-", "z+"]
        ]
        size_handles = [
            Line2D([], [], marker="o", linestyle="", markerfacecolor="white", markeredgecolor="black",
                   markersize=np.sqrt(titan_marker_area(value)), label=f"|z| = {value:g}")
            for value in TITAN_Z_LEGEND_VALUES
        ]
        absent_handle = Line2D([], [], marker="o", linestyle="", markerfacecolor="#F2F2F2",
                               markeredgecolor="#D9D9D9", markersize=6, label="No TITAN response")
        ax.set_title(f"SPIEC-EASI network: TITAN response to {variable}")
        ax.set_axis_off()
        ax.legend(
            handles=direction_handles + ([absent_handle] if absent else []) + size_handles,
            frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0),
            title="TITAN response and strength",
        )
        fig.tight_layout()
        savefig(fig, outdir / f"spieceasi_titan_{re.sub(r'[^A-Za-z0-9]+', '_', variable).strip('_').lower()}", formats)


def plot_titan_heterogeneous_overlays(
    graph_path: Path,
    response: pd.DataFrame,
    variables: list[str],
    outdir: Path,
    prefix: str,
    formats: list[str],
) -> None:
    """Overlay each TITAN gradient on the established ASV--MAG network layout."""
    if not graph_path.is_file():
        raise FileNotFoundError(f"Heterogeneous ASV--MAG GraphML not found: {graph_path}")
    graph = nx.read_graphml(graph_path)
    (outdir / "network").mkdir(parents=True, exist_ok=True)
    for variable in variables:
        data = response.loc[
            response.environmental_variable.astype(str).str.casefold().eq(variable.casefold())
        ].copy()
        if data.empty:
            continue
        data["abs_z"] = pd.to_numeric(data.get("z_score"), errors="coerce").abs()
        data = data.sort_values("abs_z", ascending=False).drop_duplicates("ASV_ID")
        titan_responses = {
            clean_asv(row.ASV_ID): {
                "response_direction": str(row.response_direction),
                "z_score": row.z_score,
            }
            for row in data.itertuples(index=False)
            if str(row.response_direction) in TITAN_DIRECTION_COLORS
        }
        plot_heterogeneous_network(
            graph=graph,
            outdir=outdir,
            prefix=prefix,
            formats=formats,
            overlay_mode="titan",
            titan_responses=titan_responses,
            titan_variable=variable,
        )


def plot_context_sheets(sample: pd.DataFrame, response: pd.DataFrame, selected: list[str],
                        context_col: str, depth_col: str, date_col: str, outdir: Path, formats: list[str]) -> None:
    for asv in selected:
        data = sample.loc[sample.ASV_ID.eq(asv)].copy()
        if data.empty:
            continue
        data[depth_col] = pd.to_numeric(data.get(depth_col), errors="coerce")
        data[date_col] = pd.to_datetime(data.get(date_col), errors="coerce")
        fig, axes = plt.subplots(2, 2, figsize=(15, 11))
        fig._aspire_compact_publication_typography = True
        size = 15 + 900 * np.sqrt(data.relative_abundance.clip(lower=0))
        axes[0, 0].scatter(data[date_col], data[depth_col], s=size, c=data.relative_abundance, cmap="viridis", edgecolors="black", linewidths=.25)
        axes[0, 0].invert_yaxis(); axes[0, 0].set_xlabel("Date"); axes[0, 0].set_ylabel("Depth (m)"); axes[0, 0].set_title("Occurrence and relative abundance", fontsize=16, pad=14)
        if context_col in data:
            shown = data.loc[~data[context_col].astype(str).str.contains(r"(?i)(?:outlier|other)", regex=True, na=False)].copy()
            summary = shown.groupby(context_col, observed=True).relative_abundance.agg(["mean", "median"]).sort_values("mean")
            labels = summary.index.astype(str).str.replace("__", "-", regex=False)
            axes[0, 1].barh(labels, summary["mean"], color="#777777")
            axes[0, 1].set_xlabel("Mean relative abundance"); axes[0, 1].set_title("Hybrid compartment", fontsize=16, pad=14)
        resp = response.loc[response.ASV_ID.eq(asv)].copy()
        resp = resp.sort_values("z_score", key=lambda x: pd.to_numeric(x, errors="coerce").abs(), ascending=False).head(12)
        colors = resp.response_direction.map({"z+": "#B2182B", "z-": "#2166AC"}).fillna("#999999")
        axes[1, 0].barh(resp.environmental_variable.astype(str), pd.to_numeric(resp.z_score, errors="coerce").fillna(0), color=colors)
        axes[1, 0].invert_yaxis(); axes[1, 0].set_xlabel("TITAN z score"); axes[1, 0].set_title("Continuous-gradient threshold responses", fontsize=16, pad=14)
        rho = pd.to_numeric(resp.spearman_rho, errors="coerce")
        axes[1, 1].axvline(0, color="#777777", lw=.8, ls="--")
        axes[1, 1].scatter(rho, np.arange(len(resp)), c=colors, edgecolors="black", linewidths=.4)
        axes[1, 1].set_yticks(np.arange(len(resp)), resp.environmental_variable.astype(str)); axes[1, 1].invert_yaxis()
        axes[1, 1].set_xlim(-1, 1); axes[1, 1].set_xlabel("Spearman ρ"); axes[1, 1].set_title("Monotonic association and TITAN direction", fontsize=16, pad=14)
        for ax in axes.flat:
            ax.tick_params(labelsize=10); ax.xaxis.label.set_fontsize(12); ax.yaxis.label.set_fontsize(12)
        rank_values = []
        for col, prefix in [("phylum", "P"), ("family", "F"), ("genus", "G")]:
            if col in data and data[col].notna().any():
                value = str(data[col].dropna().iloc[0]).replace("_", " ").strip()
                if value: rank_values.append(f"{prefix}: {value}")
        tax_label = "; ".join(rank_values) if rank_values else "taxonomy unavailable"
        fig.text(.5, .995, f"{asv} — {tax_label}", ha="center", va="top", fontsize=15)
        fig.tight_layout(rect=[0, 0, 1, .97])
        savefig(fig, outdir / f"{asv}_ecological_context", formats)


def plot_heterogeneous_subnetworks(mappings: pd.DataFrame, targets: pd.DataFrame,
                                   selected: list[str], outdir: Path, formats: list[str]) -> None:
    """Render readable ASV--genome-lineage--target-function summaries.

    Genome records remain uncollapsed in the atlas tables.  The figure collapses
    clonal/near-clonal records to species and records their count in the label.
    """
    if targets.empty:
        return
    key = "mag_match_id" if "mag_match_id" in targets and "mag_match_id" in mappings else "genome_id"
    keep = [c for c in ["ASV_ID", key, "mag_tax_species", "mag_tax_genus", "asv_link_weight"] if c in mappings]
    joined = mappings[keep].drop_duplicates().merge(targets, on=key, how="inner", suffixes=("", "_target"))
    joined = joined.loc[joined.present.astype(str).str.lower().eq("true")].copy()
    for asv in selected:
        data = joined.loc[joined.ASV_ID.eq(asv)].copy()
        if data.empty:
            continue
        species_col = "mag_tax_species" if "mag_tax_species" in data else "mag_tax_genus"
        data["lineage"] = data[species_col].fillna("").astype(str).str.strip()
        if "mag_tax_genus" in data:
            data.loc[data.lineage.eq(""), "lineage"] = data.loc[data.lineage.eq(""), "mag_tax_genus"].fillna("unclassified")
        data.loc[data.lineage.eq(""), "lineage"] = "unclassified"
        counts = data[["lineage", key]].drop_duplicates().groupby("lineage")[key].nunique()
        modules = data[["module_id", "module_name"]].drop_duplicates().sort_values("module_id")
        lineages = sorted(data.lineage.unique())
        graph = nx.Graph(); graph.add_node(asv, kind="asv")
        for lineage in lineages:
            label = f"{lineage} (n={int(counts[lineage])})"; graph.add_node(label, kind="lineage"); graph.add_edge(asv, label, weight=1)
        for row in modules.itertuples(index=False): graph.add_node(row.module_id, kind="function", label=row.module_name)
        for (lineage, module), group in data.groupby(["lineage", "module_id"]):
            label = f"{lineage} (n={int(counts[lineage])})"
            fraction = group[key].nunique() / counts[lineage]
            graph.add_edge(label, module, weight=fraction)
        def evenly_spaced(n: int) -> np.ndarray:
            return np.linspace(1.0, -1.0, n) if n > 1 else np.array([0.0])
        pos = {asv: (-2.0, 0.0)}
        for lineage, y in zip(lineages, evenly_spaced(len(lineages))):
            pos[f"{lineage} (n={int(counts[lineage])})"] = (-0.55, float(y))
        for module, y in zip(modules.module_id, evenly_spaced(len(modules))):
            pos[module] = (1.15, float(y))
        fig, (ax, legend_ax) = plt.subplots(1, 2, figsize=(15, max(8, .55 * len(modules) + 4)), gridspec_kw={"width_ratios": [1.45, 1.0]})
        widths = [.5 + 2 * graph.edges[e].get("weight", 1) for e in graph.edges]
        nx.draw_networkx_edges(graph, pos, ax=ax, edge_color="#A0A0A0", width=widths, alpha=.8)
        for kind, color, shape, size in [("asv", "#FFFFFF", "o", 850), ("lineage", "#8DA0CB", "s", 320), ("function", "#66C2A5", "D", 210)]:
            nodes = [n for n, d in graph.nodes(data=True) if d.get("kind") == kind]
            nx.draw_networkx_nodes(graph, pos, nodelist=nodes, node_color=color, node_shape=shape, node_size=size, edgecolors="black", linewidths=.8, ax=ax)
        nx.draw_networkx_labels(graph, pos, labels={asv: asv}, font_size=9, ax=ax)
        for node, attrs in graph.nodes(data=True):
            x, y = pos[node]
            if attrs.get("kind") == "lineage": ax.text(x - .08, y, node, ha="right", va="center", fontsize=7)
            elif attrs.get("kind") == "function": ax.text(x + .08, y, node, ha="left", va="center", fontsize=7)
        module_legend = "\n".join(f"{r.module_id}: {r.module_name}" for r in modules.itertuples(index=False))
        legend_ax.text(0.0, .5, module_legend, transform=legend_ax.transAxes, va="center", fontsize=7)
        legend_ax.set_axis_off()
        ax.set_title(f"{asv}: linked genome lineages and present target N/S modules")
        ax.set_axis_off(); fig.tight_layout()
        savefig(fig, outdir / f"{asv}_mag_target_function_subnetwork", formats)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--relative-abundance", required=True); p.add_argument("--metadata", required=True)
    p.add_argument("--microbial-compartments", required=True); p.add_argument("--taxonomy", required=True)
    p.add_argument("--modules", required=True); p.add_argument("--node-features", required=True); p.add_argument("--network-layout", required=True); p.add_argument("--network-graph", required=True)
    p.add_argument("--titan-results", required=True); p.add_argument("--correlations", required=True); p.add_argument("--indicator-dir", required=True)
    p.add_argument("--module-compartment-associations", required=True)
    p.add_argument("--module-measurement-associations", required=True)
    p.add_argument("--module-measurement-member-support", required=True)
    p.add_argument("--mag-mappings", required=True); p.add_argument("--functional-modules", required=True); p.add_argument("--multiomics-dir", required=True)
    p.add_argument("--target-functional-modules", default="")
    p.add_argument("--heterogeneous-graph", required=True)
    p.add_argument("--heterogeneous-prefix", default="asv_mag_network_metagenome")
    p.add_argument("--rna-dna-log2-tpm-ratio")
    p.add_argument("--rna-dna-log2-tpm-se")
    p.add_argument("--rna-dna-max-log2-se", type=float, default=0.3)
    p.add_argument(
        "--rna-dna-mag-id-mode",
        choices=["exact", "suffix_after_double_underscore"], default="exact",
    )
    p.add_argument("--outdir", required=True); p.add_argument("--sample-col", default="sampleID"); p.add_argument("--cruise-col", default="Cruise"); p.add_argument("--depth-col", default="Depth"); p.add_argument("--date-col", default="Date")
    p.add_argument("--context-cols", default="o2_compartment,gmm_component,o2_subcompartment_final,microbial_compartment,Season,renewal_phase,cruise_group")
    p.add_argument("--network-variables", default="Oxygen,H2S,Nitrate,Ammonium")
    p.add_argument("--linkage-q-threshold", type=float, default=0.05)
    p.add_argument("--selected-asvs", default=""); p.add_argument("--context-sheet-source", choices=["explicit", "mag_linked"], default="mag_linked")
    p.add_argument("--formats", default="pdf,png,svg")
    args = p.parse_args()
    if not 0 < args.linkage_q_threshold <= 1:
        p.error("--linkage-q-threshold must be in (0, 1]")
    out = Path(args.outdir); tables = out / "tables"; audit = out / "audit"; plots = out / "plots"
    ra = read_table(Path(args.relative_abundance)); metadata = read_table(Path(args.metadata)); microbial = read_table(Path(args.microbial_compartments))
    tax = taxonomy_table(Path(args.taxonomy)); network = load_network(Path(args.modules), Path(args.node_features), Path(args.network_layout))
    titan = read_table(Path(args.titan_results)); correlations = read_table(Path(args.correlations))
    module_compartment = read_table(Path(args.module_compartment_associations))
    module_measurement = read_table(Path(args.module_measurement_associations))
    member_support = read_table(Path(args.module_measurement_member_support))
    mappings = read_table(Path(args.mag_mappings)); mappings["ASV_ID"] = mappings.ASV_ID.map(clean_asv)
    functional = read_table(Path(args.functional_modules))
    sample = make_sample_context(ra, metadata, microbial, tax, network, args.sample_col)
    response = environmental_response(sample, titan, correlations, args.cruise_col, args.depth_col, args.date_col)
    linkage, linkage_supported, linkage_mag, linkage_definitions = ecological_linkage_crosswalk(
        network, tax, response, module_compartment, module_measurement,
        member_support, mappings, args.linkage_q_threshold,
    )
    module_mag_crosswalk = asv_ecological_module_mag_crosswalk(
        network, tax, mappings,
    )
    indicators = indicator_evidence(Path(args.indicator_dir))
    contexts = context_profiles(sample, split_values(args.context_cols), args.cruise_col, args.depth_col, args.date_col)
    magfun = mag_function_links(mappings, functional, network, tax)
    multiomics = multiomics_context(Path(args.multiomics_dir), sample)
    paired = paired_multiomics(
        multiomics,
        log2_ratio_path=(
            Path(args.rna_dna_log2_tpm_ratio)
            if args.rna_dna_log2_tpm_ratio else None
        ),
        log2_se_path=(
            Path(args.rna_dna_log2_tpm_se) if args.rna_dna_log2_tpm_se else None
        ),
        mag_id_mode=args.rna_dna_mag_id_mode,
        max_log2_se=args.rna_dna_max_log2_se,
    )
    index = evidence_index(sample, response, indicators, network, mappings, functional)
    write_table(sample, tables / "asv_sample_context.tsv"); write_table(response, tables / "asv_environmental_response.tsv")
    write_table(contexts, tables / "asv_context_profiles.tsv"); write_table(indicators, tables / "asv_indicator_evidence.tsv")
    write_table(magfun, tables / "asv_mag_function_links.tsv"); write_table(multiomics, tables / "asv_mag_multiomics_context.tsv")
    write_table(paired, tables / "asv_mag_paired_dna_rna_context.tsv")
    write_table(index, tables / "asv_evidence_index.tsv")
    write_table(
        module_mag_crosswalk,
        tables / "asv_ecological_module_mag_crosswalk.tsv",
    )
    write_table(linkage, tables / "asv_module_compartment_measurement_crosswalk.tsv")
    write_table(linkage_supported, tables / "asv_module_compartment_measurement_significant.tsv")
    write_table(linkage_mag, tables / "asv_module_compartment_measurement_mag_linked.tsv")
    write_table(linkage_definitions, audit / "ecological_linkage_evidence_definitions.tsv")
    manifest_rows = []
    for path in sorted(tables.glob("*.tsv")):
        frame = read_table(path); manifest_rows.append({"table": path.name, "rows": len(frame), "columns": len(frame.columns), "primary_key_description": {"asv_sample_context.tsv":"ASV_ID + sample_ID", "asv_environmental_response.tsv":"ASV_ID + environmental_variable", "asv_context_profiles.tsv":"ASV_ID + context_type + context_level", "asv_indicator_evidence.tsv":"ASV_ID + grouping_analysis", "asv_mag_function_links.tsv":"ASV_ID + mag_match_id + module_id", "asv_mag_multiomics_context.tsv":"ASV_ID + mag_match_id + sample + data_modality", "asv_mag_paired_dna_rna_context.tsv":"ASV_ID + mag_match_id + sample", "asv_evidence_index.tsv":"ASV_ID", "asv_ecological_module_mag_crosswalk.tsv":"ASV_ID + ecological_module", "asv_module_compartment_measurement_crosswalk.tsv":"ASV_ID + ecological_module + compartment_strategy + environmental_variable", "asv_module_compartment_measurement_significant.tsv":"ASV_ID + ecological_module + compartment_strategy + environmental_variable", "asv_module_compartment_measurement_mag_linked.tsv":"ASV_ID + ecological_module + compartment_strategy + environmental_variable (MAG IDs aggregated within ASV)"}.get(path.name, "see schema")})
    write_table(pd.DataFrame(manifest_rows), audit / "ecological_context_atlas_manifest.tsv")
    params = vars(args).copy(); params.update({"python_version": platform.python_version(), "pandas_version": pd.__version__, "numpy_version": np.__version__, "scipy_version": scipy.__version__, "networkx_version": nx.__version__})
    audit.mkdir(parents=True, exist_ok=True); (audit / "ecological_context_atlas_parameters.json").write_text(json.dumps(params, indent=2) + "\n")
    formats = split_values(args.formats)
    plot_network_overlays(Path(args.network_graph), network, response, split_values(args.network_variables), plots / "network_overlays", formats)
    plot_titan_heterogeneous_overlays(
        Path(args.heterogeneous_graph), response, split_values(args.network_variables),
        plots, args.heterogeneous_prefix, formats,
    )
    selected = split_values(args.selected_asvs)
    if not selected and args.context_sheet_source == "mag_linked": selected = sorted(mappings.ASV_ID.dropna().unique())
    plot_context_sheets(sample, response, selected, "o2_subcompartment_final", args.depth_col, args.date_col, plots / "asv_context_sheets", formats)
    target_path = Path(args.target_functional_modules) if args.target_functional_modules else None
    targets = read_table(target_path) if target_path and target_path.exists() else pd.DataFrame()
    plot_heterogeneous_subnetworks(mappings, targets, selected, plots / "heterogeneous_subnetworks", formats)


if __name__ == "__main__":
    main()
