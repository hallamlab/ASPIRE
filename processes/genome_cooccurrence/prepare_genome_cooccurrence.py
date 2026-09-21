#!/usr/bin/env python3
"""Prepare matched species-level raw recruitment matrices for co-occurrence analysis.

The representative selection intentionally mirrors BASIN's deterministic
species representative hierarchy without consuming BASIN output.  This keeps
ASPIRE independently reproducible while ensuring both workflows select the
same genome from the shared, versioned genome-quality table.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


def read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t" if path.suffix.lower() in {".tsv", ".tab"} else ",", low_memory=False)


def genome_key(value: object) -> str:
    text = str(value).strip().lower()
    text = re.sub(r"\.(?:fna|fa|fasta|faa|gz)$", "", text)
    return re.sub(r"[^a-z0-9]+", "_", text).strip("_")


def species_name(value: object) -> str:
    """Normalize a species-quantification label to its terminal GTDB species."""
    text = str(value).strip()
    if "|" in text:
        text = text.rsplit("|", 1)[-1].strip()
    return re.sub(r"\s+", " ", text)


def truthy(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"true", "t", "1", "yes"})


def select_representatives(qc: pd.DataFrame, expected_n: int) -> pd.DataFrame:
    required = {
        "Genome_Id", "Species", "Phylum", "Class", "Order", "Family", "Genus",
        "mimag_tier", "sum_len", "Completeness", "Contamination", "N50",
    }
    missing = required - set(qc.columns)
    if missing:
        raise ValueError(f"Genome quality table lacks columns: {sorted(missing)}")
    work = qc.copy()
    work["genome_key"] = work["Genome_Id"].map(genome_key)
    work["Species"] = work["Species"].astype(str).str.strip()
    work["mimag_tier"] = work["mimag_tier"].astype(str).str.strip().str.lower()
    work = work.loc[
        work["Species"].ne("")
        & ~work["Species"].str.lower().isin({"nan", "none", "unclassified"})
        & work["mimag_tier"].isin({"high", "medium"})
    ].copy()
    if "mag_in_final_fasta_set" in work:
        work = work.loc[truthy(work["mag_in_final_fasta_set"])]
    elif "fasta_path" in work:
        work = work.loc[work["fasta_path"].notna() & work["fasta_path"].astype(str).str.strip().ne("")]
    for column in ("sum_len", "Completeness", "Contamination", "N50"):
        work[column] = pd.to_numeric(work[column], errors="coerce")
    work["mimag_priority"] = work["mimag_tier"].map({"high": 0, "medium": 1})
    work = work.sort_values(
        ["Species", "mimag_priority", "sum_len", "Completeness", "Contamination", "N50", "Genome_Id"],
        ascending=[True, True, False, False, True, False, True],
        kind="mergesort",
    )
    work["selection_rank_within_species"] = work.groupby("Species").cumcount() + 1
    work["selected_species_representative"] = work["selection_rank_within_species"].eq(1)
    work["selection_hierarchy"] = (
        "mimag_tier(high>medium)>sum_len>Completeness>Contamination(lower)>N50>genome_id"
    )
    selected_n = int(work["selected_species_representative"].sum())
    if expected_n > 0 and selected_n != expected_n:
        raise ValueError(
            f"Expected {expected_n} species representatives but deterministic selection produced {selected_n}"
        )
    return work


def prepare_modality(
    abundance: pd.DataFrame,
    selected: pd.DataFrame,
    *,
    sample_col: str,
    genome_col: str,
    value_col: str,
    closure_total: float,
    input_scale: str,
    input_feature_level: str,
    label: str,
    outdir: Path,
    minimum_read_count: float = 0.0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    required = {sample_col, genome_col, value_col}
    missing = required - set(abundance.columns)
    if missing:
        raise ValueError(f"{label} abundance table lacks columns: {sorted(missing)}")
    work = abundance[[sample_col, genome_col, value_col]].copy()
    work[value_col] = pd.to_numeric(work[value_col], errors="coerce")
    negative_n = int(work[value_col].lt(0).fillna(False).sum())
    if negative_n:
        raise ValueError(f"{label} abundance contains {negative_n} negative values")
    if minimum_read_count < 0:
        raise ValueError("minimum_read_count cannot be negative")
    if minimum_read_count > 0 and input_scale != "raw_counts":
        raise ValueError("minimum_read_count is defined only for raw-count input")
    selected_rows = selected.loc[selected["selected_species_representative"]].copy()
    selected_rows["species_key"] = selected_rows["Species"].map(species_name)
    selected_keys = selected_rows["species_key"].tolist()
    id_lookup = selected_rows.set_index("species_key")["Genome_Id"].to_dict()
    if input_feature_level == "genome":
        # Build a complete genome-by-sample grid. Missing long-table rows are
        # true zero recruitment, not genomes to omit from the species mean.
        eligible = selected.copy()
        eligible["species_key"] = eligible["Species"].map(species_name)
        eligible = eligible.loc[eligible["species_key"].isin(selected_keys)].copy()
        genome_to_species = eligible.drop_duplicates("genome_key").set_index("genome_key")["species_key"]
        work["genome_key"] = work[genome_col].map(genome_key)
        work = work.loc[work["genome_key"].isin(genome_to_species.index)].copy()
        input_genomes = set(work["genome_key"])
        missing_genomes = sorted(set(genome_to_species.index) - input_genomes)
        if missing_genomes:
            raise ValueError(
                f"{label} lacks {len(missing_genomes)} eligible genome features: "
                + ", ".join(missing_genomes[:10])
            )
        samples = work[sample_col].dropna().astype(str).drop_duplicates().tolist()
        work[sample_col] = work[sample_col].astype(str)
        genome_sample = work.pivot_table(
            index="genome_key", columns=sample_col, values=value_col,
            aggfunc="sum", fill_value=0.0,
        ).reindex(index=genome_to_species.index.tolist(), columns=samples, fill_value=0.0)
        genome_sample = genome_sample.apply(pd.to_numeric, errors="coerce").fillna(0.0)
        genome_before_floor = genome_sample.copy()
        below_floor = genome_sample.gt(0) & genome_sample.lt(minimum_read_count)
        genome_sample = genome_sample.mask(below_floor, 0.0)
        genome_long = genome_before_floor.rename_axis("genome_key").reset_index().melt(
            id_vars="genome_key", var_name="sample_id", value_name="read_count_before_floor"
        )
        retained_long = genome_sample.rename_axis("genome_key").reset_index().melt(
            id_vars="genome_key", var_name="sample_id", value_name="read_count_after_floor"
        )
        threshold_audit = genome_long.merge(
            retained_long, on=["genome_key", "sample_id"], how="left", validate="one_to_one"
        )
        threshold_audit["species_key"] = threshold_audit["genome_key"].map(genome_to_species)
        threshold_audit["minimum_read_count"] = minimum_read_count
        threshold_audit["set_to_zero_by_read_floor"] = (
            threshold_audit["read_count_before_floor"].gt(0)
            & threshold_audit["read_count_after_floor"].eq(0)
        )
        threshold_audit.to_csv(
            outdir / f"{label}_genome_sample_read_floor_audit.tsv", sep="\t", index=False
        )
        # Arithmetic means include all eligible genomes for the species,
        # including zero-recruitment genomes in a sample.
        pivot_before_floor = genome_before_floor.groupby(genome_to_species, sort=False).mean()
        pivot = genome_sample.groupby(genome_to_species, sort=False).mean()
        pivot = pivot.reindex(selected_keys, fill_value=0.0)
        floor_cells_by_sample = below_floor.sum(axis=0)
    elif input_feature_level == "species":
        work["feature_key"] = work[genome_col].map(species_name)
        if minimum_read_count > 0:
            raise ValueError(
                "A per-genome minimum_read_count requires genome-level input; "
                "a pre-aggregated species table cannot audit the genome-level threshold"
            )
        work = work.loc[work["feature_key"].isin(selected_keys)].copy()
        pivot = work.pivot_table(
            index="feature_key", columns=sample_col, values=value_col,
            aggfunc="sum", fill_value=0.0,
        ).reindex(selected_keys, fill_value=0.0)
        pivot = pivot.apply(pd.to_numeric, errors="coerce").fillna(0.0)
        pivot_before_floor = pivot.copy()
        floor_cells_by_sample = pd.Series(0, index=pivot.columns, dtype=int)
    else:
        raise ValueError(f"Unsupported input feature level: {input_feature_level}")
    sample_totals_before_floor = pivot_before_floor.sum(axis=0)
    sample_totals = pivot.sum(axis=0)
    usable_samples = sample_totals.loc[sample_totals.gt(0)].index
    pivot = pivot.loc[:, usable_samples]
    if pivot.shape[1] < 5:
        raise ValueError(f"{label} has fewer than five usable samples across selected representatives")
    if input_scale == "raw_counts":
        prepared = pivot.astype(float)
    elif input_scale == "closed_abundance":
        prepared = pivot.divide(pivot.sum(axis=0), axis=1) * float(closure_total)
    else:
        raise ValueError(f"Unsupported input scale: {input_scale}")
    prepared.index = [id_lookup[key] for key in prepared.index]
    prepared.index.name = "Feature_ID"
    prepared.reset_index().to_csv(outdir / f"{label}_species_matrix.tsv", sep="\t", index=False)

    audit = pd.DataFrame({
        "modality": label,
        "sample_id": sample_totals.index.astype(str),
        "selected_species_abundance_sum_before_read_floor": sample_totals_before_floor.to_numpy(float),
        "genome_sample_cells_set_to_zero": floor_cells_by_sample.reindex(sample_totals.index).to_numpy(int),
        "minimum_read_count": minimum_read_count,
        "selected_species_abundance_sum_before_closure": sample_totals.to_numpy(float),
        "included_in_network": sample_totals.gt(0).to_numpy(bool),
        "exclusion_reason": np.where(sample_totals.gt(0), "", "zero abundance across selected species"),
        "input_scale": input_scale,
        "input_feature_level": input_feature_level,
        "closure_total": float(closure_total) if input_scale == "closed_abundance" else np.nan,
        "selected_species_n": pivot.shape[0],
    })
    prepared_sums = prepared.sum(axis=0).to_dict()
    audit["prepared_matrix_sum"] = audit["sample_id"].map(prepared_sums)
    return prepared, audit


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome-qc", required=True, type=Path)
    parser.add_argument("--metagenome-abundance", required=True, type=Path)
    parser.add_argument("--metatranscriptome-abundance", required=True, type=Path)
    parser.add_argument("--asv-mag-mappings", required=True, type=Path)
    parser.add_argument("--asv-modules", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--sample-col", default="sample_id")
    parser.add_argument("--genome-col", default="genome_id")
    parser.add_argument("--metagenome-value-col", default="read_count")
    parser.add_argument("--metatranscriptome-value-col", default="read_count")
    parser.add_argument("--input-scale", choices=["raw_counts", "closed_abundance"], default="raw_counts")
    parser.add_argument("--input-feature-level", choices=["genome", "species"], default="species")
    parser.add_argument("--expected-species", type=int, default=24)
    parser.add_argument("--closure-total", type=float, default=1_000_000.0)
    parser.add_argument("--minimum-read-count", type=float, default=0.0)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    selection = select_representatives(read_table(args.genome_qc), args.expected_species)
    selection.to_csv(args.outdir / "species_representative_selection_audit.tsv", sep="\t", index=False)
    selected = selection.loc[selection["selected_species_representative"]].copy()

    _, metag_audit = prepare_modality(
        read_table(args.metagenome_abundance), selection,
        sample_col=args.sample_col, genome_col=args.genome_col,
        value_col=args.metagenome_value_col, closure_total=args.closure_total, input_scale=args.input_scale,
        input_feature_level=args.input_feature_level,
        label="metagenome", outdir=args.outdir,
        minimum_read_count=args.minimum_read_count,
    )
    _, metat_audit = prepare_modality(
        read_table(args.metatranscriptome_abundance), selection,
        sample_col=args.sample_col, genome_col=args.genome_col,
        value_col=args.metatranscriptome_value_col, closure_total=args.closure_total, input_scale=args.input_scale,
        input_feature_level=args.input_feature_level,
        label="metatranscriptome", outdir=args.outdir,
        minimum_read_count=args.minimum_read_count,
    )
    pd.concat([metag_audit, metat_audit], ignore_index=True).to_csv(
        args.outdir / "sample_closure_audit.tsv", sep="\t", index=False
    )

    mappings = read_table(args.asv_mag_mappings)
    if "analysis_eligible_pair" in mappings:
        mappings = mappings.loc[truthy(mappings["analysis_eligible_pair"])].copy()
    mappings["genome_key"] = mappings["genome_id"].map(genome_key)
    mappings["ASV_ID"] = mappings["ASV_ID"].astype(str).str.replace(r";.*$", "", regex=True)
    modules = read_table(args.asv_modules)
    modules = modules[["Taxon", "module_label"]].rename(columns={"Taxon": "ASV_ID"})
    crosswalk = mappings.merge(modules, on="ASV_ID", how="left")
    crosswalk = crosswalk.loc[crosswalk["genome_key"].isin(set(selected["genome_key"]))].copy()
    crosswalk = selected[["Genome_Id", "genome_key", "Species"]].merge(
        crosswalk, on="genome_key", how="left", suffixes=("_selected", "")
    )
    crosswalk.to_csv(args.outdir / "representative_asv_crosswalk.tsv", sep="\t", index=False)

    grouped = crosswalk.groupby("genome_key", dropna=False).agg(
        linked_asvs=("ASV_ID", lambda s: ";".join(sorted(set(s.dropna().astype(str))))),
        linked_asvs_n=("ASV_ID", lambda s: s.dropna().astype(str).nunique()),
        propagated_ecological_modules=("module_label", lambda s: ";".join(sorted(set(s.dropna().astype(str))))),
        propagated_ecological_modules_n=("module_label", lambda s: s.dropna().astype(str).nunique()),
    ).reset_index()
    node_columns = [
        "Genome_Id", "genome_key", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species",
        "mimag_tier", "sum_len", "Completeness", "Contamination", "N50",
    ]
    node_metadata = selected[node_columns].merge(grouped, on="genome_key", how="left")
    node_metadata[["linked_asvs_n", "propagated_ecological_modules_n"]] = node_metadata[
        ["linked_asvs_n", "propagated_ecological_modules_n"]
    ].fillna(0).astype(int)
    node_metadata[["linked_asvs", "propagated_ecological_modules"]] = node_metadata[
        ["linked_asvs", "propagated_ecological_modules"]
    ].fillna("")
    node_metadata.to_csv(args.outdir / "representative_node_metadata.tsv", sep="\t", index=False)

    parameters = {
        "expected_species": args.expected_species,
        "selected_species": int(len(selected)),
        "closure_total": args.closure_total,
        "input_scale": args.input_scale,
        "input_feature_level": args.input_feature_level,
        "minimum_read_count_per_genome_sample": args.minimum_read_count,
        "species_sample_aggregation": (
            "arithmetic mean across all eligible genomes after applying the per-genome read-count floor; "
            "zero-recruitment genomes included in the denominator"
            if args.input_feature_level == "genome" else
            "pre-aggregated species input; duplicate species-sample rows summed"
        ),
        "zero_replacement_effective_relative_abundance": (
            1.0 / args.closure_total if args.input_scale == "closed_abundance" else None
        ),
        "zero_handling": "performed downstream by sample-wise multiplicative replacement before CLR transformation",
        "selection_hierarchy": selected["selection_hierarchy"].iloc[0],
        "network_interpretation": "signed proportionality among species-level relative recruitment profiles",
    }
    (args.outdir / "genome_cooccurrence_prepare_parameters.json").write_text(
        json.dumps(parameters, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
