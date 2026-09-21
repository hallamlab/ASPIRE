#!/usr/bin/env python3
"""Prepare a strictly ASV-derived CLR/Aitchison cohort for microbial clustering."""

from __future__ import annotations

import argparse
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd


def read_table(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix.lower() in {".tsv", ".tab", ".txt"} else ","
    return pd.read_csv(path, sep=sep, low_memory=False)


def truthy(values: pd.Series) -> pd.Series:
    return values.astype(str).str.strip().str.lower().isin({"true", "t", "1", "yes", "y"})


def normalize_asv_ids(values: pd.Series | pd.Index) -> pd.Index:
    return pd.Index(values.astype(str).str.strip().str.replace(r";.*$", "", regex=True))


def load_counts(path: Path, transpose: bool, identifier_col: str) -> pd.DataFrame:
    table = read_table(path)
    if table.empty or table.shape[1] < 3:
        raise ValueError("ASV table must contain identifiers and at least two samples/features")
    identifier = identifier_col if identifier_col in table.columns else table.columns[0]
    ids = normalize_asv_ids(table[identifier])
    values = table.drop(columns=identifier).apply(pd.to_numeric, errors="coerce")
    if values.isna().any().any():
        raise ValueError("ASV table contains missing or nonnumeric abundance values")
    if (values < 0).any().any():
        raise ValueError("ASV abundances cannot be negative")
    if transpose:
        values.index = ids
        matrix = values.T
    else:
        values.index = ids
        matrix = values
    matrix.index = matrix.index.astype(str).str.strip()
    matrix.columns = normalize_asv_ids(matrix.columns)
    if matrix.index.duplicated().any():
        raise ValueError("ASV table contains duplicate sample identifiers")
    if matrix.columns.duplicated().any():
        raise ValueError("ASV table contains duplicate ASV identifiers")
    return matrix.astype(float)


def multiplicative_replace(composition: np.ndarray, delta: float) -> tuple[np.ndarray, np.ndarray]:
    replaced = np.empty_like(composition, dtype=float)
    zero_counts = (composition == 0).sum(axis=1)
    for row_index, row in enumerate(composition):
        zero_mask = row == 0
        zeros = int(zero_mask.sum())
        if zeros == 0:
            replaced[row_index] = row / row.sum()
            continue
        replacement_mass = zeros * delta
        if replacement_mass >= 1:
            raise ValueError(
                f"Multiplicative replacement is invalid for sample row {row_index}: "
                f"{zeros} zeros × delta {delta:g} is at least one"
            )
        nonzero_total = row[~zero_mask].sum()
        if nonzero_total <= 0:
            raise ValueError(f"Sample row {row_index} is all zero")
        replaced[row_index, zero_mask] = delta
        replaced[row_index, ~zero_mask] = row[~zero_mask] * (1.0 - replacement_mass) / nonzero_total
    return replaced, zero_counts


def pseudocount_replace(composition: np.ndarray, pseudocount: float) -> tuple[np.ndarray, np.ndarray]:
    if pseudocount <= 0:
        raise ValueError("Fixed pseudocount must be positive")
    zero_counts = (composition == 0).sum(axis=1)
    replaced = composition + pseudocount
    replaced /= replaced.sum(axis=1, keepdims=True)
    return replaced, zero_counts


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--counts", type=Path, required=True)
    parser.add_argument("--spieceasi-filter-audit", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--sample-col", default="sampleID")
    parser.add_argument(
        "--stability-block-col",
        default="",
        help="Metadata column defining profiles for within-block perturbation and blocked robustness validation",
    )
    parser.add_argument(
        "--stability-stratum-col",
        default="",
        help="Optional metadata column used to balance whole-block robustness resampling",
    )
    parser.add_argument("--asv-id-col", default="ASV_ID")
    parser.add_argument("--transpose", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--minimum-prevalence", type=float, default=0.0)
    parser.add_argument("--minimum-max-relative-abundance", type=float, default=0.0)
    parser.add_argument("--zero-replacement", choices=["multiplicative", "pseudocount"], default="multiplicative")
    parser.add_argument("--multiplicative-delta", type=float, default=0.0)
    parser.add_argument("--pseudocount", type=float, default=1e-6)
    args = parser.parse_args()

    if not 0 <= args.minimum_prevalence <= 1:
        raise ValueError("--minimum-prevalence must be between 0 and 1")
    if not 0 <= args.minimum_max_relative_abundance <= 1:
        raise ValueError("--minimum-max-relative-abundance must be between 0 and 1")
    if args.multiplicative_delta < 0:
        raise ValueError("--multiplicative-delta must be zero (automatic) or positive")

    args.outdir.mkdir(parents=True, exist_ok=True)
    counts = load_counts(args.counts, args.transpose, args.asv_id_col)
    metadata = read_table(args.metadata)
    if args.sample_col not in metadata.columns:
        raise ValueError(f"Metadata lacks sample identifier column: {args.sample_col}")
    metadata_ids = metadata[args.sample_col].astype(str).str.strip()
    if metadata_ids.duplicated().any():
        duplicate = metadata_ids[metadata_ids.duplicated()].unique().tolist()
        raise ValueError(f"Metadata sample identifiers are not unique: {duplicate[:10]}")
    metadata_set = set(metadata_ids)
    common_samples = [sample for sample in counts.index if sample in metadata_set]
    if len(common_samples) < 3:
        raise ValueError("Fewer than three samples match explicitly between ASV counts and metadata")

    sample_union = list(dict.fromkeys([*counts.index.tolist(), *metadata_ids.tolist()]))
    sample_audit = pd.DataFrame({args.sample_col: sample_union})
    sample_audit["present_in_asv_counts"] = sample_audit[args.sample_col].isin(counts.index)
    sample_audit["present_in_metadata"] = sample_audit[args.sample_col].isin(metadata_set)
    sample_audit["matched_by_explicit_identifier"] = (
        sample_audit.present_in_asv_counts & sample_audit.present_in_metadata
    )

    source_audit = read_table(args.spieceasi_filter_audit)
    if not {"ASV_ID", "retained_final"}.issubset(source_audit.columns):
        raise ValueError("SPIEC-EASI audit must contain ASV_ID and retained_final")
    source_audit = source_audit.copy()
    source_audit["ASV_ID"] = normalize_asv_ids(source_audit["ASV_ID"])
    retained_final = set(source_audit.loc[truthy(source_audit.retained_final), "ASV_ID"])
    missing = sorted(retained_final.difference(counts.columns))
    if missing:
        raise ValueError(f"SPIEC-EASI retained_final ASVs are absent from counts: {missing[:10]}")
    cohort_order = [asv for asv in counts.columns if asv in retained_final]
    if len(cohort_order) < 2:
        raise ValueError("Fewer than two SPIEC-EASI retained_final ASVs are available")

    cohort_counts = counts.loc[common_samples, cohort_order]
    initial_library = cohort_counts.sum(axis=1)
    nonzero_samples = initial_library.gt(0)
    sample_audit["retained_for_clustering"] = sample_audit[args.sample_col].isin(
        cohort_counts.index[nonzero_samples]
    )
    sample_audit["exclusion_reason"] = np.select(
        [
            ~sample_audit.matched_by_explicit_identifier,
            sample_audit.matched_by_explicit_identifier & ~sample_audit.retained_for_clustering,
        ],
        ["unmatched_sample_identifier", "zero_spieceasi_cohort_library"],
        default="retained",
    )
    cohort_counts = cohort_counts.loc[nonzero_samples]
    validation_blocks = pd.DataFrame({args.sample_col: cohort_counts.index})
    if args.stability_block_col:
        if args.stability_block_col not in metadata.columns:
            raise ValueError(
                f"Metadata lacks configured stability block column: {args.stability_block_col}"
            )
        block_lookup = metadata.assign(
            **{args.sample_col: metadata_ids}
        ).set_index(args.sample_col)[args.stability_block_col]
        validation_blocks["stability_block"] = validation_blocks[args.sample_col].map(block_lookup)
        missing_blocks = validation_blocks["stability_block"].isna() | validation_blocks[
            "stability_block"
        ].astype(str).str.strip().eq("")
        if missing_blocks.any():
            missing_ids = validation_blocks.loc[missing_blocks, args.sample_col].head(10).tolist()
            raise ValueError(
                "Configured stability blocks are missing for retained samples: "
                f"{missing_ids}"
            )
        validation_blocks["stability_block"] = (
            validation_blocks["stability_block"].astype(str).str.strip()
        )
        if validation_blocks["stability_block"].nunique() < 3:
            raise ValueError("Blocked stability validation requires at least three distinct blocks")
        validation_blocks["stability_block_source"] = args.stability_block_col
    else:
        validation_blocks["stability_block"] = validation_blocks[args.sample_col]
        validation_blocks["stability_block_source"] = "sample_ID"
    if args.stability_stratum_col:
        if args.stability_stratum_col not in metadata.columns:
            raise ValueError(
                f"Metadata lacks configured stability stratum column: {args.stability_stratum_col}"
            )
        stratum_lookup = metadata.assign(
            **{args.sample_col: metadata_ids}
        ).set_index(args.sample_col)[args.stability_stratum_col]
        validation_blocks["stability_stratum"] = validation_blocks[args.sample_col].map(
            stratum_lookup
        )
        missing_strata = validation_blocks["stability_stratum"].isna() | validation_blocks[
            "stability_stratum"
        ].astype(str).str.strip().eq("")
        if missing_strata.any():
            missing_ids = validation_blocks.loc[missing_strata, args.sample_col].head(10).tolist()
            raise ValueError(
                "Configured stability strata are missing for retained samples: "
                f"{missing_ids}"
            )
        validation_blocks["stability_stratum"] = (
            validation_blocks["stability_stratum"].astype(str).str.strip()
        )
        validation_blocks["stability_stratum_source"] = args.stability_stratum_col
        block_strata = validation_blocks.groupby("stability_block")["stability_stratum"].nunique()
        if block_strata.gt(1).any():
            invalid = block_strata[block_strata.gt(1)].index.tolist()[:10]
            raise ValueError(
                "Each stability block must belong to exactly one stratum; invalid blocks: "
                f"{invalid}"
            )
    else:
        validation_blocks["stability_stratum"] = "all"
        validation_blocks["stability_stratum_source"] = "unstratified"
    relative = cohort_counts.div(cohort_counts.sum(axis=1), axis=0)

    prevalence = relative.gt(0).mean(axis=0)
    max_abundance = relative.max(axis=0)
    nonzero_variance = relative.var(axis=0).gt(0)
    keep = (
        prevalence.ge(args.minimum_prevalence)
        & max_abundance.ge(args.minimum_max_relative_abundance)
    )
    feature_audit = source_audit.copy()
    feature_audit = feature_audit.merge(
        pd.DataFrame({
            "ASV_ID": cohort_order,
            "clustering_prevalence": prevalence.reindex(cohort_order).values,
            "clustering_max_relative_abundance": max_abundance.reindex(cohort_order).values,
            "nonzero_variance": nonzero_variance.reindex(cohort_order).values,
        }),
        on="ASV_ID",
        how="left",
    )
    keep_set = set(relative.columns[keep])
    feature_audit["retained_for_microbiome_clustering"] = feature_audit.ASV_ID.isin(keep_set)

    def exclusion_reason(row: pd.Series) -> str:
        if not truthy(pd.Series([row.retained_final])).iloc[0]:
            return "not_retained_final_by_spieceasi"
        reasons = []
        if pd.notna(row.clustering_prevalence) and row.clustering_prevalence < args.minimum_prevalence:
            reasons.append("prevalence_below_configured_minimum")
        if pd.notna(row.clustering_max_relative_abundance) and row.clustering_max_relative_abundance < args.minimum_max_relative_abundance:
            reasons.append("max_relative_abundance_below_configured_minimum")
        return ";".join(reasons) if reasons else "retained"

    feature_audit["microbiome_clustering_exclusion_reason"] = feature_audit.apply(exclusion_reason, axis=1)
    relative = relative.loc[:, [asv for asv in cohort_order if asv in keep_set]]
    if relative.shape[1] < 2:
        raise ValueError("Configured filters left fewer than two ASVs for clustering")

    dimensions = relative.shape[1]
    resolved_delta = args.multiplicative_delta if args.multiplicative_delta > 0 else 1.0 / dimensions**2
    if args.zero_replacement == "multiplicative":
        replaced_values, zero_counts = multiplicative_replace(relative.to_numpy(float), resolved_delta)
        replacement_value = resolved_delta
    else:
        replaced_values, zero_counts = pseudocount_replace(relative.to_numpy(float), args.pseudocount)
        replacement_value = args.pseudocount
    if not np.allclose(replaced_values.sum(axis=1), 1.0, atol=1e-10):
        raise RuntimeError("Zero-replaced compositions do not close to one")
    if (replaced_values <= 0).any():
        raise RuntimeError("Zero replacement left nonpositive values")

    log_values = np.log(replaced_values)
    clr_values = log_values - log_values.mean(axis=1, keepdims=True)
    if not np.isfinite(clr_values).all():
        raise RuntimeError("CLR transformation produced nonfinite values")
    differences = clr_values[:, np.newaxis, :] - clr_values[np.newaxis, :, :]
    aitchison = np.sqrt(np.einsum("ijk,ijk->ij", differences, differences))
    np.fill_diagonal(aitchison, 0.0)

    replaced = pd.DataFrame(replaced_values, index=relative.index, columns=relative.columns)
    clr = pd.DataFrame(clr_values, index=relative.index, columns=relative.columns)
    distance = pd.DataFrame(aitchison, index=relative.index, columns=relative.index)
    for table in (relative, replaced, clr, distance):
        table.index.name = args.sample_col
    relative.to_csv(args.outdir / "asv_relative_abundance.tsv", sep="\t")
    replaced.to_csv(args.outdir / "asv_zero_replaced_composition.tsv", sep="\t")
    clr.to_csv(args.outdir / "asv_clr_matrix.tsv", sep="\t")
    distance.to_csv(args.outdir / "aitchison_distance_matrix.tsv", sep="\t")
    validation_blocks.to_csv(
        args.outdir / "microbial_clustering_stability_blocks.tsv", sep="\t", index=False
    )
    sample_audit.to_csv(args.outdir / "microbial_clustering_sample_audit.tsv", sep="\t", index=False)
    feature_audit.to_csv(args.outdir / "microbial_clustering_asv_filter_audit.tsv", sep="\t", index=False)
    zero_audit = pd.DataFrame({
        args.sample_col: relative.index,
        "number_of_ASVs": dimensions,
        "number_of_zeros_replaced": zero_counts,
        "zero_fraction": zero_counts / dimensions,
        "replacement_method": args.zero_replacement,
        "resolved_replacement_value": replacement_value,
        "total_replacement_mass": zero_counts * replacement_value,
    })
    zero_audit.to_csv(args.outdir / "zero_replacement_audit.tsv", sep="\t", index=False)
    config = {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()}
    config.update({
        "matched_samples": len(relative),
        "spieceasi_retained_final_asvs": len(cohort_order),
        "clustering_asvs": dimensions,
        "resolved_multiplicative_delta": resolved_delta,
        "software_versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
        "inference_variables": ["ASV relative abundance only"],
        "excluded_from_inference": ["depth", "environmental measurements", "environmental compartments", "date"],
        "validation_design_variables": (
            [value for value in [args.stability_block_col, args.stability_stratum_col] if value]
            or ["sample_ID"]
        ),
        "validation_design_note": (
            "Blocks define within-profile sample perturbation and whole-profile robustness resampling; "
            "strata balance the latter. Neither is supplied to PAM or the Aitchison distance calculation."
        ),
    })
    (args.outdir / "microbial_clustering_preparation_config.json").write_text(json.dumps(config, indent=2))


if __name__ == "__main__":
    main()
