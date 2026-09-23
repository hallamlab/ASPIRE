#!/usr/bin/env python3
"""Prepare the SPARK manuscript long/wide inputs from an ASPIRE final output.

The historical analysis intentionally used two related representations:
the positive-count long table (1,425 ASVs for SPARK) and the corrected wide
matrix (1,377 ASVs).  They must not be coerced to have the same ASV set.
"""

import argparse
from pathlib import Path

import pandas as pd


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--final-dir", required=True, type=Path)
    p.add_argument("--reference-wide", required=True, type=Path,
                   help="Deposited post-three-tier table (SD3); supplies the manuscript sample cohort")
    p.add_argument("--removed-asvs", required=True, type=Path,
                   help="Deposited three-tier audit table (SD1)")
    p.add_argument("--outdir", required=True, type=Path)
    args = p.parse_args()

    final_wide_path = args.final_dir / "ASVs" / "ASV_final.micro.tsv"
    final_long_path = args.final_dir / "metadata" / "ASV_meta_micro.tsv"
    for path in (final_wide_path, final_long_path, args.reference_wide, args.removed_asvs):
        if not path.is_file():
            raise SystemExit(f"Required input does not exist: {path}")

    deposited = pd.read_csv(args.reference_wide, sep="\t", low_memory=False)
    removed = pd.read_csv(args.removed_asvs, sep="\t", low_memory=False)
    wide = pd.read_csv(final_wide_path, sep="\t", low_memory=False)
    long = pd.read_csv(final_long_path, sep="\t", low_memory=False)
    wide = wide.rename(columns={wide.columns[0]: "ASV_ID"})
    deposited = deposited.rename(columns={deposited.columns[0]: "ASV_ID"})

    sample_ids = list(deposited.columns[1:])
    missing_samples = [x for x in sample_ids if x not in wide.columns]
    if missing_samples:
        raise SystemExit(f"Authoritative wide table lacks manuscript samples: {missing_samples}")

    removed_ids = set(removed["ASV_ID"].astype(str))
    deposited_ids = set(deposited["ASV_ID"].astype(str))
    final_ids = set(wide["ASV_ID"].astype(str))
    # Two SD1 ASVs occur only in the positive-count long table: they are rows
    # in ASPIRE's broader final matrix but are all zero in the 125-sample
    # manuscript cohort. Do not introduce those zero rows into the wide table.
    cohort_sums = wide.set_index("ASV_ID")[sample_ids].sum(axis=1)
    positive_removed_ids = removed_ids & set(cohort_sums[cohort_sums > 0].index.astype(str))
    pre_ids = deposited_ids | positive_removed_ids
    pre_wide = wide.loc[wide["ASV_ID"].astype(str).isin(pre_ids), ["ASV_ID"] + sample_ids].copy()

    expected_n = len(deposited_ids) + len(positive_removed_ids)
    if len(pre_wide) != expected_n:
        raise SystemExit(f"Could not reconstruct pre-tier wide matrix: got {len(pre_wide)}, expected {expected_n}")
    reapplied = pre_wide.loc[~pre_wide["ASV_ID"].astype(str).isin(removed_ids)].reset_index(drop=True)
    if not reapplied.equals(deposited.reset_index(drop=True)):
        raise SystemExit("Authoritative final output does not reproduce deposited SD3 after applying SD1")

    manuscript_long = long.loc[
        long["Sample"].astype(str).isin(sample_ids) & (pd.to_numeric(long["count"], errors="coerce") > 0)
    ].copy()
    if "corr_count" not in manuscript_long.columns:
        raise SystemExit("Authoritative long table lacks required corrected-count column: corr_count")
    # SPARK's downstream analyses use the corrected counts represented by SD3.
    # Keep raw-positive rows so the broader long-format taxonomy set (including
    # the two long-only SD1 ASVs) remains available to the three-tier audit.
    manuscript_long["count"] = pd.to_numeric(manuscript_long["corr_count"], errors="coerce").fillna(0).astype(int)
    # Names used by the deposited SPARK analysis drivers; retain ASPIRE names too.
    manuscript_long["lmp_id"] = manuscript_long["Sample"].astype(str)
    manuscript_long["sample"] = manuscript_long["Sample"].astype(str)
    manuscript_long["type_group"] = manuscript_long["Type_Group"]

    args.outdir.mkdir(parents=True, exist_ok=True)
    pre_wide.to_csv(args.outdir / "ASV_master_count_wide.tsv", sep="\t", index=False)
    manuscript_long.to_csv(args.outdir / "ASV_master_long.tsv", sep="\t", index=False)

    summary = [
        f"source_final_dir\t{args.final_dir.resolve()}",
        f"samples\t{len(sample_ids)}",
        f"long_asvs\t{manuscript_long['ASV_ID'].nunique()}",
        f"long_rows\t{len(manuscript_long)}",
        f"long_reads\t{int(manuscript_long['count'].sum())}",
        f"wide_asvs\t{len(pre_wide)}",
        f"wide_reads\t{int(pre_wide.iloc[:, 1:].to_numpy().sum())}",
        f"sd1_ids_in_wide\t{len(removed_ids & set(pre_wide['ASV_ID'].astype(str)))}",
    ]
    (args.outdir / "preparation_summary.tsv").write_text("metric\tvalue\n" + "\n".join(summary) + "\n")


if __name__ == "__main__":
    main()
