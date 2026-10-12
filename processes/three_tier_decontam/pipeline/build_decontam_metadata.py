#!/usr/bin/env python3
"""
build_decontam_metadata.py

Build the sample-level metadata table that run_decontam.R consumes.

Adds the following columns to the upstream Biofactorial metadata file:
  - is_negative_control (T/F): true for PBS, PBS_twz, Negative_96, Negative_man
  - is_positive_control (T/F): true for Positive_96, Positive_man
  - DNA_conc: numeric from upstream metadata (NA where missing)

Aligns sample IDs to those in ASV_counts.tsv (handles "-" vs "_" mismatches).
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


NEG_LABELS = {"PBS", "PBS_twz", "Negative-96", "Negative_96",
              "Negative-man", "Negative_man"}
POS_LABELS = {"Positive-96", "Positive_96", "Positive-man", "Positive_man"}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--counts", required=True, type=Path,
                    help="Raw ASV count table TSV (column 1 = ASV ID, remaining cols = samples)")
    ap.add_argument("--metadata-in", required=True, type=Path,
                    help="Upstream Biofactorial metadata TSV (metadata_updated_mito_raw.tsv)")
    ap.add_argument("--metadata-out", required=True, type=Path,
                    help="Output metadata TSV to feed into run_decontam.R")
    ap.add_argument("--sample-col-in", default="Sample",
                    help="Sample-ID column in upstream metadata [default: Sample]")
    ap.add_argument("--sample-col-out", default="Sample",
                    help="Sample-ID column written for downstream decontam scripts [default: Sample]")
    ap.add_argument("--negative-control-labels",
                    default=",".join(sorted(NEG_LABELS)),
                    help="Comma-separated sample IDs to flag as negative controls")
    ap.add_argument("--positive-control-labels",
                    default=",".join(sorted(POS_LABELS)),
                    help="Comma-separated sample IDs to flag as positive controls")
    ap.add_argument("--negative-control-col", default="is_negative_control",
                    help="Name of the generated negative-control column")
    ap.add_argument("--positive-control-col", default="is_positive_control",
                    help="Name of the generated positive-control column")
    ap.add_argument("--concentration-col", default="DNA_conc",
                    help="Input/output DNA-concentration column name")
    ap.add_argument("--type-col", default="Type_Group",
                    help="Sample-type column to retain in the generated metadata")
    args = ap.parse_args()

    # Read count-table sample IDs (column headers, excluding first column)
    counts_header = pd.read_csv(args.counts, sep="\t", nrows=0).columns.tolist()
    counts_samples = counts_header[1:]
    print(f"[INFO] Count table has {len(counts_samples)} sample columns", file=sys.stderr)

    # Read metadata
    meta = pd.read_csv(args.metadata_in, sep="\t", low_memory=False)
    print(f"[INFO] Metadata has {len(meta)} rows, columns: {list(meta.columns)[:8]}...",
          file=sys.stderr)
    for required_col in (args.sample_col_in, args.type_col):
        if required_col not in meta.columns:
            sys.exit(f"ERROR: required metadata column not found: {required_col}")

    # Build normalized lookup so we can match "Negative-96" <-> "Negative_96"
    def normalize(s):
        return str(s).replace("-", "_") if pd.notna(s) else s

    meta["_sample_norm"] = meta[args.sample_col_in].map(normalize)
    counts_samples_norm = [normalize(s) for s in counts_samples]

    # Reindex metadata to match count-table sample order
    meta_idx = meta.set_index("_sample_norm")
    matched = meta_idx.reindex(counts_samples_norm).reset_index(drop=True)
    matched_count = matched[args.sample_col_in].notna().sum()
    # Write the configured output sample column so it matches count-table column
    # names exactly. This may differ from the input metadata column name.
    # (handles any "-" vs "_" mismatches between metadata and count headers)
    matched[args.sample_col_out] = counts_samples

    # Diagnose
    print(f"[INFO] Matched {matched_count} of {len(counts_samples)} count-table samples to metadata",
          file=sys.stderr)

    # Build the control flags
    sample_normed = matched[args.sample_col_out].map(normalize)
    neg_labels = {x.strip() for x in args.negative_control_labels.split(",") if x.strip()}
    pos_labels = {x.strip() for x in args.positive_control_labels.split(",") if x.strip()}
    is_neg = sample_normed.isin({normalize(x) for x in neg_labels})
    is_pos = sample_normed.isin({normalize(x) for x in pos_labels})

    matched[args.negative_control_col] = is_neg.values
    matched[args.positive_control_col] = is_pos.values

    # Ensure DNA_conc is numeric (or empty for missing)
    if args.concentration_col in matched.columns:
        matched[args.concentration_col] = pd.to_numeric(
            matched[args.concentration_col], errors="coerce"
        )
    else:
        matched[args.concentration_col] = float("nan")

    print(f"[INFO] Negative controls flagged: {int(is_neg.sum())}", file=sys.stderr)
    print(f"[INFO] Positive controls flagged: {int(is_pos.sum())}", file=sys.stderr)
    print(f"[INFO] Samples with non-NA {args.concentration_col}: "
          f"{matched[args.concentration_col].notna().sum()}",
          file=sys.stderr)

    # Trim to just what run_decontam.R needs (keeps the file tidy)
    keep_cols = [args.sample_col_out, args.negative_control_col, args.positive_control_col,
                 args.concentration_col]
    for opt_col in [args.type_col, "Set", "Type", "Participant_ID", "Case"]:
        if opt_col in matched.columns:
            keep_cols.append(opt_col)
    keep_cols = list(dict.fromkeys(keep_cols))
    out = matched[keep_cols]

    args.metadata_out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.metadata_out, sep="\t", index=False)
    print(f"[INFO] Wrote {args.metadata_out} ({len(out)} rows, {len(out.columns)} cols)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
