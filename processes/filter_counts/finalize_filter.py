#!/usr/bin/env python3
"""Finalize combined filtering with nonzero ASVs, matching FASTA and an audit.

No sample-depth threshold is applied here: control decontamination has already
selected the biological cohort. Intermediate tables remain unchanged.
"""
import argparse
from pathlib import Path

import pandas as pd
from Bio import SeqIO


def finalize(input_counts, table_counts, counts_path, fasta, output_fasta):
    frames = [pd.read_csv(p, sep="\t", index_col=0) for p in
              (input_counts, table_counts, counts_path)]
    final = frames[-1].loc[frames[-1].sum(axis=1) > 0]
    frames.append(final)
    steps = ["input", "initial_abundance_nonzero", "group_nontarget_abundance_taxonomy", "final_nonzero"]
    audit = []
    removed = []
    for i, (step, frame) in enumerate(zip(steps, frames)):
        row = dict(step=step, samples=frame.shape[1], asvs=frame.shape[0],
                   reads=int(frame.to_numpy().sum()))
        audit.append(row)
        print(row)
        if i:
            removed.extend(dict(ASV_ID=asv, removal_step=step)
                           for asv in frames[i - 1].index.difference(frame.index))
    wanted = set(final.index.astype(str))
    seen = set()
    with open(output_fasta, "w") as handle:
        for record in SeqIO.parse(fasta, "fasta"):
            asv = record.id.split(";")[0]
            if asv in wanted:
                if asv in seen:
                    raise ValueError(f"Duplicate FASTA ASV: {asv}")
                seen.add(asv)
                SeqIO.write(record, handle, "fasta")
    if seen != wanted:
        raise ValueError(f"Final ASVs missing from FASTA: {sorted(wanted - seen)}")
    final.to_csv(counts_path, sep="\t")
    output_dir = Path(counts_path).parent
    pd.DataFrame(audit).to_csv(output_dir / "filter_audit.tsv", sep="\t", index=False)
    pd.DataFrame(removed, columns=["ASV_ID", "removal_step"]).to_csv(
        output_dir / "filter_removed_asvs.tsv", sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("input-counts", "table-counts", "counts", "fasta", "output-fasta"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    finalize(args.input_counts, args.table_counts, args.counts, args.fasta, args.output_fasta)
