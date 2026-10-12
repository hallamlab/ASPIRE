#!/usr/bin/env python3
"""Project the audit long table onto the validated downstream wide matrix."""

import argparse
from pathlib import Path

import pandas as pd


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--long", required=True, type=Path)
    p.add_argument("--wide", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args()
    long = pd.read_csv(args.long, sep="\t", low_memory=False)
    wide = pd.read_csv(args.wide, sep="\t", low_memory=False).rename(columns=lambda c: "ASV_ID" if c == "#OTU ID" else c)
    samples = list(wide.columns[1:])
    ids = set(wide["ASV_ID"].astype(str))
    out = long.loc[
        long["ASV_ID"].astype(str).isin(ids)
        & long["sample"].astype(str).isin(samples)
        & (pd.to_numeric(long["count"], errors="coerce") > 0)
    ].copy()
    pivot = out.pivot_table(index="ASV_ID", columns="sample", values="count", aggfunc="sum", fill_value=0)
    expected = wide.set_index("ASV_ID")
    pivot = pivot.reindex(index=expected.index, columns=expected.columns, fill_value=0).astype(int)
    if not pivot.equals(expected.astype(int)):
        raise SystemExit("Projected downstream long table does not reproduce the validated wide matrix")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
