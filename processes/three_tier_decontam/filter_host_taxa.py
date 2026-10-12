#!/usr/bin/env python3
"""Remove residual host-classified ASVs after the validated SPARK checkpoint."""

import argparse
from pathlib import Path

import pandas as pd


HOST_RULES = {
    "Domain": {"Eukaryota"},
    "Phylum": {"Chordata", "Vertebrata"},
    "Class": {"Mammalia"},
    "Species": {"Homo sapiens", "Homo_sapiens"},
}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--long", required=True, type=Path)
    p.add_argument("--wide", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    args = p.parse_args()
    long = pd.read_csv(args.long, sep="\t", low_memory=False)
    wide = pd.read_csv(args.wide, sep="\t", low_memory=False)
    wide = wide.rename(columns={wide.columns[0]: "ASV_ID"})
    taxonomy = long.drop_duplicates("ASV_ID").set_index("ASV_ID")

    reasons = {}
    for asv, row in taxonomy.iterrows():
        matched = []
        for column, values in HOST_RULES.items():
            value = str(row.get(column, "")).strip()
            if value in values:
                matched.append(f"{column}:{value}")
        lineage = str(row.get("Taxon", ""))
        if "Vertebrata" in lineage and not any(x.startswith("Phylum:") for x in matched):
            matched.append("Taxon:Vertebrata")
        if matched:
            reasons[str(asv)] = "; ".join(matched)

    ids = set(reasons) & set(wide["ASV_ID"].astype(str))
    removed_wide = wide[wide["ASV_ID"].astype(str).isin(ids)].copy()
    audit = taxonomy.loc[sorted(ids)].reset_index()
    audit["reads_removed"] = audit["ASV_ID"].map(
        removed_wide.set_index("ASV_ID").iloc[:, :].sum(axis=1).astype(int)
    )
    audit["host_filter_reason"] = audit["ASV_ID"].map(reasons)

    filtered_wide = wide[~wide["ASV_ID"].astype(str).isin(ids)].copy()
    filtered_long = long[~long["ASV_ID"].astype(str).isin(ids)].copy()
    args.outdir.mkdir(parents=True, exist_ok=True)
    audit.to_csv(args.outdir / "host_asvs_removed.tsv", sep="\t", index=False)
    filtered_wide.to_csv(args.outdir / "ASV_host_filtered.tsv", sep="\t", index=False)
    filtered_long.to_csv(args.outdir / "ASV_master_long_host_filtered.tsv", sep="\t", index=False)
    summary = pd.DataFrame([{
        "input_asvs": len(wide), "removed_host_asvs": len(ids),
        "output_asvs": len(filtered_wide),
        "input_reads": int(wide.iloc[:, 1:].to_numpy().sum()),
        "removed_host_reads": int(removed_wide.iloc[:, 1:].to_numpy().sum()),
        "output_reads": int(filtered_wide.iloc[:, 1:].to_numpy().sum()),
        "samples": len(wide.columns) - 1,
    }])
    summary.to_csv(args.outdir / "host_filter_summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
