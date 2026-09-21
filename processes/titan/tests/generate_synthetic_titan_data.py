#!/usr/bin/env python3
"""Generate a small TITAN validation dataset with known opposing transitions."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def generate(outdir: Path, n_samples: int = 24) -> None:
    if n_samples < 20:
        raise ValueError("The synthetic validation dataset requires at least 20 samples")
    outdir.mkdir(parents=True, exist_ok=True)
    samples = [f"S{i + 1:02d}" for i in range(n_samples)]
    gradient = np.linspace(0.0, 100.0, n_samples)
    midpoint = n_samples // 2
    counts = pd.DataFrame({"ASV_ID": ["ASV_declining", "ASV_increasing", "ASV_background", "ASV_rare", "ASV_filtered"]})
    for index, sample in enumerate(samples):
        counts[sample] = [
            500 if index < midpoint else 5,
            5 if index < midpoint else 500,
            100 + (index % 3),
            20 if index < 2 else 0,
            50,
        ]
    counts.to_csv(outdir / "asv_counts.tsv", sep="\t", index=False)
    pd.DataFrame({
        "ASV_ID": counts["ASV_ID"],
        "retained_final": [True, True, True, True, False],
        "retention_reason": ["standard_filter"] * 4 + ["filtered_out"],
    }).to_csv(outdir / "spieceasi_filtering_audit.csv", index=False)
    pd.DataFrame({
        "sampleID": samples,
        "synthetic_gradient": gradient,
        "invariant_gradient": np.ones(n_samples),
        "sparse_gradient": [float(i) if i < 9 else np.nan for i in range(n_samples)],
    }).to_csv(outdir / "measurement_matrix.samples_by_measurement.tsv", sep="\t", index=False)
    pd.DataFrame({
        "Feature ID": counts["ASV_ID"],
        "Taxon": [
            "d__Bacteria;p__Declinata", "d__Bacteria;p__Incremata",
            "d__Bacteria;p__Backgroundota", "d__Bacteria;p__Rareota",
            "d__Bacteria;p__Filteredota",
        ],
    }).to_csv(outdir / "taxonomy.tsv", sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--n-samples", type=int, default=24)
    args = parser.parse_args()
    generate(args.outdir, args.n_samples)


if __name__ == "__main__":
    main()
