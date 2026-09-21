#!/usr/bin/env python3
"""Generate deterministic clear, overlapping, and continuous microbiome fixtures."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def _multinomial_table(probabilities: np.ndarray, depth: int, rng: np.random.Generator) -> np.ndarray:
    return np.vstack([rng.multinomial(depth, row / row.sum()) for row in probabilities])


def _write_scenario(root: Path, counts: np.ndarray, metadata: pd.DataFrame) -> None:
    root.mkdir(parents=True, exist_ok=True)
    asv_ids = [f"ASV{index + 1}" for index in range(counts.shape[1])]
    sample_ids = metadata["sampleID"].tolist()
    count_table = pd.DataFrame(counts.T, columns=sample_ids)
    count_table.insert(0, "ASV_ID", asv_ids)
    count_table.to_csv(root / "asv_counts.tsv", sep="\t", index=False)
    pd.DataFrame({
        "ASV_ID": asv_ids,
        "retained_final": True,
        "retained_prevalence": True,
        "retained_abundance": True,
    }).to_csv(root / "spieceasi_filtering_audit.csv", index=False)
    metadata.to_csv(root / "metadata.tsv", sep="\t", index=False)
    pd.DataFrame({
        "Feature ID": asv_ids,
        "Taxon": [f"d__Bacteria; p__Synthetic_{index % 4}; g__Taxon_{index + 1}" for index in range(len(asv_ids))],
    }).to_csv(root / "taxonomy.tsv", sep="\t", index=False)


def generate(output_root: Path, seed: int = 42) -> None:
    rng = np.random.default_rng(seed)
    n_samples, n_asvs = 90, 30
    sample_ids = [f"S{index + 1:03d}" for index in range(n_samples)]
    n_profiles = n_samples // 3
    common_metadata = pd.DataFrame({
        "sampleID": sample_ids,
        "Cruise": np.repeat([f"P{index + 1:03d}" for index in range(n_profiles)], 3),
        "Depth": np.tile([10, 100, 200], n_samples // 3),
        "date": np.repeat(
            pd.date_range("2010-01-01", periods=n_profiles, freq="30D").astype(str), 3
        ),
        "Season": np.repeat(np.tile(["Winter", "Spring", "Summer", "Autumn", "Winter"], 6), 3),
    })

    # Three deliberately separated assemblages with distinct ten-ASV blocks.
    clear_labels = np.tile(["E1", "E2", "E3"], n_profiles)
    clear_prob = np.full((n_samples, n_asvs), 0.01)
    for row, label in enumerate(clear_labels):
        block = int(label[-1]) - 1
        clear_prob[row, block * 10:(block + 1) * 10] = rng.lognormal(1.8, 0.18, 10)
    clear_counts = _multinomial_table(clear_prob, 20000, rng)
    clear_metadata = common_metadata.assign(
        o2_compartment=clear_labels,
        gmm_component=np.where(clear_labels == "E1", "G0", np.where(clear_labels == "E2", "G1", "G2")),
        o2_subcompartment_final=[f"{label}_hybrid" for label in clear_labels],
    )
    _write_scenario(output_root / "clear", clear_counts, clear_metadata)

    # Partially overlapping assemblages retain a weak class signal for sensitivity checks.
    overlap_labels = np.tile(["E1", "E2", "E3"], n_profiles)
    overlap_prob = rng.lognormal(0.0, 0.55, size=(n_samples, n_asvs))
    for row, label in enumerate(overlap_labels):
        block = int(label[-1]) - 1
        overlap_prob[row, block * 10:(block + 1) * 10] *= 1.45
    overlap_counts = _multinomial_table(overlap_prob, 20000, rng)
    overlap_metadata = common_metadata.assign(
        o2_compartment=overlap_labels,
        gmm_component="G0",
        o2_subcompartment_final=[f"{label}_hybrid" for label in overlap_labels],
    )
    _write_scenario(output_root / "overlapping", overlap_counts, overlap_metadata)

    # A single unimodal logistic-normal cloud has continuous compositional variation but no classes.
    latent = rng.normal(0, 1, size=(n_samples, 5))
    loadings = rng.normal(0, 0.32, size=(5, n_asvs))
    clr_like = latent @ loadings + rng.normal(0, 0.45, size=(n_samples, n_asvs))
    continuous_prob = np.exp(clr_like - clr_like.max(axis=1, keepdims=True))
    continuous_counts = _multinomial_table(continuous_prob, 20000, rng)
    continuous_metadata = common_metadata.assign(
        o2_compartment=np.where(common_metadata.Depth.eq(10), "oxic", "low_oxygen"),
        gmm_component="G0",
        o2_subcompartment_final="continuous_reference_only",
    )
    _write_scenario(output_root / "continuous", continuous_counts, continuous_metadata)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    generate(args.outdir, args.seed)


if __name__ == "__main__":
    main()
