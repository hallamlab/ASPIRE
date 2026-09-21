#!/usr/bin/env python3
"""Generate deterministic turnover fixtures with known ecological behavior."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


DEPTHS = [10, 50, 100, 150, 200]
ASVS = [f"ASV{i}" for i in range(1, 13)]


def community(state: str, variant: int = 0) -> np.ndarray:
    values = np.full(len(ASVS), 0.005)
    if state == "upper":
        values[0:4] = [0.34, 0.26, 0.18, 0.12]
    else:
        values[4:8] = [0.31, 0.27, 0.20, 0.12]
    if variant:
        values[8:12] = [0.30, 0.25, 0.20, 0.15]
        values[0:8] *= 0.20
    return values / values.sum()


def write_scenario(root: Path, rule) -> None:
    root.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(42)
    profiles = ["P1", "P2", "P3", "P4"]
    dates = pd.to_datetime(["2020-01-01", "2020-02-05", "2020-04-20", "2020-06-01"])
    samples: list[dict[str, object]] = []
    count_columns: dict[str, np.ndarray] = {}
    for profile_index, (profile, date) in enumerate(zip(profiles, dates)):
        profile_depths = DEPTHS if profile != "P3" else [10, 50, 150, 200]
        for depth in profile_depths:
            state, variant = rule(profile_index, depth)
            sample = f"{profile}_{depth}m"
            probabilities = community(state, variant)
            count_columns[sample] = rng.multinomial(25000, probabilities)
            samples.append({
                "sampleID": sample, "Cruise": profile, "date": date.strftime("%Y-%m-%d"),
                "Depth": depth, "o2_compartment": "oxic" if state == "upper" else "suboxic",
                "gmm_component": "G0" if state == "upper" else "G1",
                "o2_subcompartment_final": f"{'oxic' if state == 'upper' else 'suboxic'}__{'gmm0' if state == 'upper' else 'gmm1'}",
            })
            if profile == "P1" and depth == 50:
                replicate = f"{profile}_{depth}m_biological_repeat"
                count_columns[replicate] = rng.multinomial(25000, probabilities)
                samples.append({**samples[-1], "sampleID": replicate})
    counts = pd.DataFrame({"ASV_ID": ASVS, **count_columns})
    counts.to_csv(root / "asv_counts.tsv", sep="\t", index=False)
    pd.DataFrame(samples).to_csv(root / "metadata.tsv", sep="\t", index=False)
    pd.DataFrame({"ASV_ID": ASVS, "retained_final": True}).to_csv(
        root / "spieceasi_filtering_audit.csv", index=False
    )


def generate(root: Path) -> None:
    # Stable, sharp transition between 50 and 100 m.
    write_scenario(root / "sharp_depth", lambda _time, depth: ("upper" if depth < 100 else "lower", 0))
    # The same upper/lower assemblages move downward through successive profiles.
    boundaries = [100, 150, 200, 150]
    write_scenario(
        root / "moving_boundary",
        lambda time, depth: ("upper" if depth < boundaries[time] else "lower", 0),
    )
    # Compartment identity remains depth-stable while its microbial composition changes after time 2.
    write_scenario(
        root / "within_compartment_change",
        lambda time, depth: ("upper" if depth < 100 else "lower", 1 if time >= 2 else 0),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    generate(args.outdir)


if __name__ == "__main__":
    main()
