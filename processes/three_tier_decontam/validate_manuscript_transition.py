#!/usr/bin/env python3
"""Fail unless recomputed three-tier outputs reproduce deposited SD1 and SD3."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--observed-sd1", required=True, type=Path)
    p.add_argument("--observed-wide", required=True, type=Path)
    p.add_argument("--expected-sd1", required=True, type=Path)
    p.add_argument("--expected-wide", required=True, type=Path)
    p.add_argument("--report", required=True, type=Path)
    args = p.parse_args()

    obs1 = pd.read_csv(args.observed_sd1, sep="\t", low_memory=False).sort_values("ASV_ID").reset_index(drop=True)
    exp1 = pd.read_csv(args.expected_sd1, sep="\t", low_memory=False).sort_values("ASV_ID").reset_index(drop=True)
    obsw = pd.read_csv(args.observed_wide, sep="\t", low_memory=False)
    expw = pd.read_csv(args.expected_wide, sep="\t", low_memory=False)

    flag_cols = ["ASV_ID", "flagged_prevalence_pooled", "flagged_frequency_within_type",
                 "flagged_biological_plausibility", "removal_reason"]
    score_cols = ["decontam_prev_p", "decontam_pooled_freq_p",
                  "decontam_pooled_combined_p", "decontam_within_type_combined_p"]
    ids_equal = obs1["ASV_ID"].tolist() == exp1["ASV_ID"].tolist()
    flags_equal = ids_equal and obs1[flag_cols].equals(exp1[flag_cols])
    scores_equal = ids_equal and all(np.allclose(obs1[c], exp1[c], rtol=0, atol=1e-12, equal_nan=True) for c in score_cols)
    wide_equal = obsw.equals(expw)
    report = {
        "sd1_ids_equal": ids_equal,
        "sd1_flags_equal": flags_equal,
        "sd1_scores_equal_at_1e-12": scores_equal,
        "sd3_table_equal": wide_equal,
        "observed_removed_asvs": len(obs1),
        "observed_wide_asvs": len(obsw),
        "observed_samples": len(obsw.columns) - 1,
        "observed_reads": int(obsw.iloc[:, 1:].to_numpy().sum()),
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    if not all((ids_equal, flags_equal, scores_equal, wide_equal)):
        raise SystemExit(f"SPARK manuscript transition validation failed; see {args.report}")


if __name__ == "__main__":
    main()
