#!/usr/bin/env python3
"""Apply minimal compatibility fixes to the deposited SPARK analysis copy."""

import argparse
from pathlib import Path


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--code-dir", required=True, type=Path)
    args = p.parse_args()
    target = args.code_dir / "plot_lung_status_analysis.py"
    text = target.read_text()
    replacements = {
        "pcoa_result.proportion_explained[0]": "pcoa_result.proportion_explained.iloc[0]",
        "pcoa_result.proportion_explained[1]": "pcoa_result.proportion_explained.iloc[1]",
    }
    for old, new in replacements.items():
        if old not in text and new not in text:
            raise SystemExit(f"Expected deposited-code expression not found: {old}")
        text = text.replace(old, new)
    target.write_text(text)


if __name__ == "__main__":
    main()
