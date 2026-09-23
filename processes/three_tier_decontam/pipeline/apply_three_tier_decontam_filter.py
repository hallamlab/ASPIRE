#!/usr/bin/env python3
"""
apply_three_tier_decontam_filter.py

Apply the three-tier contamination filter to the analyzed ASV tables and
produce a decontam-filtered version of ASV_master_long_filtered.tsv and
ASV_master_count_wide_filtered.tsv ready for re-running run_main_analysis_pipeline.sh.

Three filters (union; an ASV is removed if flagged by ANY of them):

  1. decontam PREVALENCE method (pooled, vs 4 negative extraction controls)
     - applied at p < 0.1
     - input: decontam_per_asv_scores.tsv from the pooled decontam run

  2. decontam FREQUENCY method, within sample type (BAL, Bronchial Brush,
     Oral Rinse, each with all 4 negative controls)
     - applied at p < 0.1 (combined min across the three within-type runs)
     - input: decontam_by_sample_type_per_asv.tsv from the within-type run

  3. BIOLOGICAL-PLAUSIBILITY screen — removes ASVs assigned to thermophilic
     anaerobic / methanogenic environmental taxa whose growth optima
     (typically 55–70 C, anaerobic-digester / hot-spring environments) are
     incompatible with human airway physiology (~37 C). This class of
     contamination is detectable by neither decontam method (absent from
     extraction blanks; contamination route is likely run-level index hopping
     or barcode bleed from co-multiplexed environmental libraries) and must
     be removed on biological grounds rather than statistical evidence.

     Removed taxa (by Phylum or by Genus name pattern, whichever applies):
       Phyla:  Coprothermobacterota, Thermotogota, Caldatribacteriota,
               Aquificota
       Archaea: Methanothermobacteriaceae and other methanogenic Euryarchaeota;
                Thermoplasmatota Marine_Group_II
       Genera with thermophilic / pyrophilic name patterns:
         Thermo*, Pyro*, Coprothermo*, Defluviitoga, Fervidobacterium,
         Thermovirga, Caldicellulosiruptor, Methanothermo*, Thermogutta

Outputs (written to --outdir):
  - flagged_asvs_combined.tsv        — every flagged ASV with each filter's call,
                                        taxonomy, and final include/exclude decision
  - removed_asvs_with_taxonomy.tsv   — final list to remove, with taxonomy + rationale
  - kept_asvs.txt                    — final list to retain
  - ASV_master_long_filtered.tsv     — long table with flagged ASVs removed
  - ASV_master_count_wide_filtered.tsv  — wide table with flagged ASVs removed
  - filter_summary.txt               — counts at each filter step
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd


# --------------------------------------------------------------------------- #
# Biological-plausibility screen definitions
# --------------------------------------------------------------------------- #

THERMOPHILE_PHYLA = {
    "Coprothermobacterota",  # thermophilic anaerobes
    "Thermotogota",          # thermophilic anaerobes
    "Caldatribacteriota",    # thermophilic
    "Aquificota",            # hyperthermophilic
}

# Archaeal lineages that are methanogenic or marine-environmental and
# implausible as human airway residents.
ARCHAEAL_PHYLA_TO_REMOVE = {
    "Euryarchaeota",     # most methanogens — case-by-case OK below
    "Crenarchaeota",     # hyperthermophiles + some mesophilic — see GENUS_PATTERNS
    "Thermoplasmatota",  # marine + soil environmental
}

# Genus-name regex patterns (case-insensitive). If a genus matches, the ASV
# is flagged as a biological-plausibility removal regardless of phylum.
THERMOPHILE_GENUS_PATTERNS = [
    r"^Thermo",          # Thermovirga, Thermomonas, Thermogutta, etc.
    r"^Pyro",            # Pyrobaculum, Pyrococcus
    r"^Coprothermo",     # Coprothermobacter
    r"^Defluviitoga",
    r"^Fervidobacterium",
    r"^Caldicellulosiruptor",
    r"^Methanothermo",   # Methanothermobacter, Methanothermus
    r"^Caldisericum",
    r"^Bathyarchae",     # uncultured sediment archaea
]

# ASVs whose genus matches the thermophile pattern but are documented as
# mesophilic (would override the screen). None known for this dataset.
THERMOPHILE_GENUS_OVERRIDES = set()  # e.g., {"Thermomonas hydrothermalis subgenus"}


def is_thermophile(phylum, genus):
    """Return True if this ASV's taxonomy matches the biological-plausibility removal criteria."""
    if pd.isna(phylum):
        phylum = ""
    if pd.isna(genus):
        genus = ""

    # Phylum-level removal
    if phylum in THERMOPHILE_PHYLA:
        return True

    # Archaeal phyla — remove unless explicitly documented as mesophilic
    if phylum in ARCHAEAL_PHYLA_TO_REMOVE:
        # Don't remove the Bathyarchaeia genus from Crenarchaeota by default
        # since some are mesophilic — but Crenarchaeota in airway is implausible,
        # so we include it for removal.
        return True

    # Genus pattern matching
    for pattern in THERMOPHILE_GENUS_PATTERNS:
        if re.match(pattern, str(genus), re.IGNORECASE):
            if genus not in THERMOPHILE_GENUS_OVERRIDES:
                return True

    return False


# --------------------------------------------------------------------------- #
# I/O
# --------------------------------------------------------------------------- #

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pooled-scores", required=True, type=Path,
                    help="decontam_per_asv_scores.tsv from pooled run")
    ap.add_argument("--within-type-scores", required=True, type=Path,
                    help="decontam_by_sample_type_per_asv.tsv from within-type run")
    ap.add_argument("--analyzed-long", required=True, type=Path,
                    help="Current ASV_master_long_filtered.tsv")
    ap.add_argument("--analyzed-wide", required=True, type=Path,
                    help="Current ASV_master_count_wide_filtered.tsv")
    ap.add_argument("--outdir", required=True, type=Path,
                    help="Output directory for filtered tables")
    ap.add_argument("--prev-threshold", type=float, default=0.1,
                    help="Decontam prevalence p-value threshold [default: 0.1]")
    ap.add_argument("--freq-threshold", type=float, default=0.1,
                    help="Decontam frequency (within-type) p-value threshold [default: 0.1]")
    ap.add_argument("--disable-biological-plausibility", action="store_true",
                    help="Disable the third-tier thermophile/environmental-taxon screen")
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    log_lines = []
    def log(msg):
        print(msg, file=sys.stderr)
        log_lines.append(msg)

    log("=" * 70)
    log("Three-tier contamination filter")
    log("=" * 70)
    log(f"Pooled decontam scores : {args.pooled_scores}")
    log(f"Within-type decontam   : {args.within_type_scores}")
    log(f"Analyzed long table    : {args.analyzed_long}")
    log(f"Analyzed wide table    : {args.analyzed_wide}")
    log(f"Output directory       : {args.outdir}")
    log(f"Prevalence threshold   : p < {args.prev_threshold}")
    log(f"Frequency threshold    : p < {args.freq_threshold}")
    log("")

    # ----- Load decontam outputs ------------------------------------------
    pooled = pd.read_csv(args.pooled_scores, sep="\t", low_memory=False)
    within = pd.read_csv(args.within_type_scores, sep="\t", low_memory=False)

    # ----- Filter 1: decontam prevalence (pooled) -------------------------
    # Use the `prev_p` column from the pooled run, not the combined min p
    if "prev_p" not in pooled.columns:
        sys.exit("ERROR: pooled scores file missing 'prev_p' column")
    prev_flagged = pooled.loc[pooled["prev_p"] < args.prev_threshold, "ASV_ID"]
    prev_flagged = set(prev_flagged.dropna().astype(str))
    log(f"Filter 1 (prevalence pooled, p < {args.prev_threshold}): "
        f"{len(prev_flagged)} ASVs flagged")

    # ----- Filter 2: decontam frequency (within-type) ---------------------
    # Use the `within_type_combined_p` column from the within-type run
    if "within_type_combined_p" not in within.columns:
        sys.exit("ERROR: within-type scores file missing 'within_type_combined_p' column")
    freq_flagged = within.loc[within["within_type_combined_p"] < args.freq_threshold, "ASV_ID"]
    freq_flagged = set(freq_flagged.dropna().astype(str))
    log(f"Filter 2 (frequency within-type, p < {args.freq_threshold}): "
        f"{len(freq_flagged)} ASVs flagged")

    # ----- Load analyzed long for taxonomy lookup --------------------------
    long_df = pd.read_csv(args.analyzed_long, sep="\t", low_memory=False)
    log(f"Analyzed long table: {len(long_df):,} rows, "
        f"{long_df['ASV_ID'].nunique()} unique ASVs")
    tax = (long_df[["ASV_ID", "Domain", "Phylum", "Class", "Order",
                    "Family", "Genus", "Species", "Taxon"]]
           .drop_duplicates(subset=["ASV_ID"]))

    # ----- Filter 3: biological-plausibility (thermophiles) ---------------
    thermo_mask = (
        pd.Series(False, index=tax.index)
        if args.disable_biological_plausibility
        else tax.apply(lambda r: is_thermophile(r["Phylum"], r["Genus"]), axis=1)
    )
    thermo_flagged = set(tax.loc[thermo_mask, "ASV_ID"].astype(str))
    log(f"Filter 3 (biological plausibility / thermophiles): "
        f"{len(thermo_flagged)} ASVs flagged")
    log("  Thermophile taxa flagged:")
    for _, row in tax.loc[thermo_mask].sort_values(["Phylum", "Genus"]).iterrows():
        log(f"    {row['ASV_ID']:10s}  {row['Phylum']:25s}  "
            f"{row['Family']:30s}  {row['Genus']}")

    # ----- Combine: union of all three filters ----------------------------
    all_flagged = prev_flagged | freq_flagged | thermo_flagged
    log("")
    log(f"Total flagged (union of 3 filters): {len(all_flagged)} ASVs")
    log(f"  Flagged by prevalence only       : "
        f"{len(prev_flagged - freq_flagged - thermo_flagged)}")
    log(f"  Flagged by within-freq only      : "
        f"{len(freq_flagged - prev_flagged - thermo_flagged)}")
    log(f"  Flagged by thermophile only      : "
        f"{len(thermo_flagged - prev_flagged - freq_flagged)}")
    log(f"  Flagged by prev + within-freq    : "
        f"{len(prev_flagged & freq_flagged - thermo_flagged)}")
    log(f"  Flagged by prev + thermophile    : "
        f"{len(prev_flagged & thermo_flagged - freq_flagged)}")
    log(f"  Flagged by within-freq + thermo  : "
        f"{len(freq_flagged & thermo_flagged - prev_flagged)}")
    log(f"  Flagged by all three             : "
        f"{len(prev_flagged & freq_flagged & thermo_flagged)}")

    # ----- Restrict to ASVs actually in the analyzed table ----------------
    analyzed_set = set(tax["ASV_ID"].astype(str))
    to_remove = all_flagged & analyzed_set
    log("")
    log(f"Of {len(all_flagged)} flagged ASVs, {len(to_remove)} are in the analyzed table.")
    log(f"  (The other {len(all_flagged - analyzed_set)} were already filtered upstream "
        f"by abundance/prevalence/mito filters.)")

    # ----- Build combined report ------------------------------------------
    flagged_df = pd.DataFrame({"ASV_ID": sorted(all_flagged)})
    flagged_df["flagged_prevalence_pooled"] = flagged_df["ASV_ID"].isin(prev_flagged)
    flagged_df["flagged_frequency_within_type"] = flagged_df["ASV_ID"].isin(freq_flagged)
    flagged_df["flagged_biological_plausibility"] = flagged_df["ASV_ID"].isin(thermo_flagged)
    flagged_df["in_analyzed_table"] = flagged_df["ASV_ID"].isin(analyzed_set)
    flagged_df = flagged_df.merge(tax, on="ASV_ID", how="left")

    # Add p-values for transparency
    pooled_pvals = pooled.set_index("ASV_ID")[["prev_p", "freq_p", "combined_min_p"]]
    pooled_pvals = pooled_pvals.rename(columns={
        "prev_p": "decontam_prev_p",
        "freq_p": "decontam_pooled_freq_p",
        "combined_min_p": "decontam_pooled_combined_p",
    })
    within_pvals = within.set_index("ASV_ID")[["within_type_combined_p"]]
    within_pvals.columns = ["decontam_within_type_combined_p"]
    flagged_df = flagged_df.merge(pooled_pvals, left_on="ASV_ID", right_index=True, how="left")
    flagged_df = flagged_df.merge(within_pvals, left_on="ASV_ID", right_index=True, how="left")

    # Reason column
    def reason(row):
        reasons = []
        if row["flagged_prevalence_pooled"]:
            reasons.append("prevalence (vs negs)")
        if row["flagged_frequency_within_type"]:
            reasons.append("frequency (within-type)")
        if row["flagged_biological_plausibility"]:
            reasons.append("biological plausibility")
        return "; ".join(reasons)
    flagged_df["removal_reason"] = flagged_df.apply(reason, axis=1)

    flagged_df.to_csv(args.outdir / "flagged_asvs_combined.tsv",
                      sep="\t", index=False)
    log(f"Wrote flagged_asvs_combined.tsv ({len(flagged_df)} rows)")

    removed_in_analyzed = flagged_df[flagged_df["in_analyzed_table"]]
    removed_in_analyzed.to_csv(args.outdir / "removed_asvs_with_taxonomy.tsv",
                               sep="\t", index=False)
    log(f"Wrote removed_asvs_with_taxonomy.tsv ({len(removed_in_analyzed)} rows)")

    kept = sorted(analyzed_set - to_remove)
    Path(args.outdir / "kept_asvs.txt").write_text("\n".join(kept) + "\n")
    log(f"Wrote kept_asvs.txt ({len(kept)} ASVs kept of {len(analyzed_set)} analyzed)")

    # ----- Filter the long table ------------------------------------------
    log("")
    log("Filtering analyzed long table...")
    long_filtered = long_df[~long_df["ASV_ID"].astype(str).isin(to_remove)]
    out_long = args.outdir / "ASV_master_long_filtered.tsv"
    long_filtered.to_csv(out_long, sep="\t", index=False)
    log(f"  Original rows : {len(long_df):,}")
    log(f"  Filtered rows : {len(long_filtered):,}")
    log(f"  Removed       : {len(long_df) - len(long_filtered):,} "
        f"({(len(long_df) - len(long_filtered)) * 100 / len(long_df):.2f}%)")
    log(f"  Wrote {out_long}")

    # ----- Read-loss accounting -------------------------------------------
    if "count" in long_df.columns:
        total_reads = long_df["count"].sum()
        kept_reads = long_filtered["count"].sum()
        log(f"  Reads in original table : {total_reads:,}")
        log(f"  Reads in filtered table : {kept_reads:,}")
        log(f"  Reads removed           : {total_reads - kept_reads:,} "
            f"({(total_reads - kept_reads) * 100 / total_reads:.3f}%)")
        # By sample type
        if "type_group" in long_df.columns:
            log("")
            log("  Per-sample-type read loss:")
            for tg in long_df["type_group"].dropna().unique():
                orig_st = long_df.loc[long_df["type_group"] == tg, "count"].sum()
                kept_st = long_filtered.loc[long_filtered["type_group"] == tg, "count"].sum()
                pct = (orig_st - kept_st) * 100 / orig_st if orig_st else 0
                log(f"    {tg:20s}  removed {orig_st - kept_st:>10,} of {orig_st:>10,} "
                    f"reads ({pct:.2f}%)")

    # ----- Filter the wide table ------------------------------------------
    log("")
    log("Filtering analyzed wide table...")
    wide_df = pd.read_csv(args.analyzed_wide, sep="\t", low_memory=False)
    first_col = wide_df.columns[0]
    wide_filtered = wide_df[~wide_df[first_col].astype(str).isin(to_remove)]
    out_wide = args.outdir / "ASV_master_count_wide_filtered.tsv"
    wide_filtered.to_csv(out_wide, sep="\t", index=False)
    log(f"  Original ASVs : {len(wide_df):,}")
    log(f"  Filtered ASVs : {len(wide_filtered):,}")
    log(f"  Wrote {out_wide}")

    # ----- Final summary --------------------------------------------------
    summary_path = args.outdir / "filter_summary.txt"
    summary_path.write_text("\n".join(log_lines) + "\n")
    log("")
    log(f"Summary written to {summary_path}")
    log("Done.")


if __name__ == "__main__":
    main()
