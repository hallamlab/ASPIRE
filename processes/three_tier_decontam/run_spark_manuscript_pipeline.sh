#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: run_spark_manuscript_pipeline.sh \
  --final-dir DIR --metadata TSV --submission-package DIR --output-dir DIR \
  [--rscript PATH] [--downstream-rscript PATH] [--downstream-python PATH]

Recomputes the deposited SPARK three-tier decontamination from an authoritative
ASPIRE final output, validates it against SD1/SD3, then rebuilds the deposited
main downstream analysis in a separate output directory.
EOF
}

FINAL_DIR=""
PACKAGE_DIR=""
OUTPUT_DIR=""
METADATA=""
RSCRIPT="Rscript"
DOWNSTREAM_RSCRIPT="Rscript"
DOWNSTREAM_PYTHON="python3"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --final-dir) FINAL_DIR="$2"; shift 2 ;;
    --submission-package) PACKAGE_DIR="$2"; shift 2 ;;
    --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
    --metadata) METADATA="$2"; shift 2 ;;
    --rscript) RSCRIPT="$2"; shift 2 ;;
    --downstream-rscript) DOWNSTREAM_RSCRIPT="$2"; shift 2 ;;
    --downstream-python) DOWNSTREAM_PYTHON="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done
[[ -n "$FINAL_DIR" && -n "$METADATA" && -n "$PACKAGE_DIR" && -n "$OUTPUT_DIR" ]] || { usage >&2; exit 2; }
[[ -d "$FINAL_DIR" ]] || { echo "Missing final directory: $FINAL_DIR" >&2; exit 2; }
[[ -f "$METADATA" ]] || { echo "Missing metadata table: $METADATA" >&2; exit 2; }
[[ -f "$PACKAGE_DIR/SPARK_analysis_code.zip" ]] || { echo "Missing analysis archive in: $PACKAGE_DIR" >&2; exit 2; }
[[ ! -e "$OUTPUT_DIR" ]] || { echo "Refusing to overwrite existing output: $OUTPUT_DIR" >&2; exit 2; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ZENODO="$PACKAGE_DIR/ZENODO_UPLOAD"
mkdir -p "$OUTPUT_DIR" "$OUTPUT_DIR/_runtime/code" "$OUTPUT_DIR/pre_decontam" \
  "$OUTPUT_DIR/decontam/pooled" "$OUTPUT_DIR/decontam/within_type" "$OUTPUT_DIR/decontam/filtered"
unzip -q "$PACKAGE_DIR/SPARK_analysis_code.zip" -d "$OUTPUT_DIR/_runtime/code"
SPARK_CODE="$OUTPUT_DIR/_runtime/code/SPARK"
python3 "$SCRIPT_DIR/patch_deposited_code.py" --code-dir "$SPARK_CODE"

python3 "$SCRIPT_DIR/prepare_authoritative_tables.py" \
  --final-dir "$FINAL_DIR" \
  --reference-wide "$ZENODO/DatasetSD3_ASV_count_table.tsv" \
  --removed-asvs "$ZENODO/DatasetSD1_removed_ASVs.tsv" \
  --outdir "$OUTPUT_DIR/pre_decontam"

python3 "$SPARK_CODE/decontam/build_decontam_metadata.py" \
  --counts "$FINAL_DIR/ASVs/ASV_counts.tsv" \
  --metadata-in "$METADATA" \
  --metadata-out "$OUTPUT_DIR/decontam/decontam_metadata.tsv"

"$RSCRIPT" "$SPARK_CODE/decontam/run_decontam.R" \
  --counts "$FINAL_DIR/ASVs/ASV_counts.tsv" \
  --metadata "$OUTPUT_DIR/decontam/decontam_metadata.tsv" \
  --outdir "$OUTPUT_DIR/decontam/pooled"

"$RSCRIPT" "$SPARK_CODE/decontam/run_decontam_by_sample_type.R" \
  --counts "$FINAL_DIR/ASVs/ASV_counts.tsv" \
  --metadata "$OUTPUT_DIR/decontam/decontam_metadata.tsv" \
  --pooled "$OUTPUT_DIR/decontam/pooled/decontam_per_asv_scores.tsv" \
  --outdir "$OUTPUT_DIR/decontam/within_type"

python3 "$SPARK_CODE/decontam/apply_three_tier_decontam_filter.py" \
  --pooled-scores "$OUTPUT_DIR/decontam/pooled/decontam_per_asv_scores.tsv" \
  --within-type-scores "$OUTPUT_DIR/decontam/within_type/decontam_by_sample_type_per_asv.tsv" \
  --analyzed-long "$OUTPUT_DIR/pre_decontam/ASV_master_long.tsv" \
  --analyzed-wide "$OUTPUT_DIR/pre_decontam/ASV_master_count_wide.tsv" \
  --outdir "$OUTPUT_DIR/decontam/filtered" \
  --prev-threshold 0.1 --freq-threshold 0.1

python3 "$SCRIPT_DIR/validate_manuscript_transition.py" \
  --observed-sd1 "$OUTPUT_DIR/decontam/filtered/removed_asvs_with_taxonomy.tsv" \
  --observed-wide "$OUTPUT_DIR/decontam/filtered/ASV_master_count_wide_filtered.tsv" \
  --expected-sd1 "$ZENODO/DatasetSD1_removed_ASVs.tsv" \
  --expected-wide "$ZENODO/DatasetSD3_ASV_count_table.tsv" \
  --report "$OUTPUT_DIR/decontam/validation.json"

python3 "$SCRIPT_DIR/build_downstream_long.py" \
  --long "$OUTPUT_DIR/decontam/filtered/ASV_master_long_filtered.tsv" \
  --wide "$OUTPUT_DIR/decontam/filtered/ASV_master_count_wide_filtered.tsv" \
  --output "$OUTPUT_DIR/decontam/filtered/ASV_master_long_downstream.tsv"

# Only validated tables are allowed into downstream manuscript analyses.
# The 125-sample manuscript cohort was already selected upstream. Reapplying
# a read threshold here would incorrectly discard low-biomass brush samples.
mkdir -p "$OUTPUT_DIR/_runtime/tool_bin" "$OUTPUT_DIR/_runtime/matplotlib"
ln -s "$DOWNSTREAM_RSCRIPT" "$OUTPUT_DIR/_runtime/tool_bin/Rscript"
ln -s "$DOWNSTREAM_PYTHON" "$OUTPUT_DIR/_runtime/tool_bin/python3"
MPLCONFIGDIR="$OUTPUT_DIR/_runtime/matplotlib" PATH="$OUTPUT_DIR/_runtime/tool_bin:$PATH" \
  bash "$SPARK_CODE/run_main_analysis_pipeline.sh" \
  --data-long "$OUTPUT_DIR/decontam/filtered/ASV_master_long_downstream.tsv" \
  --data-wide "$OUTPUT_DIR/decontam/filtered/ASV_master_count_wide_filtered.tsv" \
  --outdir "$OUTPUT_DIR/downstream" \
  --sample-col sample --patient-col Participant_ID --case-col Case --type-col type_group \
  --min-sample-reads 0

printf 'SPARK manuscript pipeline completed and validated: %s\n' "$OUTPUT_DIR"
