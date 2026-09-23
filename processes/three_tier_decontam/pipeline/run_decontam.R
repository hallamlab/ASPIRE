#!/usr/bin/env Rscript
# run_decontam.R
# Apply decontam (Davis et al. 2018) to a raw ASV count table using both
# the prevalence method (vs negative controls) and the frequency method
# (vs per-sample DNA concentration), and add an abundance-ratio sanity
# check for each flagged ASV.
#
# Outputs (written to --outdir):
#   - decontam_per_asv_scores.tsv  : per-ASV p-values, flags, and ratios
#   - decontam_flagged_summary.tsv : count summary by threshold and method
#   - decontam_flagged_at_0.1.txt  : flagged ASV IDs at the default threshold
#   - decontam_flagged_at_0.5.txt  : flagged ASV IDs at the aggressive threshold
#   - decontam_run_log.txt         : run-time summary

suppressPackageStartupMessages({
  library(optparse)
  library(decontam)
})

option_list <- list(
  make_option("--counts", type="character",
              help="Raw ASV count table TSV (rows=ASVs, columns=samples, including negative controls)"),
  make_option("--metadata", type="character",
              help="Sample metadata TSV. Must contain a sample-id column and an is_negative_control logical column."),
  make_option("--sample-col", type="character", default="Sample",
              help="Metadata column matching count-table column names [default: %default]"),
  make_option("--neg-col", type="character", default="is_negative_control",
              help="Metadata logical column flagging negative controls [default: %default]"),
  make_option("--conc-col", type="character", default="DNA_conc",
              help="Metadata column with per-sample DNA concentration (numeric, NA where missing) [default: %default]"),
  make_option("--exclude-positives-col", type="character", default="is_positive_control",
              help="Metadata logical column flagging positive controls to exclude from analysis [default: %default]"),
  make_option("--threshold-default", type="numeric", default=0.1,
              help="Decontam p-value threshold for primary contaminant flag [default: %default]"),
  make_option("--threshold-aggressive", type="numeric", default=0.5,
              help="Decontam p-value threshold for sensitivity (aggressive) filter [default: %default]"),
  make_option("--outdir", type="character", default="decontam_output",
              help="Output directory [default: %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$counts) || is.null(opt$metadata)) {
  stop("Both --counts and --metadata are required.")
}

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(opt$outdir, "decontam_run_log.txt")
log_lines <- c()
log <- function(msg) {
  cat(msg, "\n")
  log_lines[length(log_lines)+1] <<- msg
}

log(paste("=== run_decontam.R ==="))
log(paste("counts:", opt$counts))
log(paste("metadata:", opt$metadata))
log(paste("sample-col:", opt$`sample-col`))
log(paste("neg-col:", opt$`neg-col`))
log(paste("conc-col:", opt$`conc-col`))
log(paste("primary threshold:", opt$`threshold-default`))
log(paste("aggressive threshold:", opt$`threshold-aggressive`))
log("")

# --- Load inputs ----------------------------------------------------------
log("Loading count table...")
counts <- read.table(opt$counts, sep="\t", header=TRUE, row.names=1,
                     check.names=FALSE, comment.char="", quote="")
log(paste("  ASVs:", nrow(counts), " Samples:", ncol(counts)))

log("Loading metadata...")
meta <- read.table(opt$metadata, sep="\t", header=TRUE, check.names=FALSE,
                   comment.char="", quote="")
log(paste("  Rows:", nrow(meta)))

if (!(opt$`sample-col` %in% colnames(meta))) {
  stop(paste("Sample column", opt$`sample-col`, "not in metadata."))
}
if (!(opt$`neg-col` %in% colnames(meta))) {
  stop(paste("Neg column", opt$`neg-col`, "not in metadata."))
}

# --- Align metadata to count-table column order --------------------------
sample_ids <- colnames(counts)
meta_match <- meta[match(sample_ids, meta[[opt$`sample-col`]]), , drop=FALSE]
if (any(is.na(meta_match[[opt$`sample-col`]]))) {
  unmatched <- sample_ids[is.na(meta_match[[opt$`sample-col`]])]
  log(paste("WARN: count-table columns missing from metadata (n=", length(unmatched), "):", sep=""))
  log(paste("  ", paste(head(unmatched, 10), collapse=", ")))
  log("  These will be dropped from decontam analysis.")
  keep <- !is.na(meta_match[[opt$`sample-col`]])
  counts <- counts[, keep, drop=FALSE]
  meta_match <- meta_match[keep, , drop=FALSE]
}

# Exclude positive controls
if (opt$`exclude-positives-col` %in% colnames(meta_match)) {
  pos <- as.logical(meta_match[[opt$`exclude-positives-col`]])
  pos[is.na(pos)] <- FALSE
  if (any(pos)) {
    log(paste("Excluding", sum(pos), "positive-control samples from decontam analysis."))
    counts <- counts[, !pos, drop=FALSE]
    meta_match <- meta_match[!pos, , drop=FALSE]
  }
}

# Parse neg-control logical
neg <- as.logical(meta_match[[opt$`neg-col`]])
neg[is.na(neg)] <- FALSE
log(paste("Negative-control samples:", sum(neg), "of", length(neg), "total."))
if (sum(neg) < 2) {
  stop("Need at least 2 negative controls for prevalence-based decontam.")
}

# Parse DNA conc (allows NA)
has_conc <- opt$`conc-col` %in% colnames(meta_match)
if (has_conc) {
  conc <- suppressWarnings(as.numeric(meta_match[[opt$`conc-col`]]))
  n_with_conc <- sum(!is.na(conc) & !neg)  # non-NA biological samples
  log(paste("Samples with non-NA DNA concentration (excluding negatives):", n_with_conc))
} else {
  conc <- NULL
  n_with_conc <- 0
  log(paste("DNA concentration column not found; frequency method will be skipped."))
}

# Counts must be transposed for decontam (samples x ASVs)
mat <- t(as.matrix(counts))
mode(mat) <- "integer"

# --- Run decontam prevalence (uses ALL negatives) ------------------------
log("")
log("Running decontam prevalence method...")
prev <- isContaminant(mat, method="prevalence", neg=neg,
                      threshold=opt$`threshold-default`)
log(paste("  ASVs flagged at p<", opt$`threshold-default`, ":",
          sum(prev$contaminant, na.rm=TRUE)))
log(paste("  ASVs flagged at p<", opt$`threshold-aggressive`, ":",
          sum(prev$p < opt$`threshold-aggressive`, na.rm=TRUE)))

# --- Run decontam frequency on the subset with DNA_conc ------------------
freq <- NULL
if (has_conc && n_with_conc >= 3) {
  log("")
  log("Running decontam frequency method on samples with DNA_conc...")
  # Subset to non-neg samples with non-NA conc
  keep_freq <- !neg & !is.na(conc) & conc > 0
  if (sum(keep_freq) < 3) {
    log("  Too few samples with non-NA DNA_conc; skipping frequency method.")
  } else {
    mat_freq <- mat[keep_freq, , drop=FALSE]
    conc_freq <- conc[keep_freq]
    log(paste("  Samples used for frequency test:", nrow(mat_freq)))
    freq <- isContaminant(mat_freq, method="frequency", conc=conc_freq,
                          threshold=opt$`threshold-default`)
    log(paste("  ASVs flagged at p<", opt$`threshold-default`, ":",
              sum(freq$contaminant, na.rm=TRUE)))
    log(paste("  ASVs flagged at p<", opt$`threshold-aggressive`, ":",
              sum(freq$p < opt$`threshold-aggressive`, na.rm=TRUE)))
  }
} else {
  log("Skipping frequency method (no DNA_conc data).")
}

# --- Build per-ASV combined output ----------------------------------------
asv_ids <- rownames(counts)
out <- data.frame(
  ASV_ID = asv_ids,
  prev_p = prev$p[match(asv_ids, rownames(prev))],
  prev_contam_at_default = prev$contaminant[match(asv_ids, rownames(prev))],
  stringsAsFactors = FALSE
)
if (!is.null(freq)) {
  out$freq_p <- freq$p[match(asv_ids, rownames(freq))]
  out$freq_contam_at_default <- freq$contaminant[match(asv_ids, rownames(freq))]
  # Combined: minimum of the two p-values (na.rm to handle ASVs absent in freq arm)
  out$combined_min_p <- pmin(out$prev_p, out$freq_p, na.rm=TRUE)
} else {
  out$freq_p <- NA_real_
  out$freq_contam_at_default <- NA
  out$combined_min_p <- out$prev_p
}

out$flag_at_0.1  <- out$combined_min_p < 0.1
out$flag_at_0.5  <- out$combined_min_p < 0.5
out$flag_at_0.1[is.na(out$flag_at_0.1)] <- FALSE
out$flag_at_0.5[is.na(out$flag_at_0.5)] <- FALSE

# --- Abundance-ratio sanity check ----------------------------------------
# For each ASV, compute median reads in biological samples vs negs
log("")
log("Computing abundance ratios (biological samples vs negative controls)...")
bio <- !neg
neg_means <- rowMeans(counts[, neg, drop=FALSE])
bio_means <- rowMeans(counts[, bio, drop=FALSE])
neg_max   <- apply(counts[, neg, drop=FALSE], 1, max)
bio_max   <- apply(counts[, bio, drop=FALSE], 1, max)

out$neg_mean_reads <- neg_means[match(asv_ids, names(neg_means))]
out$bio_mean_reads <- bio_means[match(asv_ids, names(bio_means))]
out$neg_max_reads  <- neg_max[match(asv_ids, names(neg_max))]
out$bio_max_reads  <- bio_max[match(asv_ids, names(bio_max))]
# Use mean for the ratio; add 1 to denominator to avoid div by zero
out$bio_to_neg_mean_ratio <- (out$bio_mean_reads + 0.1) / (out$neg_mean_reads + 0.1)

# Confidence category for flagged ASVs
classify <- function(flag, ratio) {
  ifelse(!flag, "not_flagged",
         ifelse(is.na(ratio), "flagged_no_neg_signal",
                ifelse(ratio < 5,  "high_confidence_contam",
                       ifelse(ratio < 50, "borderline_likely_cross_contam",
                              "likely_false_positive_flag"))))
}
out$confidence_at_0.1 <- classify(out$flag_at_0.1, out$bio_to_neg_mean_ratio)
out$confidence_at_0.5 <- classify(out$flag_at_0.5, out$bio_to_neg_mean_ratio)

# Sort by combined p, then by ratio (most contaminant-like first)
out <- out[order(out$combined_min_p, -out$bio_to_neg_mean_ratio), ]

# --- Write outputs --------------------------------------------------------
write.table(out, file.path(opt$outdir, "decontam_per_asv_scores.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)
writeLines(out$ASV_ID[out$flag_at_0.1],
           file.path(opt$outdir, "decontam_flagged_at_0.1.txt"))
writeLines(out$ASV_ID[out$flag_at_0.5],
           file.path(opt$outdir, "decontam_flagged_at_0.5.txt"))

# Summary
summ <- data.frame(
  threshold = c("0.1 (default)", "0.5 (aggressive)"),
  n_flagged = c(sum(out$flag_at_0.1), sum(out$flag_at_0.5)),
  high_confidence_contam = c(
    sum(out$confidence_at_0.1 == "high_confidence_contam"),
    sum(out$confidence_at_0.5 == "high_confidence_contam")
  ),
  borderline = c(
    sum(out$confidence_at_0.1 == "borderline_likely_cross_contam"),
    sum(out$confidence_at_0.5 == "borderline_likely_cross_contam")
  ),
  likely_false_positive = c(
    sum(out$confidence_at_0.1 == "likely_false_positive_flag"),
    sum(out$confidence_at_0.5 == "likely_false_positive_flag")
  ),
  flagged_no_neg_signal = c(
    sum(out$confidence_at_0.1 == "flagged_no_neg_signal"),
    sum(out$confidence_at_0.5 == "flagged_no_neg_signal")
  )
)
write.table(summ, file.path(opt$outdir, "decontam_flagged_summary.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

log("")
log("=== Summary ===")
log(paste("Total ASVs:", nrow(out)))
log(paste("Flagged at p<0.1:", sum(out$flag_at_0.1)))
log(paste("  high-confidence contaminants (ratio<5x):",
          sum(out$confidence_at_0.1 == "high_confidence_contam")))
log(paste("  borderline (ratio 5-50x):",
          sum(out$confidence_at_0.1 == "borderline_likely_cross_contam")))
log(paste("  likely false-positive flags (ratio>50x):",
          sum(out$confidence_at_0.1 == "likely_false_positive_flag")))
log(paste("Flagged at p<0.5:", sum(out$flag_at_0.5)))
log("")
log("Outputs:")
log(paste("  ", file.path(opt$outdir, "decontam_per_asv_scores.tsv")))
log(paste("  ", file.path(opt$outdir, "decontam_flagged_summary.tsv")))
log(paste("  ", file.path(opt$outdir, "decontam_flagged_at_0.1.txt")))
log(paste("  ", file.path(opt$outdir, "decontam_flagged_at_0.5.txt")))

writeLines(log_lines, log_file)
