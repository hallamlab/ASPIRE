#!/usr/bin/env Rscript
# run_decontam_by_sample_type.R
# Run decontam separately within each sample type (with all negative controls
# included in each run), then combine results across sample types and compare
# to the pooled-decontam output.
#
# Rationale: pooled frequency-method decontam can be confounded with sample
# type when DNA concentration is structurally lower in some sample types
# (e.g., Bronchial Brush). Running within each sample type removes this
# cross-type confound, isolating true contamination signal (which should
# appear within at least one sample type if real).
#
# Inputs:
#   --counts:    raw ASV count table TSV (rows=ASVs, columns=samples, includes
#                negative controls)
#   --metadata:  same metadata TSV used by run_decontam.R
#                (Sample, is_negative_control, is_positive_control, DNA_conc,
#                Type_Group)
#   --pooled:    optional path to decontam_per_asv_scores.tsv from a pooled run,
#                for side-by-side comparison
#
# Outputs (in --outdir):
#   - decontam_by_sample_type_per_asv.tsv : per-ASV per-sample-type flags
#   - decontam_by_sample_type_summary.tsv : count summary
#   - decontam_within_type_flagged_at_0.1.txt
#   - decontam_within_type_flagged_at_0.5.txt
#   - decontam_pooled_vs_within_type_comparison.tsv (if --pooled provided)
#   - decontam_by_sample_type_log.txt

suppressPackageStartupMessages({
  library(optparse)
  library(decontam)
})

option_list <- list(
  make_option("--counts", type="character", help="Raw ASV count table TSV"),
  make_option("--metadata", type="character", help="Sample metadata TSV"),
  make_option("--sample-col", type="character", default="Sample",
              help="Metadata sample-id col [default: %default]"),
  make_option("--neg-col", type="character", default="is_negative_control",
              help="Metadata neg-control logical col [default: %default]"),
  make_option("--conc-col", type="character", default="DNA_conc",
              help="Metadata DNA-concentration col [default: %default]"),
  make_option("--pos-col", type="character", default="is_positive_control",
              help="Metadata positive-control logical col [default: %default]"),
  make_option("--type-col", type="character", default="Type_Group",
              help="Metadata sample-type col [default: %default]"),
  make_option("--sample-types", type="character", default="BAL,Bronchial Brush,Oral Rinse",
              help="Comma-separated sample types to test [default: %default]"),
  make_option("--threshold-default", type="numeric", default=0.1,
              help="Decontam p threshold for primary flag [default: %default]"),
  make_option("--threshold-aggressive", type="numeric", default=0.5,
              help="Decontam p threshold for aggressive flag [default: %default]"),
  make_option("--pooled", type="character", default=NULL,
              help="Optional: pooled decontam_per_asv_scores.tsv for comparison"),
  make_option("--outdir", type="character", default="decontam_by_type_output",
              help="Output directory [default: %default]"),
  make_option("--combine-mode", type="character", default="min",
              help="How to combine per-type p-values: 'min' (flag if any type) or 'fisher' [default: %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$counts) || is.null(opt$metadata)) {
  stop("--counts and --metadata are required")
}

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)
log_lines <- c()
log <- function(msg) { cat(msg, "\n"); log_lines[length(log_lines)+1] <<- msg }

log("=== run_decontam_by_sample_type.R ===")
log(paste("counts:", opt$counts))
log(paste("metadata:", opt$metadata))
log(paste("sample-types:", opt$`sample-types`))
log(paste("combine-mode:", opt$`combine-mode`))
log("")

# --- Load data ----------------------------------------------------------
log("Loading inputs...")
counts <- read.table(opt$counts, sep="\t", header=TRUE, row.names=1,
                     check.names=FALSE, comment.char="", quote="")
meta <- read.table(opt$metadata, sep="\t", header=TRUE, check.names=FALSE,
                   comment.char="", quote="")
log(paste("  Counts:", nrow(counts), "ASVs x", ncol(counts), "samples"))
log(paste("  Metadata:", nrow(meta), "rows"))

# Align metadata to count table
sample_ids <- colnames(counts)
meta_match <- meta[match(sample_ids, meta[[opt$`sample-col`]]), , drop=FALSE]

# Strip positive controls from analysis
pos <- as.logical(meta_match[[opt$`pos-col`]])
pos[is.na(pos)] <- FALSE
if (any(pos)) {
  log(paste("Excluding", sum(pos), "positive controls."))
  counts <- counts[, !pos, drop=FALSE]
  meta_match <- meta_match[!pos, , drop=FALSE]
  sample_ids <- colnames(counts)
}

neg <- as.logical(meta_match[[opt$`neg-col`]])
neg[is.na(neg)] <- FALSE
log(paste("Negative controls:", sum(neg)))
if (sum(neg) < 2) stop("Need at least 2 negative controls.")

conc <- suppressWarnings(as.numeric(meta_match[[opt$`conc-col`]]))
sample_type <- as.character(meta_match[[opt$`type-col`]])

# --- Per-sample-type decontam -------------------------------------------
sample_types <- strsplit(opt$`sample-types`, ",", fixed=TRUE)[[1]]
sample_types <- trimws(sample_types)

# Accumulator: per-ASV per-type p-values (prevalence + frequency)
asv_ids <- rownames(counts)
all_results <- list()

for (st in sample_types) {
  log("")
  log(paste("=== Sample type:", st, "==="))

  # Subset: this sample type + all negative controls
  in_type <- (sample_type == st) & !is.na(sample_type)
  keep <- in_type | neg
  if (sum(in_type) < 5) {
    log(paste("  WARN: only", sum(in_type), "samples of type '", st, "'; skipping."))
    next
  }
  if (sum(neg) < 2) {
    log(paste("  WARN: too few negatives; skipping."))
    next
  }
  log(paste("  Biological samples:", sum(in_type)))
  log(paste("  Negative controls used:", sum(neg)))

  sub_counts <- counts[, keep, drop=FALSE]
  sub_neg <- neg[keep]
  sub_conc <- conc[keep]

  # Drop ASVs that are all-zero in this subset
  asv_keep <- rowSums(sub_counts) > 0
  sub_counts_nz <- sub_counts[asv_keep, , drop=FALSE]
  log(paste("  ASVs with non-zero reads in this subset:", nrow(sub_counts_nz)))

  mat <- t(as.matrix(sub_counts_nz))
  mode(mat) <- "integer"

  # Prevalence (this sample type vs negatives)
  prev <- isContaminant(mat, method="prevalence", neg=sub_neg,
                        threshold=opt$`threshold-default`)
  log(paste("  Prevalence flagged at p<", opt$`threshold-default`, ":",
            sum(prev$contaminant, na.rm=TRUE)))

  # Frequency (only on biological samples of this type, w/ DNA_conc)
  freq <- NULL
  bio_keep <- in_type[keep] & !is.na(sub_conc) & sub_conc > 0
  n_freq <- sum(bio_keep)
  if (n_freq >= 5) {
    mat_freq <- mat[bio_keep, , drop=FALSE]
    conc_freq <- sub_conc[bio_keep]
    # drop ASVs all-zero in this subset
    asv_keep_freq <- colSums(mat_freq) > 0
    if (sum(asv_keep_freq) > 0) {
      mat_freq <- mat_freq[, asv_keep_freq, drop=FALSE]
      freq <- isContaminant(mat_freq, method="frequency", conc=conc_freq,
                            threshold=opt$`threshold-default`)
      log(paste("  Frequency flagged at p<", opt$`threshold-default`, ":",
                sum(freq$contaminant, na.rm=TRUE),
                "(n biological samples with DNA_conc:", n_freq, ")"))
    } else {
      log("  No ASVs left for frequency test.")
    }
  } else {
    log(paste("  Too few samples with DNA_conc for frequency test (n=", n_freq, ")", sep=""))
  }

  # Combine prev + freq for this sample type, per ASV
  per_asv <- data.frame(
    ASV_ID = asv_ids,
    prev_p = NA_real_,
    prev_flag = NA,
    freq_p = NA_real_,
    freq_flag = NA,
    combined_min_p = NA_real_,
    flag = NA,
    stringsAsFactors = FALSE
  )
  rownames(per_asv) <- asv_ids
  per_asv[rownames(prev), "prev_p"] <- prev$p
  per_asv[rownames(prev), "prev_flag"] <- prev$contaminant
  if (!is.null(freq)) {
    per_asv[rownames(freq), "freq_p"] <- freq$p
    per_asv[rownames(freq), "freq_flag"] <- freq$contaminant
  }
  per_asv$combined_min_p <- pmin(per_asv$prev_p, per_asv$freq_p, na.rm=TRUE)
  per_asv$flag <- per_asv$combined_min_p < opt$`threshold-default`
  per_asv$flag[is.na(per_asv$flag)] <- FALSE

  colnames(per_asv)[-1] <- paste0(st, "_", colnames(per_asv)[-1])
  colnames(per_asv) <- gsub(" ", "_", colnames(per_asv))
  all_results[[st]] <- per_asv
}

# --- Combine per-sample-type results into single per-ASV table ---------
log("")
log("=== Combining per-type results ===")

# Start with all ASVs
out <- data.frame(ASV_ID = asv_ids, stringsAsFactors=FALSE)
for (st in names(all_results)) {
  out <- merge(out, all_results[[st]], by="ASV_ID", all.x=TRUE, sort=FALSE)
}

# For each ASV, take the minimum combined p across sample types
p_cols <- grep("_combined_min_p$", colnames(out), value=TRUE)
out$within_type_min_p <- apply(out[, p_cols, drop=FALSE], 1, function(x) {
  x_clean <- x[!is.na(x)]
  if (length(x_clean) == 0) NA_real_ else min(x_clean)
})

# Fisher's combined p (alternative method)
if (opt$`combine-mode` == "fisher") {
  out$fisher_p <- apply(out[, p_cols, drop=FALSE], 1, function(x) {
    x_clean <- x[!is.na(x) & x > 0]
    if (length(x_clean) == 0) return(NA_real_)
    if (length(x_clean) == 1) return(x_clean)
    chi2 <- -2 * sum(log(x_clean))
    df <- 2 * length(x_clean)
    pchisq(chi2, df = df, lower.tail = FALSE)
  })
  out$within_type_combined_p <- out$fisher_p
} else {
  out$within_type_combined_p <- out$within_type_min_p
}

# Count how many sample types flagged each ASV at default threshold
flag_cols <- grep("_flag$", colnames(out), value=TRUE)
flag_cols <- flag_cols[!grepl("prev|freq", flag_cols)]  # just the per-type combined flags
out$n_types_flagged_at_0.1 <- apply(out[, flag_cols, drop=FALSE], 1, function(x) {
  sum(x == TRUE, na.rm=TRUE)
})

# Final flags
out$flag_within_type_at_0.1 <- out$within_type_combined_p < 0.1
out$flag_within_type_at_0.1[is.na(out$flag_within_type_at_0.1)] <- FALSE
out$flag_within_type_at_0.5 <- out$within_type_combined_p < 0.5
out$flag_within_type_at_0.5[is.na(out$flag_within_type_at_0.5)] <- FALSE

# Sort by within_type_combined_p
out <- out[order(out$within_type_combined_p), ]

# --- Write outputs --------------------------------------------------------
write.table(out, file.path(opt$outdir, "decontam_by_sample_type_per_asv.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)
writeLines(out$ASV_ID[out$flag_within_type_at_0.1],
           file.path(opt$outdir, "decontam_within_type_flagged_at_0.1.txt"))
writeLines(out$ASV_ID[out$flag_within_type_at_0.5],
           file.path(opt$outdir, "decontam_within_type_flagged_at_0.5.txt"))

# Summary
summ <- data.frame(
  threshold = c("0.1 (default)", "0.5 (aggressive)"),
  n_flagged = c(sum(out$flag_within_type_at_0.1), sum(out$flag_within_type_at_0.5)),
  flagged_in_multiple_types = c(
    sum(out$n_types_flagged_at_0.1 >= 2 & out$flag_within_type_at_0.1),
    sum(out$n_types_flagged_at_0.1 >= 2 & out$flag_within_type_at_0.5)
  ),
  flagged_in_single_type = c(
    sum(out$n_types_flagged_at_0.1 == 1 & out$flag_within_type_at_0.1),
    sum(out$n_types_flagged_at_0.1 == 1 & out$flag_within_type_at_0.5)
  )
)
write.table(summ, file.path(opt$outdir, "decontam_by_sample_type_summary.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

log(paste("Total ASVs flagged within-type at p<0.1:", sum(out$flag_within_type_at_0.1)))
log(paste("  flagged in multiple sample types:",
          sum(out$n_types_flagged_at_0.1 >= 2 & out$flag_within_type_at_0.1)))
log(paste("  flagged in only one sample type:",
          sum(out$n_types_flagged_at_0.1 == 1 & out$flag_within_type_at_0.1)))
log(paste("Total ASVs flagged within-type at p<0.5:", sum(out$flag_within_type_at_0.5)))

# --- Pooled comparison (if provided) ------------------------------------
if (!is.null(opt$pooled)) {
  log("")
  log("=== Comparing within-type to pooled decontam ===")
  pooled <- read.table(opt$pooled, sep="\t", header=TRUE, check.names=FALSE)
  comp <- merge(
    pooled[, c("ASV_ID", "combined_min_p", "flag_at_0.1", "flag_at_0.5")],
    out[, c("ASV_ID", "within_type_combined_p",
            "n_types_flagged_at_0.1",
            "flag_within_type_at_0.1", "flag_within_type_at_0.5")],
    by="ASV_ID", all=TRUE
  )
  colnames(comp) <- c("ASV_ID", "pooled_p", "pooled_flag_0.1", "pooled_flag_0.5",
                      "within_p", "n_types_flagged", "within_flag_0.1", "within_flag_0.5")
  # 2x2 table at p<0.1
  comp$pooled_flag_0.1[is.na(comp$pooled_flag_0.1)] <- FALSE
  comp$within_flag_0.1[is.na(comp$within_flag_0.1)] <- FALSE
  log("")
  log("Pooled-vs-within at p<0.1 (rows=pooled, cols=within):")
  log(paste(capture.output(print(table(
    pooled=ifelse(comp$pooled_flag_0.1, "flagged", "not"),
    within=ifelse(comp$within_flag_0.1, "flagged", "not")
  ))), collapse="\n"))
  # Sort comparison by min p
  comp$any_p <- pmin(comp$pooled_p, comp$within_p, na.rm=TRUE)
  comp <- comp[order(comp$any_p), ]
  comp$any_p <- NULL
  write.table(comp, file.path(opt$outdir, "decontam_pooled_vs_within_type_comparison.tsv"),
              sep="\t", row.names=FALSE, quote=FALSE)
  log("")
  log(paste("Comparison table written to:",
            file.path(opt$outdir, "decontam_pooled_vs_within_type_comparison.tsv")))
}

writeLines(log_lines, file.path(opt$outdir, "decontam_by_sample_type_log.txt"))
log("")
log("Done.")
