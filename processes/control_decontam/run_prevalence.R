#!/usr/bin/env Rscript
# One arm only: cohort-level prevalence against technical OR biological/sample controls.
suppressPackageStartupMessages({library(decontam); library(optparse)})
opt <- parse_args(OptionParser(option_list=list(
  make_option("--counts", type="character"),
  make_option("--metadata", type="character"),
  make_option("--threshold", type="double", default=0.1),
  make_option("--output", type="character")
)))
counts <- read.delim(opt$counts, row.names=1, check.names=FALSE)
meta <- read.delim(opt$metadata, check.names=FALSE)
stopifnot(identical(colnames(counts), meta$Sample), opt$threshold > 0, opt$threshold < 1)
neg <- as.logical(meta$is_negative)
stopifnot(!anyNA(neg), any(neg), any(!neg))
# Biological QC and exclusion of zero-depth samples occur upstream.
mat <- t(as.matrix(counts))
stopifnot(all(rowSums(mat) > 0))
present <- colSums(mat) > 0
out <- data.frame(ASV_ID=colnames(mat), score=NA_real_, contaminant=FALSE,
                  status="absent_in_arm")
if (any(present)) {
  tested <- mat[, present, drop=FALSE]
  # decontam 1.22 drops matrix dimensions internally with one feature.
  # A zero-only placeholder preserves dimensions without altering prevalence.
  if (ncol(tested) == 1) {
    placeholder <- make.unique(c(colnames(mat), "__absent_placeholder__"))
    placeholder <- tail(placeholder, 1)
    tested <- cbind(tested, 0)
    colnames(tested)[2] <- placeholder
  }
  result <- isContaminant(tested, method="prevalence",
                         neg=neg, threshold=opt$threshold)
  result <- result[rownames(result) %in% out$ASV_ID, , drop=FALSE]
  positions <- match(rownames(result), out$ASV_ID)
  out$score[positions] <- result$p
  out$contaminant[positions] <- !is.na(result$contaminant) & result$contaminant
  out$status[positions] <- ifelse(is.na(result$p), "not_estimable", "tested")
}
write.table(out, opt$output, sep="\t", row.names=FALSE, quote=FALSE, na="NA")
cat(sprintf("samples=%d biological=%d negatives=%d ASVs=%d tested=%d removed=%d threshold=%g\n",
            nrow(mat), sum(!neg), sum(neg), ncol(mat), sum(present),
            sum(out$contaminant), opt$threshold))
