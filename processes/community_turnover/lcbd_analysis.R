#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(adespatial)
  library(optparse)
})

options(stringsAsFactors = FALSE)
option_list <- list(
  make_option("--prepared-dir", type = "character"),
  make_option("--metadata", type = "character"),
  make_option("--outdir", type = "character"),
  make_option("--sample-col", type = "character", default = "sampleID"),
  make_option("--profile-col", type = "character", default = "Cruise"),
  make_option("--date-col", type = "character", default = "date"),
  make_option("--depth-col", type = "character", default = "Depth"),
  make_option("--environmental-compartment-cols", type = "character", default = ""),
  make_option("--permutations", type = "integer", default = 999),
  make_option("--seed", type = "integer", default = 42),
  make_option("--within-profile-enabled", type = "logical", default = TRUE),
  make_option("--within-profile-min-samples", type = "integer", default = 3)
)
args <- parse_args(OptionParser(option_list = option_list), convert_hyphens_to_underscores = TRUE)
if (is.null(args$prepared_dir) || is.null(args$metadata) || is.null(args$outdir)) stop("Required paths are missing", call. = FALSE)
if (args$permutations < 1 || args$within_profile_min_samples < 3) stop("Invalid permutation or within-profile minimum", call. = FALSE)
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
read_tsv <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

clr_table <- read_tsv(file.path(args$prepared_dir, "asv_clr_matrix.tsv"))
sample_ids <- trimws(as.character(clr_table[[1]]))
clr <- as.matrix(clr_table[, -1, drop = FALSE])
storage.mode(clr) <- "double"
rownames(clr) <- sample_ids
metadata <- read_tsv(args$metadata)
required <- c(args$sample_col, args$profile_col, args$date_col, args$depth_col)
missing <- setdiff(required, names(metadata))
if (length(missing)) stop(paste("Metadata lacks required columns:", paste(missing, collapse = ", ")), call. = FALSE)
metadata[[args$sample_col]] <- trimws(as.character(metadata[[args$sample_col]]))
if (anyDuplicated(metadata[[args$sample_col]])) stop("Metadata sample identifiers are not unique", call. = FALSE)
metadata <- metadata[match(sample_ids, metadata[[args$sample_col]]), , drop = FALSE]
valid <- !is.na(metadata[[args$sample_col]]) & !is.na(metadata[[args$profile_col]]) &
  nzchar(trimws(as.character(metadata[[args$profile_col]]))) &
  !is.na(as.Date(metadata[[args$date_col]])) &
  is.finite(suppressWarnings(as.numeric(metadata[[args$depth_col]])))
metadata <- metadata[valid, , drop = FALSE]
clr <- clr[valid, , drop = FALSE]
if (nrow(clr) < 3) stop("Fewer than three samples are eligible for LCBD", call. = FALSE)

# beta.div rejects negative inputs before dispatching to the Euclidean method,
# whereas valid CLR coordinates necessarily contain negative entries. A single
# uniform translation preserves all Euclidean distances. It also preserves the
# package's column-wise permutation null because the same constant remains in
# every entry before and after each within-column permutation.
clr_translation <- if (min(clr) < 0) -min(clr) else 0
clr_for_adespatial <- clr + clr_translation
set.seed(args$seed)
global <- adespatial::beta.div(
  clr_for_adespatial, method = "euclidean", sqrt.D = FALSE,
  nperm = args$permutations, adj = FALSE
)
direct_centered <- sweep(clr, 2, colMeans(clr), "-")
direct_squared <- rowSums(direct_centered^2)
direct_lcbd <- direct_squared / sum(direct_squared)
lcbd_max_abs_difference <- max(abs(as.numeric(global$LCBD) - direct_lcbd))
if (!is.finite(lcbd_max_abs_difference) || lcbd_max_abs_difference > 1e-10) {
  stop("Translated adespatial LCBD does not reproduce direct CLR LCBD", call. = FALSE)
}
global_rows <- data.frame(
  sample_ID = rownames(clr), reference_population = "global", reference_group = "all_samples",
  reference_population_n = nrow(clr), LCBD = as.numeric(global$LCBD),
  LCBD_pvalue = as.numeric(global$p.LCBD),
  adjusted_pvalue = p.adjust(as.numeric(global$p.LCBD), method = "BH"),
  permutation_tested = TRUE, stringsAsFactors = FALSE
)

within_rows <- list()
within_audit <- list()
if (isTRUE(args$within_profile_enabled)) {
  profiles <- unique(as.character(metadata[[args$profile_col]]))
  for (profile in profiles) {
    selected <- which(as.character(metadata[[args$profile_col]]) == profile)
    if (length(selected) < args$within_profile_min_samples) {
      within_audit[[length(within_audit) + 1]] <- data.frame(
        profile_ID = profile, number_of_samples = length(selected), status = "below_minimum_samples"
      )
      next
    }
    centered <- sweep(clr[selected, , drop = FALSE], 2, colMeans(clr[selected, , drop = FALSE]), "-")
    squared <- rowSums(centered^2)
    total <- sum(squared)
    lcbd <- if (total > 0) squared / total else rep(NA_real_, length(selected))
    within_rows[[length(within_rows) + 1]] <- data.frame(
      sample_ID = rownames(clr)[selected], reference_population = "within_profile",
      reference_group = profile, reference_population_n = length(selected), LCBD = lcbd,
      LCBD_pvalue = NA_real_, adjusted_pvalue = NA_real_, permutation_tested = FALSE,
      stringsAsFactors = FALSE
    )
    within_audit[[length(within_audit) + 1]] <- data.frame(
      profile_ID = profile, number_of_samples = length(selected), status = ifelse(total > 0, "completed", "zero_total_variance")
    )
  }
}
results <- if (length(within_rows)) rbind(global_rows, do.call(rbind, within_rows)) else global_rows
metadata_columns <- unique(c(args$sample_col, args$profile_col, args$date_col, args$depth_col,
  Filter(nzchar, strsplit(args$environmental_compartment_cols, "[|,]")[[1]])))
metadata_columns <- intersect(metadata_columns, names(metadata))
metadata_join <- metadata[, metadata_columns, drop = FALSE]
names(metadata_join)[names(metadata_join) == args$sample_col] <- "sample_ID"
results <- merge(results, metadata_join, by = "sample_ID", all.x = TRUE, sort = FALSE)
write_tsv(results, file.path(args$outdir, "lcbd.tsv"))
audit <- if (length(within_audit)) do.call(rbind, within_audit) else data.frame(
  profile_ID = character(), number_of_samples = integer(), status = character()
)
write_tsv(audit, file.path(args$outdir, "lcbd_within_profile_audit.tsv"))
writeLines(capture.output(sessionInfo()), file.path(args$outdir, "lcbd_session_info.txt"))
write_tsv(data.frame(
  method = "adespatial::beta.div on CLR coordinates", distance_geometry = "Aitchison",
  adespatial_input_translation = clr_translation,
  translation_effect = "uniform_translation_preserves_euclidean_distances_and_column_permutation_null",
  direct_CLR_LCBD_max_abs_difference = lcbd_max_abs_difference,
  global_permutations = args$permutations, p_adjustment = "Benjamini-Hochberg",
  seed = args$seed, within_profile_permutation_tested = FALSE,
  within_profile_min_samples = args$within_profile_min_samples,
  adespatial_version = as.character(packageVersion("adespatial")), R_version = R.version.string
), file.path(args$outdir, "lcbd_parameters.tsv"))
