#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(TITAN2)
  library(jsonlite)
})

options(stringsAsFactors = FALSE)
RNGkind("L'Ecuyer-CMRG")

option_list <- list(
  make_option("--prepared-dir", dest = "prepared_dir", type = "character"),
  make_option("--outdir", type = "character"),
  make_option("--min-split", dest = "min_split", type = "integer", default = 5),
  make_option("--permutations", type = "integer", default = 250),
  make_option("--bootstrap-count", dest = "bootstrap_count", type = "integer", default = 500),
  make_option("--seed", type = "integer", default = 42),
  make_option("--imax", type = "logical", default = FALSE),
  make_option("--iv-total", dest = "iv_total", type = "logical", default = FALSE),
  make_option("--purity-cutoff", dest = "purity_cutoff", type = "double", default = 0.95),
  make_option("--reliability-cutoff", dest = "reliability_cutoff", type = "double", default = 0.95),
  make_option("--ncpus", type = "integer", default = 1),
  make_option("--memory", type = "logical", default = FALSE)
)
args <- parse_args(OptionParser(option_list = option_list))
if (is.null(args$prepared_dir) || is.null(args$outdir)) {
  stop("--prepared-dir and --outdir are required", call. = FALSE)
}
if (args$min_split < 3 || args$permutations < 1 || args$bootstrap_count < 1 ||
    args$ncpus < 1 || args$purity_cutoff < 0 || args$purity_cutoff > 1 ||
    args$reliability_cutoff < 0 || args$reliability_cutoff > 1) {
  stop("Invalid TITAN parameter value", call. = FALSE)
}

prepared_dir <- normalizePath(args$prepared_dir)
outdir <- args$outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "variables"), recursive = TRUE, showWarnings = FALSE)
manifest <- read_tsv(file.path(prepared_dir, "titan_input_manifest.tsv"), show_col_types = FALSE)

clean_names <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("%", "pct", x, fixed = TRUE)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  tolower(x)
}

value_or_na <- function(row, candidates) {
  for (candidate in candidates) {
    if (candidate %in% names(row)) return(as.numeric(row[[candidate]]))
  }
  NA_real_
}

matrix_column_or_na <- function(value, column, length_out) {
  if (is.null(value)) return(rep(NA_real_, length_out))
  matrix_value <- as.matrix(value)
  if (nrow(matrix_value) != length_out || ncol(matrix_value) < column) {
    return(rep(NA_real_, length_out))
  }
  as.numeric(matrix_value[, column])
}

status_rows <- list()
for (i in seq_len(nrow(manifest))) {
  entry <- manifest[i, ]
  variable <- as.character(entry$environmental_variable)
  slug <- as.character(entry$variable_slug)
  variable_out <- file.path(outdir, "variables", slug)
  dir.create(variable_out, recursive = TRUE, showWarnings = FALSE)
  if (as.character(entry$status) != "ready") {
    status_rows[[length(status_rows) + 1]] <- data.frame(
      environmental_variable = variable, variable_slug = slug,
      status = "skipped", status_reason = as.character(entry$status_reason),
      number_of_samples = as.integer(entry$number_of_samples),
      number_of_ASVs_tested = as.integer(entry$number_of_ASVs_tested),
      titan2_version = as.character(packageVersion("TITAN2")),
      stringsAsFactors = FALSE
    )
    next
  }

  tryCatch({
    variable_in <- file.path(prepared_dir, "variables", slug)
    env_table <- read_tsv(file.path(variable_in, "environment.tsv"), show_col_types = FALSE)
    taxa_table <- read_tsv(file.path(variable_in, "taxa_relative_abundance.tsv"), show_col_types = FALSE)
    taxon_metadata <- read_tsv(file.path(variable_in, "taxon_metadata.tsv"), show_col_types = FALSE)
    sample_col <- names(env_table)[1]
    if (names(taxa_table)[1] != sample_col) stop("Environment and taxa sample-ID columns differ")
    if (anyDuplicated(env_table[[sample_col]]) || anyDuplicated(taxa_table[[sample_col]])) {
      stop("Duplicate sample identifiers in prepared TITAN input")
    }
    if (!identical(as.character(env_table[[sample_col]]), as.character(taxa_table[[sample_col]]))) {
      stop("Environment and taxa sample order differs after explicit identifier matching")
    }
    env <- as.numeric(env_table$environmental_value)
    taxa <- as.data.frame(taxa_table[, -1, drop = FALSE], check.names = FALSE)
    taxa[] <- lapply(taxa, as.numeric)
    if (any(!is.finite(env)) || any(!is.finite(as.matrix(taxa)))) {
      stop("Non-finite values in prepared TITAN input")
    }

    set.seed(args$seed)
    result <- TITAN2::titan(
      env = env,
      txa = taxa,
      minSplt = args$min_split,
      numPerm = args$permutations,
      boot = TRUE,
      nBoot = args$bootstrap_count,
      imax = args$imax,
      ivTot = args$iv_total,
      pur.cut = args$purity_cutoff,
      rel.cut = args$reliability_cutoff,
      ncpus = args$ncpus,
      memory = args$memory,
      messaging = TRUE
    )

    spp <- as.data.frame(result$sppmax, check.names = FALSE)
    spp$ASV_ID <- rownames(spp)
    names(spp) <- clean_names(names(spp))
    metadata_index <- match(spp$asv_id, taxon_metadata$ASV_ID)
    response_code <- as.integer(spp$maxgrp)
    observed_cp <- if (isTRUE(args$imax)) spp$ienv_cp else spp$zenv_cp
    taxon_out <- data.frame(
      environmental_variable = variable,
      ASV_ID = spp$asv_id,
      taxonomy = taxon_metadata$taxonomy[metadata_index],
      change_point = as.numeric(observed_cp),
      response_direction = ifelse(response_code == 1, "z-", ifelse(response_code == 2, "z+", NA_character_)),
      indicator_score = as.numeric(spp$indval),
      z_score = as.numeric(spp$zscore),
      purity = as.numeric(spp$purity),
      reliability = as.numeric(spp$reliability),
      bootstrap_cp_q05 = value_or_na(spp, c("5pct", "pct5", "q05")),
      bootstrap_cp_q10 = value_or_na(spp, c("10pct", "pct10", "q10")),
      bootstrap_cp_q50 = value_or_na(spp, c("50pct", "pct50", "q50")),
      bootstrap_cp_q90 = value_or_na(spp, c("90pct", "pct90", "q90")),
      bootstrap_cp_q95 = value_or_na(spp, c("95pct", "pct95", "q95")),
      number_of_samples_used = length(env),
      prevalence = as.numeric(taxon_metadata$prevalence[metadata_index]),
      mean_relative_abundance = as.numeric(taxon_metadata$mean_relative_abundance[metadata_index]),
      occurrence_frequency = as.integer(spp$freq),
      response_group_code = response_code,
      indval_change_point = as.numeric(spp$ienv_cp),
      zscore_change_point = as.numeric(spp$zenv_cp),
      observed_indval_probability = as.numeric(spp$obsiv_prob),
      median_bootstrap_z_score = as.numeric(spp$z_median),
      titan_filter_code = as.integer(spp$filter),
      passes_purity = as.numeric(spp$purity) >= args$purity_cutoff,
      passes_reliability = as.numeric(spp$reliability) >= args$reliability_cutoff,
      passes_purity_and_reliability = as.numeric(spp$purity) >= args$purity_cutoff &
        as.numeric(spp$reliability) >= args$reliability_cutoff,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    write_tsv(taxon_out, file.path(variable_out, "titan_taxon_results.tsv"), na = "")

    cp <- as.data.frame(result$sumz.cp, check.names = FALSE)
    cp$threshold_type <- rownames(cp)
    names(cp) <- clean_names(names(cp))
    cp_index <- function(name, column) {
      row <- cp[cp$threshold_type == name, , drop = FALSE]
      if (!nrow(row) || !column %in% names(row)) return(NA_real_)
      as.numeric(row[[column]][1])
    }
    community_out <- data.frame(
      environmental_variable = variable,
      z_minus_threshold = cp_index("sumz-", "cp"),
      z_plus_threshold = cp_index("sumz+", "cp"),
      z_minus_bootstrap_q05 = cp_index("sumz-", "0_05"),
      z_minus_bootstrap_q10 = cp_index("sumz-", "0_10"),
      z_minus_bootstrap_q50 = cp_index("sumz-", "0_50"),
      z_minus_bootstrap_q90 = cp_index("sumz-", "0_90"),
      z_minus_bootstrap_q95 = cp_index("sumz-", "0_95"),
      z_plus_bootstrap_q05 = cp_index("sumz+", "0_05"),
      z_plus_bootstrap_q10 = cp_index("sumz+", "0_10"),
      z_plus_bootstrap_q50 = cp_index("sumz+", "0_50"),
      z_plus_bootstrap_q90 = cp_index("sumz+", "0_90"),
      z_plus_bootstrap_q95 = cp_index("sumz+", "0_95"),
      filtered_z_minus_threshold = cp_index("fsumz-", "cp"),
      filtered_z_plus_threshold = cp_index("fsumz+", "cp"),
      filtered_z_minus_bootstrap_q05 = cp_index("fsumz-", "0_05"),
      filtered_z_minus_bootstrap_q95 = cp_index("fsumz-", "0_95"),
      filtered_z_plus_bootstrap_q05 = cp_index("fsumz+", "0_05"),
      filtered_z_plus_bootstrap_q95 = cp_index("fsumz+", "0_95"),
      number_of_samples = length(env),
      number_of_ASVs_tested = ncol(taxa),
      number_of_pure_reliable_z_minus = sum(taxon_out$response_direction == "z-" & taxon_out$passes_purity_and_reliability, na.rm = TRUE),
      number_of_pure_reliable_z_plus = sum(taxon_out$response_direction == "z+" & taxon_out$passes_purity_and_reliability, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    write_tsv(community_out, file.path(variable_out, "titan_community_thresholds.tsv"), na = "")

    ivz <- as.matrix(result$ivz)
    curve_length <- length(result$envcls)
    curve_out <- data.frame(
      environmental_variable = variable,
      environmental_value = as.numeric(result$envcls),
      sum_z_minus = matrix_column_or_na(ivz, 1, curve_length),
      sum_z_plus = matrix_column_or_na(ivz, 2, curve_length),
      filtered_sum_z_minus = matrix_column_or_na(result$ivz.f, 1, curve_length),
      filtered_sum_z_plus = matrix_column_or_na(result$ivz.f, 2, curve_length)
    )
    write_tsv(curve_out, file.path(variable_out, "titan_community_response_curve.tsv"), na = "")

    bootstrap_length <- nrow(as.matrix(result$maxSumz))
    bootstrap_out <- data.frame(
      environmental_variable = variable,
      bootstrap_replicate = seq_len(bootstrap_length),
      z_minus_threshold = matrix_column_or_na(result$maxSumz, 1, bootstrap_length),
      z_plus_threshold = matrix_column_or_na(result$maxSumz, 2, bootstrap_length),
      filtered_z_minus_threshold = matrix_column_or_na(result$maxFsumz, 1, bootstrap_length),
      filtered_z_plus_threshold = matrix_column_or_na(result$maxFsumz, 2, bootstrap_length)
    )
    write_tsv(bootstrap_out, file.path(variable_out, "titan_community_bootstrap_thresholds.tsv"), na = "")
    saveRDS(result, file.path(variable_out, "titan_result.rds"))

    status_rows[[length(status_rows) + 1]] <- data.frame(
      environmental_variable = variable, variable_slug = slug,
      status = "completed", status_reason = "",
      number_of_samples = length(env), number_of_ASVs_tested = ncol(taxa),
      titan2_version = as.character(packageVersion("TITAN2")),
      stringsAsFactors = FALSE
    )
  }, error = function(error) {
    writeLines(conditionMessage(error), file.path(variable_out, "titan_error.txt"))
    status_rows[[length(status_rows) + 1]] <<- data.frame(
      environmental_variable = variable, variable_slug = slug,
      status = "failed", status_reason = conditionMessage(error),
      number_of_samples = as.integer(entry$number_of_samples),
      number_of_ASVs_tested = as.integer(entry$number_of_ASVs_tested),
      titan2_version = as.character(packageVersion("TITAN2")),
      stringsAsFactors = FALSE
    )
  })
}

status <- do.call(rbind, status_rows)
write_tsv(status, file.path(outdir, "titan_variable_status.tsv"), na = "")
analysis_config <- list(
  min_split = args$min_split,
  permutations = args$permutations,
  bootstrap_count = args$bootstrap_count,
  seed = args$seed,
  imax = args$imax,
  iv_total = args$iv_total,
  purity_cutoff = args$purity_cutoff,
  reliability_cutoff = args$reliability_cutoff,
  ncpus = args$ncpus,
  memory = args$memory,
  bootstrap_unit = "sample",
  TITAN2_version = as.character(packageVersion("TITAN2")),
  R_version = R.version.string
)
write_json(analysis_config, file.path(outdir, "titan_analysis_config.json"), pretty = TRUE, auto_unbox = TRUE)
capture.output(sessionInfo(), file = file.path(outdir, "titan_session_info.txt"))
