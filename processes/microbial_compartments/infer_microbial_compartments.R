#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(cluster)
  library(optparse)
})

options(stringsAsFactors = FALSE)
RNGkind("L'Ecuyer-CMRG")

option_list <- list(
  make_option("--prepared-dir", type = "character"),
  make_option("--outdir", type = "character"),
  make_option("--sample-col", type = "character", default = "sampleID"),
  make_option("--k-min", type = "integer", default = 2),
  make_option("--k-max", type = "integer", default = 10),
  make_option("--min-cluster-size", type = "integer", default = 5),
  make_option("--min-cluster-fraction", type = "double", default = 0.02),
  # Retained for command-line compatibility and reporting only. Silhouette is
  # descriptive because overlapping ecological states are expected.
  make_option("--min-mean-silhouette", type = "double", default = 0.25),
  make_option("--min-stability-ari", type = "double", default = 0.70),
  make_option("--min-cluster-jaccard", type = "double", default = 0.75),
  make_option("--stability-replicates", type = "integer", default = 200),
  make_option("--stability-sample-fraction", type = "double", default = 0.80),
  make_option("--stability-primary-quantile", type = "double", default = 0.25),
  make_option("--stability-near-tie-tolerance", type = "double", default = 0.02),
  make_option("--blocked-robustness-enabled", type = "logical", default = TRUE),
  make_option("--prediction-strength-enabled", type = "logical", default = TRUE),
  make_option("--hierarchical-enabled", type = "logical", default = TRUE),
  make_option("--seed", type = "integer", default = 42),
  make_option("--ncpus", type = "integer", default = 1)
)
args <- parse_args(
  OptionParser(option_list = option_list),
  convert_hyphens_to_underscores = TRUE
)

read_tsv <- function(path, ...) {
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, ...)
}

write_tsv <- function(table, path, na = "") {
  write.table(table, path, sep = "\t", quote = FALSE, row.names = FALSE, na = na)
}

json_escape <- function(value) {
  value <- gsub("\\\\", "\\\\\\\\", as.character(value), fixed = TRUE)
  gsub('"', '\\\\"', value, fixed = TRUE)
}

json_value <- function(value) {
  if (length(value) > 1) {
    return(paste0("[", paste(vapply(value, json_value, character(1)), collapse = ", "), "]"))
  }
  if (is.logical(value)) return(tolower(as.character(value)))
  if (is.numeric(value)) return(as.character(value))
  paste0('"', json_escape(value), '"')
}

write_json_object <- function(values, path) {
  entries <- vapply(names(values), function(name) {
    paste0('  "', json_escape(name), '": ', json_value(values[[name]]))
  }, character(1))
  writeLines(c("{", paste0(entries, collapse = ",\n"), "}"), path)
}
if (is.null(args$prepared_dir) || is.null(args$outdir)) {
  stop("--prepared-dir and --outdir are required", call. = FALSE)
}
if (args$k_min < 2 || args$k_max < args$k_min) stop("Invalid K range", call. = FALSE)
if (args$min_cluster_size < 2) stop("--min-cluster-size must be at least 2", call. = FALSE)
if (args$min_cluster_fraction < 0 || args$min_cluster_fraction >= 1) stop("Invalid minimum cluster fraction", call. = FALSE)
if (args$min_mean_silhouette < -1 || args$min_mean_silhouette > 1) stop("Invalid silhouette threshold", call. = FALSE)
if (args$min_stability_ari < -1 || args$min_stability_ari > 1) stop("Invalid ARI threshold", call. = FALSE)
if (args$min_cluster_jaccard < 0 || args$min_cluster_jaccard > 1) stop("Invalid cluster Jaccard threshold", call. = FALSE)
if (args$stability_replicates < 1 || args$ncpus < 1) stop("Replicates and CPUs must be positive", call. = FALSE)
if (args$stability_sample_fraction <= 0 || args$stability_sample_fraction > 1) stop("Invalid subsampling fraction", call. = FALSE)
if (args$stability_primary_quantile <= 0 || args$stability_primary_quantile >= 0.5) stop("Stability primary quantile must be between zero and 0.5", call. = FALSE)
if (args$stability_near_tie_tolerance < 0 || args$stability_near_tie_tolerance > 1) stop("Invalid stability near-tie tolerance", call. = FALSE)

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
distance_table <- read_tsv(
  file.path(args$prepared_dir, "aitchison_distance_matrix.tsv")
)
if (ncol(distance_table) < 4) stop("Aitchison matrix contains fewer than three samples", call. = FALSE)
sample_ids <- as.character(distance_table[[1]])
distance_matrix <- as.matrix(distance_table[, -1, drop = FALSE])
storage.mode(distance_matrix) <- "double"
rownames(distance_matrix) <- sample_ids
if (!identical(colnames(distance_matrix), sample_ids)) stop("Aitchison columns do not match explicit sample-ID order", call. = FALSE)
if (anyDuplicated(sample_ids)) stop("Duplicate sample identifiers in Aitchison matrix", call. = FALSE)
if (any(!is.finite(distance_matrix)) || any(distance_matrix < 0)) stop("Invalid Aitchison distances", call. = FALSE)
if (!isTRUE(all.equal(distance_matrix, t(distance_matrix), tolerance = 1e-10))) stop("Aitchison matrix is not symmetric", call. = FALSE)
if (max(abs(diag(distance_matrix))) > 1e-10) stop("Aitchison matrix diagonal is not zero", call. = FALSE)

block_path <- file.path(args$prepared_dir, "microbial_clustering_stability_blocks.tsv")
if (!file.exists(block_path)) stop("Prepared stability-block table is missing", call. = FALSE)
block_table <- read_tsv(block_path)
if (!all(c(args$sample_col, "stability_block", "stability_stratum") %in% names(block_table))) {
  stop("Stability-block table lacks required columns", call. = FALSE)
}
if (anyDuplicated(block_table[[args$sample_col]])) stop("Duplicate sample IDs in stability-block table", call. = FALSE)
block_index <- match(sample_ids, as.character(block_table[[args$sample_col]]))
if (anyNA(block_index)) stop("Stability-block table does not cover every clustering sample", call. = FALSE)
stability_blocks <- as.character(block_table$stability_block[block_index])
stability_strata <- as.character(block_table$stability_stratum[block_index])
if (anyNA(stability_blocks) || any(!nzchar(stability_blocks))) stop("Invalid stability block identifiers", call. = FALSE)
if (anyNA(stability_strata) || any(!nzchar(stability_strata))) stop("Invalid stability stratum identifiers", call. = FALSE)
unique_blocks <- sort(unique(stability_blocks))
if (length(unique_blocks) < 3) stop("Stability validation requires at least three resampling units", call. = FALSE)
block_stratum_counts <- tapply(stability_strata, stability_blocks, function(values) length(unique(values)))
if (any(block_stratum_counts != 1)) stop("Each stability block must belong to one robustness stratum", call. = FALSE)
block_strata <- vapply(unique_blocks, function(block) unique(stability_strata[stability_blocks == block])[1], character(1))
names(block_strata) <- unique_blocks

n_samples <- length(sample_ids)
k_values <- seq.int(args$k_min, min(args$k_max, n_samples - 1L))
if (!length(k_values)) stop("No candidate K is valid for the sample count", call. = FALSE)
minimum_cluster_n <- max(args$min_cluster_size, ceiling(args$min_cluster_fraction * n_samples))

adjusted_rand <- function(labels_a, labels_b) {
  table_ab <- table(labels_a, labels_b)
  choose2 <- function(x) x * (x - 1) / 2
  sum_cells <- sum(choose2(table_ab))
  sum_rows <- sum(choose2(rowSums(table_ab)))
  sum_cols <- sum(choose2(colSums(table_ab)))
  total_pairs <- choose2(sum(table_ab))
  if (total_pairs <= 0) return(NA_real_)
  expected <- sum_rows * sum_cols / total_pairs
  maximum <- 0.5 * (sum_rows + sum_cols)
  denominator <- maximum - expected
  if (abs(denominator) < .Machine$double.eps) return(ifelse(sum_cells == maximum, 1, 0))
  (sum_cells - expected) / denominator
}

canonicalize_pam <- function(fit, ids) {
  medoid_ids <- ids[fit$id.med]
  old_labels <- seq_along(medoid_ids)
  ordered_old <- old_labels[order(medoid_ids)]
  map <- setNames(paste0("MC", seq_along(ordered_old)), ordered_old)
  list(
    labels = unname(map[as.character(fit$clustering)]),
    medoids = data.frame(
      microbial_compartment = paste0("MC", seq_along(ordered_old)),
      medoid_sample_ID = medoid_ids[ordered_old],
      stringsAsFactors = FALSE
    )
  )
}

canonicalize_hierarchical <- function(labels, ids) {
  groups <- sort(unique(labels))
  keys <- vapply(groups, function(group) min(ids[labels == group]), character(1))
  ordered <- groups[order(keys)]
  map <- setNames(paste0("HC", seq_along(ordered)), ordered)
  unname(map[as.character(labels)])
}

silhouette_values <- function(labels, matrix) {
  if (length(unique(labels)) < 2 || length(unique(labels)) >= length(labels)) return(rep(NA_real_, length(labels)))
  as.numeric(cluster::silhouette(as.integer(factor(labels)), as.dist(matrix))[, "sil_width"])
}

prediction_strength <- function(train_fit, train_indices, test_indices, k) {
  if (!isTRUE(args$prediction_strength_enabled) || length(test_indices) <= k) return(NA_real_)
  test_distance <- distance_matrix[test_indices, test_indices, drop = FALSE]
  test_fit <- tryCatch(cluster::pam(as.dist(test_distance), k = k, diss = TRUE), error = function(error) NULL)
  if (is.null(test_fit)) return(NA_real_)
  medoid_global <- train_indices[train_fit$id.med]
  medoid_clusters <- train_fit$clustering[train_fit$id.med]
  test_to_medoids <- distance_matrix[test_indices, medoid_global, drop = FALSE]
  nearest_medoid <- apply(test_to_medoids, 1, which.min)
  predicted <- medoid_clusters[nearest_medoid]
  independent <- test_fit$clustering
  scores <- vapply(seq_len(k), function(group) {
    members <- which(predicted == group)
    if (length(members) < 2) return(0)
    pairs <- utils::combn(members, 2)
    mean(independent[pairs[1, ]] == independent[pairs[2, ]])
  }, numeric(1))
  min(scores)
}

cluster_recovery <- function(full_labels, subset_labels, selected, medoid_ids, k, replicate_index) {
  full_groups <- sort(unique(full_labels))
  sub_groups <- sort(unique(subset_labels))
  rows <- lapply(full_groups, function(group) {
    full_members <- selected[full_labels[selected] == group]
    jaccards <- vapply(sub_groups, function(sub_group) {
      sub_members <- selected[subset_labels == sub_group]
      union_n <- length(union(full_members, sub_members))
      if (union_n == 0) return(0)
      length(intersect(full_members, sub_members)) / union_n
    }, numeric(1))
    best_index <- if (length(jaccards)) which.max(jaccards) else NA_integer_
    medoid_id <- medoid_ids$medoid_sample_ID[medoid_ids$microbial_compartment == group][1]
    medoid_global <- match(medoid_id, sample_ids)
    medoid_subset_position <- match(medoid_global, selected)
    medoid_recovered <- if (is.na(medoid_subset_position) || is.na(best_index)) {
      NA
    } else {
      subset_labels[medoid_subset_position] == sub_groups[best_index]
    }
    data.frame(
      K = k,
      replicate = replicate_index,
      microbial_compartment = group,
      selected_members = length(full_members),
      best_jaccard = if (length(jaccards)) max(jaccards) else 0,
      medoid_sample_ID = medoid_id,
      medoid_included = !is.na(medoid_subset_position),
      medoid_recovered = medoid_recovered,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

primary_stability_for_k <- function(k, full_labels, medoid_ids) {
  worker <- function(replicate_index) {
    set.seed(args$seed + k * 100000L + replicate_index)
    selected <- sort(unlist(lapply(unique_blocks, function(block) {
      members <- which(stability_blocks == block)
      retain_n <- max(1L, floor(args$stability_sample_fraction * length(members)))
      sort(sample(members, min(retain_n, length(members)), replace = FALSE))
    }), use.names = FALSE))
    held_out <- setdiff(seq_len(n_samples), selected)
    if (length(selected) <= k) return(NULL)
    sub_distance <- distance_matrix[selected, selected, drop = FALSE]
    fit <- tryCatch(cluster::pam(as.dist(sub_distance), k = k, diss = TRUE), error = function(error) NULL)
    if (is.null(fit)) return(NULL)
    recovery <- cluster_recovery(full_labels, fit$clustering, selected, medoid_ids, k, replicate_index)
    data.frame_row <- data.frame(
      K = k,
      replicate = replicate_index,
      validation_role = "primary_selection",
      resampling_unit = "samples_within_every_block",
      selected_blocks = length(unique(stability_blocks[selected])),
      held_out_blocks = 0L,
      selected_samples = length(selected),
      held_out_samples = length(held_out),
      adjusted_rand_index = adjusted_rand(full_labels[selected], fit$clustering),
      minimum_cluster_jaccard = min(recovery$best_jaccard),
      mean_cluster_jaccard = mean(recovery$best_jaccard),
      prediction_strength = NA_real_,
      stringsAsFactors = FALSE
    )
    list(replicate = data.frame_row, clusters = recovery)
  }
  results <- if (args$ncpus > 1 && .Platform$OS.type != "windows") {
    parallel::mclapply(
      seq_len(args$stability_replicates), worker,
      mc.cores = args$ncpus, mc.preschedule = TRUE, mc.set.seed = FALSE
    )
  } else {
    lapply(seq_len(args$stability_replicates), worker)
  }
  results <- Filter(Negate(is.null), results)
  list(
    replicates = if (length(results)) do.call(rbind, lapply(results, `[[`, "replicate")) else data.frame(),
    clusters = if (length(results)) do.call(rbind, lapply(results, `[[`, "clusters")) else data.frame()
  )
}

season_balanced_block_stability_for_k <- function(k, full_labels) {
  empty_replicates <- data.frame(
    K = integer(), replicate = integer(), validation_role = character(),
    resampling_unit = character(), selected_blocks = integer(), held_out_blocks = integer(),
    selected_samples = integer(), held_out_samples = integer(),
    adjusted_rand_index = double(), prediction_strength = double(),
    stringsAsFactors = FALSE
  )
  empty_strata <- data.frame(
    K = integer(), replicate = integer(), stability_stratum = character(),
    available_blocks = integer(), selected_blocks = integer(), held_out_blocks = integer(),
    selected_samples = integer(), held_out_samples = integer(),
    stringsAsFactors = FALSE
  )
  if (!isTRUE(args$blocked_robustness_enabled)) {
    return(list(replicates = empty_replicates, strata = empty_strata))
  }
  worker <- function(replicate_index) {
    set.seed(args$seed + 50000000L + k * 100000L + replicate_index)
    strata_rows <- list()
    selected_blocks <- character(0)
    for (stratum in sort(unique(block_strata))) {
      available <- names(block_strata)[block_strata == stratum]
      if (length(available) == 1L) {
        retained <- available
      } else {
        retain_n <- max(1L, floor(args$stability_sample_fraction * length(available)))
        retain_n <- min(retain_n, length(available) - 1L)
        retained <- sort(sample(available, retain_n, replace = FALSE))
      }
      held_out_in_stratum <- setdiff(available, retained)
      selected_blocks <- c(selected_blocks, retained)
      strata_rows[[length(strata_rows) + 1]] <- data.frame(
        K = k, replicate = replicate_index, stability_stratum = stratum,
        available_blocks = length(available), selected_blocks = length(retained),
        held_out_blocks = length(held_out_in_stratum),
        selected_samples = sum(stability_blocks %in% retained),
        held_out_samples = sum(stability_blocks %in% held_out_in_stratum),
        stringsAsFactors = FALSE
      )
    }
    selected_blocks <- sort(unique(selected_blocks))
    selected <- which(stability_blocks %in% selected_blocks)
    held_out <- setdiff(seq_len(n_samples), selected)
    if (length(selected) <= k || length(held_out) <= k) return(NULL)
    fit <- tryCatch(
      cluster::pam(as.dist(distance_matrix[selected, selected, drop = FALSE]), k = k, diss = TRUE),
      error = function(error) NULL
    )
    if (is.null(fit)) return(NULL)
    replicate_row <- data.frame(
      K = k, replicate = replicate_index,
      validation_role = "secondary_robustness",
      resampling_unit = "whole_blocks_balanced_within_stratum",
      selected_blocks = length(selected_blocks),
      held_out_blocks = length(unique_blocks) - length(selected_blocks),
      selected_samples = length(selected), held_out_samples = length(held_out),
      adjusted_rand_index = adjusted_rand(full_labels[selected], fit$clustering),
      prediction_strength = prediction_strength(fit, selected, held_out, k),
      stringsAsFactors = FALSE
    )
    list(replicate = replicate_row, strata = do.call(rbind, strata_rows))
  }
  results <- if (args$ncpus > 1 && .Platform$OS.type != "windows") {
    parallel::mclapply(
      seq_len(args$stability_replicates), worker,
      mc.cores = args$ncpus, mc.preschedule = TRUE, mc.set.seed = FALSE
    )
  } else {
    lapply(seq_len(args$stability_replicates), worker)
  }
  results <- Filter(Negate(is.null), results)
  list(
    replicates = if (length(results)) do.call(rbind, lapply(results, `[[`, "replicate")) else empty_replicates,
    strata = if (length(results)) do.call(rbind, lapply(results, `[[`, "strata")) else empty_strata
  )
}

validation_rows <- list()
candidate_rows <- list()
medoid_rows <- list()
stability_replicate_rows <- list()
stability_cluster_rows <- list()
blocked_replicate_rows <- list()
blocked_strata_rows <- list()
pam_models <- list()
for (k in k_values) {
  fit <- tryCatch(cluster::pam(as.dist(distance_matrix), k = k, diss = TRUE), error = function(error) error)
  if (inherits(fit, "error")) {
    validation_rows[[length(validation_rows) + 1]] <- data.frame(
      cluster_method = "PAM", K = k, mean_silhouette = NA_real_, median_silhouette = NA_real_,
      cluster_sizes = "", minimum_cluster_size = NA_integer_, minimum_cluster_fraction = NA_real_,
      median_subsampling_ARI = NA_real_, subsampling_ARI_q25 = NA_real_, subsampling_ARI_q75 = NA_real_,
      primary_stability_quantile = NA_real_, minimum_cluster_median_jaccard = NA_real_,
      minimum_cluster_jaccard_q25 = NA_real_, median_prediction_strength = NA_real_,
      prediction_strength_q25 = NA_real_, median_season_balanced_block_ARI = NA_real_,
      season_balanced_block_ARI_q25 = NA_real_, negative_silhouette_fraction = NA_real_,
      agreement_with_PAM_ARI = NA_real_,
      successful_stability_replicates = 0L, feasible_cluster_sizes = FALSE,
      candidate_status = "failed", status_reason = conditionMessage(fit), stringsAsFactors = FALSE
    )
    next
  }
  canonical <- canonicalize_pam(fit, sample_ids)
  labels <- canonical$labels
  sil <- silhouette_values(labels, distance_matrix)
  sizes <- sort(table(labels))
  stability <- primary_stability_for_k(k, labels, canonical$medoids)
  blocked_stability <- season_balanced_block_stability_for_k(k, labels)
  replicate_stability <- stability$replicates
  cluster_stability <- stability$clusters
  if (nrow(replicate_stability)) stability_replicate_rows[[length(stability_replicate_rows) + 1]] <- replicate_stability
  if (nrow(cluster_stability)) stability_cluster_rows[[length(stability_cluster_rows) + 1]] <- cluster_stability
  if (nrow(blocked_stability$replicates)) blocked_replicate_rows[[length(blocked_replicate_rows) + 1]] <- blocked_stability$replicates
  if (nrow(blocked_stability$strata)) blocked_strata_rows[[length(blocked_strata_rows) + 1]] <- blocked_stability$strata
  successful <- replicate_stability$adjusted_rand_index[is.finite(replicate_stability$adjusted_rand_index)]
  blocked_ari <- blocked_stability$replicates$adjusted_rand_index[is.finite(blocked_stability$replicates$adjusted_rand_index)]
  blocked_prediction <- blocked_stability$replicates$prediction_strength[is.finite(blocked_stability$replicates$prediction_strength)]
  cluster_summary <- if (nrow(cluster_stability)) {
    do.call(rbind, lapply(split(cluster_stability, cluster_stability$microbial_compartment), function(values) {
      data.frame(
        microbial_compartment = values$microbial_compartment[1],
        median_jaccard = median(values$best_jaccard, na.rm = TRUE),
        jaccard_q25 = unname(quantile(values$best_jaccard, args$stability_primary_quantile, na.rm = TRUE)),
        stringsAsFactors = FALSE
      )
    }))
  } else data.frame()
  feasible <- min(sizes) >= minimum_cluster_n
  validation_rows[[length(validation_rows) + 1]] <- data.frame(
    cluster_method = "PAM", K = k, mean_silhouette = mean(sil, na.rm = TRUE),
    median_silhouette = median(sil, na.rm = TRUE),
    cluster_sizes = paste(paste0(names(sizes), ":", as.integer(sizes)), collapse = ";"),
    minimum_cluster_size = min(sizes), minimum_cluster_fraction = min(sizes) / n_samples,
    median_subsampling_ARI = if (length(successful)) median(successful) else NA_real_,
    subsampling_ARI_q25 = if (length(successful)) unname(quantile(successful, 0.25)) else NA_real_,
    subsampling_ARI_q75 = if (length(successful)) unname(quantile(successful, 0.75)) else NA_real_,
    primary_stability_quantile = if (length(successful)) unname(quantile(successful, args$stability_primary_quantile)) else NA_real_,
    minimum_cluster_median_jaccard = if (nrow(cluster_summary)) min(cluster_summary$median_jaccard) else NA_real_,
    minimum_cluster_jaccard_q25 = if (nrow(cluster_summary)) min(cluster_summary$jaccard_q25) else NA_real_,
    median_prediction_strength = if (length(blocked_prediction)) median(blocked_prediction) else NA_real_,
    prediction_strength_q25 = if (length(blocked_prediction)) unname(quantile(blocked_prediction, 0.25)) else NA_real_,
    median_season_balanced_block_ARI = if (length(blocked_ari)) median(blocked_ari) else NA_real_,
    season_balanced_block_ARI_q25 = if (length(blocked_ari)) unname(quantile(blocked_ari, args$stability_primary_quantile)) else NA_real_,
    negative_silhouette_fraction = mean(sil < 0, na.rm = TRUE),
    agreement_with_PAM_ARI = NA_real_,
    successful_stability_replicates = length(successful), feasible_cluster_sizes = feasible,
    candidate_status = "completed", status_reason = ifelse(feasible, "", "minimum_cluster_size_or_fraction_failed"),
    stringsAsFactors = FALSE
  )
  candidate_rows[[length(candidate_rows) + 1]] <- data.frame(
    sample_ID = sample_ids, cluster_method = "PAM", K = k,
    candidate_compartment = labels, silhouette_width_for_sample = sil,
    stringsAsFactors = FALSE
  )
  medoids <- canonical$medoids
  medoids$K <- k
  medoid_rows[[length(medoid_rows) + 1]] <- medoids[, c("K", "microbial_compartment", "medoid_sample_ID")]
  pam_models[[as.character(k)]] <- list(labels = labels, silhouette = sil)
}

if (isTRUE(args$hierarchical_enabled)) {
  hierarchy <- hclust(as.dist(distance_matrix), method = "average")
  for (k in k_values) {
    raw_labels <- cutree(hierarchy, k = k)
    labels <- canonicalize_hierarchical(raw_labels, sample_ids)
    sil <- silhouette_values(labels, distance_matrix)
    sizes <- sort(table(labels))
    pam_model <- pam_models[[as.character(k)]]
    pam_labels <- if (!is.null(pam_model)) pam_model$labels else NULL
    validation_rows[[length(validation_rows) + 1]] <- data.frame(
      cluster_method = "hierarchical_average", K = k,
      mean_silhouette = mean(sil, na.rm = TRUE), median_silhouette = median(sil, na.rm = TRUE),
      cluster_sizes = paste(paste0(names(sizes), ":", as.integer(sizes)), collapse = ";"),
      minimum_cluster_size = min(sizes), minimum_cluster_fraction = min(sizes) / n_samples,
      median_subsampling_ARI = NA_real_,
      subsampling_ARI_q25 = NA_real_, subsampling_ARI_q75 = NA_real_,
      primary_stability_quantile = NA_real_, minimum_cluster_median_jaccard = NA_real_,
      minimum_cluster_jaccard_q25 = NA_real_, median_prediction_strength = NA_real_,
      prediction_strength_q25 = NA_real_, median_season_balanced_block_ARI = NA_real_,
      season_balanced_block_ARI_q25 = NA_real_, negative_silhouette_fraction = mean(sil < 0, na.rm = TRUE),
      agreement_with_PAM_ARI = if (!is.null(pam_labels)) adjusted_rand(pam_labels, labels) else NA_real_,
      successful_stability_replicates = NA_integer_, feasible_cluster_sizes = min(sizes) >= minimum_cluster_n,
      candidate_status = "completed",
      status_reason = "hierarchical_sensitivity_only",
      stringsAsFactors = FALSE
    )
    candidate_rows[[length(candidate_rows) + 1]] <- data.frame(
      sample_ID = sample_ids, cluster_method = "hierarchical_average", K = k,
      candidate_compartment = labels, silhouette_width_for_sample = sil,
      stringsAsFactors = FALSE
    )
  }
}

validation <- do.call(rbind, validation_rows)
pam_validation <- validation[
  validation$cluster_method == "PAM" & validation$candidate_status == "completed" & validation$feasible_cluster_sizes,
  , drop = FALSE
]
rankable_pam <- pam_validation[is.finite(pam_validation$primary_stability_quantile), , drop = FALSE]
if (nrow(rankable_pam)) {
  best_score <- max(rankable_pam$primary_stability_quantile)
  best_rows <- rankable_pam[
    abs(rankable_pam$primary_stability_quantile - best_score) < .Machine$double.eps^0.5,
    , drop = FALSE
  ]
  best <- best_rows[which.min(best_rows$K), , drop = FALSE]
  best_k <- as.integer(best$K)
  supported_pool <- pam_validation[
    is.finite(pam_validation$primary_stability_quantile) &
      pam_validation$primary_stability_quantile >= args$min_stability_ari &
      is.finite(pam_validation$minimum_cluster_median_jaccard) &
      pam_validation$minimum_cluster_median_jaccard >= args$min_cluster_jaccard,
    , drop = FALSE
  ]
  if (nrow(supported_pool)) {
    selected_score <- max(supported_pool$primary_stability_quantile)
    selected_rows <- supported_pool[
      supported_pool$primary_stability_quantile >= selected_score - args$stability_near_tie_tolerance,
      , drop = FALSE
    ]
    selected <- selected_rows[which.min(selected_rows$K), , drop = FALSE]
    selected_k <- as.integer(selected$K)
    supported <- TRUE
    support_reasons <- character(0)
  } else {
    selected <- NULL
    selected_k <- NA_integer_
    supported <- FALSE
    support_reasons <- "no_candidate_met_within_block_stability_and_cluster_recovery_requirements"
    if (!is.finite(best$primary_stability_quantile) || best$primary_stability_quantile < args$min_stability_ari) {
      support_reasons <- c(support_reasons, "best_candidate_below_primary_stability_quantile_threshold")
    }
    if (!is.finite(best$minimum_cluster_median_jaccard) || best$minimum_cluster_median_jaccard < args$min_cluster_jaccard) {
      support_reasons <- c(support_reasons, "best_candidate_below_cluster_jaccard_threshold")
    }
  }
} else {
  best <- NULL
  best_k <- NA_integer_
  selected <- NULL
  selected_k <- NA_integer_
  support_reasons <- if (nrow(pam_validation)) {
    "no_feasible_candidate_completed_stability_validation"
  } else {
    "no_candidate_satisfied_cluster_size_constraints"
  }
  supported <- FALSE
}

validation$selected_candidate <- FALSE
if (is.finite(selected_k)) {
  validation$selected_candidate <- validation$cluster_method == "PAM" & validation$K == selected_k
}
validation$meets_primary_support_rule <- FALSE
primary_rows <- validation$cluster_method == "PAM" &
  validation$candidate_status == "completed" &
  !is.na(validation$feasible_cluster_sizes) & validation$feasible_cluster_sizes &
  is.finite(validation$primary_stability_quantile) &
  is.finite(validation$minimum_cluster_median_jaccard)
validation$meets_primary_support_rule[primary_rows] <-
  validation$primary_stability_quantile[primary_rows] >= args$min_stability_ari &
  validation$minimum_cluster_median_jaccard[primary_rows] >= args$min_cluster_jaccard
write_tsv(validation, file.path(args$outdir, "cluster_validation.tsv"), na = "")
stability_replicates <- if (length(stability_replicate_rows)) do.call(rbind, stability_replicate_rows) else data.frame()
stability_cluster_replicates <- if (length(stability_cluster_rows)) do.call(rbind, stability_cluster_rows) else data.frame()
write_tsv(stability_replicates, file.path(args$outdir, "stability_resampling_replicates.tsv"), na = "")
write_tsv(stability_cluster_replicates, file.path(args$outdir, "stability_cluster_recovery_replicates.tsv"), na = "")
blocked_replicates <- if (length(blocked_replicate_rows)) do.call(rbind, blocked_replicate_rows) else data.frame(
  K = integer(), replicate = integer(), validation_role = character(), resampling_unit = character(),
  selected_blocks = integer(), held_out_blocks = integer(), selected_samples = integer(),
  held_out_samples = integer(), adjusted_rand_index = double(), prediction_strength = double()
)
blocked_strata_audit <- if (length(blocked_strata_rows)) do.call(rbind, blocked_strata_rows) else data.frame(
  K = integer(), replicate = integer(), stability_stratum = character(), available_blocks = integer(),
  selected_blocks = integer(), held_out_blocks = integer(), selected_samples = integer(), held_out_samples = integer()
)
write_tsv(blocked_replicates, file.path(args$outdir, "season_balanced_block_replicates.tsv"), na = "")
write_tsv(blocked_strata_audit, file.path(args$outdir, "season_balanced_block_strata_audit.tsv"), na = "")
cluster_recovery_summary <- if (nrow(stability_cluster_replicates)) {
  do.call(rbind, lapply(
    split(stability_cluster_replicates, interaction(stability_cluster_replicates$K, stability_cluster_replicates$microbial_compartment, drop = TRUE)),
    function(values) data.frame(
      K = values$K[1],
      microbial_compartment = values$microbial_compartment[1],
      full_data_cluster_size = sum(
        pam_models[[as.character(values$K[1])]]$labels == values$microbial_compartment[1]
      ),
      medoid_sample_ID = values$medoid_sample_ID[1],
      median_best_jaccard = median(values$best_jaccard, na.rm = TRUE),
      best_jaccard_primary_quantile = unname(quantile(values$best_jaccard, args$stability_primary_quantile, na.rm = TRUE)),
      medoid_inclusion_replicates = sum(values$medoid_included, na.rm = TRUE),
      medoid_recovery_fraction_when_included = mean(values$medoid_recovered[values$medoid_included], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  ))
} else data.frame()
write_tsv(cluster_recovery_summary, file.path(args$outdir, "stability_cluster_recovery_summary.tsv"), na = "")
candidate_assignments <- do.call(rbind, candidate_rows)
write_tsv(candidate_assignments, file.path(args$outdir, "microbial_candidate_assignments.tsv"), na = "")
write_tsv(if (length(medoid_rows)) do.call(rbind, medoid_rows) else data.frame(
  K = integer(), microbial_compartment = character(), medoid_sample_ID = character()
), file.path(args$outdir, "pam_medoids.tsv"), na = "")

if (is.finite(best_k)) {
  best_model <- pam_models[[as.character(best_k)]]
  best_candidate <- data.frame(
    sample_ID = sample_ids, best_candidate_K = best_k,
    best_candidate_compartment = best_model$labels,
    best_candidate_silhouette_width = best_model$silhouette,
    support_status = ifelse(supported, "supported", "weak_unresolved"), stringsAsFactors = FALSE
  )
} else {
  best_candidate <- data.frame(
    sample_ID = sample_ids, best_candidate_K = NA_integer_,
    best_candidate_compartment = NA_character_, best_candidate_silhouette_width = NA_real_,
    support_status = "weak_unresolved", stringsAsFactors = FALSE
  )
}
write_tsv(best_candidate, file.path(args$outdir, "microbial_best_candidate_assignments.tsv"), na = "")

primary <- data.frame(
  sample_ID = sample_ids,
  microbial_compartment = if (supported) pam_models[[as.character(selected_k)]]$labels else "unresolved",
  cluster_method = "PAM",
  selected_K = if (supported) selected_k else NA_integer_,
  silhouette_width_for_sample = if (supported) pam_models[[as.character(selected_k)]]$silhouette else NA_real_,
  clustering_support_status = ifelse(supported, "supported", "weak_unresolved"),
  stringsAsFactors = FALSE
)
names(primary)[1] <- args$sample_col
write_tsv(primary, file.path(args$outdir, "microbial_compartments.tsv"), na = "")

if (supported) {
  summary_base <- do.call(rbind, lapply(sort(unique(primary$microbial_compartment)), function(group) {
    values <- primary$silhouette_width_for_sample[primary$microbial_compartment == group]
    data.frame(
      microbial_compartment = group, number_of_samples = length(values),
      mean_silhouette = mean(values), median_silhouette = median(values), stringsAsFactors = FALSE
    )
  }))
} else {
  summary_base <- data.frame(
    microbial_compartment = character(), number_of_samples = integer(),
    mean_silhouette = double(), median_silhouette = double()
  )
}
write_tsv(summary_base, file.path(args$outdir, "microbial_cluster_summary_base.tsv"), na = "")

decision <- data.frame(
  support_status = ifelse(supported, "supported", "weak_unresolved"),
  selected_K = if (supported) selected_k else NA_integer_,
  selected_mean_silhouette = if (!is.null(selected)) selected$mean_silhouette else NA_real_,
  selected_median_subsampling_ARI = if (!is.null(selected)) selected$median_subsampling_ARI else NA_real_,
  selected_primary_stability_quantile = if (!is.null(selected)) selected$primary_stability_quantile else NA_real_,
  selected_minimum_cluster_median_jaccard = if (!is.null(selected)) selected$minimum_cluster_median_jaccard else NA_real_,
  selected_median_season_balanced_block_ARI = if (!is.null(selected)) selected$median_season_balanced_block_ARI else NA_real_,
  selected_median_prediction_strength = if (!is.null(selected)) selected$median_prediction_strength else NA_real_,
  best_candidate_K = best_k,
  best_candidate_mean_silhouette = if (!is.null(best)) best$mean_silhouette else NA_real_,
  best_candidate_median_subsampling_ARI = if (!is.null(best)) best$median_subsampling_ARI else NA_real_,
  best_candidate_primary_stability_quantile = if (!is.null(best)) best$primary_stability_quantile else NA_real_,
  best_candidate_minimum_cluster_median_jaccard = if (!is.null(best)) best$minimum_cluster_median_jaccard else NA_real_,
  silhouette_selection_role = "descriptive_only",
  configured_legacy_mean_silhouette_threshold_not_applied = args$min_mean_silhouette,
  required_primary_stability_quantile = args$min_stability_ari,
  required_minimum_cluster_median_jaccard = args$min_cluster_jaccard,
  required_minimum_cluster_n = minimum_cluster_n,
  stability_near_tie_tolerance = args$stability_near_tie_tolerance,
  selection_rule = "among candidates meeting size, within-block lower-tail ARI, and per-state Jaccard requirements, choose the smallest K within the configured tolerance of the maximum lower-tail ARI",
  support_failure_reasons = paste(support_reasons, collapse = ";"),
  stringsAsFactors = FALSE
)
write_tsv(decision, file.path(args$outdir, "microbial_cluster_selection_decision.tsv"), na = "")

config <- list(
  k_min = args$k_min, k_max = args$k_max,
  min_cluster_size = args$min_cluster_size,
  min_cluster_fraction = args$min_cluster_fraction,
  resolved_minimum_cluster_n = minimum_cluster_n,
  min_mean_silhouette = args$min_mean_silhouette,
  silhouette_selection_role = "descriptive_only",
  min_stability_ari = args$min_stability_ari,
  min_cluster_jaccard = args$min_cluster_jaccard,
  stability_replicates = args$stability_replicates,
  stability_sample_fraction = args$stability_sample_fraction,
  primary_stability_resampling_unit = "samples_within_every_block",
  stability_primary_quantile = args$stability_primary_quantile,
  stability_near_tie_tolerance = args$stability_near_tie_tolerance,
  number_of_stability_blocks = length(unique_blocks),
  blocked_robustness_enabled = args$blocked_robustness_enabled,
  blocked_robustness_resampling_unit = "whole_blocks_balanced_within_stratum",
  number_of_blocked_robustness_strata = length(unique(block_strata)),
  prediction_strength_enabled = args$prediction_strength_enabled,
  selection_rule = "size feasibility plus within-block lower-tail ARI plus minimum per-cluster median Jaccard; smallest K within configured tolerance of maximum qualifying lower-tail ARI",
  hierarchical_enabled = args$hierarchical_enabled,
  hierarchical_linkage = "average",
  seed = args$seed, ncpus = args$ncpus,
  clustering_fit_inputs = c("Aitchison distances", "sample identifiers"),
  selection_validation_inputs = c("configured block identifiers for within-block sample perturbation only"),
  secondary_robustness_inputs = c("configured block and stratum identifiers for season-balanced whole-block withholding"),
  prohibited_clustering_fit_inputs = c("depth", "date", "environmental measurements", "environmental compartments", "stability block identifiers", "stability stratum identifiers"),
  R_version = R.version.string,
  cluster_version = as.character(packageVersion("cluster"))
)
write_json_object(config, file.path(args$outdir, "microbial_clustering_inference_config.json"))
capture.output(sessionInfo(), file = file.path(args$outdir, "microbial_clustering_session_info.txt"))
