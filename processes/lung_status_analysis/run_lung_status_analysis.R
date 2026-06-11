#!/usr/bin/env Rscript
################################################################################
# Planned contrast analysis: TumorSide vs Contralateral vs Healthy
#
# Three planned contrasts:
#   A: TumorSide vs Contralateral (paired within cancer patients)
#   B: Contralateral vs Healthy (between-patient, patient-level)
#   C: TumorSide vs Healthy (between-patient, patient-level)
#
# Primary analysis: Lung Brush (n=8 paired cancer patients)
# Exploratory: BAL (n=4 paired cancer patients)
#
# FDR correction across the 3 contrasts per sample type.
################################################################################

library(vegan)
library(permute)  # For how() function
library(dplyr)
library(tidyr)
library(tibble)
library(readr)

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript run_lung_status_analysis.R <metadata.tsv> <asv_table.tsv> <outdir>")
}

metadata_file <- args[1]
asv_file <- args[2]
outdir <- args[3]

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("Loading data...\n")

if (!file.exists(metadata_file)) {
  stop(sprintf("Metadata file not found: %s", metadata_file))
}
if (!file.exists(asv_file)) {
  stop(sprintf("ASV table file not found: %s", asv_file))
}

metadata <- read_tsv(metadata_file, show_col_types = FALSE)
asv_table <- read_tsv(asv_file, show_col_types = FALSE) %>%
  column_to_rownames("sample")

# Align samples
shared_samples <- intersect(metadata$sample, rownames(asv_table))

if (length(shared_samples) == 0) {
  stop("No overlapping samples between metadata and ASV table")
}

metadata <- metadata %>% filter(sample %in% shared_samples)
asv_table <- asv_table[shared_samples, ]
asv_table <- as.matrix(asv_table)
storage.mode(asv_table) <- "numeric"

cat(sprintf("  %d samples in analysis\n", nrow(metadata)))

# Patient column (try multiple names)
patient_col <- intersect(c("Participant_ID", "patient_code"), colnames(metadata))[1]
if (is.na(patient_col)) stop("No patient ID column found")

# Convert to relative abundance
rs <- rowSums(asv_table)
rs[rs == 0] <- 1
asv_rel <- sweep(asv_table, 1, rs, "/")

# Bray-Curtis distance
bray_dist <- vegdist(asv_rel, method = "bray")

# Shannon diversity
shannon <- diversity(asv_table, index = "shannon")
metadata$shannon <- shannon[metadata$sample]

################################################################################
# CONTRAST A: TumorSide vs Contralateral (paired, cancer patients only)
################################################################################
cat("\n=== CONTRAST A: TumorSide vs Contralateral (paired) ===\n")

meta_paired <- metadata %>%
  filter(lung_status %in% c("TumorSide", "Contralateral"))

# Check pairing
paired_check <- meta_paired %>%
  group_by(!!sym(patient_col)) %>%
  summarise(
    n = n(),
    has_both = all(c("TumorSide", "Contralateral") %in% lung_status),
    .groups = "drop"
  ) %>%
  filter(has_both)

cat(sprintf("  %d cancer patients with both TumorSide and Contralateral samples\n", nrow(paired_check)))

if (nrow(paired_check) < 3) {
  cat("  WARNING: Too few paired patients for reliable analysis\n")
}

# Subset to only fully paired patients
paired_patients <- paired_check[[patient_col]]
meta_paired <- meta_paired %>%
  filter(!!sym(patient_col) %in% paired_patients)

# Build one profile per patient x lung_status (robust to duplicate/replicate samples)
paired_level <- meta_paired %>%
  group_by(!!sym(patient_col), lung_status) %>%
  summarise(
    samples = list(sample),
    shannon = mean(shannon, na.rm = TRUE),
    .groups = "drop"
  )

asv_paired <- matrix(0, nrow = nrow(paired_level), ncol = ncol(asv_rel))
rownames(asv_paired) <- paste(paired_level[[patient_col]], paired_level$lung_status, sep = "_")
colnames(asv_paired) <- colnames(asv_rel)

for (i in seq_len(nrow(paired_level))) {
  s <- paired_level$samples[[i]]
  if (length(s) == 1) {
    asv_paired[i, ] <- as.numeric(asv_rel[s, , drop = FALSE])
  } else {
    asv_paired[i, ] <- colMeans(asv_rel[s, , drop = FALSE], na.rm = TRUE)
  }
}

# Beta diversity (PERMANOVA with patient blocking)
dist_paired <- vegdist(asv_paired, method = "bray")

set.seed(12345)
perm_a <- vegan::adonis2(
  dist_paired ~ lung_status,
  data = paired_level,
  permutations = how(nperm = 9999, blocks = paired_level[[patient_col]]),
  method = "bray"
)

cat("\nPERMANOVA (patient-blocked):\n")
print(perm_a)

# PERMDISP check
disp_a <- betadisper(dist_paired, paired_level$lung_status)
perm_disp_a <- permutest(disp_a, pairwise = FALSE)
cat("\nPERMDISP:\n")
print(perm_disp_a)

# Alpha diversity (paired Wilcoxon signed-rank test)
alpha_wide_a <- paired_level %>%
  select(!!sym(patient_col), lung_status, shannon) %>%
  pivot_wider(names_from = lung_status, values_from = shannon)

wilcox_a <- wilcox.test(
  alpha_wide_a$TumorSide,
  alpha_wide_a$Contralateral,
  paired = TRUE,
  exact = FALSE
)

cat("\nAlpha diversity (paired Wilcoxon):\n")
cat(sprintf("  V = %.2f, p = %.4g\n", wilcox_a$statistic, wilcox_a$p.value))
cat(sprintf("  Median TumorSide: %.3f\n", median(alpha_wide_a$TumorSide, na.rm = TRUE)))
cat(sprintf("  Median Contralateral: %.3f\n", median(alpha_wide_a$Contralateral, na.rm = TRUE)))

# Save results
results_a <- tibble(
  contrast = "A_TumorSide_vs_Contralateral",
  comparison_type = "paired",
  n_patients = nrow(paired_check),
  n_samples = nrow(paired_level),
  permanova_R2 = perm_a$R2[1],
  permanova_F = perm_a$F[1],
  permanova_p = perm_a$`Pr(>F)`[1],
  permdisp_F = perm_disp_a$tab$F[1],
  permdisp_p = perm_disp_a$tab$`Pr(>F)`[1],
  alpha_statistic = wilcox_a$statistic,
  alpha_p = wilcox_a$p.value,
  alpha_median_group1 = median(alpha_wide_a$TumorSide, na.rm = TRUE),
  alpha_median_group2 = median(alpha_wide_a$Contralateral, na.rm = TRUE)
)

# Save pairwise distances for visualization
dist_mat <- as.matrix(dist_paired)
pairdist_a <- tibble()
for (pat in paired_patients) {
  s <- paired_level %>% filter(!!sym(patient_col) == pat)
  if (nrow(s) == 2 && all(c("TumorSide", "Contralateral") %in% s$lung_status)) {
    key_t <- paste(pat, "TumorSide", sep = "_")
    key_c <- paste(pat, "Contralateral", sep = "_")
    pairdist_a <- bind_rows(
      pairdist_a,
      tibble(
        patient = pat,
        sample1 = key_t,
        sample2 = key_c,
        bray = dist_mat[key_t, key_c],
        comparison = "TumorSide_vs_Contralateral"
      )
    )
  }
}

write_tsv(pairdist_a, file.path(outdir, "contrast_A_pairwise_distances.tsv"))

################################################################################
# CONTRAST B & C: Patient-level comparisons
# Average control samples (left/right) to create one profile per patient
################################################################################
cat("\n=== Preparing patient-level data for contrasts B & C ===\n")

# Function to average ASV profiles for a patient
average_patient_profile <- function(samples, asv_mat) {
  if (length(samples) == 1) {
    return(as.numeric(asv_mat[samples, , drop = FALSE]))
  } else {
    return(colMeans(asv_mat[samples, , drop = FALSE], na.rm = TRUE))
  }
}

# Create patient-level ASV table
patient_level <- metadata %>%
  group_by(!!sym(patient_col), lung_status) %>%
  summarise(
    n_samples = n(),
    samples = list(sample),
    shannon_mean = mean(shannon, na.rm = TRUE),
    .groups = "drop"
  )

cat(sprintf("  %d patient-level profiles\n", nrow(patient_level)))

# Build patient-level ASV matrix
asv_patient <- matrix(0, nrow = nrow(patient_level), ncol = ncol(asv_rel))
rownames(asv_patient) <- paste(patient_level[[patient_col]], patient_level$lung_status, sep = "_")
colnames(asv_patient) <- colnames(asv_rel)

for (i in seq_len(nrow(patient_level))) {
  samples <- patient_level$samples[[i]]
  asv_patient[i, ] <- average_patient_profile(samples, asv_rel)
}

# Recalculate distances on patient-level profiles
bray_patient <- vegdist(asv_patient, method = "bray")

################################################################################
# CONTRAST B: Contralateral vs Healthy
################################################################################
cat("\n=== CONTRAST B: Contralateral vs Healthy ===\n")

meta_b <- patient_level %>%
  filter(lung_status %in% c("Contralateral", "Healthy"))

cat(sprintf("  %d patients (Contralateral: %d, Healthy: %d)\n",
            nrow(meta_b),
            sum(meta_b$lung_status == "Contralateral"),
            sum(meta_b$lung_status == "Healthy")))

dist_b <- as.dist(as.matrix(bray_patient)[
  paste(meta_b[[patient_col]], meta_b$lung_status, sep = "_"),
  paste(meta_b[[patient_col]], meta_b$lung_status, sep = "_")
])

set.seed(12345)
perm_b <- adonis2(
  dist_b ~ lung_status,
  data = meta_b,
  permutations = 9999,
  method = "bray"
)

cat("\nPERMANOVA:\n")
print(perm_b)

disp_b <- betadisper(dist_b, meta_b$lung_status)
perm_disp_b <- permutest(disp_b, pairwise = FALSE)
cat("\nPERMDISP:\n")
print(perm_disp_b)

# Alpha diversity (Wilcoxon rank-sum test)
alpha_b <- meta_b %>% select(lung_status, shannon_mean)
wilcox_b <- wilcox.test(
  shannon_mean ~ lung_status,
  data = alpha_b,
  exact = FALSE
)

cat("\nAlpha diversity (Wilcoxon rank-sum):\n")
cat(sprintf("  W = %.2f, p = %.4g\n", wilcox_b$statistic, wilcox_b$p.value))
cat(sprintf("  Median Contralateral: %.3f\n", median(alpha_b$shannon_mean[alpha_b$lung_status == "Contralateral"], na.rm = TRUE)))
cat(sprintf("  Median Healthy: %.3f\n", median(alpha_b$shannon_mean[alpha_b$lung_status == "Healthy"], na.rm = TRUE)))

results_b <- tibble(
  contrast = "B_Contralateral_vs_Healthy",
  comparison_type = "between_patient",
  n_patients = nrow(meta_b),
  n_samples = nrow(meta_b),
  permanova_R2 = perm_b$R2[1],
  permanova_F = perm_b$F[1],
  permanova_p = perm_b$`Pr(>F)`[1],
  permdisp_F = perm_disp_b$tab$F[1],
  permdisp_p = perm_disp_b$tab$`Pr(>F)`[1],
  alpha_statistic = wilcox_b$statistic,
  alpha_p = wilcox_b$p.value,
  alpha_median_group1 = median(alpha_b$shannon_mean[alpha_b$lung_status == "Contralateral"], na.rm = TRUE),
  alpha_median_group2 = median(alpha_b$shannon_mean[alpha_b$lung_status == "Healthy"], na.rm = TRUE)
)

################################################################################
# CONTRAST C: TumorSide vs Healthy
################################################################################
cat("\n=== CONTRAST C: TumorSide vs Healthy ===\n")

meta_c <- patient_level %>%
  filter(lung_status %in% c("TumorSide", "Healthy"))

cat(sprintf("  %d patients (TumorSide: %d, Healthy: %d)\n",
            nrow(meta_c),
            sum(meta_c$lung_status == "TumorSide"),
            sum(meta_c$lung_status == "Healthy")))

dist_c <- as.dist(as.matrix(bray_patient)[
  paste(meta_c[[patient_col]], meta_c$lung_status, sep = "_"),
  paste(meta_c[[patient_col]], meta_c$lung_status, sep = "_")
])

set.seed(12345)
perm_c <- adonis2(
  dist_c ~ lung_status,
  data = meta_c,
  permutations = 9999,
  method = "bray"
)

cat("\nPERMANOVA:\n")
print(perm_c)

disp_c <- betadisper(dist_c, meta_c$lung_status)
perm_disp_c <- permutest(disp_c, pairwise = FALSE)
cat("\nPERMDISP:\n")
print(perm_disp_c)

# Alpha diversity
alpha_c <- meta_c %>% select(lung_status, shannon_mean)
wilcox_c <- wilcox.test(
  shannon_mean ~ lung_status,
  data = alpha_c,
  exact = FALSE
)

cat("\nAlpha diversity (Wilcoxon rank-sum):\n")
cat(sprintf("  W = %.2f, p = %.4g\n", wilcox_c$statistic, wilcox_c$p.value))
cat(sprintf("  Median TumorSide: %.3f\n", median(alpha_c$shannon_mean[alpha_c$lung_status == "TumorSide"], na.rm = TRUE)))
cat(sprintf("  Median Healthy: %.3f\n", median(alpha_c$shannon_mean[alpha_c$lung_status == "Healthy"], na.rm = TRUE)))

results_c <- tibble(
  contrast = "C_TumorSide_vs_Healthy",
  comparison_type = "between_patient",
  n_patients = nrow(meta_c),
  n_samples = nrow(meta_c),
  permanova_R2 = perm_c$R2[1],
  permanova_F = perm_c$F[1],
  permanova_p = perm_c$`Pr(>F)`[1],
  permdisp_F = perm_disp_c$tab$F[1],
  permdisp_p = perm_disp_c$tab$`Pr(>F)`[1],
  alpha_statistic = wilcox_c$statistic,
  alpha_p = wilcox_c$p.value,
  alpha_median_group1 = median(alpha_c$shannon_mean[alpha_c$lung_status == "TumorSide"], na.rm = TRUE),
  alpha_median_group2 = median(alpha_c$shannon_mean[alpha_c$lung_status == "Healthy"], na.rm = TRUE)
)

################################################################################
# Combine results and apply FDR correction
################################################################################
all_results <- bind_rows(results_a, results_b, results_c)

# FDR correction across the 3 planned contrasts
all_results <- all_results %>%
  mutate(
    permanova_q = p.adjust(permanova_p, method = "fdr"),
    alpha_q = p.adjust(alpha_p, method = "fdr")
  )

cat("\n=== Summary of all contrasts (FDR-corrected) ===\n")
print(all_results)

write_tsv(all_results, file.path(outdir, "lung_status_contrasts_summary.tsv"))

# Save patient-level metadata for visualization
write_tsv(patient_level, file.path(outdir, "patient_level_metadata.tsv"))

# Save patient-level distance matrix
write.table(
  as.matrix(bray_patient),
  file.path(outdir, "patient_level_bray_distances.tsv"),
  sep = "\t",
  quote = FALSE
)

cat(sprintf("\n=== Analysis complete. Results saved to %s ===\n", outdir))
