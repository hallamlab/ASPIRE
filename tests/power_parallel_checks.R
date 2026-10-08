# Standalone scientific parity check; optional second argument is a serial baseline.
argv <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(argv[1])
suppressPackageStartupMessages({library(indicspecies); library(permute)})
load_functions <- function(path, env) {
  for (expr in parse(path)) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        is.call(expr[[3]]) && identical(expr[[3]][[1]], as.name("function"))) eval(expr, env)
  }
}
script_path <- file.path(root, "processes/power_analysis_pipeline/power_indicspecies.R")
source(file.path(dirname(script_path), "power_parallel.R"))
load_functions(script_path, environment())
checkpoint_dir <- tempfile("isa-parallel-")
dir.create(checkpoint_dir)
args <- list(workers=1L, `checkpoint-dir`=checkpoint_dir, outdir=checkpoint_dir)
set.seed(78)
counts <- matrix(rpois(24*6, 20), 24, 6)
patients <- rep(paste0("P", 1:12), each=2)
case <- factor(rep(rep(c("Cancer", "Control"), each=6), each=2))
types <- factor(rep(c("A", "B"), 12))
counts[case == "Cancer", 1] <- counts[case == "Cancer", 1] * 5
for (blocking in c(FALSE, TRUE)) for (null in c(FALSE, TRUE)) {
  input <- list(counts=counts, patient_ids=patients, grouping=if (blocking) types else case,
                spike_asv_idx=1L, spike_fc=2, n_size=4, n_simulations=8L,
                n_perm=19L, seed=14, use_blocking=blocking, use_true_null=null, alpha=.2)
  args$workers <- 1L
  serial <- do.call(run_isa_power_generic, input)
  if (length(argv) > 1) {
    old <- new.env(parent=globalenv())
    load_functions(argv[2], old)
    reference <- do.call(old$run_isa_power_generic, input)
    stopifnot(identical(serial, reference))
  }
  files <- list.files(checkpoint_dir, full.names=TRUE, pattern="\\.rds$")
  # Use an empty directory to exercise actual workers, not cache reuse.
  unlink(files)
  args$workers <- 2L
  parallel <- do.call(run_isa_power_generic, input)
  stopifnot(identical(serial, parallel))
  checkpoint <- list.files(checkpoint_dir, full.names=TRUE, pattern="\\.rds$")
  stopifnot(length(checkpoint) == 1)
  saveRDS(readRDS(checkpoint)[1:3], checkpoint)
  args$workers <- 3L
  resumed <- do.call(run_isa_power_generic, input)
  stopifnot(identical(serial, resumed), length(readRDS(checkpoint)) == 8)
  unlink(checkpoint)
}
unlink(checkpoint_dir, recursive=TRUE)
cat("ISA serial / parallel / resumed results match across blocked and null designs.\n")
