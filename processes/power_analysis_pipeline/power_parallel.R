# Indexed seeds are set inside each replicate by the existing bootstrap code.
# Parent-only atomic checkpoints and ordered reduction preserve resume behavior.
power_replicates <- function(fun, count, inputs, workers, checkpoint_dir, outdir, script_path) {
  if (count < 1) stop("simulation count must be positive")
  if (workers > 1 && .Platform$OS.type == "windows") stop("Parallel power requires POSIX; use workers=1")
  if (is.null(checkpoint_dir)) checkpoint_dir <- file.path(outdir, "checkpoints")
  dir.create(checkpoint_dir, recursive=TRUE, showWarnings=FALSE)
  sources <- sort(list.files(dirname(script_path), pattern="\\.R$", full.names=TRUE))
  versions <- vapply(c("indicspecies", "permute", "dplyr"), function(p) as.character(packageVersion(p)), "")
  keyfile <- tempfile(tmpdir=checkpoint_dir)
  on.exit(unlink(keyfile), add=TRUE)
  saveRDS(list(inputs=inputs, source=tools::md5sum(sources), versions=versions,
               R=R.version.string), keyfile, version=2)
  key <- unname(tools::md5sum(keyfile))
  checkpoint <- file.path(checkpoint_dir, paste0(key, ".rds"))
  completed <- if (file.exists(checkpoint)) readRDS(checkpoint) else list()
  if (!is.list(completed) || length(completed) > count ||
      any(!vapply(completed, function(x) is.numeric(x) && length(x) == 5L, TRUE))) {
    stop("Invalid ISA checkpoint: ", checkpoint)
  }
  cat(sprintf("    ISA: %d/%d cached; workers=%d\n", length(completed), count, workers))
  start <- proc.time()[[3]]
  while (length(completed) < count) {
    indices <- seq.int(length(completed)+1L, min(count, length(completed)+max(25L, workers)))
    results <- if (workers == 1L) lapply(indices, fun) else
      parallel::mclapply(indices, fun, mc.cores=min(workers, length(indices)),
                        mc.set.seed=FALSE, mc.preschedule=TRUE)
    if (any(vapply(results, inherits, TRUE, what="try-error")) ||
        any(!vapply(results, function(x) is.numeric(x) && length(x) == 5L, TRUE))) {
      stop("ISA worker failed; preserving the last completed checkpoint")
    }
    completed <- c(completed, results)
    temporary <- paste0(checkpoint, ".", Sys.getpid(), ".tmp")
    saveRDS(completed, temporary)
    if (!file.rename(temporary, checkpoint)) stop("Could not save ISA checkpoint")
    cat(sprintf("    ISA: %d/%d; elapsed %.1fs\n", length(completed), count, proc.time()[[3]]-start))
    flush.console()
  }
  completed
}
