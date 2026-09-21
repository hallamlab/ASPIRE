#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

options_list <- list(
  make_option("--library", type = "character"),
  make_option("--repository", type = "character", default = "dkahle/TITAN2"),
  make_option("--revision", type = "character"),
  make_option("--version", type = "character", default = "2.4.4"),
  make_option("--manifest", type = "character", default = "titan2_installation.tsv")
)
args <- parse_args(OptionParser(option_list = options_list))

if (is.null(args$library) || !nzchar(args$library)) {
  stop("--library is required", call. = FALSE)
}
if (is.null(args$revision) || !grepl("^[0-9a-f]{40}$", args$revision)) {
  stop("--revision must be a full 40-character Git commit SHA", call. = FALSE)
}

dir.create(args$library, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(normalizePath(args$library), .libPaths()))

remote <- paste0(args$repository, "@", args$revision)
remotes::install_github(
  remote,
  lib = args$library,
  dependencies = FALSE,
  upgrade = "never",
  build_vignettes = FALSE,
  quiet = FALSE
)

description <- utils::packageDescription("TITAN2", lib.loc = args$library)
installed_version <- as.character(description$Version)
installed_sha <- if (is.null(description$RemoteSha)) "" else as.character(description$RemoteSha)
if (!identical(installed_version, args$version)) {
  stop(
    sprintf("Installed TITAN2 version %s; expected %s", installed_version, args$version),
    call. = FALSE
  )
}
if (!identical(tolower(installed_sha), tolower(args$revision))) {
  stop(
    sprintf("Installed TITAN2 SHA %s; expected %s", installed_sha, args$revision),
    call. = FALSE
  )
}

manifest <- data.frame(
  package = "TITAN2",
  version = installed_version,
  repository = paste0("https://github.com/", args$repository),
  revision = installed_sha,
  installation_method = "remotes::install_github",
  stringsAsFactors = FALSE
)
utils::write.table(
  manifest,
  file = args$manifest,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
