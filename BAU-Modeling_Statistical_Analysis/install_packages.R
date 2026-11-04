#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
library_dir <- Sys.getenv("BAU_R_LIBRARY", "/home/dnanexus/bau-work/r-library")
dir.create(library_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(library_dir, .libPaths()))

versions <- c(
  arrow = "24.0.0", data.table = "1.18.4", digest = "0.6.39", dplyr = "1.2.1",
  ggplot2 = "4.0.3", haven = "2.5.5", igraph = "2.3.3", jsonlite = "2.0.0",
  ggrepel = "0.9.8", matrixStats = "1.5.0", patchwork = "1.3.2", readr = "2.2.0",
  sandwich = "3.1.1", survival = "3.8.6", svglite = "2.2.2", tidyr = "1.3.2",
  fgsea = "1.38.0"
)
installed_version <- function(package) {
  if (requireNamespace(package, quietly = TRUE)) as.character(packageVersion(package)) else NA_character_
}
actual <- vapply(names(versions), installed_version, character(1))
needed <- names(versions)[is.na(actual) | actual != versions]
if (length(needed)) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org", lib = library_dir)
  }
  for (package in needed) {
    if (package == "fgsea") {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org", lib = library_dir)
      }
      options(repos = BiocManager::repositories(version = "3.23"))
      remotes::install_url(
        "https://bioconductor.org/packages/3.23/bioc/src/contrib/fgsea_1.38.0.tar.gz",
        lib = library_dir, dependencies = NA, upgrade = "never"
      )
    } else {
      remotes::install_version(package,
        version = versions[[package]], lib = library_dir,
        repos = "https://cloud.r-project.org", upgrade = "never"
      )
    }
  }
}
if ("--gene-sets" %in% args && !requireNamespace("msigdbr", quietly = TRUE)) {
  install.packages("msigdbr", repos = "https://cloud.r-project.org", lib = library_dir)
}
