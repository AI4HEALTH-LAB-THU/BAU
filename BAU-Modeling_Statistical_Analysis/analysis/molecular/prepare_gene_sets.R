#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(msigdbr)
})

args_all <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
default_root <- normalizePath(file.path(dirname(script_file), "../.."), mustWork = TRUE)
code_root <- Sys.getenv("BAU_CODE_ROOT", default_root)
data_root <- Sys.getenv("BAU_DATA_ROOT", "/home/dnanexus/bau-work/data")
reference_dir <- Sys.getenv("GENE_SET_ROOT", file.path(data_root, "references/msigdb"))
dir.create(reference_dir, recursive = TRUE, showWarnings = FALSE)
specifications <- list(
  hallmark = c("H", NA_character_),
  reactome = c("C2", "CP:REACTOME"),
  go_bp = c("C5", "GO:BP"),
  go_cc = c("C5", "GO:CC"),
  go_mf = c("C5", "GO:MF")
)
for (name in names(specifications)) {
  spec <- specifications[[name]]
  arguments <- list(species = "Homo sapiens", collection = spec[[1]])
  if (!is.na(spec[[2]])) arguments$subcollection <- spec[[2]]
  raw <- as.data.table(do.call(msigdbr, arguments))
  if (!identical(unique(raw$db_version), "2026.1.Hs")) {
    stop("MSigDB 2026.1.Hs is required for this analysis")
  }
  reference <- unique(raw[, .(
    db_version,
    collection = gs_collection,
    subcollection = gs_subcollection,
    gene_set = gs_name,
    description = gs_description,
    gene_symbol
  )])
  setorder(reference, gene_set, gene_symbol)
  fwrite(
    reference,
    file.path(reference_dir, paste0("MSigDB_2026.1_Hs_", name, ".tsv.gz")),
    sep = "\t", compress = "gzip"
  )
}
