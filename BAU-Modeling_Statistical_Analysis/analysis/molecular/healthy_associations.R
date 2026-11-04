#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
default_root <- normalizePath(file.path(dirname(script_file), "../.."), mustWork = TRUE)
code_root <- normalizePath(
  Sys.getenv("BAU_CODE_ROOT", default_root),
  mustWork = FALSE
)
ukb_root <- normalizePath(Sys.getenv("UKB_DATA_ROOT", "/mnt/project/UKB"), mustWork = FALSE)
results_root <- Sys.getenv("BAU_RESULTS_ROOT", "/home/dnanexus/bau-work/results")
source(file.path(code_root, "R/molecular.R"))
result_dir <- file.path(results_root, "molecular/healthy")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
all_cohorts <- list()

for (modality in c("nmr", "protein")) {
  analysis_id <- if (modality == "nmr") "nmr_to_protein" else "protein_to_noninput_protein"
  dataset <- load_protein_targets(
    code_root, ukb_root, modality,
    exclude_inputs = modality == "protein"
  )
  discovery_rows <- which(dataset$meta$OuterFold %in% 1:3)
  validation_rows <- which(dataset$meta$OuterFold %in% 4:5)
  transformed <- transform_molecular_targets(dataset, discovery_rows)
  design <- make_molecular_design(dataset$meta, discovery_rows)
  healthy_discovery <- discovery_rows[
    dataset$meta$Baseline_unhealthy[discovery_rows] == 0L
  ]
  healthy_validation <- validation_rows[
    dataset$meta$Baseline_unhealthy[validation_rows] == 0L
  ]
  discovery <- fit_molecular_split(
    transformed$y, healthy_discovery, design, dataset$meta$AssessmentCenter
  )
  validation <- fit_molecular_split(
    transformed$y, healthy_validation, design, dataset$meta$AssessmentCenter
  )
  association <- combine_molecular_splits(discovery, validation)
  association[, `:=`(analysis = analysis_id, age_model = modality)]
  fwrite(
    association,
    file.path(result_dir, paste0(analysis_id, "_associations.tsv")),
    sep = "\t", na = "NA"
  )
  all_cohorts[[modality]] <- data.table(
    age_model = modality,
    analysis = analysis_id,
    discovery_n = length(healthy_discovery),
    validation_n = length(healthy_validation),
    discovery_centres = uniqueN(dataset$meta$AssessmentCenter[healthy_discovery]),
    validation_centres = uniqueN(dataset$meta$AssessmentCenter[healthy_validation]),
    targets = length(dataset$features)
  )
  rm(dataset, transformed, design, discovery, validation, association)
  gc()
}
fwrite(rbindlist(all_cohorts), file.path(result_dir, "cohort_summary.tsv"), sep = "\t")
