#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
value_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[[1]]) else default
}
args_all <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
default_root <- normalizePath(file.path(dirname(script_file), "../.."), mustWork = TRUE)
code_root <- normalizePath(value_arg(
  "code-root", Sys.getenv("BAU_CODE_ROOT", default_root)
), mustWork = FALSE)
ukb_root <- normalizePath(value_arg(
  "ukb-root", Sys.getenv("UKB_DATA_ROOT", "/mnt/project/UKB")
), mustWork = FALSE)
results_root <- Sys.getenv("BAU_RESULTS_ROOT", "/home/dnanexus/bau-work/results")
source(file.path(code_root, "R/molecular.R"))

result_dir <- file.path(results_root, "molecular/cross_omic")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

analysis_specs <- list(
  nmr_to_protein = list(
    model = "nmr",
    target = "protein",
    expected_features = 2919L,
    loader = function() load_protein_targets(code_root, ukb_root, "nmr")
  ),
  protein_to_nmr = list(
    model = "protein",
    target = "nmr",
    expected_features = 168L,
    loader = function() load_nmr_targets(code_root, ukb_root, "protein")
  ),
  protein_to_noninput_protein = list(
    model = "protein",
    target = "non-input protein",
    expected_features = 2715L,
    loader = function() {
      load_protein_targets(
        code_root, ukb_root, "protein",
        exclude_inputs = TRUE
      )
    }
  )
)

association_tables <- list()
transform_tables <- list()
cohort_tables <- list()

for (analysis_id in names(analysis_specs)) {
  spec <- analysis_specs[[analysis_id]]
  message("Preparing ", analysis_id)
  dataset <- spec$loader()
  if (length(dataset$features) != spec$expected_features) {
    stop(paste("Unexpected target count for", analysis_id))
  }
  discovery_rows <- which(dataset$meta$OuterFold %in% 1:3)
  validation_rows <- which(dataset$meta$OuterFold %in% 4:5)
  discovery_centres <- unique(as.character(
    dataset$meta$AssessmentCenter[discovery_rows]
  ))
  validation_centres <- unique(as.character(
    dataset$meta$AssessmentCenter[validation_rows]
  ))
  if (length(intersect(discovery_centres, validation_centres))) {
    stop(paste("Assessment centres cross splits for", analysis_id))
  }

  transformed <- transform_molecular_targets(dataset, discovery_rows)
  design <- make_molecular_design(dataset$meta, discovery_rows, "linear")
  message("Fitting discovery centres for ", analysis_id)
  discovery <- fit_molecular_split(
    transformed$y, discovery_rows, design, dataset$meta$AssessmentCenter
  )
  message("Fitting validation centres for ", analysis_id)
  validation <- fit_molecular_split(
    transformed$y, validation_rows, design, dataset$meta$AssessmentCenter
  )
  association <- combine_molecular_splits(discovery, validation)
  association[, `:=`(
    analysis = analysis_id,
    age_model = spec$model,
    target_layer = spec$target
  )]
  setcolorder(association, c(
    "analysis", "age_model", "target_layer", "feature",
    setdiff(names(association), c("analysis", "age_model", "target_layer", "feature"))
  ))
  transformed$parameters[, `:=`(
    analysis = analysis_id,
    age_model = spec$model,
    target_layer = spec$target
  )]
  cohort <- data.table(
    analysis = analysis_id,
    age_model = spec$model,
    target_layer = spec$target,
    source_prediction_n = dataset$source_n,
    linked_n = nrow(dataset$meta),
    discovery_n = length(discovery_rows),
    validation_n = length(validation_rows),
    healthy_discovery_n = sum(
      dataset$meta$Baseline_unhealthy[discovery_rows] == 0L
    ),
    discovery_centres = length(discovery_centres),
    validation_centres = length(validation_centres),
    targets = length(dataset$features)
  )

  fwrite(
    association,
    file.path(result_dir, paste0(analysis_id, "_associations.tsv")),
    sep = "\t", na = "NA"
  )
  fwrite(
    transformed$parameters,
    file.path(result_dir, paste0(analysis_id, "_transformations.tsv")),
    sep = "\t", na = "NA"
  )
  association_tables[[analysis_id]] <- association
  transform_tables[[analysis_id]] <- transformed$parameters
  cohort_tables[[analysis_id]] <- cohort
  rm(dataset, transformed, design, discovery, validation, association)
  gc()
}

fwrite(
  rbindlist(association_tables, fill = TRUE),
  file.path(result_dir, "all_associations.tsv"),
  sep = "\t", na = "NA"
)
fwrite(
  rbindlist(transform_tables, fill = TRUE),
  file.path(result_dir, "all_transformations.tsv"),
  sep = "\t", na = "NA"
)
fwrite(
  rbindlist(cohort_tables),
  file.path(result_dir, "cohort_summary.tsv"),
  sep = "\t", na = "NA"
)
