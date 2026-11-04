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
result_dir <- file.path(results_root, "molecular/sensitivity")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

specifications <- c(
  "absolute_bag", "bag_spline", "central_95", "baseline_healthy", "low_missingness"
)
all_coefficients <- list()
all_summaries <- list()
all_cohorts <- list()

for (analysis_id in c("nmr_to_protein", "protein_to_nmr")) {
  dataset <- if (analysis_id == "nmr_to_protein") {
    load_protein_targets(code_root, ukb_root, "nmr")
  } else {
    load_nmr_targets(code_root, ukb_root, "protein")
  }
  discovery_rows <- which(dataset$meta$OuterFold %in% 1:3)
  transformed <- transform_molecular_targets(dataset, discovery_rows)
  designs <- list(
    linear = make_molecular_design(dataset$meta, discovery_rows, "linear"),
    absolute = make_molecular_design(dataset$meta, discovery_rows, "absolute"),
    spline = make_molecular_design(dataset$meta, discovery_rows, "spline")
  )
  bag_limits <- quantile(dataset$meta$BAG[discovery_rows], c(0.025, 0.975), names = FALSE)
  rows <- list(
    absolute_bag = discovery_rows,
    bag_spline = discovery_rows,
    central_95 = discovery_rows[
      dataset$meta$BAG[discovery_rows] >= bag_limits[[1]] &
        dataset$meta$BAG[discovery_rows] <= bag_limits[[2]]
    ],
    baseline_healthy = discovery_rows[
      dataset$meta$Baseline_unhealthy[discovery_rows] == 0L
    ],
    low_missingness = discovery_rows[
      dataset$meta$target_missing_fraction[discovery_rows] <= 0.01
    ]
  )
  primary <- fread(file.path(
    results_root, "molecular/cross_omic",
    paste0(analysis_id, "_associations.tsv")
  ))[, .(
    feature,
    primary_beta = beta_BAU_discovery,
    primary_replicated = replicated_BAU
  )]
  recalculated <- fit_single_term_split(
    transformed$y, discovery_rows, designs$linear,
    dataset$meta$AssessmentCenter, "BAU"
  )
  check <- merge(primary, recalculated, by = "feature")
  if (max(abs(check$primary_beta - check$beta)) > 1e-8) {
    stop(paste("Primary coefficients differ for", analysis_id))
  }

  for (specification in specifications) {
    message("Fitting ", analysis_id, " ", specification)
    design <- switch(specification,
      absolute_bag = designs$absolute,
      bag_spline = designs$spline,
      designs$linear
    )
    estimate <- fit_single_term_split(
      transformed$y, rows[[specification]], design,
      dataset$meta$AssessmentCenter, "BAU"
    )
    estimate <- merge(estimate, primary, by = "feature", sort = FALSE)
    estimate[, `:=`(analysis = analysis_id, specification = specification)]
    key <- paste(analysis_id, specification, sep = "_")
    replicated <- estimate[primary_replicated == TRUE]
    all_coefficients[[key]] <- estimate
    all_summaries[[key]] <- data.table(
      analysis = analysis_id,
      specification = specification,
      molecules = nrow(estimate),
      pearson_r = cor(estimate$primary_beta, estimate$beta),
      zero_intercept_slope = sum(estimate$primary_beta * estimate$beta) /
        sum(estimate$primary_beta^2),
      primary_replicated_n = nrow(replicated),
      same_direction_n = sum(sign(replicated$primary_beta) == sign(replicated$beta)),
      direction_agreement = mean(sign(replicated$primary_beta) == sign(replicated$beta))
    )
    all_cohorts[[key]] <- data.table(
      analysis = analysis_id,
      specification = specification,
      nominal_n = length(rows[[specification]]),
      minimum_model_n = min(estimate$n),
      maximum_model_n = max(estimate$n),
      centres = uniqueN(dataset$meta$AssessmentCenter[rows[[specification]]]),
      bag_lower = bag_limits[[1]],
      bag_upper = bag_limits[[2]]
    )
  }
  rm(dataset, transformed, designs)
  gc()
}

fwrite(rbindlist(all_coefficients), file.path(result_dir, "bau_coefficients.tsv"), sep = "\t")
fwrite(rbindlist(all_summaries), file.path(result_dir, "coefficient_agreement.tsv"), sep = "\t")
fwrite(rbindlist(all_cohorts), file.path(result_dir, "cohort_summary.tsv"), sep = "\t")
