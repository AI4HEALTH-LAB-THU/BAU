#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
code_arg <- grep("^--code-root=", args, value = TRUE)
args_all <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
default_root <- normalizePath(file.path(dirname(script_file), "../.."), mustWork = TRUE)
code_root <- if (length(code_arg)) {
  sub("^--code-root=", "", code_arg[[1]])
} else {
  Sys.getenv("BAU_CODE_ROOT", default_root)
}
code_root <- normalizePath(code_root, mustWork = FALSE)
data_root <- Sys.getenv("BAU_DATA_ROOT", "/home/dnanexus/bau-work/data")
results_root <- Sys.getenv("BAU_RESULTS_ROOT", "/home/dnanexus/bau-work/results")
result_dir <- file.path(results_root, "phenotypes/model_diagnostics")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

config <- fromJSON(file.path(
  code_root, "config/ukb.json"
))
modality_spec <- data.table(
  modality = c("nmr", "protein"),
  label = c("NMR", "Proteomic"),
  expected_n = c(
    config$modalities$nmr$expected_rows,
    config$modalities$protein$expected_rows
  )
)

read_oof <- function(modality, expected_n) {
  path <- file.path(
    results_root, "models", modality, "ngboost",
    "all_participant_predictions.csv"
  )
  needed <- c(
    "eid_int", "Age", "Sex_male", "Baseline_unhealthy",
    "AssessmentCenter", "OuterFold", "SplitRole", "BA_Calibrated",
    "Var_Calibrated", "BAG", "BAU", "Z_AgeGap_Calibrated",
    "Z_Uncertainty_Calibrated"
  )
  x <- fread(path, select = needed, showProgress = FALSE)
  if (nrow(x) != expected_n || uniqueN(x$eid_int) != expected_n) {
    stop(paste("Unexpected participant count for", modality))
  }
  if (!all(x$SplitRole == "outer_test") || !setequal(x$OuterFold, 1:5)) {
    stop(paste("Outer-fold assignment is incomplete for", modality))
  }
  if (!isTRUE(all.equal(x$BAG, x$Z_AgeGap_Calibrated)) ||
    !isTRUE(all.equal(x$BAU, x$Z_Uncertainty_Calibrated))) {
    stop(paste("Phenotype aliases differ for", modality))
  }
  x[, modality := modality]
  x
}

predictions <- rbindlist(lapply(seq_len(nrow(modality_spec)), function(i) {
  read_oof(modality_spec$modality[[i]], modality_spec$expected_n[[i]])
}))

population_summary <- predictions[, .(
  n = .N,
  healthy_n = sum(Baseline_unhealthy == 0L),
  female_n = sum(Sex_male == 0L),
  male_n = sum(Sex_male == 1L),
  age_min = min(Age),
  age_q1 = quantile(Age, 0.25),
  age_median = median(Age),
  age_q3 = quantile(Age, 0.75),
  age_max = max(Age),
  centres = uniqueN(AssessmentCenter)
), by = modality]

all_population <- copy(predictions)
all_population[, population := "all eligible"]
healthy_population <- copy(predictions[Baseline_unhealthy == 0L])
healthy_population[, population := "healthy reference"]
age_sex_histogram <- rbindlist(list(
  all_population,
  healthy_population
))[, .(n = .N), by = .(
  modality, population,
  age_year = floor(Age),
  sex = fifelse(Sex_male == 1L, "male", "female")
)][, proportion := n / sum(n), by = .(modality, population)]

fold_roles <- rbindlist(lapply(modality_spec$modality, function(modality) {
  candidate_path <- file.path(
    data_root, sub("^data/", "", config$modalities[[modality]]$candidate)
  )
  candidate <- as.data.table(read_feather(
    candidate_path,
    col_select = c("eid", "Age", "Sex_male", "AssessmentCenter")
  ))
  candidate[, eid_int := as.integer(eid)]
  rbindlist(lapply(1:5, function(fold) {
    role_path <- file.path(
      results_root, "models", modality, "ngboost",
      paste0("outer_fold_", fold), "development_roles.csv"
    )
    roles <- fread(role_path)
    roles[, eid_int := as.integer(eid_int)]
    roles <- merge(
      roles, candidate[, .(eid_int, Age, Sex_male, AssessmentCenter)],
      by = "eid_int", all.x = TRUE, sort = FALSE
    )
    roles[, .(
      n = .N,
      female_n = sum(Sex_male == 0L),
      age_q1 = quantile(Age, 0.25),
      age_median = median(Age),
      age_q3 = quantile(Age, 0.75),
      centres = uniqueN(AssessmentCenter)
    ), by = .(role)][, `:=`(modality = modality, outer_fold = fold)]
  }))
}))

held_roles <- predictions[, .(
  n = .N,
  female_n = sum(Sex_male == 0L),
  age_q1 = quantile(Age, 0.25),
  age_median = median(Age),
  age_q3 = quantile(Age, 0.75),
  centres = uniqueN(AssessmentCenter)
), by = .(modality, outer_fold = OuterFold)][, role := "held_centre_test"]
fold_roles <- rbindlist(list(fold_roles, held_roles), fill = TRUE)

performance <- predictions[,
  {
    error <- BA_Calibrated - Age
    fit <- lm(BA_Calibrated ~ Age)
    .(
      n = .N,
      mae_years = mean(abs(error)),
      rmse_years = sqrt(mean(error^2)),
      pearson_r = cor(Age, BA_Calibrated),
      r_squared = 1 - sum(error^2) / sum((Age - mean(Age))^2),
      calibration_intercept = coef(fit)[[1]],
      calibration_slope = coef(fit)[[2]]
    )
  },
  by = .(modality, outer_fold = fifelse(OuterFold > 0L, OuterFold, NA_integer_))
]
pooled_performance <- predictions[,
  {
    error <- BA_Calibrated - Age
    fit <- lm(BA_Calibrated ~ Age)
    .(
      n = .N,
      mae_years = mean(abs(error)),
      rmse_years = sqrt(mean(error^2)),
      pearson_r = cor(Age, BA_Calibrated),
      r_squared = 1 - sum(error^2) / sum((Age - mean(Age))^2),
      calibration_intercept = coef(fit)[[1]],
      calibration_slope = coef(fit)[[2]]
    )
  },
  by = modality
][, outer_fold := NA_integer_]
performance <- rbindlist(list(performance, pooled_performance), use.names = TRUE)

prediction_density <- predictions[, .(
  n = .N
), by = .(
  modality,
  chronological_age = floor(Age * 2) / 2 + 0.25,
  predicted_age = floor(BA_Calibrated * 2) / 2 + 0.25
)]

relationship <- predictions[,
  {
    r <- cor(BAG, BAU)
    fit <- lm(BAU ~ BAG)
    z <- atanh(r)
    se <- 1 / sqrt(.N - 3)
    interval <- tanh(z + c(-1, 1) * qnorm(0.975) * se)
    .(
      n = .N,
      pearson_r = r,
      pearson_ci_lower = interval[[1]],
      pearson_ci_upper = interval[[2]],
      spearman_rho = cor(BAG, BAU, method = "spearman"),
      linear_r_squared = summary(fit)$r.squared,
      intercept = coef(fit)[[1]],
      slope = coef(fit)[[2]]
    )
  },
  by = modality
]

display_points <- rbindlist(lapply(seq_len(nrow(modality_spec)), function(i) {
  modality_name <- modality_spec$modality[[i]]
  x <- predictions[modality == modality_name, .(eid_int, BAG, BAU)]
  set.seed(202608230L + i)
  x[sample.int(nrow(x), 6000L)][, modality := modality_name]
}))

fwrite(population_summary, file.path(result_dir, "population_summary.tsv"), sep = "\t")
fwrite(age_sex_histogram, file.path(result_dir, "age_sex_histogram.tsv"), sep = "\t")
fwrite(fold_roles, file.path(result_dir, "fold_roles.tsv"), sep = "\t")
fwrite(performance, file.path(result_dir, "age_prediction_performance.tsv"), sep = "\t")
fwrite(prediction_density, file.path(result_dir, "age_prediction_density.tsv"), sep = "\t")
fwrite(relationship, file.path(result_dir, "bag_bau_relationship.tsv"), sep = "\t")
fwrite(display_points, file.path(result_dir, "bag_bau_display_sample.tsv"), sep = "\t")
