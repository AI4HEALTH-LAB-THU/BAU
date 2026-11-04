#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(splines)
})
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- Sys.getenv("BAU_CODE_ROOT", normalizePath(file.path(dirname(script), "../..")))
data_root <- Sys.getenv("BAU_DATA_ROOT", "/home/dnanexus/bau-work/data")
results_root <- Sys.getenv("BAU_RESULTS_ROOT", "/home/dnanexus/bau-work/results")
age_grid <- seq(40, 70, by = 0.25)

fit_support <- function(modality) {
  output <- file.path(results_root, "phenotypes/molecular", modality)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  models <- file.path(results_root, "models", modality, "ngboost")
  population <- fread(file.path(models, "all_participant_predictions.csv"))
  setnames(population, "eid_int", "eid")
  setorder(population, eid)
  features <- readLines(file.path(root, "config/features", paste0(modality, ".txt")))
  features <- features[nzchar(features)]
  candidate <- as.data.table(read_feather(file.path(data_root, "ukb/modeling", paste0(modality, ".feather"))))
  stopifnot(!anyDuplicated(candidate$eid), setequal(population$eid, candidate$eid))
  candidate <- candidate[match(population$eid, eid)]
  x <- as.matrix(candidate[, ..features])
  storage.mode(x) <- "double"
  x[!is.finite(x)] <- NA_real_
  population <- population[, .(eid, Age, Sex_male, AssessmentCenter, Baseline_unhealthy, OuterFold, BAG, BAU)]
  population[, `:=`(n_observed = rowSums(is.finite(x)), input_missing_fraction = rowMeans(!is.finite(x)))]
  write_feather(population, file.path(output, "population.feather"))
  write_feather(cbind(population[, .(eid)], as.data.table(x)), file.path(output, "measurements.feather"))
  parameters <- list()
  for (fold in 1:5) {
    roles <- fread(file.path(models, paste0("outer_fold_", fold), "development_roles.csv"))
    stopifnot(!anyDuplicated(roles$eid_int))
    role <- roles$role[match(population$eid, roles$eid_int)]
    fit <- which(role == "model_fit")
    calibration <- which(role == "calibration_reference")
    test <- which(population$OuterFold == fold)
    stopifnot(
      all(population$Baseline_unhealthy[c(fit, calibration)] == 0),
      !length(intersect(population$AssessmentCenter[c(fit, calibration)], population$AssessmentCenter[test]))
    )
    for (sex in 0:1) {
      fi <- fit[population$Sex_male[fit] == sex]
      ca <- calibration[population$Sex_male[calibration] == sex]
      knots <- as.numeric(quantile(population$Age[fi], c(0.25, 0.5, 0.75)))
      design <- function(age) cbind(1, ns(age, knots = knots, Boundary.knots = c(40, 70)))
      fit_design <- design(population$Age[fi])
      cal_design <- design(population$Age[ca])
      grid_design <- design(age_grid)
      means <- matrix(NA_real_, length(age_grid), length(features))
      scales <- shifts <- numeric(length(features))
      coefficients <- matrix(NA_real_, ncol(fit_design), length(features))
      for (j in seq_along(features)) {
        observed <- is.finite(x[fi, j])
        beta <- lm.fit(fit_design[observed, , drop = FALSE], x[fi, j][observed])$coefficients
        residual <- x[ca, j] - drop(cal_design %*% beta)
        shifts[j] <- mean(residual, na.rm = TRUE)
        scales[j] <- sd(residual, na.rm = TRUE)
        coefficients[, j] <- beta
        stopifnot(all(is.finite(beta)), is.finite(shifts[j]), is.finite(scales[j]), scales[j] > 0)
        means[, j] <- drop(grid_design %*% beta) + shifts[j]
      }
      key <- paste(fold, sex, sep = "_")
      fwrite(data.table(feature = features, sd = scales, shift = shifts),
        file.path(output, paste0("parameters_", key, ".tsv")),
        sep = "\t"
      )
      colnames(means) <- features
      write_feather(
        data.table(age = age_grid, as.data.table(means)),
        file.path(output, paste0("mean_grid_", key, ".feather"))
      )
      parameters[[key]] <- list(
        knots = knots, boundary = c(40, 70), coefficients = coefficients,
        shift = shifts, sd = scales, features = features
      )
    }
    message(modality, ": completed fold ", fold)
  }
  saveRDS(parameters, file.path(output, "support_models.rds"))
}

for (modality in c("protein", "nmr")) fit_support(modality)
