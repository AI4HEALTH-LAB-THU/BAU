#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(jsonlite)
  library(matrixStats)
  library(splines)
})

args_all <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
default_root <- normalizePath(file.path(dirname(script_file), "../.."), mustWork = TRUE)
project_root <- Sys.getenv("BAU_CODE_ROOT", default_root)
data_root <- Sys.getenv("BAU_DATA_ROOT", "/home/dnanexus/bau-work/data")
results_root <- Sys.getenv("BAU_RESULTS_ROOT", "/home/dnanexus/bau-work/results")
config <- jsonlite::read_json(
  file.path(project_root, "config/organ_clocks.json"),
  simplifyVector = TRUE
)

phenotypes_dir <- file.path(results_root, "phenotypes")
output_dir <- file.path(phenotypes_dir, "domain_scores")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

model_root <- file.path(results_root, "models")
organ_root <- file.path(data_root, "ukb/domain_clocks")
role_root <- file.path(results_root, "phenotypes/domain_reference/protein")
protein_pred <- fread(file.path(model_root, "protein/ngboost/all_participant_predictions.csv"))
setnames(protein_pred, "eid_int", "eid")
setorder(protein_pred, eid)
stopifnot(!anyDuplicated(protein_pred$eid), nrow(protein_pred) == 43867L)
protein_roles <- lapply(seq_len(config$outer_folds), function(fold) {
  x <- fread(file.path(role_root, paste0("fold_", fold, ".csv")))
  setnames(x, c("eid_int", "role"), c("eid", "DevelopmentRole"))
  stopifnot(!anyDuplicated(x$eid))
  x
})
protein_domain_features <- setNames(lapply(config$protein_domains, function(organ) {
  readLines(file.path(project_root, "config/organ_features", paste0(organ, ".txt")))
}), config$protein_domains)
protein_domain_map <- rbindlist(lapply(config$protein_domains, function(organ) {
  data.table(domain = organ, feature = protein_domain_features[[organ]])
}))
fwrite(protein_domain_map, file.path(phenotypes_dir, "protein_domain_map.tsv"), sep = "\t")

impute_and_scale <- function(x_train, x_new) {
  med <- matrixStats::colMedians(x_train, na.rm = TRUE)
  med[!is.finite(med)] <- 0
  for (j in seq_len(ncol(x_train))) {
    miss_train <- is.na(x_train[, j])
    if (any(miss_train)) x_train[miss_train, j] <- med[[j]]
    miss_new <- is.na(x_new[, j])
    if (any(miss_new)) x_new[miss_new, j] <- med[[j]]
  }
  centre <- colMeans(x_train)
  scale <- matrixStats::colSds(x_train)
  scale[!is.finite(scale) | scale < 1e-8] <- 1
  x_train <- sweep(sweep(x_train, 2L, centre, "-"), 2L, scale, "/")
  x_new <- sweep(sweep(x_new, 2L, centre, "-"), 2L, scale, "/")
  list(train = x_train, new = x_new, median = med, centre = centre, scale = scale)
}

ridge_path_predict <- function(x_train, y_train, x_valid, alphas) {
  pp <- impute_and_scale(x_train, x_valid)
  x_train <- pp$train
  x_valid <- pp$new
  y_mean <- mean(y_train)
  y_centered <- y_train - y_mean
  n <- nrow(x_train)
  gram <- crossprod(x_train) / n
  rhs <- as.numeric(crossprod(x_train, y_centered) / n)
  eig <- eigen(gram, symmetric = TRUE)
  projected <- as.numeric(crossprod(eig$vectors, rhs))
  predictions <- vapply(alphas, function(alpha) {
    beta <- eig$vectors %*% (projected / (pmax(eig$values, 0) + alpha))
    as.numeric(y_mean + x_valid %*% beta)
  }, numeric(nrow(x_valid)))
  if (length(alphas) == 1L) predictions <- matrix(predictions, ncol = 1L)
  predictions
}

select_ridge_alpha <- function(x, y, inner_fold, alphas) {
  folds <- sort(unique(inner_fold))
  rmse <- matrix(NA_real_, nrow = length(folds), ncol = length(alphas))
  rownames(rmse) <- folds
  colnames(rmse) <- format(alphas, scientific = TRUE)
  for (i in seq_along(folds)) {
    fold <- folds[[i]]
    train <- which(inner_fold != fold)
    valid <- which(inner_fold == fold)
    pred <- ridge_path_predict(x[train, , drop = FALSE], y[train], x[valid, , drop = FALSE], alphas)
    error <- sweep(pred, 1L, y[valid], "-")
    rmse[i, ] <- sqrt(colMeans(error^2))
  }
  mean_rmse <- colMeans(rmse)
  se_rmse <- apply(rmse, 2L, sd) / sqrt(nrow(rmse))
  best <- which.min(mean_rmse)
  threshold <- mean_rmse[[best]] + se_rmse[[best]]
  eligible <- which(mean_rmse <= threshold + 1e-12)
  selected <- eligible[[which.max(alphas[eligible])]]
  list(
    alpha = alphas[[selected]],
    selected_index = selected,
    best_index = best,
    threshold = threshold,
    mean_rmse = mean_rmse,
    se_rmse = se_rmse,
    fold_rmse = rmse
  )
}

fit_ridge_predict <- function(x_train, y_train, x_new, alpha) {
  as.numeric(ridge_path_predict(x_train, y_train, x_new, alpha)[, 1L])
}

calibrate_domain_age <- function(cal_meta, cal_raw, test_meta, test_raw, spline_df) {
  cal_out <- copy(cal_meta)
  test_out <- copy(test_meta)
  cal_out[, predicted_age_raw := cal_raw]
  test_out[, predicted_age_raw := test_raw]
  cal_out[, `:=`(predicted_age = NA_real_, domain_BAG = NA_real_)]
  test_out[, `:=`(predicted_age = NA_real_, domain_BAG = NA_real_)]
  calibration_records <- list()

  for (sex_value in c(0L, 1L)) {
    cal_index <- which(cal_out$Sex_male == sex_value)
    test_index <- which(test_out$Sex_male == sex_value)
    stopifnot(length(cal_index) >= 250L)

    affine_data <- data.frame(
      Age = cal_out$Age[cal_index],
      predicted_age_raw = cal_out$predicted_age_raw[cal_index]
    )
    affine <- lm(Age ~ predicted_age_raw, data = affine_data)
    slope <- unname(coef(affine)[["predicted_age_raw"]])
    stopifnot(is.finite(slope), slope > config$minimum_calibration_slope)
    cal_out$predicted_age[cal_index] <- predict(affine, newdata = affine_data)
    test_out$predicted_age[test_index] <- predict(
      affine,
      newdata = data.frame(predicted_age_raw = test_out$predicted_age_raw[test_index])
    )

    cal_gap <- cal_out$predicted_age[cal_index] - cal_out$Age[cal_index]
    bias_data <- data.frame(Age = cal_out$Age[cal_index], age_gap = cal_gap)
    bias_model <- lm(age_gap ~ ns(Age, df = spline_df), data = bias_data)
    cal_bias <- predict(bias_model, newdata = bias_data)
    test_bias <- predict(
      bias_model,
      newdata = data.frame(Age = test_out$Age[test_index])
    )
    cal_residual <- cal_gap - cal_bias
    residual_sd <- sd(cal_residual)
    stopifnot(is.finite(residual_sd), residual_sd > 1e-6)
    cal_out$domain_BAG[cal_index] <- cal_residual / residual_sd
    test_gap <- test_out$predicted_age[test_index] - test_out$Age[test_index]
    test_out$domain_BAG[test_index] <- (test_gap - test_bias) / residual_sd

    calibration_records[[length(calibration_records) + 1L]] <- data.table(
      Sex_male = sex_value,
      n_calibration = length(cal_index),
      affine_intercept = unname(coef(affine)[["(Intercept)"]]),
      affine_slope = slope,
      domain_gap_reference_sd = residual_sd
    )
  }

  list(
    calibration = cal_out,
    test = test_out,
    records = rbindlist(calibration_records)
  )
}

fit_domain_across_outer_folds <- function(
  candidate, features, formal_meta, roles, score_ids, modality, domain
) {
  candidate[, eid := as.integer(eid)]
  stopifnot(!anyDuplicated(candidate$eid))
  keep_columns <- c("eid", features)
  xdt <- candidate[, ..keep_columns]
  dat <- merge(
    formal_meta[, .(eid, Age, Sex_male, Baseline_unhealthy, AssessmentCenter, OuterFold)],
    xdt,
    by = "eid",
    all = FALSE,
    sort = FALSE
  )
  setorder(dat, eid)
  x <- as.matrix(dat[, ..features])
  storage.mode(x) <- "double"
  alphas <- as.numeric(config$ridge_alpha_grid)

  test_scores <- list()
  calibration_scores <- list()
  tuning_records <- list()
  calibration_records <- list()

  for (outer_fold in seq_len(config$outer_folds)) {
    message("Fitting ", modality, " / ", domain, " / fold ", outer_fold)
    role <- roles[[outer_fold]]$DevelopmentRole[match(dat$eid, roles[[outer_fold]]$eid)]
    model_fit <- which(role == "model_fit")
    calibration <- which(role == "calibration_reference")
    outer_test <- which(dat$OuterFold == outer_fold)
    stopifnot(
      length(model_fit) > 1000L, length(calibration) > 250L, length(outer_test) > 250L,
      all(dat$Baseline_unhealthy[c(model_fit, calibration)] == 0L),
      !length(intersect(
        unique(as.character(dat$AssessmentCenter[outer_test])),
        unique(as.character(dat$AssessmentCenter[c(model_fit, calibration)]))
      ))
    )

    inner_fold <- dat$OuterFold[model_fit]
    stopifnot(length(unique(inner_fold)) == config$outer_folds - 1L)
    tuning <- select_ridge_alpha(
      x[model_fit, , drop = FALSE],
      dat$Age[model_fit],
      inner_fold,
      alphas
    )
    final_new <- rbind(x[calibration, , drop = FALSE], x[outer_test, , drop = FALSE])
    raw_prediction <- fit_ridge_predict(
      x[model_fit, , drop = FALSE],
      dat$Age[model_fit],
      final_new,
      tuning$alpha
    )
    n_cal <- length(calibration)
    raw_cal <- raw_prediction[seq_len(n_cal)]
    raw_test <- raw_prediction[n_cal + seq_along(outer_test)]
    calibrated <- calibrate_domain_age(
      dat[calibration, .(eid, Age, Sex_male, AssessmentCenter, OuterFold)],
      raw_cal,
      dat[outer_test, .(eid, Age, Sex_male, Baseline_unhealthy, AssessmentCenter, OuterFold)],
      raw_test,
      config$age_spline_df
    )

    calibrated$calibration[, `:=`(modality = modality, domain = domain, role = "calibration_reference")]
    calibrated$test[, `:=`(modality = modality, domain = domain, role = "outer_test")]
    calibration_scores[[outer_fold]] <- calibrated$calibration[, .(
      eid, Age, Sex_male, AssessmentCenter, OuterFold,
      ModelOuterFold = outer_fold,
      modality, domain, role, predicted_age, domain_BAG
    )]
    test_scores[[outer_fold]] <- calibrated$test[eid %in% score_ids, .(
      eid, Age, Sex_male, Baseline_unhealthy, AssessmentCenter, OuterFold,
      ModelOuterFold = outer_fold,
      modality, domain, role, predicted_age, domain_BAG
    )]

    tuning_records[[outer_fold]] <- data.table(
      modality = modality,
      domain = domain,
      outer_fold = outer_fold,
      alpha = alphas,
      mean_inner_rmse = tuning$mean_rmse,
      se_inner_rmse = tuning$se_rmse,
      selected = seq_along(alphas) == tuning$selected_index,
      minimum_rmse_alpha = seq_along(alphas) == tuning$best_index
    )
    cal_record <- copy(calibrated$records)
    cal_record[, `:=`(modality = modality, domain = domain, outer_fold = outer_fold)]
    calibration_records[[outer_fold]] <- cal_record
  }

  list(
    test = rbindlist(test_scores),
    calibration = rbindlist(calibration_scores),
    tuning = rbindlist(tuning_records),
    calibration_parameters = rbindlist(calibration_records)
  )
}

formal_meta_protein <- protein_pred[, .(
  eid, Age, Sex_male, Baseline_unhealthy, AssessmentCenter, OuterFold
)]
protein_results <- setNames(vector("list", length(config$protein_domains)), config$protein_domains)
missingness <- list()
for (organ in config$protein_domains) {
  candidate <- as.data.table(read_feather(file.path(organ_root, paste0(organ, ".feather"))))
  features <- protein_domain_features[[organ]]
  protein_results[[organ]] <- fit_domain_across_outer_folds(
    candidate, features, formal_meta_protein, protein_roles,
    protein_pred$eid, "protein", organ
  )
  missingness[[organ]] <- data.table(
    eid = candidate$eid,
    missing_n = rowSums(is.na(as.matrix(candidate[, ..features]))),
    feature_n = length(features)
  )
}
scores <- rbindlist(lapply(protein_results, `[[`, "test"))
stopifnot(!anyDuplicated(scores[, .(eid, domain)]), all(is.finite(scores$domain_BAG)))
write_feather(scores, file.path(output_dir, "protein_domain_test_scores.feather"))
write_feather(
  rbindlist(lapply(protein_results, `[[`, "calibration")),
  file.path(output_dir, "protein_domain_calibration_scores.feather")
)
fwrite(rbindlist(lapply(protein_results, `[[`, "tuning")),
  file.path(phenotypes_dir, "domain_clock_tuning.tsv"),
  sep = "\t"
)
fwrite(rbindlist(lapply(protein_results, `[[`, "calibration_parameters")),
  file.path(phenotypes_dir, "domain_clock_calibration_parameters.tsv"),
  sep = "\t"
)
performance <- scores[, .(
  n = .N, pearson_r = cor(Age, predicted_age),
  mae_years = mean(abs(Age - predicted_age))
), by = domain]
fwrite(performance, file.path(phenotypes_dir, "domain_clock_performance.tsv"), sep = "\t")
wide <- dcast(scores, eid ~ domain, value.var = "domain_BAG")
wide <- wide[complete.cases(wide)]
organs <- config$protein_domains
values <- as.matrix(wide[, ..organs])
wide[, `:=`(
  organ_mean = rowMeans(values),
  organ_disagreement = sqrt(rowMeans((values - rowMeans(values))^2))
)]
miss <- rbindlist(missingness)[, .(organ_missing = sum(missing_n) / sum(feature_n)), by = eid]
wide <- merge(wide, miss, by = "eid")
wide <- merge(wide, protein_pred[, .(
  eid, Age, Sex_male, AssessmentCenter, OuterFold,
  Baseline_unhealthy, BAG, BAU
)], by = "eid")
stopifnot(nrow(wide) == 35059L)
write_feather(wide, file.path(phenotypes_dir, "organ_disagreement.feather"))
fwrite(data.table(n_all = uniqueN(scores$eid), n_complete = nrow(wide), n_organs = length(organs)),
  file.path(phenotypes_dir, "organ_coverage.tsv"),
  sep = "\t"
)
