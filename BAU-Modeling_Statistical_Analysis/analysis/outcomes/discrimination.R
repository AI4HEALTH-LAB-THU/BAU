#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

args_all <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
default_root <- normalizePath(file.path(dirname(script_file), "../.."), mustWork = TRUE)
code_root <- Sys.getenv("BAU_CODE_ROOT", default_root)
results_root <- Sys.getenv("BAU_RESULTS_ROOT", "/home/dnanexus/bau-work/results")
source(file.path(code_root, "R/survival.R"))

horizon <- 10
models <- c("clinical", "clinical_bag", "clinical_bau", "clinical_bag_bau")
outcomes <- c(
  "all-cause mortality", "loss of healthspan", "major disease multimorbidity",
  "coronary heart disease", "stroke", "chronic obstructive pulmonary disease",
  "liver diseases", "diabetes mellitus", "renal failure", "dementia",
  "digestive system cancer", "respiratory system cancer"
)

model_matrix <- function(risk, model) {
  columns <- c(
    "Age", "Sex_male", "current_smoking", "diabetes_history",
    "systolic_blood_pressure", "total_cholesterol", "hdl_cholesterol",
    "antihypertensive_treatment"
  )
  if (model %in% c("clinical_bag", "clinical_bag_bau")) columns <- c(columns, "BAG")
  if (model %in% c("clinical_bau", "clinical_bag_bau")) columns <- c(columns, "BAU")
  as.matrix(risk[, ..columns])
}

fit_fold <- function(x_train, time_train, event_train, x_test) {
  variable <- apply(x_train, 2L, function(x) is.finite(sd(x)) && sd(x) > 1e-10)
  train <- x_train[, variable, drop = FALSE]
  test <- x_test[, variable, drop = FALSE]
  fit <- coxph(
    Surv(time_train, event_train) ~ train,
    ties = "efron",
    singular.ok = FALSE,
    control = coxph.control(iter.max = 40L, eps = 1e-8)
  )
  beta <- unname(coef(fit))
  baseline <- basehaz(fit, centered = TRUE)
  eligible <- baseline$time <= horizon
  hazard <- if (any(eligible)) tail(baseline$hazard[eligible], 1L) else 0
  offset <- sum(fit$means * beta)
  linear_predictor <- drop(test %*% beta) - offset
  -expm1(-hazard * exp(linear_predictor))
}

analyse_outcome <- function(risk, modality, outcome) {
  scores <- matrix(NA_real_, nrow(risk), length(models), dimnames = list(NULL, models))
  for (fold in 1:5) {
    train <- risk$OuterFold != fold
    test <- !train
    for (model in models) {
      x <- model_matrix(risk, model)
      scores[test, model] <- fit_fold(
        x[train, , drop = FALSE],
        risk$time_years[train],
        risk$event[train],
        x[test, , drop = FALSE]
      )
    }
  }
  if (!all(is.finite(scores))) stop(paste("Incomplete held-fold scores for", outcome))

  score_data <- data.frame(
    time = risk$time_years,
    event = risk$event,
    centre = risk$centre,
    scores,
    check.names = FALSE
  )
  concordance_fit <- concordance(
    Surv(time, event) ~ clinical + clinical_bag + clinical_bau + clinical_bag_bau,
    data = score_data,
    cluster = centre,
    ymax = horizon,
    timewt = "n/G2",
    reverse = TRUE,
    influence = 0
  )
  estimates <- unname(coef(concordance_fit))
  names(estimates) <- models
  covariance <- vcov(concordance_fit)
  rownames(covariance) <- colnames(covariance) <- models
  se <- sqrt(diag(covariance))
  absolute <- data.table(
    modality = modality,
    outcome = outcome,
    model = models,
    estimate = estimates,
    se = se,
    ci_lower = estimates - qnorm(0.975) * se,
    ci_upper = estimates + qnorm(0.975) * se,
    n_risk = nrow(risk),
    n_events = sum(risk$event)
  )

  increments <- rbindlist(lapply(models[-1L], function(model) {
    contrast <- setNames(rep(0, length(models)), models)
    contrast[[model]] <- 1
    contrast[["clinical"]] <- -1
    delta <- sum(contrast * estimates)
    delta_se <- sqrt(as.numeric(t(contrast) %*% covariance %*% contrast))
    data.table(
      modality = modality,
      outcome = outcome,
      model = model,
      reference = "clinical",
      estimate_x1000 = 1000 * delta,
      se_x1000 = 1000 * delta_se,
      ci_lower_x1000 = 1000 * (delta - qnorm(0.975) * delta_se),
      ci_upper_x1000 = 1000 * (delta + qnorm(0.975) * delta_se),
      n_risk = nrow(risk),
      n_events = sum(risk$event)
    )
  }))

  folds <- rbindlist(lapply(1:5, function(fold) {
    keep <- risk$OuterFold == fold
    rbindlist(lapply(models, function(model) {
      held <- data.frame(
        time = risk$time_years[keep], event = risk$event[keep],
        score = scores[keep, model]
      )
      fit <- concordance(
        Surv(time, event) ~ score,
        data = held,
        ymax = horizon,
        timewt = "n/G2",
        reverse = TRUE
      )
      data.table(
        modality = modality,
        outcome = outcome,
        fold = fold,
        model = model,
        c_index = unname(coef(fit)),
        n_test = sum(keep),
        events_test = sum(risk$event[keep])
      )
    }))
  }))
  list(absolute = absolute, increments = increments, folds = folds)
}

support <- read_ukb_support(code_root)
absolute_rows <- list()
increment_rows <- list()
fold_rows <- list()
for (modality in c("nmr", "protein")) {
  message("Preparing discrimination analysis for ", modality)
  analysis <- prepare_survival_data(code_root, modality, "ngboost", support)
  analysis[, `:=`(
    current_smoking = as.numeric(as.character(`Current smoker`)),
    diabetes_history = as.numeric(as.character(`Diabetes history`)),
    systolic_blood_pressure = as.numeric(`Systolic blood pressure`),
    total_cholesterol = as.numeric(Cholesterol),
    hdl_cholesterol = as.numeric(`HDL cholesterol`),
    antihypertensive_treatment = as.numeric(as.character(`Hypotensive drugs`))
  )]
  for (outcome in outcomes) {
    message("Scoring ", modality, " / ", outcome)
    risk <- make_risk_set(analysis, outcome)
    result <- analyse_outcome(risk, modality, outcome)
    absolute_rows[[length(absolute_rows) + 1L]] <- result$absolute
    increment_rows[[length(increment_rows) + 1L]] <- result$increments
    fold_rows[[length(fold_rows) + 1L]] <- result$folds
  }
}

output_dir <- file.path(results_root, "outcomes/discrimination")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(rbindlist(absolute_rows), file.path(output_dir, "absolute_c_index.tsv"), sep = "\t")
fwrite(rbindlist(increment_rows), file.path(output_dir, "incremental_c_index.tsv"), sep = "\t")
fwrite(rbindlist(fold_rows), file.path(output_dir, "fold_c_index.tsv"), sep = "\t")

message("Discrimination analysis complete")
