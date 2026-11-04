#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(survival)
})

args_all <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
default_root <- normalizePath(file.path(dirname(script_file), "../.."), mustWork = TRUE)
root <- Sys.getenv("BAU_CODE_ROOT", default_root)
data_root <- Sys.getenv("BAU_DATA_ROOT", "/home/dnanexus/bau-work/data")
results_root <- Sys.getenv("BAU_RESULTS_ROOT", "/home/dnanexus/bau-work/results")
external_root <- file.path(data_root, "external")
prediction_root <- file.path(results_root, "models", "external")
output_dir <- file.path(results_root, "outcomes", "external")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

normalise_id <- function(x) sub("\\.0$", "", trimws(as.character(x)))
as_factor <- function(x) {
  value <- trimws(as.character(x))
  value[is.na(x) | value %in% c("", "NA", "NaN")] <- "Missing"
  factor(value)
}
mode_value <- function(x, missing = c("", "NA", "NaN", "unknown")) {
  value <- trimws(as.character(x))
  bad <- is.na(x) | value %in% missing
  observed <- value[!bad]
  value[bad] <- names(sort(table(observed), decreasing = TRUE))[1]
  value
}

read_predictions <- function(cohort) {
  x <- as.data.table(read_feather(
    file.path(prediction_root, cohort, "predictions.feather")
  ))
  x[, participant_id := normalise_id(eid)]
  x[, .(
    participant_id,
    model_age = as.numeric(Age),
    model_sex = tolower(as.character(Sex)),
    BAG = as.numeric(BAG),
    BAU = as.numeric(BAU)
  )]
}

read_cohort <- function(cohort) {
  if (cohort == "charls") {
    covariates <- as.data.table(read_feather(
      file.path(external_root, "CHARLS", "candidate_qc.feather")
    ))
    survival <- as.data.table(read_feather(
      file.path(external_root, "CHARLS", "survival_outcomes.feather")
    ))
    covariates[, participant_id := normalise_id(participant_id)]
    survival[, participant_id := normalise_id(participant_id)]
    source <- merge(covariates, survival, by = "participant_id", all.x = TRUE)
    source <- source[
      age >= 40 & age <= 70 & sex %chin% c("female", "male") &
        blood_eligible_initial_lt20 == TRUE
    ]
  } else {
    source <- as.data.table(read_feather(file.path(
      external_root, toupper(cohort), "analysis_wide.feather"
    )))
    source[, participant_id := normalise_id(participant_id)]
    source <- source[age >= 40 & age <= 70]
  }
  predictions <- read_predictions(cohort)
  if (anyDuplicated(source$participant_id) || anyDuplicated(predictions$participant_id)) {
    stop("Participant identifiers must be unique")
  }
  x <- merge(predictions, source, by = "participant_id", all.x = TRUE)
  if (nrow(x) != nrow(predictions)) stop("Cohort join is incomplete")
  if (any(abs(x$model_age - x$age) > 1e-8) || any(x$model_sex != x$sex)) {
    stop("Age or sex differs between the prediction and cohort tables")
  }
  x[survival_eligible == TRUE & is.finite(time_midpoint_years) & time_midpoint_years > 0]
}

prepare_covariates <- function(x, cohort) {
  x[, sex_f := as_factor(sex)]
  x[, health_f := as_factor(baseline_health_status)]
  x[, smoking_f := as_factor(smoking_status)]
  if (cohort == "nhanes") {
    x[, cycle_f := as_factor(cycle)]
    x[, race_f := as_factor(race_ethnicity_raw)]
    x[, citizenship_f := factor(fcase(
      citizenship_raw == 1, "citizen",
      citizenship_raw == 2, "not_citizen",
      default = "unknown"
    ))]
    x[, education_f := as_factor(education_raw)]
    x[, income_f := factor(fcase(
      is.na(income_poverty_ratio), "missing",
      income_poverty_ratio < 1.3, "lt_1.3",
      income_poverty_ratio < 3.5, "1.3_to_3.5",
      default = "ge_3.5"
    ))]
    x[, alcohol_f := factor(fcase(
      is.na(alcohol_drinks_day), "unknown",
      alcohol_drinks_day <= 0, "none",
      alcohol_drinks_day <= 1, "up_to_1_per_day",
      default = "gt_1_per_day"
    ))]
  }
  if (cohort == "charls") {
    education_text <- as.character(x$education_raw)
    x[, education_f := factor(fcase(
      is.na(education_raw), "missing",
      grepl("^(1|2|3) ", education_text), "none_or_incomplete_primary",
      grepl("^4 ", education_text), "elementary",
      grepl("^5 ", education_text), "middle_school",
      default = "secondary_or_higher"
    ))]
    x[, alcohol_f := as_factor(alcohol_ever_status)]
  }
  if (cohort == "chns") {
    ethnicity <- mode_value(x$ethnicity_raw, c("", "NA", "NaN"))
    urban <- mode_value(x$urban_rural_raw, c("", "NA", "NaN"))
    working <- mode_value(x$working_raw, c("", "NA", "NaN"))
    insurance <- mode_value(x$medical_insurance_raw, c("", "NA", "NaN"))
    education <- as.numeric(mode_value(x$education_raw, c("", "NA", "NaN")))
    smoking <- mode_value(x$smoking_status)
    x[, ethnicity_f := factor(ifelse(as.numeric(ethnicity) == 1, "han", "other"))]
    x[, urban_f := factor(urban)]
    x[, working_f := factor(working)]
    x[, insurance_f := factor(insurance)]
    x[, education_f := factor(fcase(
      education <= 0, "none",
      education == 1, "primary",
      education == 2, "middle",
      default = "secondary_or_higher"
    ))]
    x[, smoking_f := factor(smoking)]
    x[, alcohol_f := as_factor(alcohol_frequency)]
    x[, health_f := factor(ifelse(
      health_status_excluding_cancer == "disease", "disease", "not_confirmed"
    ))]
  }
  x
}

covariates_for <- function(cohort) {
  switch(cohort,
    nhanes = paste(
      "age + sex_f + health_f + race_f + citizenship_f + education_f +",
      "income_f + smoking_f + alcohol_f + strata(cycle_f)"
    ),
    charls = "age + sex_f + education_f + smoking_f + alcohol_f + health_f",
    chns = paste(
      "age + sex_f + ethnicity_f + urban_f + working_f + education_f +",
      "insurance_f + smoking_f + alcohol_f + health_f"
    )
  )
}

fit_cox <- function(x, cohort, time_rule) {
  time <- x$time_midpoint_years
  event <- as.logical(x$event)
  if (time_rule == "earliest") {
    time[event] <- pmax(x$time_lower_years[event] + 1 / 365.25, 1 / 365.25)
  }
  if (time_rule == "upper") time[event] <- x$time_upper_bound_years[event]
  x[, analysis_time := time]
  formula <- as.formula(paste(
    "Surv(analysis_time, event) ~ BAG + BAU +", covariates_for(cohort)
  ))
  fit <- coxph(formula, data = x, ties = "efron", x = TRUE, y = TRUE)
  summary_fit <- summary(fit)
  coefficients <- as.data.table(summary_fit$coefficients, keep.rownames = "term")
  se_name <- if ("robust se" %chin% names(coefficients)) "robust se" else "se(coef)"
  output <- coefficients[term %chin% c("BAG", "BAU")]
  output[, .(
    cohort = toupper(cohort),
    time_rule,
    term,
    n = fit$n,
    events = fit$nevent,
    beta = coef,
    se = get(se_name),
    hr = `exp(coef)`,
    ci_lower = exp(coef - qnorm(0.975) * get(se_name)),
    ci_upper = exp(coef + qnorm(0.975) * get(se_name)),
    p = `Pr(>|z|)`
  )]
}

expected <- data.table(
  cohort = c("nhanes", "charls", "chns"),
  n = c(24225L, 9604L, 5137L),
  events = c(3003L, 806L, 142L)
)
results <- list()
index <- 1L
for (cohort_name in expected$cohort) {
  x <- prepare_covariates(read_cohort(cohort_name), cohort_name)
  row <- expected[cohort == cohort_name]
  if (nrow(x) != row$n || sum(x$event) != row$events) {
    stop(sprintf("Unexpected %s risk set", cohort_name))
  }
  rules <- if (cohort_name == "nhanes") "midpoint" else c("midpoint", "earliest", "upper")
  for (rule in rules) {
    results[[index]] <- fit_cox(x, cohort_name, rule)
    index <- index + 1L
  }
}

estimates <- rbindlist(results)
estimates[, primary := time_rule == "midpoint"]
estimates[primary == TRUE, q_primary := p.adjust(p, method = "BH")]
fwrite(estimates, file.path(output_dir, "mortality_estimates.tsv"), sep = "\t")
fwrite(
  estimates[primary == TRUE],
  file.path(output_dir, "primary_mortality_estimates.tsv"),
  sep = "\t"
)
