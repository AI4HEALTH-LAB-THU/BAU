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

outcomes <- c(
  "all-cause mortality", "loss of healthspan", "major disease multimorbidity"
)
modifiers <- list(
  age = list(
    levels = c("40-49", "50-59", "60+"),
    make = function(x) fifelse(x$Age < 50, "40-49", fifelse(x$Age < 60, "50-59", "60+"))
  ),
  sex = list(
    levels = c("female", "male"),
    make = function(x) fifelse(x$Sex_male == 1L, "male", "female")
  ),
  baseline_health = list(
    levels = c("no_history", "disease_history"),
    make = function(x) fifelse(x$Baseline_unhealthy == 0L, "no_history", "disease_history")
  ),
  bmi = list(
    levels = c("<25", "25-<30", "30+"),
    make = function(x) fifelse(x$BMI < 25, "<25", fifelse(x$BMI < 30, "25-<30", "30+"))
  ),
  current_smoking = list(
    levels = c("not_current", "current"),
    make = function(x) fifelse(x$smoker == 1, "current", "not_current")
  ),
  physical_activity = list(
    levels = c("<30_min_day", "30+_min_day"),
    make = function(x) {
      fifelse(
        as.numeric(as.character(x[["Daily moderate to vigorous physical activity"]])) < 30,
        "<30_min_day",
        "30+_min_day"
      )
    }
  )
)

interaction_name <- function(names, exposure, level) {
  candidates <- c(
    paste0(exposure, ":subgroup", level),
    paste0("subgroup", level, ":", exposure)
  )
  found <- intersect(candidates, names)
  if (length(found) != 1L) stop("Interaction coefficient was not found")
  found[[1L]]
}

wald_test <- function(beta, variance) {
  statistic <- as.numeric(crossprod(beta, solve(variance, beta)))
  data.table(df = length(beta), chisq = statistic, p_value = pchisq(statistic, length(beta), lower.tail = FALSE))
}

fit_subgroup <- function(risk, modality, outcome, modifier, specification) {
  risk <- copy(risk)
  risk[, subgroup := factor(specification$make(risk), levels = specification$levels)]
  covariates <- c(
    "age_z", "Sex_male", "smoker", "daily_alcohol", "mvpa_z",
    "education", "income", "employment"
  )
  if (modifier == "sex") covariates <- setdiff(covariates, "Sex_male")
  if (modifier == "current_smoking") covariates <- setdiff(covariates, "smoker")
  terms <- c(
    "BAG", "BAU", "subgroup", "BAG:subgroup", "BAU:subgroup", covariates,
    "strata(centre)",
    if (outcome == "major disease multimorbidity") "strata(baseline_h7)" else character(),
    "cluster(centre)"
  )
  formula <- as.formula(paste("Surv(time_years, event) ~", paste(terms, collapse = " + ")))
  fit <- coxph(formula, data = risk, ties = "efron", x = TRUE)
  coefficient <- coef(fit)
  variance <- fit$var
  names(coefficient) <- rownames(variance) <- colnames(variance) <- names(coef(fit))
  effect_rows <- list()
  interaction_rows <- list()

  for (exposure in c("BAG", "BAU")) {
    interaction_terms <- vapply(
      specification$levels[-1L],
      function(level) interaction_name(names(coefficient), exposure, level),
      character(1)
    )
    test <- wald_test(
      coefficient[interaction_terms],
      variance[interaction_terms, interaction_terms, drop = FALSE]
    )
    interaction_rows[[exposure]] <- data.table(
      modality = modality,
      outcome = outcome,
      modifier = modifier,
      exposure = exposure,
      df = test$df,
      chisq = test$chisq,
      p_value = test$p_value
    )

    for (level in specification$levels) {
      contrast <- setNames(rep(0, length(coefficient)), names(coefficient))
      contrast[[exposure]] <- 1
      if (level != specification$levels[[1L]]) {
        contrast[[interaction_name(names(coefficient), exposure, level)]] <- 1
      }
      beta <- sum(contrast * coefficient)
      se <- sqrt(as.numeric(t(contrast) %*% variance %*% contrast))
      count <- risk[subgroup == level, .(n = .N, events = sum(event))]
      effect_rows[[paste(exposure, level)]] <- data.table(
        modality = modality,
        outcome = outcome,
        modifier = modifier,
        level = level,
        exposure = exposure,
        n_risk = count$n,
        n_events = count$events,
        beta = beta,
        se = se,
        hr = exp(beta),
        ci_lower = exp(beta - qnorm(0.975) * se),
        ci_upper = exp(beta + qnorm(0.975) * se)
      )
    }
  }
  list(
    effects = rbindlist(effect_rows),
    interactions = rbindlist(interaction_rows),
    model = data.table(
      modality = modality,
      outcome = outcome,
      modifier = modifier,
      n_risk = nrow(risk),
      n_events = sum(risk$event),
      coefficients = sum(is.finite(coefficient))
    )
  )
}

support <- read_ukb_support(code_root)
effect_rows <- list()
interaction_rows <- list()
model_rows <- list()

for (modality in c("nmr", "protein")) {
  message("Preparing subgroup analyses for ", modality)
  analysis <- prepare_survival_data(code_root, modality, "ngboost", support)
  for (outcome in outcomes) {
    risk <- make_risk_set(analysis, outcome)
    for (modifier in names(modifiers)) {
      fitted <- fit_subgroup(risk, modality, outcome, modifier, modifiers[[modifier]])
      effect_rows[[length(effect_rows) + 1L]] <- fitted$effects
      interaction_rows[[length(interaction_rows) + 1L]] <- fitted$interactions
      model_rows[[length(model_rows) + 1L]] <- fitted$model
    }
  }
}

effects <- rbindlist(effect_rows)
interactions <- rbindlist(interaction_rows)
interactions[, q_value := p.adjust(p_value, method = "BH")]
models <- rbindlist(model_rows)
output_dir <- file.path(results_root, "outcomes/subgroups")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(effects, file.path(output_dir, "subgroup_estimates.tsv"), sep = "\t")
fwrite(interactions, file.path(output_dir, "interaction_tests.tsv"), sep = "\t")
fwrite(models, file.path(output_dir, "model_summary.tsv"), sep = "\t")

message("Subgroup analysis complete")
