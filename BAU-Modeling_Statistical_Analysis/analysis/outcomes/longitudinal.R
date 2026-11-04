#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(survival)
})

args_all <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
default_root <- normalizePath(file.path(dirname(script_file), "../.."), mustWork = TRUE)
code_root <- Sys.getenv("BAU_CODE_ROOT", default_root)
results_root <- Sys.getenv("BAU_RESULTS_ROOT", "/home/dnanexus/bau-work/results")
source(file.path(code_root, "R/survival.R"))

outcomes <- data.table(
  outcome = c(
    "all-cause mortality", "loss of healthspan", "major disease multimorbidity",
    "coronary heart disease", "stroke", "chronic obstructive pulmonary disease",
    "liver diseases", "diabetes mellitus", "renal failure", "dementia",
    "digestive system cancer", "respiratory system cancer"
  ),
  type = c("mortality", "healthspan", "multimorbidity", rep("disease", 9L))
)

pairs <- as.data.table(read_feather(
  file.path(results_root, "phenotypes/repeat_pairs.feather")
))
support <- read_ukb_support(code_root)
base <- merge(pairs, support$covariates, by = "eid", all.x = TRUE, sort = FALSE)
base <- merge(base, support$outcomes, by = "eid", all.x = TRUE, sort = FALSE)
base <- merge(base, support$h7, by = "eid", all.x = TRUE, sort = FALSE)
if (nrow(base) != 13620L || anyDuplicated(base$eid)) stop("Incomplete repeat cohort linkage")

base[, `:=`(
  repeat_age_z = standardize(repeat_age),
  mvpa_z = standardize(as.numeric(as.character(`Daily moderate to vigorous physical activity`))),
  smoker = as.numeric(as.character(`Current smoker`)),
  daily_alcohol = as.numeric(as.character(`Daily alcohol intake`)),
  education = factor(as.character(Education)),
  income = factor(as.character(Income)),
  employment = factor(as.character(Employment)),
  centre = factor(as.character(baseline_center))
)]

base <- prepare_longitudinal_change(base)


formula <- Surv(time_days, event) ~
  baseline_bag + baseline_bau + change_bag + change_bau +
  repeat_age_z + sex_male + smoker + daily_alcohol + mvpa_z +
  education + income + employment + interval_z + strata(centre)

effect_rows <- list()
risk_rows <- list()
ph_rows <- list()
for (i in seq_len(nrow(outcomes))) {
  outcome <- outcomes$outcome[[i]]
  risk <- make_longitudinal_risk_set(base, outcome, outcomes$type[[i]])
  fit <- coxph(formula,
    data = risk, ties = "efron", x = TRUE, y = TRUE,
    control = coxph.control(iter.max = 60L, eps = 1e-9)
  )
  table <- summary(fit)$coefficients
  effect_rows[[i]] <- rbindlist(lapply(c("change_bag", "change_bau"), function(term) {
    beta <- table[term, "coef"]
    se <- table[term, "se(coef)"]
    data.table(
      outcome = outcome,
      exposure = if (term == "change_bag") "BAG change" else "BAU change",
      n_risk = nrow(risk),
      n_events = sum(risk$event),
      beta = beta,
      se = se,
      hr = exp(beta),
      ci_lower = exp(beta - qnorm(0.975) * se),
      ci_upper = exp(beta + qnorm(0.975) * se),
      p_value = 2 * pnorm(-abs(beta / se))
    )
  }))
  risk_rows[[i]] <- data.table(
    outcome = outcome,
    n_risk = nrow(risk),
    n_events = sum(risk$event),
    follow_up_q1_years = as.numeric(quantile(risk$time_days / 365.25, 0.25)),
    median_follow_up_years = median(risk$time_days / 365.25),
    follow_up_q3_years = as.numeric(quantile(risk$time_days / 365.25, 0.75))
  )
  zph <- cox.zph(fit, transform = "km")
  ph_rows[[i]] <- as.data.table(as.data.frame(zph$table), keep.rownames = "term")[
    , outcome := outcome
  ]
}

effects <- rbindlist(effect_rows)
effects[, q_value := p.adjust(p_value, method = "BH")]
change_summary <- rbindlist(lapply(c("bag", "bau"), function(phenotype) {
  difference <- base[[paste0("repeat_", phenotype)]] - base[[paste0("baseline_", phenotype)]]
  data.table(
    phenotype = toupper(phenotype), n = nrow(base),
    difference_mean = mean(difference), difference_sd = sd(difference)
  )
}))
output_dir <- file.path(results_root, "outcomes/longitudinal")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(effects, file.path(output_dir, "simple_change_estimates.tsv"), sep = "\t")
fwrite(rbindlist(risk_rows), file.path(output_dir, "landmark_risk_sets.tsv"), sep = "\t")
fwrite(rbindlist(ph_rows, fill = TRUE), file.path(output_dir, "proportional_hazards.tsv"), sep = "\t")
fwrite(change_summary, file.path(output_dir, "change_standardization.tsv"), sep = "\t")

message("Longitudinal analysis complete")
