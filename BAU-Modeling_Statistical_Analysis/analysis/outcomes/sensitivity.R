#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

args_all <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
default_root <- normalizePath(file.path(dirname(script_file), "../.."), mustWork = TRUE)
code_root <- Sys.getenv("BAU_CODE_ROOT", default_root)
ukb_root <- Sys.getenv("UKB_DATA_ROOT", "/mnt/project/UKB")
results_root <- Sys.getenv("BAU_RESULTS_ROOT", "/home/dnanexus/bau-work/results")
source(file.path(code_root, "R/survival.R"))

outcomes <- c(
  "all-cause mortality", "loss of healthspan", "major disease multimorbidity",
  "coronary heart disease", "stroke", "chronic obstructive pulmonary disease",
  "liver diseases", "diabetes mellitus", "renal failure", "dementia",
  "digestive system cancer", "respiratory system cancer"
)
competing_outcomes <- outcomes[-c(1L, 2L)]
output_dir <- file.path(results_root, "outcomes/sensitivity")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

effect_table <- function(fit, modality, outcome, analysis, period = NA_character_,
                         n_risk = fit$n, n_events = fit$nevent) {
  table <- summary(fit)$coefficients
  se_name <- if ("robust se" %in% colnames(table)) "robust se" else "se(coef)"
  rbindlist(lapply(c("BAG", "BAU"), function(term) {
    beta <- table[term, "coef"]
    se <- table[term, se_name]
    data.table(
      analysis = analysis,
      period = period,
      modality = modality,
      outcome = outcome,
      exposure = term,
      n_risk = n_risk,
      n_events = n_events,
      beta = beta,
      se = se,
      hr = exp(beta),
      ci_lower = exp(beta - qnorm(0.975) * se),
      ci_upper = exp(beta + qnorm(0.975) * se),
      p_value = 2 * pnorm(-abs(beta / se))
    )
  }))
}

fit_with_formula <- function(risk, formula, modality, outcome, analysis, period = NA_character_) {
  fit <- coxph(formula, data = risk, ties = "efron", robust = TRUE)
  effect_table(fit, modality, outcome, analysis, period)
}

support <- read_ukb_support(code_root)
history_columns <- setdiff(names(support$history), c("eid", "Baseline_unhealthy"))

drug_path <- file.path(
  ukb_root, "RawData/02-AssessmentCenter/02-Touchscreen/drug_0.rds"
)
drug <- as.data.table(readRDS(drug_path))
drug_columns <- setdiff(names(drug), "eid")
drug_values <- as.matrix(drug[, ..drug_columns])
known <- rowSums(matrix(
  drug_values %in% c(-7, 1, 2, 3, 4, 5),
  nrow = nrow(drug_values)
), na.rm = TRUE) > 0L
medication <- data.table(
  eid = as.integer(drug$eid),
  lipid_lowering = factor(
    fifelse(!known, "Unknown", fifelse(rowSums(drug_values == 1, na.rm = TRUE) > 0, "Yes", "No")),
    levels = c("No", "Yes", "Unknown")
  )
)
rm(drug, drug_values, known)

extended_formula <- function(multimorbidity) {
  terms <- c(
    "BAG", "BAU", "age_z", "Sex_male", "smoker", "daily_alcohol",
    "mvpa_z", "education", "income", "employment", "healthy_food",
    "processed_meat", "red_meat", "frequent_salt", "family_cvd",
    "family_diabetes", "hypotensive_drugs", "insulin", "lipid_lowering",
    "bmi_z", "systolic_bp_z", "baseline_disease_burden", "strata(centre)",
    if (multimorbidity) "strata(baseline_h7)" else character(),
    "cluster(centre)"
  )
  as.formula(paste("Surv(time_years, event) ~", paste(terms, collapse = " + ")))
}

extended_rows <- list()
landmark_rows <- list()
piecewise_rows <- list()
competing_rows <- list()

for (modality in c("nmr", "protein")) {
  message("Preparing sensitivity analyses for ", modality)
  analysis <- prepare_survival_data(code_root, modality, "ngboost", support)
  analysis <- merge(analysis, medication, by = "eid", all.x = TRUE, sort = FALSE)
  analysis[, `:=`(
    healthy_food = as.numeric(as.character(`Daily healthy food`)),
    processed_meat = as.numeric(as.character(`Frequently processed meat`)),
    red_meat = as.numeric(as.character(`Frequently red meat`)),
    frequent_salt = as.numeric(as.character(`Frequently salt intake`)),
    family_cvd = as.numeric(as.character(`Family cardiovascular disease history`)),
    family_diabetes = as.numeric(as.character(`Family diabetes history`)),
    hypotensive_drugs = as.numeric(as.character(`Hypotensive drugs`)),
    insulin = as.numeric(as.character(Insulin)),
    bmi_z = standardize(BMI),
    systolic_bp_z = standardize(`Systolic blood pressure`),
    baseline_disease_count = rowSums(.SD)
  ), .SDcols = history_columns]
  analysis[, baseline_disease_burden := factor(
    fifelse(
      baseline_disease_count == 0L, "0",
      fifelse(baseline_disease_count == 1L, "1", "2+")
    ),
    levels = c("0", "1", "2+")
  )]

  for (outcome in outcomes) {
    risk <- make_risk_set(analysis, outcome)
    extended_rows[[length(extended_rows) + 1L]] <- fit_with_formula(
      risk,
      extended_formula(outcome == "major disease multimorbidity"),
      modality,
      outcome,
      "extended_covariates"
    )

    landmark <- copy(risk[time_years > 1])
    landmark[, time_years := time_years - 1]
    landmark_rows[[length(landmark_rows) + 1L]] <- fit_with_formula(
      landmark,
      joint_cox_formula(outcome == "major disease multimorbidity"),
      modality,
      outcome,
      "one_year_landmark"
    )

    windows <- data.table(
      period = c("0-5 years", "5-10 years", ">10 years"),
      start = c(0, 5, 10),
      end = c(5, 10, Inf)
    )
    for (i in seq_len(nrow(windows))) {
      window <- windows[i]
      period_risk <- copy(risk[time_years > window$start])
      period_risk[, event := as.integer(event == 1L & time_years <= window$end)]
      period_risk[, time_years := pmin(time_years, window$end) - window$start]
      piecewise_rows[[length(piecewise_rows) + 1L]] <- fit_with_formula(
        period_risk,
        joint_cox_formula(outcome == "major disease multimorbidity"),
        modality,
        outcome,
        "follow_up_period",
        window$period
      )
    }

    if (outcome %in% competing_outcomes) {
      cox <- fit_joint_cox(risk, outcome, modality)$estimates
      cox[, analysis := "cause_specific_cox"]
      competing_rows[[length(competing_rows) + 1L]] <- cox

      fg_input <- copy(risk)
      fg_input[, status := factor(
        fifelse(event == 1L, "target", fifelse(competing_event == 1L, "death", "censor")),
        levels = c("censor", "target", "death")
      )]
      fg_input[, `:=`(row_id = .I, centre_keep = centre)]
      if (outcome == "major disease multimorbidity") {
        fg_input[, `:=`(
          baseline_keep = baseline_h7,
          weight_stratum = interaction(centre, baseline_h7, drop = TRUE)
        )]
        carry <- paste(
          "row_id + BAG + BAU + age_z + Sex_male + smoker + daily_alcohol +",
          "mvpa_z + education + income + employment + centre_keep + baseline_keep +",
          "strata(weight_stratum)"
        )
      } else {
        fg_input[, weight_stratum := centre]
        carry <- paste(
          "row_id + BAG + BAU + age_z + Sex_male + smoker + daily_alcohol +",
          "mvpa_z + education + income + employment + centre_keep +",
          "strata(weight_stratum)"
        )
      }
      expanded <- finegray(
        as.formula(paste("Surv(time_years, status) ~", carry)),
        data = fg_input,
        etype = "target",
        prefix = "fg"
      )
      rhs <- c(
        "BAG", "BAU", "age_z", "Sex_male", "smoker", "daily_alcohol",
        "mvpa_z", "education", "income", "employment", "strata(centre_keep)",
        if (outcome == "major disease multimorbidity") "strata(baseline_keep)" else character(),
        "cluster(centre_keep)"
      )
      fg_fit <- coxph(
        as.formula(paste(
          "Surv(fgstart, fgstop, fgstatus) ~", paste(rhs, collapse = " + ")
        )),
        data = expanded,
        weights = fgwt,
        ties = "efron",
        robust = TRUE
      )
      competing_rows[[length(competing_rows) + 1L]] <- effect_table(
        fg_fit, modality, outcome, "fine_gray",
        n_risk = nrow(risk), n_events = sum(risk$event)
      )
    }
    rm(risk)
    gc(FALSE)
  }
  rm(analysis)
  gc(FALSE)
}

fwrite(rbindlist(extended_rows), file.path(output_dir, "extended_covariates.tsv"), sep = "\t")
fwrite(rbindlist(landmark_rows), file.path(output_dir, "one_year_landmark.tsv"), sep = "\t")
fwrite(rbindlist(piecewise_rows), file.path(output_dir, "follow_up_periods.tsv"), sep = "\t")
fwrite(rbindlist(competing_rows, fill = TRUE), file.path(output_dir, "competing_risk.tsv"), sep = "\t")

message("Outcome sensitivity analyses complete")
