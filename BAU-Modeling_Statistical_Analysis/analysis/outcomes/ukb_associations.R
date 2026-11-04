#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args_all <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
default_root <- normalizePath(file.path(dirname(script_file), "../.."), mustWork = TRUE)
code_root <- Sys.getenv("BAU_CODE_ROOT", default_root)
results_root <- Sys.getenv("BAU_RESULTS_ROOT", "/home/dnanexus/bau-work/results")
source(file.path(code_root, "R/survival.R"))

arg_value <- function(name, default) {
  item <- grep(paste0("^--", name, "="), commandArgs(trailingOnly = TRUE), value = TRUE)
  if (length(item)) sub(paste0("^--", name, "="), "", item[[1L]]) else default
}

model <- arg_value("model", "ngboost")
modalities <- strsplit(arg_value(
  "modalities",
  if (model == "ngboost") "nmr,protein,blood" else "nmr,protein"
), ",", fixed = TRUE)[[1L]]
if (!model %in% c("ngboost", "lightgbmlss", "xgboostlss")) stop("Unknown model")
if (!all(modalities %in% c("nmr", "protein", "blood"))) stop("Unknown modality")

primary_outcomes <- c(
  "all-cause mortality",
  "loss of healthspan",
  "major disease multimorbidity"
)
disease_outcomes <- c(
  "coronary heart disease",
  "stroke",
  "chronic obstructive pulmonary disease",
  "liver diseases",
  "diabetes mellitus",
  "renal failure",
  "dementia",
  "digestive system cancer",
  "respiratory system cancer"
)
outcomes <- c(primary_outcomes, disease_outcomes)

add_primary_outcome_fdr <- function(estimates, primary_outcomes) {
  x <- copy(estimates)
  x[, `:=`(q_primary_omics_12 = NA_real_, q_primary_blood_6 = NA_real_)]
  families <- list(
    list(modalities = c("nmr", "protein"), n = 12L, column = "q_primary_omics_12"),
    list(modalities = "blood", n = 6L, column = "q_primary_blood_6")
  )
  for (family in families) {
    rows <- which(x$model == "ngboost" & x$outcome %chin% primary_outcomes &
      x$modality %chin% family$modalities)
    if (length(rows) != family$n) next
    if (uniqueN(x[rows, .(modality, outcome, exposure)]) != family$n) {
      stop("Repeated estimates in a primary-outcome FDR family")
    }
    x[rows, (family$column) := p.adjust(p_value, method = "BH", n = family$n)]
  }
  x
}

output_dir <- file.path(results_root, "outcomes", model)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

support <- read_ukb_support(code_root)
estimate_rows <- list()
ph_rows <- list()
curve_rows <- list()
risk_rows <- list()
state_rows <- list()

for (modality in modalities) {
  message("Preparing ", modality, " data")
  analysis <- prepare_survival_data(code_root, modality, model, support)
  analysis <- assign_bag_bau_state(analysis)
  curve_outcomes <- if (model != "ngboost") {
    character()
  } else if (modality == "blood") {
    primary_outcomes
  } else {
    outcomes
  }

  for (outcome in outcomes) {
    message("Fitting ", modality, " / ", outcome)
    risk <- make_risk_set(analysis, outcome)
    fitted <- fit_joint_cox(risk, outcome, modality, model)
    estimate_rows[[length(estimate_rows) + 1L]] <- fitted$estimates
    if (!is.null(fitted$ph)) {
      ph <- as.data.table(as.data.frame(fitted$ph$table), keep.rownames = "term")
      ph <- ph[term %in% c("BAG", "BAU", "GLOBAL")]
      ph[, `:=`(modality = modality, model = model, outcome = outcome)]
      ph_rows[[length(ph_rows) + 1L]] <- ph
    }

    if (outcome %in% curve_outcomes) {
      curve <- if (outcome == "major disease multimorbidity") {
        extract_aj_curve(risk, outcome)
      } else {
        extract_km_curve(risk, outcome)
      }
      curve <- curve[time_years <= 15]
      curve[, `:=`(modality = modality, model = model)]
      curve_rows[[length(curve_rows) + 1L]] <- curve
      risk_rows[[length(risk_rows) + 1L]] <- risk_table(risk, outcome)[
        , `:=`(modality = modality, model = model)
      ]
      state_rows[[length(state_rows) + 1L]] <- risk[, .(
        n = .N,
        events = sum(event),
        competing_deaths = sum(competing_event)
      ), by = state][, `:=`(modality = modality, model = model, outcome = outcome)]
    }
    rm(risk, fitted)
    gc(FALSE)
  }
  rm(analysis)
  gc(FALSE)
}

estimates <- rbindlist(estimate_rows, use.names = TRUE, fill = TRUE)
estimates <- add_primary_outcome_fdr(estimates, primary_outcomes)
fwrite(estimates, file.path(output_dir, "joint_cox_estimates.tsv"), sep = "\t")
if (length(ph_rows)) {
  fwrite(rbindlist(ph_rows, fill = TRUE), file.path(output_dir, "proportional_hazards.tsv"), sep = "\t")
}
if (length(curve_rows)) {
  fwrite(rbindlist(curve_rows), file.path(output_dir, "cumulative_event_curves.tsv"), sep = "\t")
  fwrite(rbindlist(risk_rows), file.path(output_dir, "numbers_at_risk.tsv"), sep = "\t")
  fwrite(rbindlist(state_rows), file.path(output_dir, "state_counts.tsv"), sep = "\t")
}

message("UK Biobank outcome analysis complete")
