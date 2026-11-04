#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(splines)
})

args <- commandArgs(trailingOnly = TRUE)
value_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[[1]]) else default
}
args_all <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
default_root <- normalizePath(file.path(dirname(script_file), "../.."), mustWork = TRUE)
code_root <- normalizePath(value_arg(
  "code-root", Sys.getenv("BAU_CODE_ROOT", default_root)
), mustWork = FALSE)
ukb_root <- normalizePath(value_arg(
  "ukb-root", Sys.getenv("UKB_DATA_ROOT", "/mnt/project/UKB")
), mustWork = FALSE)
data_root <- Sys.getenv("BAU_DATA_ROOT", "/home/dnanexus/bau-work/data")
results_root <- Sys.getenv("BAU_RESULTS_ROOT", "/home/dnanexus/bau-work/results")
source(file.path(code_root, "R/molecular.R"))

result_dir <- file.path(results_root, "molecular/physiology")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
definitions <- fread(file.path(
  code_root, "config/physiology.tsv"
))

source_features <- unique(unlist(strsplit(
  definitions$source_feature, " | ",
  fixed = TRUE
)))
clinical <- as.data.table(read_feather(
  file.path(data_root, "ukb/modeling/routine_blood.feather"),
  col_select = c(
    "eid", "Date_birth", "Date_attend",
    source_features
  )
))
clinical[, eid := as.integer(eid)]

covariate_names <- c(
  "Current smoker", "Daily alcohol intake",
  "Daily moderate to vigorous physical activity", "Education", "Income", "Employment"
)
covariates <- as.data.table(read_feather(
  file.path(data_root, "ukb/support/covariates.feather"),
  col_select = c("eid", covariate_names)
))
covariates[, eid := as.integer(eid)]
setnames(covariates, covariate_names, c(
  "smoker", "daily_alcohol", "physical_activity", "education", "income", "employment"
))
ethnicity <- read_ethnicity(ukb_root)

make_outcomes <- function(linked, discovery_rows, healthy_rows) {
  y <- matrix(NA_real_, nrow(linked), nrow(definitions))
  colnames(y) <- definitions$outcome_id
  records <- vector("list", nrow(definitions))
  for (j in seq_len(nrow(definitions))) {
    definition <- definitions[j]
    raw <- if (definition$derivation == "ratio_positive_denominator") {
      numerator <- as.numeric(linked[["Neutrophill count"]])
      denominator <- as.numeric(linked[["Lymphocyte count"]])
      fifelse(
        is.finite(numerator) & is.finite(denominator) & denominator > 0,
        numerator / denominator, NA_real_
      )
    } else {
      as.numeric(linked[[definition$source_feature]])
    }
    limits <- quantile(
      raw[discovery_rows], c(0.005, 0.995),
      na.rm = TRUE, type = 8,
      names = FALSE
    )
    transformed <- pmax(limits[[1]], pmin(limits[[2]], raw))
    if (definition$transform == "log") transformed <- log(transformed)
    reference_mean <- mean(transformed[healthy_rows], na.rm = TRUE)
    reference_sd <- sd(transformed[healthy_rows], na.rm = TRUE)
    y[, j] <- (transformed - reference_mean) / reference_sd
    records[[j]] <- data.table(
      outcome_id = definition$outcome_id,
      winsor_lower = limits[[1]],
      winsor_upper = limits[[2]],
      healthy_reference_n = sum(is.finite(transformed[healthy_rows])),
      healthy_reference_mean = reference_mean,
      healthy_reference_sd = reference_sd,
      discovery_missing_fraction = mean(!is.finite(raw[discovery_rows]))
    )
  }
  list(y = y, parameters = rbindlist(records))
}

make_design <- function(linked, discovery_rows) {
  attendance <- as.numeric(as.IDate(linked$Date_attend))
  day <- as.integer(format(as.IDate(linked$Date_attend), "%j"))
  age_fit <- ns(linked$Age_exact[discovery_rows], df = 4, Boundary.knots = c(40, 70))
  age_basis <- predict(age_fit, newx = linked$Age_exact)
  colnames(age_basis) <- paste0("age_spline_", seq_len(ncol(age_basis)))
  activity <- as.numeric(linked$physical_activity)
  activity_z <- (
    activity - mean(activity[discovery_rows])
  ) / sd(activity[discovery_rows])
  design <- data.frame(
    BAG = linked$BAG,
    BAU = linked$BAU,
    Sex_male = linked$Sex_male,
    ethnicity_group = linked$ethnicity_group,
    smoker = as.numeric(as.character(linked$smoker)),
    daily_alcohol = as.numeric(as.character(linked$daily_alcohol)),
    physical_activity_z = activity_z,
    education = factor(as.character(linked$education)),
    income = factor(as.character(linked$income)),
    employment = factor(as.character(linked$employment)),
    AssessmentCenter = linked$AssessmentCenter,
    attendance_z = (
      attendance - mean(attendance[discovery_rows])
    ) / sd(attendance[discovery_rows]),
    season_sin = sin(2 * pi * day / 365.2425),
    season_cos = cos(2 * pi * day / 365.2425),
    model_input_missing_fraction = linked$model_input_missing_fraction,
    age_basis,
    check.names = FALSE
  )
  formula <- as.formula(paste(
    "~ BAG + BAU + Sex_male + ethnicity_group + smoker + daily_alcohol +",
    "physical_activity_z + education + income + employment + AssessmentCenter +",
    "attendance_z + season_sin + season_cos + model_input_missing_fraction +",
    paste(colnames(age_basis), collapse = " + ")
  ))
  model.matrix(formula, design)
}

all_associations <- list()
all_transformations <- list()
all_cohorts <- list()

for (modality in c("nmr", "protein")) {
  predictions <- read_molecular_predictions(code_root, modality)
  model_features <- read_feature_order(file.path(
    code_root, "config/features",
    paste0(modality, ".txt")
  ))
  model_data <- as.data.table(read_feather(
    file.path(data_root, "ukb/modeling", paste0(modality, ".feather")),
    col_select = c("eid", model_features)
  ))
  model_data[, eid := as.integer(eid)]
  model_data[, model_input_missing_fraction := rowMeans(is.na(.SD)), .SDcols = model_features]
  model_data <- model_data[, .(eid, model_input_missing_fraction)]

  linked <- merge(predictions, clinical, by = "eid", all = FALSE, sort = FALSE)
  linked <- merge(linked, model_data, by = "eid", all.x = TRUE, sort = FALSE)
  linked <- merge(linked, covariates, by = "eid", all.x = TRUE, sort = FALSE)
  linked <- merge(linked, ethnicity, by = "eid", all.x = TRUE, sort = FALSE)
  linked[, `:=`(
    Age_exact = as.numeric(as.IDate(Date_attend) - as.IDate(Date_birth)) / 365.2425,
    Sex_male = as.integer(Sex_male),
    AssessmentCenter = factor(as.character(AssessmentCenter)),
    ethnicity_group = factor(
      collapse_ethnicity(as.integer(ethnicity_code)),
      levels = c("White", "Mixed", "Asian", "Black", "Chinese", "Other", "Unknown")
    )
  )]
  complete_columns <- c(
    "BAG", "BAU", "Age_exact", "Sex_male", "AssessmentCenter",
    "smoker", "daily_alcohol", "physical_activity", "education", "income",
    "employment", "model_input_missing_fraction", "Date_attend"
  )
  linked <- linked[complete.cases(linked[, ..complete_columns])]
  discovery_rows <- which(linked$OuterFold %in% 1:3)
  validation_rows <- which(linked$OuterFold %in% 4:5)
  healthy_rows <- discovery_rows[linked$Baseline_unhealthy[discovery_rows] == 0L]

  transformed <- make_outcomes(linked, discovery_rows, healthy_rows)
  design <- make_design(linked, discovery_rows)
  discovery <- fit_molecular_split(
    transformed$y, discovery_rows, design, linked$AssessmentCenter
  )
  validation <- fit_molecular_split(
    transformed$y, validation_rows, design, linked$AssessmentCenter
  )
  association <- combine_molecular_splits(discovery, validation)
  setnames(association, "feature", "outcome_id")
  association <- merge(definitions, association, by = "outcome_id", sort = FALSE)
  setorder(association, outcome_order)
  association[, modality := modality]
  transformed$parameters[, modality := modality]
  all_associations[[modality]] <- association
  all_transformations[[modality]] <- transformed$parameters
  all_cohorts[[modality]] <- data.table(
    modality = modality,
    linked_n = nrow(linked),
    discovery_n = length(discovery_rows),
    validation_n = length(validation_rows),
    healthy_discovery_n = length(healthy_rows),
    discovery_centres = uniqueN(linked$AssessmentCenter[discovery_rows]),
    validation_centres = uniqueN(linked$AssessmentCenter[validation_rows])
  )
}

fwrite(
  rbindlist(all_associations, fill = TRUE),
  file.path(result_dir, "physiological_associations.tsv"),
  sep = "\t", na = "NA"
)
fwrite(
  rbindlist(all_transformations, fill = TRUE),
  file.path(result_dir, "outcome_transformations.tsv"),
  sep = "\t", na = "NA"
)
fwrite(
  rbindlist(all_cohorts),
  file.path(result_dir, "cohort_summary.tsv"),
  sep = "\t", na = "NA"
)
