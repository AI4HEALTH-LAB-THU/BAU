#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

full_args <- commandArgs(trailingOnly = FALSE)
script_arg <- full_args[grep("^--file=", full_args)][1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
default_root <- normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
project_root <- Sys.getenv("BAU_CODE_ROOT", default_root)
results_root <- Sys.getenv("BAU_RESULTS_ROOT", "/home/dnanexus/bau-work/results")
source_dir <- file.path(results_root, "phenotypes/matched_examples")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

prediction_path <- file.path(
  results_root, "models/nmr/ngboost",
  "all_participant_predictions.csv"
)
expected_n <- 272672L
dat <- fread(
  prediction_path,
  select = c(
    "eid_int", "Age", "Sex_male", "BA_Calibrated", "Sigma_Calibrated",
    "AgeGap_Calibrated", "BAG", "BAU", "Z_AgeGap_Calibrated",
    "Z_Uncertainty_Calibrated", "OuterFold", "AssessmentCenter",
    "SplitRole", "Model"
  ),
  showProgress = FALSE
)
if (nrow(dat) != expected_n || uniqueN(dat$eid_int) != expected_n) {
  stop("Expected 272672 unique NMR participants", call. = FALSE)
}
if (!identical(dat$BAG, dat$Z_AgeGap_Calibrated) ||
  !identical(dat$BAU, dat$Z_Uncertainty_Calibrated)) {
  stop("BAG and BAU aliases differ.", call. = FALSE)
}
required_numeric <- c(
  "Age", "BA_Calibrated", "Sigma_Calibrated", "AgeGap_Calibrated", "BAG", "BAU"
)
if (!all(vapply(dat[, ..required_numeric], function(x) all(is.finite(x)), logical(1))) ||
  !all(dat$Sigma_Calibrated > 0)) {
  stop("Required figure variables contain invalid values.", call. = FALSE)
}
if (!identical(sort(unique(dat$OuterFold)), 1:5) ||
  !identical(unique(dat$SplitRole), "outer_test") ||
  !identical(unique(dat$Model), "ngboost")) {
  stop("Unexpected prediction metadata", call. = FALSE)
}
if (!all(dat[, uniqueN(OuterFold), by = AssessmentCenter]$V1 == 1L)) {
  stop("Assessment centres cross outer-fold boundaries.", call. = FALSE)
}

quantile_type <- 7L
bau_cut_low <- as.numeric(quantile(dat$BAU, 0.05, type = quantile_type))
bau_cut_high <- as.numeric(quantile(dat$BAU, 0.95, type = quantile_type))

pair_band_levels <- c("Younger age gap", "Older age gap")
uncertainty_levels <- c("Lower uncertainty", "Higher uncertainty")
dat[, pair_band := fifelse(
  BAG >= -0.75 & BAG <= -0.25,
  pair_band_levels[[1]],
  fifelse(BAG >= 0.25 & BAG <= 0.75, pair_band_levels[[2]], NA_character_)
)]
dat[, uncertainty_class := fifelse(
  BAU <= bau_cut_low,
  uncertainty_levels[[1]],
  fifelse(BAU >= bau_cut_high, uncertainty_levels[[2]], NA_character_)
)]

upper_stratum_base <- dat[, .(
  stratum_n = .N,
  median_calibrated_age_gap_years = median(AgeGap_Calibrated)
), by = .(Age, Sex_male)]
upper_cell_counts <- dat[
  !is.na(pair_band) & !is.na(uncertainty_class),
  .N,
  by = .(Age, Sex_male, pair_band, uncertainty_class)
]
upper_strata <- dcast(
  upper_cell_counts,
  Age + Sex_male ~ pair_band + uncertainty_class,
  value.var = "N",
  fill = 0L
)
required_upper_cells <- c(
  "Younger age gap_Lower uncertainty",
  "Younger age gap_Higher uncertainty",
  "Older age gap_Lower uncertainty",
  "Older age gap_Higher uncertainty"
)
if (!all(required_upper_cells %in% names(upper_strata))) {
  stop("Missing BAG/BAU groups for the upper-row examples", call. = FALSE)
}
setnames(
  upper_strata,
  required_upper_cells,
  c("younger_low", "younger_high", "older_low", "older_high")
)
upper_strata <- merge(upper_stratum_base, upper_strata, by = c("Age", "Sex_male"))
upper_strata <- upper_strata[
  abs(median_calibrated_age_gap_years) <= 0.10 &
    younger_low > 0 & younger_high > 0 & older_low > 0 & older_high > 0
]
upper_strata[, absolute_median_calibrated_age_gap_years := abs(
  median_calibrated_age_gap_years
)]
setorder(
  upper_strata,
  -stratum_n,
  absolute_median_calibrated_age_gap_years,
  Age,
  Sex_male
)
if (nrow(upper_strata) == 0L) {
  stop("Upper-row matching found no eligible age-sex group", call. = FALSE)
}
upper_selected_stratum <- upper_strata[1]
upper_age <- upper_selected_stratum$Age[[1]]
upper_sex_male <- upper_selected_stratum$Sex_male[[1]]
upper_sex_label <- if (upper_sex_male == 1) "men" else "women"
upper_stratum_population <- dat[Age == upper_age & Sex_male == upper_sex_male]

select_upper_pair <- function(pair_band_value) {
  low <- copy(upper_stratum_population[
    pair_band == pair_band_value & uncertainty_class == uncertainty_levels[[1]]
  ])
  high <- copy(upper_stratum_population[
    pair_band == pair_band_value & uncertainty_class == uncertainty_levels[[2]]
  ])
  if (nrow(low) == 0L || nrow(high) == 0L) {
    stop("Empty uncertainty group in upper-row matching", call. = FALSE)
  }
  low[, join_key := 1L]
  high[, join_key := 1L]
  pairs <- merge(
    low,
    high,
    by = "join_key",
    suffixes = c("_low", "_high"),
    allow.cartesian = TRUE
  )
  pairs[, point_prediction_difference_years := abs(
    BA_Calibrated_low - BA_Calibrated_high
  )]
  pairs[, uncertainty_target_distance :=
    abs(BAU_low - bau_cut_low) + abs(BAU_high - bau_cut_high)]
  setorder(
    pairs,
    point_prediction_difference_years,
    uncertainty_target_distance,
    eid_int_low,
    eid_int_high
  )
  winner <- pairs[1]
  selected <- rbindlist(list(
    low[eid_int == winner$eid_int_low[[1]]],
    high[eid_int == winner$eid_int_high[[1]]]
  ), use.names = TRUE)
  selected[, `:=`(
    point_prediction_difference_years = winner$point_prediction_difference_years[[1]],
    predictive_sd_ratio =
      winner$Sigma_Calibrated_high[[1]] / winner$Sigma_Calibrated_low[[1]],
    join_key = NULL
  )]
  selected
}

upper_examples <- rbindlist(
  lapply(pair_band_levels, select_upper_pair),
  use.names = TRUE
)
upper_examples[, pair_band := factor(pair_band, levels = pair_band_levels)]
upper_examples[, uncertainty_class := factor(
  uncertainty_class,
  levels = uncertainty_levels
)]
setorder(upper_examples, pair_band, uncertainty_class)
upper_examples[, participant_label := paste("Participant", seq_len(.N))]
upper_examples[, anonymous_example := paste0("P", seq_len(.N))]
if (nrow(upper_examples) != 4L || uniqueN(upper_examples$eid_int) != 4L ||
  uniqueN(upper_examples$Age) != 1L || uniqueN(upper_examples$Sex_male) != 1L ||
  any(abs(upper_examples$BAG) >= 1)) {
  stop("Expected four unique upper-row participants of the same age and sex, with abs(BAG) < 1", call. = FALSE)
}
upper_uncertainty_check <- dcast(
  upper_examples[, .(pair_band, uncertainty_class, Sigma_Calibrated)],
  pair_band ~ uncertainty_class,
  value.var = "Sigma_Calibrated"
)
if (!all(upper_uncertainty_check[[uncertainty_levels[[2]]]] >
  upper_uncertainty_check[[uncertainty_levels[[1]]]])) {
  stop("Upper higher-uncertainty distributions are not wider.", call. = FALSE)
}

anonymous_upper_source <- upper_examples[, .(
  anonymous_example,
  participant_label,
  chronological_age_years = Age,
  sex = upper_sex_label,
  age_gap_group = as.character(pair_band),
  uncertainty_group = as.character(uncertainty_class),
  biological_age_prediction_years = BA_Calibrated,
  predicted_age_sd_years = Sigma_Calibrated,
  calibrated_age_gap_years = AgeGap_Calibrated,
  biological_age_gap_sd = BAG,
  biological_age_uncertainty_sd = BAU,
  matched_point_prediction_difference_years = point_prediction_difference_years,
  higher_to_lower_predictive_sd_ratio = predictive_sd_ratio
)]
fwrite(anonymous_upper_source, file.path(source_dir, "non_extreme_examples.tsv"), sep = "\t")
