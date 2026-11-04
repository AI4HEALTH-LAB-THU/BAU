library(data.table)
library(sandwich)
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
root <- normalizePath(file.path(dirname(script), ".."))
source(file.path(root, "R/molecular.R"))
source(file.path(root, "R/survival.R"))

set.seed(31)
n <- 400L
meta <- data.frame(
  BAG = rnorm(n), BAU = rnorm(n), age = rnorm(n),
  centre = factor(rep(1:10, each = 40))
)
design <- model.matrix(~ BAG + BAU + age + centre, meta)
y <- cbind(
  protein_a = 0.3 * meta$BAG - 0.6 * meta$BAU + rnorm(n),
  protein_b = 0.8 * meta$BAU + rnorm(n)
)
y[seq(3, n, 11), 2] <- NA_real_
result <- fit_molecular_split(y, seq_len(n), design, meta$centre)
for (j in seq_len(ncol(y))) {
  observed <- is.finite(y[, j])
  fit <- lm(y[observed, j] ~ design[observed, ] - 1)
  vc <- vcovCL(fit, cluster = meta$centre[observed], type = "HC1", cadjust = TRUE)
  target <- match(c("BAG", "BAU"), colnames(design))
  expected <- unname(coef(fit)[target])
  actual <- unlist(result[j, .(beta_BAG, beta_BAU)], use.names = FALSE)
  stopifnot(max(abs(expected - actual)) < 1e-10)
  stopifnot(abs(result$se_BAG[j] - sqrt(vc[target[1], target[1]])) < 1e-10)
  stopifnot(abs(result$se_BAU[j] - sqrt(vc[target[2], target[2]])) < 1e-10)
  difference_se <- sqrt(vc[target[1], target[1]] + vc[target[2], target[2]] -
    2 * vc[target[1], target[2]])
  stopifnot(abs(result$se_difference[j] - difference_se) < 1e-10)
}
single <- fit_single_term_split(y, seq_len(n), design, meta$centre, "BAU")
stopifnot(max(abs(single$se - result$se_BAU)) < 1e-10)

x <- data.table(
  eid = 1:5, `death duration` = c(400, 500, 600, 800, 900),
  `death diagnose` = c(1, 0, 0, 1, 0),
  h7_multimorbidity_duration = c(700, 500, -1, 300, 700),
  h7_multimorbidity_diagnose = c(0, 1, 1, 1, 0),
  baseline_h7_domain_count = c(0, 1, 2, 0, 1),
  `stroke duration` = c(700, 500, -1, 300, 700),
  `stroke diagnose` = c(0, 1, 1, 1, 0),
  centre = factor(c(1, 2, 3, 1, 2)), education = factor("college"),
  income = factor("middle"), employment = factor("paid")
)
for (outcome in c("major disease multimorbidity", "stroke")) {
  risk <- make_risk_set(x, outcome)
  stopifnot(identical(risk$eid, c(1L, 2L, 4L, 5L)))
  stopifnot(identical(risk$time_days, c(400, 500, 300, 700)))
  stopifnot(identical(risk$event, c(0L, 1L, 1L, 0L)))
  stopifnot(identical(risk$competing_event, c(1L, 0L, 0L, 0L)))
}
stopifnot(identical(x$`stroke duration`, c(700, 500, -1, 300, 700)))
message("Statistical checks passed")

# Change is scaled in the paired cohort before outcome-specific exclusions.
pairs <- data.table(
  baseline_bag = c(-2, -1, 0, 1, 2),
  repeat_bag = c(-1, 2, 2, 0, 6),
  baseline_bau = c(0, 1, -1, 2, 3),
  repeat_bau = c(1, 1, 2, 1, 5),
  visit_interval_years = c(3, 4, 6, 4, 5)
)
changes <- prepare_longitudinal_change(pairs)
expected_bag <- c(1, 3, 2, -1, 4)
expected_bau <- c(1, 0, 3, -1, 2)
stopifnot(isTRUE(all.equal(
  changes$change_bag,
  (expected_bag - mean(expected_bag)) / sd(expected_bag)
)))
stopifnot(isTRUE(all.equal(
  changes$change_bau,
  (expected_bau - mean(expected_bau)) / sd(expected_bau)
)))
stopifnot(!"change_bag" %in% names(pairs))
message("Longitudinal change checks passed")


# Reference transforms are learned in discovery and shared with held-out centres.
dataset <- list(
  y = matrix(c(1, 3, 5, 7, 100, 200, 2, 4, 6, 8, 300, 400), ncol = 2),
  meta = data.table(Baseline_unhealthy = c(0L, 0L, 1L, 0L, 0L, 0L)),
  features = c("protein_a", "protein_b"),
  feature_missingness = c(0, 0)
)
reference <- transform_molecular_targets(dataset, 1:4, tail_fraction = 0.25)
dataset$y[5:6, ] <- -1000
changed_validation <- transform_molecular_targets(dataset, 1:4, tail_fraction = 0.25)
stopifnot(identical(reference$parameters, changed_validation$parameters))
stopifnot(identical(reference$y[1:4, ], changed_validation$y[1:4, ]))
stopifnot(all(reference$y[5:6, ] >= changed_validation$y[5:6, ]))
message("Discovery reference checks passed")

# Replication uses discovery selection, direction and its selected test family.
discovery <- data.table(feature = letters[1:4], cluster_df = 19L)
validation <- data.table(feature = letters[1:4], cluster_df = 9L)
for (term in c("BAG", "BAU", "difference")) {
  discovery[, (paste0("beta_", term)) := 1]
  discovery[, (paste0("p_", term)) := c(0.001, 0.001, 0.5, 0.001)]
  validation[, (paste0("beta_", term)) := c(5, -5, 100, 0.05)]
  validation[, (paste0("se_", term)) := 1]
}
replication <- combine_molecular_splits(discovery, validation)
expected_q <- p.adjust(pt(c(5, -5, 0.05), 9, lower.tail = FALSE), "BH")
for (term in c("BAG", "BAU", "difference")) {
  stopifnot(identical(replication[[paste0("replicated_", term)]], c(TRUE, FALSE, FALSE, FALSE)))
  observed_q <- replication[[paste0("validation_directional_q_", term)]]
  stopifnot(is.na(observed_q[3]), max(abs(observed_q[c(1, 2, 4)] - expected_q)) < 1e-14)
}
message("Directional replication checks passed")

landmark <- data.table(
  eid = 1:7, visit_interval_years = 1,
  `death duration` = c(900, 400, 900, 900, 800, 1000, 1000),
  `death diagnose` = c(0L, 1L, 0L, 0L, 0L, 0L, 0L),
  `stroke duration` = c(700, 500, 200, 365.25, 700, 700, 700),
  `stroke diagnose` = c(1L, 0L, 1L, 1L, 0L, 1L, 1L),
  h7_multimorbidity_duration = c(700, 500, 200, 365.25, 700, 700, 700),
  h7_multimorbidity_diagnose = c(1L, 0L, 1L, 1L, 0L, 1L, 1L),
  baseline_h7_domain_count = c(0L, 1L, 0L, 0L, 2L, 0L, 0L),
  centre = factor(1L), income = factor("middle"),
  education = factor(c(rep("college", 5), "unknown", "college")),
  employment = factor(c(rep("paid", 6), "unknown"))
)
disease_risk <- make_longitudinal_risk_set(landmark, "stroke", "disease")
stopifnot(identical(disease_risk$eid, c(1L, 2L, 5L)))
stopifnot(identical(disease_risk$event, c(1L, 0L, 0L)))
stopifnot(identical(disease_risk$competing_death, c(0L, 1L, 0L)))
stopifnot(max(abs(disease_risk$time_days - c(334.75, 34.75, 334.75))) < 1e-12)
multimorbidity_risk <- make_longitudinal_risk_set(
  landmark, "major disease multimorbidity", "multimorbidity"
)
stopifnot(identical(multimorbidity_risk$eid, c(1L, 2L)))
stopifnot(!"landmark_days" %in% names(landmark))
message("Reassessment landmark checks passed")


# The primary figures use separate 12-test omics and six-test routine-blood families.
outcome_expressions <- parse(file.path(root, "analysis/outcomes/ukb_associations.R"))
function_expression <- Filter(function(expression) {
  is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
    identical(expression[[2L]], as.name("add_primary_outcome_fdr"))
}, as.list(outcome_expressions))
stopifnot(length(function_expression) == 1L)
eval(function_expression[[1L]])
primary <- c("all-cause mortality", "loss of healthspan", "major disease multimorbidity")
fdr_input <- CJ(
  modality = c("nmr", "protein", "blood"),
  outcome = primary, exposure = c("BAG", "BAU")
)
fdr_input[, `:=`(model = "ngboost", p_value = seq_len(.N) / 1000)]
fdr_input <- rbind(fdr_input, data.table(
  modality = "nmr", outcome = "stroke", exposure = c("BAG", "BAU"),
  model = "ngboost", p_value = c(1e-8, 0.9)
))
fdr_result <- add_primary_outcome_fdr(fdr_input, primary)
omics <- which(fdr_result$modality %chin% c("nmr", "protein") &
  fdr_result$outcome %chin% primary)
blood <- which(fdr_result$modality == "blood")
stopifnot(isTRUE(all.equal(
  fdr_result$q_primary_omics_12[omics],
  p.adjust(fdr_input$p_value[omics], "BH", n = 12)
)))
stopifnot(isTRUE(all.equal(
  fdr_result$q_primary_blood_6[blood],
  p.adjust(fdr_input$p_value[blood], "BH", n = 6)
)))
stopifnot(
  all(is.na(fdr_result[outcome == "stroke", q_primary_omics_12])),
  all(is.na(fdr_result[modality != "blood", q_primary_blood_6])),
  identical(fdr_result$p_value, fdr_input$p_value),
  !"q_primary_omics_12" %in% names(fdr_input)
)
subset_result <- add_primary_outcome_fdr(fdr_input[modality == "nmr"], primary)
stopifnot(all(is.na(subset_result$q_primary_omics_12)))
alternative_input <- copy(fdr_input)[, model := "lightgbmlss"]
alternative_result <- add_primary_outcome_fdr(alternative_input, primary)
stopifnot(
  all(is.na(alternative_result$q_primary_omics_12)),
  all(is.na(alternative_result$q_primary_blood_6))
)
message("Primary-outcome FDR-family checks passed")
