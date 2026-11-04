#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(patchwork)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])
default_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
root <- Sys.getenv("BAU_CODE_ROOT", default_root)
results_root <- Sys.getenv("BAU_RESULTS_ROOT", "/home/dnanexus/bau-work/results")
source(file.path(root, "R/theme.R"))
theme_set(paper_theme())
input <- file.path(results_root, "outcomes")
output <- file.path(results_root, "figures")
read_result <- function(path) {
  x <- fread(file.path(input, path))
  if ("outcome" %in% names(x)) x[, outcome := factor(outcome, levels = outcome_order)]
  x
}
save <- function(plot, name, width = 180, height = 90) {
  save_panel(plot, file.path(output, name), width, height)
}
primary <- c("all-cause mortality", "loss of healthspan", "major disease multimorbidity")
outcome_order <- c(
  primary, "coronary heart disease", "stroke", "chronic obstructive pulmonary disease",
  "liver diseases", "diabetes mellitus", "renal failure", "dementia",
  "digestive system cancer", "respiratory system cancer"
)
phenotype_colours <- c(
  BAG = bau_colours[["bag"]], BAU = bau_colours[["bau"]],
  `BAG change` = bau_colours[["bag"]], `BAU change` = bau_colours[["bau"]]
)
states <- c(LL = "#4B657A", HL = "#E69F00", LH = "#009E73", HH = "#CC79A7")
forest <- function(x, rows = "outcome", facets = ~modality) {
  x <- copy(x)
  row_order <- if (rows == "outcome") outcome_order[outcome_order %in% x$outcome] else unique(x[[rows]])
  x[, row_label := factor(get(rows), levels = rev(row_order))]
  facets_plot <- if (length(all.vars(facets)) == 1L) {
    facet_wrap(facets, scales = "free_y")
  } else {
    facet_grid(facets, scales = "free_y")
  }
  ggplot(x, aes(hr, row_label, colour = exposure, shape = exposure, group = exposure)) +
    geom_vline(xintercept = 1, linewidth = 0.3, colour = "grey70") +
    geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper),
      orientation = "y",
      width = 0.3, position = position_dodge(0.5), linewidth = 0.4
    ) +
    geom_point(position = position_dodge(0.5), size = 1.7, fill = "white") +
    scale_shape_manual(values = c(BAG = 21, BAU = 23, `BAG change` = 21, `BAU change` = 23)) +
    scale_colour_manual(values = phenotype_colours) +
    scale_x_log10() +
    facets_plot +
    labs(x = "Hazard ratio", y = NULL) +
    theme(legend.position = "bottom")
}

curves <- read_result("ngboost/cumulative_event_curves.tsv")
at_risk <- read_result("ngboost/numbers_at_risk.tsv")
effects <- read_result("ngboost/joint_cox_estimates.tsv")
curve_plot <- function(x) {
  annotation <- effects[modality %in% unique(x$modality) & outcome %in% unique(x$outcome)]
  annotation <- annotation[, .(label = paste(sprintf(
    "%s  %.2f (%.2f-%.2f)",
    exposure, hr, ci_lower, ci_upper
  ), collapse = "\n")), by = outcome]
  ggplot(x, aes(time_years, estimate, colour = state)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = state), alpha = 0.1, colour = NA) +
    geom_step(linewidth = 0.45) +
    geom_text(
      data = annotation, aes(x = 0.4, y = Inf, label = label),
      inherit.aes = FALSE, hjust = 0, vjust = 1.2, size = 2.1
    ) +
    scale_colour_manual(values = states) +
    scale_fill_manual(values = states, guide = "none") +
    scale_x_continuous(breaks = c(0, 5, 10, 15), limits = c(0, 15)) +
    facet_wrap(~outcome, scales = "free_y", nrow = 1) +
    labs(x = "Follow-up (years)", y = "Cumulative probability") +
    theme(legend.position = "bottom")
}
for (modality_name in c("nmr", "protein")) {
  panel <- if (modality_name == "nmr") "a" else "c"
  risk_table <- ggplot(
    at_risk[modality == modality_name & outcome %in% primary],
    aes(time_years, state, label = n_risk, colour = state)
  ) +
    geom_text(size = 2) +
    scale_colour_manual(values = states, guide = "none") +
    scale_x_continuous(breaks = c(0, 5, 10, 15), limits = c(0, 15)) +
    facet_wrap(~outcome, nrow = 1) +
    labs(x = NULL, y = "At risk") +
    theme(
      axis.line = element_blank(), axis.ticks = element_blank(),
      strip.text = element_blank()
    )
  save(curve_plot(curves[modality == modality_name & outcome %in% primary]) / risk_table +
    plot_layout(heights = c(3, 1)), paste0("figure_2", panel, "_", modality_name), height = 95)
}
save(forest(effects[modality == "nmr" & outcome %in% primary]), "figure_2b", 110, 65)
save(forest(effects[modality == "protein" & outcome %in% primary]), "figure_2d", 110, 65)
save(curve_plot(curves[modality == "blood" & outcome %in% primary]), "figure_2g", height = 80)
for (modality_name in c("nmr", "protein")) {
  p <- curve_plot(curves[modality == modality_name & !outcome %in% primary]) +
    facet_wrap(~outcome, ncol = 3, scales = "free_y")
  panel <- if (modality_name == "nmr") "a" else "b"
  save(p, paste0("supplementary_figure_3", panel, "_", modality_name), height = 150)
}
save(
  forest(effects[modality == "blood" & !outcome %in% primary]),
  "supplementary_figure_3f_routine_blood", 110, 115
)

increment <- read_result("discrimination/incremental_c_index.tsv")
increment[, model := factor(model,
  levels = c("clinical_bag", "clinical_bau", "clinical_bag_bau"),
  labels = c("+ BAG", "+ BAU", "+ BAG and BAU")
)]
increment_plot <- function(x) {
  ggplot(x, aes(estimate_x1000, outcome, colour = model)) +
    geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
    geom_errorbar(aes(xmin = ci_lower_x1000, xmax = ci_upper_x1000),
      orientation = "y",
      width = 0.25, position = position_dodge(0.6), linewidth = 0.4
    ) +
    geom_point(position = position_dodge(0.6), size = 1.6) +
    facet_wrap(~modality) +
    scale_colour_manual(values = c(bau_colours[["bag"]], bau_colours[["bau"]], "#714A70")) +
    labs(x = "Increment in 10-year C-index (x 1,000)", y = NULL) +
    theme(legend.position = "bottom")
}
absolute <- read_result("discrimination/absolute_c_index.tsv")
absolute[, model := factor(model,
  levels = c("clinical", "clinical_bag", "clinical_bau", "clinical_bag_bau"),
  labels = c("Clinical", "+ BAG", "+ BAU", "+ BAG and BAU")
)]
absolute_plot <- function(x) {
  ggplot(x, aes(estimate, outcome, colour = model)) +
    geom_point(position = position_dodge(0.6), size = 1.4) +
    geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper),
      orientation = "y",
      position = position_dodge(0.6), width = 0.25
    ) +
    facet_wrap(~modality) +
    scale_colour_manual(values = c(
      "#84959B", bau_colours[["bag"]],
      bau_colours[["bau"]], "#714A70"
    )) +
    labs(x = "10-year Uno C-index", y = NULL) +
    theme(legend.position = "bottom")
}
save(absolute_plot(absolute[outcome %in% primary]) |
  increment_plot(increment[outcome %in% primary]), "figure_2e", width = 240, height = 80)
for (modality_name in c("nmr", "protein")) {
  panel <- if (modality_name == "nmr") "c" else "d"
  save(
    absolute_plot(absolute[modality == modality_name & !outcome %in% primary]) |
      increment_plot(increment[modality == modality_name & !outcome %in% primary]),
    paste0("supplementary_figure_3", panel),
    width = 190, height = 115
  )
}

change <- read_result("longitudinal/simple_change_estimates.tsv")
change[, modality := "Repeat NMR"]
save(forest(change[outcome %in% primary]), "figure_2f", 110, 65)
save(forest(change[!outcome %in% primary]), "supplementary_figure_3e", 110, 115)
external <- read_result("external/primary_mortality_estimates.tsv")
setnames(external, "term", "exposure")
external[, modality := "External cohorts"]
save(forest(external, rows = "cohort"), "figure_2h", 110, 65)

for (spec in c("extended_covariates", "one_year_landmark")) {
  x <- read_result(paste0("sensitivity/", spec, ".tsv"))
  panel <- if (spec == "extended_covariates") "a" else "b"
  save(forest(x), paste0("supplementary_figure_4", panel), height = 145)
}
period <- read_result("sensitivity/follow_up_periods.tsv")
save(forest(period, facets = modality ~ period), "supplementary_figure_4c", 240, 160)
competing <- read_result("sensitivity/competing_risk.tsv")
competing[, method := factor(analysis,
  levels = c("cause_specific_cox", "fine_gray"),
  labels = c("Cause-specific Cox", "Fine–Gray")
)]
competing[, fill_key := paste(exposure, method, sep = " / ")]
competing[, row_label := factor(outcome, levels = rev(outcome_order[outcome_order %in% outcome]))]
competing_plot <- ggplot(competing, aes(hr, row_label,
  colour = exposure,
  group = interaction(exposure, method)
)) +
  geom_vline(xintercept = 1, colour = "grey70", linewidth = 0.3) +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper),
    orientation = "y",
    position = position_dodge(0.65), width = 0.3, linewidth = 0.4
  ) +
  geom_point(aes(fill = fill_key),
    shape = 21, size = 1.7,
    position = position_dodge(0.65)
  ) +
  scale_fill_manual(values = c(
    `BAG / Cause-specific Cox` = "white", `BAU / Cause-specific Cox` = "white",
    `BAG / Fine–Gray` = bau_colours[["bag"]], `BAU / Fine–Gray` = bau_colours[["bau"]]
  )) +
  scale_colour_manual(values = phenotype_colours) +
  scale_x_log10() +
  facet_wrap(~modality) +
  labs(x = "Hazard ratio", y = NULL, fill = NULL) +
  theme(legend.position = "bottom")
save(competing_plot, "supplementary_figure_4d", 190, 145)
algorithms <- rbindlist(lapply(c("ngboost", "lightgbmlss", "xgboostlss"), function(model_name) {
  read_result(paste0(model_name, "/joint_cox_estimates.tsv"))[
    modality %in% c("nmr", "protein")
  ]
}))
for (model_name in c("xgboostlss", "lightgbmlss")) {
  panel <- if (model_name == "xgboostlss") "e" else "f"
  save(forest(algorithms[model == model_name]), paste0("supplementary_figure_4", panel), 180, 145)
}

subgroups <- read_result("subgroups/subgroup_estimates.tsv")
interactions <- read_result("subgroups/interaction_tests.tsv")
modifiers <- c("age", "sex", "current_smoking", "bmi", "baseline_health", "physical_activity")
for (i in seq_along(modifiers)) {
  modifier_name <- modifiers[[i]]
  estimates <- copy(subgroups[modifier == modifier_name])
  overall <- copy(effects[modality %in% c("nmr", "protein") & outcome %in% primary])
  overall[, level := "overall"]
  estimates <- rbindlist(list(overall, estimates), fill = TRUE)
  level_labels <- c(
    overall = "Overall", `40-49` = "40–49 years", `50-59` = "50–59 years", `60+` = "≥60 years",
    female = "Female", male = "Male", not_current = "Not current", current = "Current",
    `<25` = "BMI <25", `25-<30` = "BMI 25–<30", `30+` = "BMI ≥30",
    no_history = "No baseline disease", disease_history = "≥1 baseline disease",
    `<30_min_day` = "<30 min/day", `30+_min_day` = "≥30 min/day"
  )
  estimates[, level := unname(level_labels[level])]
  tests <- interactions[modifier == modifier_name]
  tests[, annotation := paste(sprintf(
    "%s P interaction = %s", exposure,
    format.pval(p_value, digits = 2, eps = 0.001)
  ), collapse = "\n"),
  by = .(modality, outcome)
  ]
  annotation <- merge(unique(tests[, .(modality, outcome, annotation)]),
    estimates[, .(label_x = min(ci_lower)), by = .(modality, outcome)],
    by = c("modality", "outcome")
  )
  plot <- forest(estimates, rows = "level", facets = outcome ~ modality) +
    geom_text(
      data = annotation, aes(x = label_x, y = Inf, label = annotation),
      inherit.aes = FALSE, hjust = 0, vjust = 1.15, size = 1.9
    ) +
    scale_y_discrete(expand = expansion(add = c(0.6, 1.1)))
  save(plot, paste0("supplementary_figure_5", letters[[i]]), 230, 100)
}
