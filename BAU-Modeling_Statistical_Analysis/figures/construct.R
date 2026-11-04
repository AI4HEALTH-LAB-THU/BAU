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
input <- file.path(results_root, "phenotypes")
output <- file.path(results_root, "figures")
read_result <- function(path) fread(file.path(input, path))
save <- function(plot, name, width = 150, height = 82) {
  save_panel(plot, file.path(output, name), width, height)
}
modality_colours <- c(nmr = "#2F7F78", protein = "#80669D")
phenotype_colours <- c(BAG = bau_colours[["bag"]], BAU = bau_colours[["bau"]])

examples <- read_result("matched_examples/non_extreme_examples.tsv")
densities <- rbindlist(lapply(seq_len(nrow(examples)), function(i) {
  x <- examples[i]
  age <- seq(x$biological_age_prediction_years - 4 * x$predicted_age_sd_years,
    x$biological_age_prediction_years + 4 * x$predicted_age_sd_years,
    length.out = 301
  )
  data.table(age,
    density = dnorm(age, x$biological_age_prediction_years, x$predicted_age_sd_years),
    age_gap_group = x$age_gap_group, uncertainty_group = x$uncertainty_group,
    participant = x$participant_label
  )
}))
densities[, uncertainty_level := fifelse(
  grepl("lower|low", uncertainty_group, ignore.case = TRUE),
  "Lower BAU", "Higher BAU"
)]
upper <- ggplot(densities, aes(age, density, colour = uncertainty_level)) +
  geom_line(aes(linetype = uncertainty_level), linewidth = 0.5) +
  facet_wrap(~age_gap_group) +
  geom_vline(xintercept = unique(examples$chronological_age_years), linetype = 2, linewidth = 0.3) +
  scale_colour_manual(values = c("Lower BAU" = "#D79A89", "Higher BAU" = "#C1644A")) +
  scale_linetype_manual(values = c("Lower BAU" = "dashed", "Higher BAU" = "solid")) +
  labs(x = "Predicted age (years)", y = "Density")
save(upper + theme(legend.position = "bottom"), "figure_1c", 160, 80)

points <- read_result("model_diagnostics/bag_bau_display_sample.tsv")
relationship <- read_result("model_diagnostics/bag_bau_relationship.tsv")
relationship[, label := sprintf("Pearson r = %.3f", pearson_r)]
save(ggplot(points, aes(BAG, BAU)) +
  geom_point(size = 0.45, alpha = 0.18, colour = "#687176") +
  geom_abline(
    data = relationship, aes(intercept = intercept, slope = slope),
    colour = "#B65F52", linewidth = 0.5
  ) +
  geom_text(
    data = relationship, aes(x = -Inf, y = Inf, label = label),
    hjust = -0.1, vjust = 1.2, size = 2.5, inherit.aes = FALSE
  ) +
  facet_wrap(~modality) +
  labs(x = "BAG", y = "BAU"), "figure_1d")
protein_examples <- read_result("molecular/protein/protein_examples.tsv")
example_metadata <- unique(protein_examples[, .(example, age, BAG, BAU, mean_disagreement)])
example_metadata[, label := sprintf(
  "%s\nAge %.1f; BAG %.2f; BAU %.2f; mean %.3f",
  example, age, BAG, BAU, mean_disagreement
)]
example_metadata[, display_row := rev(seq_len(.N))]
protein_examples[, display_row := example_metadata$display_row[match(example, example_metadata$example)]]
protein_examples[, bin := floor(disagreement / 0.05)]
protein_examples[, offset := seq_len(.N) - (.N + 1) / 2,
  by = .(example, bin)
]
protein_examples[, offset := 0.25 * offset / max(max(abs(offset)), 1), by = example]
protein_examples[, colour := fifelse(
  grepl("lower", example, ignore.case = TRUE), "#A4A5A8", "#80669D"
)]
save(ggplot(protein_examples, aes(disagreement, display_row + offset)) +
  geom_point(aes(colour = colour), size = 0.7, alpha = 0.55) +
  scale_colour_identity() +
  geom_point(
    data = example_metadata, aes(x = mean_disagreement, y = display_row),
    shape = 23, size = 2.2, colour = "#39434B", fill = "white", inherit.aes = FALSE
  ) +
  scale_y_continuous(breaks = example_metadata$display_row, labels = example_metadata$label) +
  labs(x = "Single-protein age-signal disagreement", y = NULL), "figure_1e", 183, 105)
for (name in c("protein", "nmr", "organ")) {
  curve <- read_result(paste0("disagreement/", name, "_adjusted_curves.tsv"))
  colour <- c(protein = "#80669D", nmr = "#2F7F78", organ = "#C1644A")[[name]]
  panel <- c(protein = "f", nmr = "g", organ = "h")[[name]]
  label <- c(
    protein = "Protein molecular disagreement", nmr = "NMR molecular disagreement",
    organ = "Organ-age disagreement"
  )[[name]]
  histogram <- read_result(paste0("disagreement/", name, "_histogram.tsv"))
  limits <- range(curve$H)
  association <- ggplot(curve, aes(H, adjusted_mean)) +
    geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = colour, alpha = 0.2) +
    geom_line(colour = colour, linewidth = 0.6) +
    coord_cartesian(xlim = limits) +
    labs(x = NULL, y = "Adjusted mean BAU") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  distribution <- ggplot(histogram, aes(xmin = bin_left, xmax = bin_right, ymin = 0, ymax = n)) +
    geom_rect(fill = colour, alpha = 0.55) +
    coord_cartesian(xlim = limits) +
    labs(x = label, y = NULL) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())
  save(association / distribution + plot_layout(heights = c(5, 1)), paste0("figure_1", panel), 110, 90)
}

population <- read_result("model_diagnostics/age_sex_histogram.tsv")
population[, mirrored_percent := 100 * proportion * fifelse(sex == "female", 1, -1)]
age_plot <- ggplot(population, aes(age_year, mirrored_percent, fill = sex)) +
  geom_col(width = 1) +
  facet_grid(population ~ modality) +
  scale_fill_manual(values = c(female = "#AEB7BC", male = "#536B78")) +
  labs(x = "Age (years)", y = "Participants per age bin (%)")
roles <- read_result("model_diagnostics/fold_roles.tsv")
role_plot <- ggplot(roles[role != "held_centre_test"], aes(factor(outer_fold), n, fill = role)) +
  geom_col(width = 0.65) +
  facet_wrap(~modality) +
  scale_fill_manual(values = c(model_fit = "#59666D", calibration_reference = "#B7BFC3")) +
  labs(x = "Outer fold", y = "Healthy development participants")
held_plot <- ggplot(roles[role == "held_centre_test"], aes(factor(outer_fold), n, fill = modality)) +
  geom_col(width = 0.65) +
  facet_wrap(~modality) +
  scale_fill_manual(values = modality_colours) +
  labs(x = "Outer fold", y = "Held-centre participants")
save(age_plot / (role_plot + held_plot), "supplementary_figure_1b", 183, 150)
density <- read_result("model_diagnostics/age_prediction_density.tsv")
performance <- read_result("model_diagnostics/age_prediction_performance.tsv")
pooled <- performance[is.na(outer_fold)]
pooled[, label := sprintf("Pearson r = %.3f\nMAE = %.2f years", pearson_r, mae_years)]
save(
  ggplot(density, aes(chronological_age, predicted_age, fill = n)) +
    geom_tile(width = 0.5, height = 0.5) +
    geom_abline(slope = 1, intercept = 0, colour = "grey35", linewidth = 0.3) +
    geom_text(
      data = pooled, aes(x = -Inf, y = Inf, label = label),
      inherit.aes = FALSE, hjust = -0.1, vjust = 1.2, size = 2.5
    ) +
    facet_wrap(~modality) +
    scale_fill_gradient(low = "#E6EEF0", high = "#285F8F", trans = "log10") +
    coord_equal() +
    labs(x = "Chronological age (years)", y = "Predicted age (years)"),
  "supplementary_figure_1c"
)
performance <- melt(performance,
  id.vars = c("modality", "outer_fold"),
  measure.vars = c("mae_years", "pearson_r")
)
save(ggplot(performance[!is.na(outer_fold)], aes(outer_fold, value, colour = modality)) +
  geom_hline(
    data = performance[is.na(outer_fold)], aes(yintercept = value),
    linetype = 2, linewidth = 0.35
  ) +
  geom_point(size = 1.5) +
  facet_wrap(~variable, scales = "free_y") +
  scale_colour_manual(values = modality_colours) +
  labs(x = "Outer fold", y = NULL), "supplementary_figure_1d", 110, 70)

persistence <- read_result("repeatability/persistence_metrics.tsv")
paired <- read_result("repeatability/paired_density.tsv")
persistence[, label := sprintf("Spearman r = %.3f", estimate)]
paired_plot <- ggplot(paired) +
  geom_rect(aes(
    xmin = baseline_low, xmax = baseline_high,
    ymin = repeat_low, ymax = repeat_high, fill = count
  )) +
  geom_abline(slope = 1, intercept = 0, colour = "grey50", linewidth = 0.3, linetype = 2) +
  geom_text(
    data = persistence, aes(x = -Inf, y = Inf, label = label),
    hjust = -0.1, vjust = 1.2, size = 2.5, inherit.aes = FALSE
  ) +
  facet_wrap(~phenotype) +
  coord_equal() +
  scale_fill_gradient(low = "#E6EEF0", high = "#285F8F", trans = "log10") +
  labs(x = "Baseline score", y = "Repeat score")
save(paired_plot, "supplementary_figure_2d", 150, 75)

explain <- read_result("explainability/total_oof_r2.tsv")
identity <- read_result("explainability/residual_identity.tsv")
setorder(explain, block_order)
explain[, block := factor(block, levels = rev(unique(block)))]
reference <- explain[block_order == 0]
explained_plot <- ggplot(explain, aes(100 * oof_r2, block, colour = modality)) +
  geom_vline(
    data = reference, aes(xintercept = 100 * oof_r2, colour = modality),
    linetype = 2, linewidth = 0.35
  ) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(xmin = 100 * ci_low, xmax = 100 * ci_high), orientation = "y", width = 0.2) +
  facet_wrap(~modality) +
  scale_colour_manual(values = modality_colours) +
  labs(x = "Total out-of-fold R² (%)", y = NULL) +
  theme(legend.position = "none")
identity_plot <- ggplot(identity, aes(modality, estimate, colour = modality)) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15) +
  scale_colour_manual(values = modality_colours) +
  labs(x = NULL, y = "Spearman r before and after adjustment") +
  theme(legend.position = "none")
save(explained_plot + identity_plot + plot_layout(widths = c(3, 1)), "supplementary_figure_2e", 183, 86)

clocks <- read_result("domain_clock_performance.tsv")
clocks[, domain := factor(domain, levels = rev(domain))]
clock_plot <- function(variable, label) {
  ggplot(clocks, aes(.data[[variable]], domain)) +
    geom_segment(aes(x = 0, xend = .data[[variable]], yend = domain), colour = "#C1644A", linewidth = 0.4) +
    geom_point(colour = "#C1644A", size = 1.7) +
    labs(x = label, y = NULL)
}
save(
  clock_plot("pearson_r", "Pearson r") + clock_plot("mae_years", "MAE (years)"),
  "supplementary_figure_2c", 135, 88
)

agreement <- read_result("algorithm_robustness/bau_rank_agreement.tsv")
algorithms <- c("ngboost", "lightgbmlss", "xgboostlss")
agreement[, algorithm_1 := factor(algorithm_1, levels = algorithms)]
agreement[, algorithm_2 := factor(algorithm_2, levels = rev(algorithms))]
heatmap <- ggplot(agreement, aes(algorithm_1, algorithm_2, fill = spearman_r)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", spearman_r)), size = 2.5) +
  facet_wrap(~modality) +
  scale_fill_gradient(low = "white", high = "#70949D", limits = c(0, 1)) +
  labs(x = NULL, y = NULL, fill = "Spearman r") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save(heatmap, "supplementary_figure_2f", 150, 80)

# Illustrative profiles for the disagreement definitions.
age <- seq(40, 70, by = 0.25)
example_curves <- rbindlist(lapply(c("Lower disagreement", "Higher disagreement"), function(group) {
  peaks <- if (group == "Lower disagreement") rep(55, 4) else c(46, 52, 58, 64)
  rbindlist(lapply(seq_along(peaks), function(j) {
    data.table(group,
      biomarker = paste("Biomarker", j), age,
      log_likelihood = -0.5 * ((age - peaks[j]) / 4)^2
    )
  }))
}))
common_points <- example_curves[age == 55]
peak_points <- example_curves[, .SD[which.max(log_likelihood)], by = .(group, biomarker)]
losses <- merge(common_points[, .(group, biomarker, common = log_likelihood)],
  peak_points[, .(group, biomarker, peak = log_likelihood)],
  by = c("group", "biomarker")
)
losses[, disagreement := peak - common]
projections <- merge(peak_points, losses[, .(group, biomarker, common)], by = c("group", "biomarker"))
projections <- projections[group == "Higher disagreement"]
curves_plot <- ggplot(example_curves, aes(age, log_likelihood, colour = biomarker)) +
  geom_line(linewidth = 0.5) +
  geom_vline(xintercept = 55, linetype = 2, colour = "grey50", linewidth = 0.3) +
  geom_segment(
    data = projections, aes(xend = age, yend = common),
    arrow = grid::arrow(ends = "both", length = grid::unit(1.3, "mm")), linewidth = 0.3
  ) +
  geom_segment(
    data = projections, aes(xend = 55, y = common, yend = common),
    linetype = 3, linewidth = 0.3
  ) +
  geom_point(data = peak_points, aes(shape = "Individual peak"), fill = "white", size = 1.8) +
  geom_point(data = common_points, aes(shape = "Common support age"), size = 1.8) +
  scale_shape_manual(values = c("Individual peak" = 21, "Common support age" = 16), name = NULL) +
  facet_wrap(~group) +
  labs(x = "Age (years)", y = "Relative log-likelihood", colour = NULL)
loss_plot <- ggplot(losses, aes(disagreement, biomarker)) +
  geom_segment(aes(x = 0, xend = disagreement, yend = biomarker), colour = "#80669D", linewidth = 0.5) +
  geom_point(colour = "#80669D", size = 1.7) +
  geom_point(
    data = losses[, .(disagreement = mean(disagreement), biomarker = "Mean disagreement"), by = group],
    shape = 18, size = 2
  ) +
  facet_wrap(~group) +
  labs(x = "Disagreement", y = NULL)
save(curves_plot / loss_plot + plot_layout(heights = c(2, 1)), "supplementary_figure_2a", 183, 120)
organs <- c("Artery", "Brain", "Heart", "Immune", "Intestine", "Kidney", "Liver", "Lung", "Muscle", "Pancreas")
pattern <- c(-1.4, 0.8, -0.2, 1.2, -0.8, 0.4, -1, 1.6, -0.4, -0.2)
pattern <- (pattern - mean(pattern)) / sqrt(mean((pattern - mean(pattern))^2))
organ_examples <- rbindlist(lapply(c(0.25, 1), function(spread) {
  data.table(
    organ = factor(organs, levels = rev(organs)), BAG = 0.5 + spread * pattern,
    group = if (spread == 0.25) "Lower disagreement (SD = 0.25)" else "Higher disagreement (SD = 1.00)"
  )
}))
organ_examples[, group := factor(group, levels = unique(group))]
save(ggplot(organ_examples, aes(BAG, organ, colour = group)) +
  geom_vline(xintercept = 0.5, linetype = 2, linewidth = 0.35) +
  geom_segment(aes(x = 0.5, xend = BAG, yend = organ), linewidth = 0.4) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = c("#89969E", "#C1644A")) +
  theme(legend.position = "none") +
  facet_wrap(~group, ncol = 1) +
  labs(x = "Organ BAG", y = NULL), "supplementary_figure_2b", 125, 140)
