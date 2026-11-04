#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(ggrepel)
library(patchwork)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])
default_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
root <- Sys.getenv("BAU_CODE_ROOT", default_root)
results_root <- Sys.getenv("BAU_RESULTS_ROOT", "/home/dnanexus/bau-work/results")
source(file.path(root, "R/theme.R"))
theme_set(paper_theme())
input <- file.path(results_root, "molecular")
output <- file.path(results_root, "figures")
save <- function(plot, name, width = 140, height = 82) {
  save_panel(plot, file.path(output, name), width, height)
}

cross_omic <- rbindlist(lapply(c("nmr_to_protein", "protein_to_nmr"), function(name) {
  fread(file.path(input, "cross_omic", paste0(name, "_associations.tsv")))
}))
molecule_labels <- fread(file.path(root, "config/displayed_molecules.tsv"))
status_colours <- c(`Other` = "#B9BDC0", `Discovery` = "#DDA58F", `Replicated` = bau_colours[["bau"]])
for (i in 1:2) {
  name <- c("nmr_to_protein", "protein_to_nmr")[[i]]
  panel_name <- paste0("figure_3", letters[[i]])
  x <- copy(cross_omic[analysis == name])
  x[, status := factor(
    fifelse(
      replicated_BAU, "Replicated",
      fifelse(discovery_hit_BAU, "Discovery", "Other")
    ),
    levels = names(status_colours)
  )]
  labels <- merge(x, molecule_labels[panel == panel_name], by = c("analysis", "feature"))
  volcano <- ggplot(x, aes(beta_BAU_discovery, -log10(p_BAU_discovery), colour = status)) +
    geom_point(size = 1.2, alpha = 0.8) +
    scale_colour_manual(values = status_colours) +
    geom_text_repel(
      data = labels, aes(label = display_label), colour = "#252525",
      size = 2.2, seed = 20260824, max.overlaps = Inf, show.legend = FALSE
    ) +
    labs(x = "BAU association coefficient", y = "-log10 P", title = "Discovery")
  validation_data <- x[discovery_hit_BAU == TRUE]
  validation_data[, validation_status := factor(
    fifelse(replicated_BAU, "Replicated", "Not replicated"),
    levels = c("Replicated", "Not replicated")
  )]
  validation <- ggplot(
    validation_data,
    aes(beta_BAU_discovery, beta_BAU_validation)
  ) +
    geom_abline(slope = 1, intercept = 0, colour = "grey65", linewidth = 0.3) +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_point(aes(fill = validation_status),
      shape = 21,
      colour = bau_colours[["bau"]], size = 1.3, stroke = 0.3
    ) +
    scale_fill_manual(values = c(`Replicated` = bau_colours[["bau"]], `Not replicated` = "white")) +
    geom_text_repel(
      data = labels[validation_label == TRUE], aes(label = display_label),
      size = 2.2, seed = 20260824, max.overlaps = Inf, show.legend = FALSE
    ) +
    annotate("text",
      x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, size = 2.4,
      label = sprintf("Pearson r = %.3f", cor(
        validation_data$beta_BAU_discovery,
        validation_data$beta_BAU_validation
      ))
    ) +
    labs(x = "Discovery coefficient", y = "Validation coefficient", title = "Held-centre validation")
  save(volcano + validation + plot_layout(guides = "collect") &
    theme(legend.position = "bottom"), panel_name)
}

coefficient_comparison <- function(analysis_id) {
  x <- copy(cross_omic[analysis == analysis_id])
  x[, absolute_coefficient_difference := abs(beta_BAU_discovery) - abs(beta_BAG_discovery)]
  x[, draw_order := abs(absolute_coefficient_difference)]
  setorder(x, draw_order)
  labels <- merge(x, molecule_labels[panel == "figure_3c"], by = c("analysis", "feature"))
  protein <- analysis_id == "nmr_to_protein"
  limit <- if (protein) 0.4 else 0.2
  ggplot(x, aes(beta_BAG_discovery, beta_BAU_discovery)) +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_abline(
      slope = c(-1, 1), intercept = 0, linetype = "dashed",
      colour = "grey65", linewidth = 0.3
    ) +
    geom_point(aes(colour = absolute_coefficient_difference),
      size = if (protein) 0.72 else 1.02, stroke = 0
    ) +
    geom_text_repel(
      data = labels, aes(label = display_label),
      colour = "#252525", family = "Helvetica", size = 2.1,
      nudge_x = labels$nudge_x, nudge_y = labels$nudge_y,
      min.segment.length = 0, segment.size = 0.23,
      seed = if (protein) 20260831 else 20260901,
      max.overlaps = Inf, show.legend = FALSE
    ) +
    annotate("text",
      x = -limit * 0.94, y = limit * 0.94, hjust = 0, vjust = 1,
      size = 2.2, label = sprintf(
        "Coefficient correlation: r = %.3f",
        cor(x$beta_BAG_discovery, x$beta_BAU_discovery)
      )
    ) +
    scale_colour_gradient2(
      low = "#4C78A8", mid = "#D9DDDC", high = "#B65F6D",
      midpoint = 0, limits = c(-0.26, 0.26), breaks = c(-0.2, 0, 0.2),
      name = expression(abs(beta[BAU]) - abs(beta[BAG])),
      guide = guide_colourbar(
        title.position = "left",
        barwidth = grid::unit(28, "mm"), barheight = grid::unit(2, "mm")
      )
    ) +
    coord_equal(xlim = c(-limit, limit), ylim = c(-limit, limit), clip = "off") +
    labs(
      x = "BAG association coefficient", y = "BAU association coefficient",
      title = if (protein) "NMR-derived BAU and BAG" else "Proteomic-derived BAU and BAG"
    ) +
    theme(legend.position = "bottom", legend.title = element_text(size = 7))
}
save(coefficient_comparison("nmr_to_protein") + coefficient_comparison("protein_to_nmr") +
  plot_layout(guides = "collect") & theme(legend.position = "bottom"), "figure_3c")

gsea_plot <- function(x) {
  x <- copy(x)
  x <- x[gene_set %in% x[replicated == TRUE, unique(gene_set)]]
  setorder(x, display_order)
  x[, label := factor(label, levels = rev(unique(label)))]
  x[, estimand := factor(estimand,
    levels = c("BAG", "BAU", "difference"),
    labels = c("BAG", "BAU", "BAU - BAG")
  )]
  ggplot(x, aes(estimand, label)) +
    geom_blank() +
    geom_point(
      data = x[replicated == TRUE],
      aes(size = pmin(-log10(padj_discovery), 20), fill = NES_discovery),
      shape = 21, colour = "#25282A", stroke = 0.35
    ) +
    scale_fill_gradient2(
      low = "#285F8F", mid = "white", high = "#B65F52",
      midpoint = 0, name = "NES"
    ) +
    scale_size_continuous(
      range = c(1, 4), name = "-log10 FDR",
      breaks = function(limits) pretty(limits, n = 3)
    ) +
    guides(
      size = guide_legend(order = 1, nrow = 1),
      fill = guide_colourbar(
        order = 2, barwidth = grid::unit(24, "mm"),
        barheight = grid::unit(2.5, "mm")
      )
    ) +
    facet_grid(collection_label ~ ., scales = "free_y", space = "free_y") +
    labs(x = NULL, y = NULL) +
    theme(
      legend.position = "bottom", legend.box = "vertical",
      legend.title = element_text(size = 7), axis.line = element_blank()
    )
}
for (population in c("primary", "healthy")) {
  for (modality in c("nmr", "protein")) {
    x <- fread(file.path(input, "enrichment", population, paste0(modality, "_gsea_displayed.tsv")))
    panel <- if (population == "primary") {
      paste0("figure_3d_", modality)
    } else {
      paste0("supplementary_figure_6", if (modality == "nmr") "d" else "e")
    }
    save(gsea_plot(x), panel, 120, 112)
  }
}
for (modality in c("nmr", "protein")) {
  x <- fread(file.path(input, "enrichment/primary", paste0(modality, "_go_displayed.tsv")))
  setorder(x, display_order)
  x[, direction := factor(direction, levels = c(
    "BAG positive", "BAG negative",
    "BAU positive", "BAU negative", "BAU > BAG", "BAU < BAG"
  ))]
  x[, label := factor(label, levels = rev(unique(label)))]
  x <- x[fdr < 0.05 & overlap >= 3 & fold_enrichment > 1]
  p <- ggplot(x, aes(direction, label, size = pmin(-log10(fdr), 20), fill = log2_enrichment)) +
    geom_point(shape = 21, colour = "grey35", stroke = 0.25) +
    scale_fill_gradient(low = "#F1DED8", high = "#B65F52", name = "log2 enrichment") +
    scale_size_continuous(range = c(1, 4), name = "-log10 FDR") +
    facet_grid(ontology ~ ., scales = "free_y", space = "free_y") +
    labs(x = NULL, y = NULL) +
    theme(legend.position = "bottom", axis.line = element_blank())
  save(p, paste0("supplementary_figure_6", if (modality == "nmr") "a" else "b"), 170, 115)
}

physiology <- fread(file.path(input, "physiology/physiological_associations.tsv"))
setnames(physiology, "label", "outcome_label")
long <- rbindlist(lapply(c("BAG", "BAU"), function(term) {
  physiology[, .(modality, outcome_id, outcome_order, domain, outcome_label,
    exposure = term, beta = get(paste0("beta_", term, "_discovery")),
    lower = get(paste0("ci_lower_", term, "_discovery")),
    upper = get(paste0("ci_upper_", term, "_discovery")),
    replicated = get(paste0("replicated_", term))
  )]
}))
setorder(long, outcome_order)
long[, `:=`(
  outcome_label = factor(outcome_label, levels = rev(unique(outcome_label))),
  domain = factor(domain, levels = unique(domain)),
  modality = factor(modality, levels = c("nmr", "protein"), labels = c("NMR-derived", "Proteomic-derived"))
)]
save(
  ggplot(long, aes(beta, outcome_label, colour = exposure, shape = exposure, group = exposure)) +
    geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.3) +
    geom_errorbar(aes(xmin = lower, xmax = upper),
      orientation = "y", width = 0.3,
      position = position_dodge(0.5), linewidth = 0.4
    ) +
    geom_point(position = position_dodge(0.5), size = 1.6, fill = "white") +
    facet_grid(domain ~ modality, scales = "free_y", space = "free_y") +
    scale_colour_manual(values = c(BAG = bau_colours[["bag"]], BAU = bau_colours[["bau"]])) +
    scale_shape_manual(values = c(BAG = 21, BAU = 23)) +
    labs(x = "Association coefficient", y = NULL) +
    theme(legend.position = "bottom"),
  "figure_3e", 168, 220
)

sensitivity <- fread(file.path(input, "sensitivity/bau_coefficients.tsv"))
agreement <- fread(file.path(input, "sensitivity/coefficient_agreement.tsv"))
specification_levels <- c("absolute_bag", "bag_spline", "central_95", "baseline_healthy", "low_missingness")
specification_labels <- c("Add absolute BAG", "BAG spline", "Central 95% BAG", "Baseline healthy", "Missingness <=1%")
for (x in list(sensitivity, agreement)) {
  x[, specification := factor(specification, levels = specification_levels, labels = specification_labels)]
  x[, analysis := factor(analysis,
    levels = c("nmr_to_protein", "protein_to_nmr"),
    labels = c("NMR-derived BAU", "Proteomic-derived BAU")
  )]
}
agreement[, annotation := sprintf(
  "r = %.3f\nSlope = %.3f\nSame direction = %.1f%%",
  pearson_r, zero_intercept_slope, 100 * direction_agreement
)]
save(ggplot(sensitivity, aes(primary_beta, beta, colour = primary_replicated)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_point(size = 0.7, alpha = 0.65) +
  facet_grid(analysis ~ specification) +
  geom_text(
    data = agreement, aes(x = -Inf, y = Inf, label = annotation),
    inherit.aes = FALSE, hjust = -0.06, vjust = 1.1, size = 2.1
  ) +
  scale_colour_manual(
    values = c(`FALSE` = "#BCC1C3", `TRUE` = bau_colours[["bau"]]),
    labels = c(`FALSE` = "Other", `TRUE` = "Primary replicated")
  ) +
  labs(x = "Primary BAU coefficient", y = "Sensitivity BAU coefficient") +
  theme(legend.position = "bottom"), "supplementary_figure_6c", 240, 95)
