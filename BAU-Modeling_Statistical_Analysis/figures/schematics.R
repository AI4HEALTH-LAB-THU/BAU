#!/usr/bin/env Rscript

library(grid)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])
default_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
root <- Sys.getenv("BAU_CODE_ROOT", default_root)
results_root <- Sys.getenv("BAU_RESULTS_ROOT", "/home/dnanexus/bau-work/results")
output <- file.path(results_root, "figures")
dir.create(output, recursive = TRUE, showWarnings = FALSE)

save_diagram <- function(draw, name, width, height) {
  grDevices::cairo_pdf(file.path(output, paste0(name, ".pdf")),
    width = width / 25.4, height = height / 25.4, family = "Helvetica"
  )
  draw()
  dev.off()
  svglite::svglite(file.path(output, paste0(name, ".svg")),
    width = width / 25.4, height = height / 25.4, bg = "white"
  )
  draw()
  dev.off()
}

local({
  W <- 183
  H <- 92
  text <- function(label, x, y, size = 7, bold = FALSE, colour = "#202426") {
    grid.text(label,
      x = unit(x, "mm"), y = unit(y, "mm"),
      gp = gpar(
        fontfamily = "Helvetica", fontsize = size, col = colour,
        fontface = if (bold) "bold" else "plain", lineheight = 1.3
      )
    )
  }
  box <- function(x1, x2, y1, y2, title, body, accent = "#687176") {
    grid.rect(
      x = unit((x1 + x2) / 2, "mm"), y = unit((y1 + y2) / 2, "mm"),
      width = unit(x2 - x1, "mm"), height = unit(y2 - y1, "mm"),
      gp = gpar(fill = "white", col = "#D8DDDF", lwd = 0.6)
    )
    grid.lines(unit(c(x1, x2), "mm"), unit(c(y2, y2), "mm"), gp = gpar(col = accent, lwd = 1.4))
    text(title, (x1 + x2) / 2, y2 - 4, 7, TRUE)
    text(body, (x1 + x2) / 2, (y1 + y2) / 2 - 2, 6.3)
  }
  arrow <- function(x1, x2, y1, y2) {
    grid.lines(unit(c(x1, x2), "mm"), unit(c(y1, y2), "mm"),
      arrow = grid::arrow(length = unit(1, "mm"), type = "open"),
      gp = gpar(col = "#687176", lwd = 0.65)
    )
  }
  draw <- function() {
    grid.newpage()
    grid.rect(gp = gpar(fill = "white", col = NA))
    text("Study overview", W / 2, 88, 8, TRUE)
    box(
      4, 101, 62, 82, "UK Biobank",
      "NMR metabolomics  n = 272,672\nPlasma proteomics  n = 43,867\nRoutine blood  n = 453,163"
    )
    box(
      105, 179, 62, 82, "External validation",
      "NHANES  n = 24,225\nCHARLS  n = 9,604\nCHNS  n = 5,137"
    )
    arrow(52.5, 52.5, 61, 57)
    arrow(142, 142, 61, 57)
    box(
      4, 179, 36, 56, "Probabilistic age models and healthy-reference phenotypes",
      "Healthy-reference model fitting  →  Independent calibration  →  Centre-held-out application\nBAG  ·  biological age gap                 BAU  ·  biological age uncertainty", "#80669D"
    )
    arrow(52.5, 52.5, 35, 30)
    arrow(137, 137, 35, 30)
    box(
      4, 91, 4, 29, "Health outcomes and longevity",
      "Mortality, healthspan and multimorbidity\nIncident diseases and risk prediction\nLongitudinal change and external validation", "#C1644A"
    )
    box(
      95, 179, 4, 29, "Molecular and physiological associations",
      "Molecular and organ-age disagreement\nProteome-wide and metabolome-wide signatures\nPathway enrichment and physiological traits", "#2F7F78"
    )
  }
  save_diagram(draw, "figure_1a", W, H)
})

# BAG and BAU definition with illustrative Normal densities.
local({
  W <- 102
  H <- 41
  mu <- 61
  sigma <- 5
  chronological <- 48
  age_min <- 42
  age_max <- 78
  baseline <- 9.5
  ink <- "#000000"
  grey <- "#777A7C"
  blue <- "#285F8F"
  clay <- "#C1644A"
  pale <- "#D79A89"
  family <- "Helvetica"
  u <- function(x) unit(x, "native")
  txt <- function(label, x, y, size = 6.2, col = ink, face = "plain", rot = 0) {
    grid.text(label, u(x), u(y), rot = rot, gp = gpar(
      fontfamily = family, fontsize = size,
      col = col, fontface = face, lineheight = 0.95
    ))
  }
  line <- function(x, y, col = grey, lwd = 0.75, lty = 1, arrow = NULL) {
    grid.lines(u(x), u(y),
      arrow = arrow,
      gp = gpar(col = col, fill = col, lwd = lwd, lty = lty, lineend = "butt")
    )
  }
  view <- function(left, right, comparison = FALSE) {
    xx <- function(age) left + (age - age_min) / (age_max - age_min) * (right - left)
    # Use the same density scale for both curves.
    scale <- if (comparison) 155 else 22 / dnorm(mu, mu, sigma)
    yy <- function(age, sd) baseline + dnorm(age, mu, sd) * scale
    curve <- function(sd, col, alpha, lty = 1) {
      ages <- seq(age_min, age_max, length.out = 1601)
      grid.polygon(u(c(xx(ages), right, left)), u(c(yy(ages, sd), baseline, baseline)),
        gp = gpar(fill = adjustcolor(col, alpha.f = alpha), col = NA)
      )
      line(xx(ages), yy(ages, sd), col = col, lwd = 1.2, lty = lty)
    }
    line(c(left, right), rep(baseline, 2), col = "#2B2D2F", lwd = 0.7)
    line(rep(left, 2), c(baseline, 32.5), col = "#2B2D2F", lwd = 0.7)
    txt("Probability density", left - 3.5, 21, size = 6.1, rot = 90)
    if (comparison) {
      curve(6.3, clay, 0.17)
      curve(2.8, pale, 0.11, lty = 2)
    } else {
      curve(sigma, clay, 0.17)
    }
    line(rep(xx(chronological), 2), c(5.4, 29), col = grey, lwd = 0.7)
    line(rep(xx(mu), 2), c(5.4, baseline), col = blue, lwd = 0.7)
    if (comparison) {
      line(rep(xx(mu), 2), c(baseline, yy(mu, 2.8)), col = clay, lwd = 0.7, lty = 2)
    } else {
      line(rep(xx(mu), 2), c(baseline, 20.8), col = clay, lwd = 0.7, lty = 2)
      line(rep(xx(mu), 2), c(23.8, yy(mu, sigma)), col = clay, lwd = 0.7, lty = 2)
    }
    grid.circle(u(xx(chronological)), u(baseline),
      r = unit(0.53, "mm"),
      gp = gpar(col = grey, fill = "white", lwd = 0.8)
    )
    grid.circle(u(xx(mu)), u(baseline),
      r = unit(0.53, "mm"),
      gp = gpar(col = blue, fill = blue, lwd = 0.8)
    )
    line(c(xx(chronological) + 0.15, xx(mu) - 0.15), c(5.6, 5.6),
      col = blue, lwd = 0.9,
      arrow = arrow(angle = 25, length = unit(0.85, "mm"), type = "closed")
    )
    txt(if (comparison) "Same BAG" else "BAG",
      mean(c(xx(chronological), xx(mu))), 7.5,
      size = if (comparison) 5.5 else 5.9, col = blue, face = "bold"
    )
    txt("Chronological\nage", xx(chronological), 1.9, size = 5.8)
    txt("Predicted\nmean age", xx(mu), 1.9, size = 5.8)
    if (!comparison) {
      width_y <- 19.7
      hw <- sigma * sqrt(-2 * log((width_y - baseline) / 22))
      line(xx(mu + c(-hw, hw)), rep(width_y, 2),
        col = clay, lwd = 0.9,
        arrow = arrow(angle = 25, length = unit(0.8, "mm"), ends = "both", type = "closed")
      )
      txt("BAU", xx(mu), 22.3, size = 5.8, col = clay, face = "bold")
    } else {
      txt("Low BAU", left + 9.3, 31.8, size = 6.2, col = clay)
      line(c(left + 13, xx(58.9)), c(30, yy(58.9, 2.8) + 0.3), col = pale, lwd = 0.75)
      txt("High BAU", right - 9, 24.8, size = 6.2, col = clay)
      line(c(right - 10.5, xx(68.5)), c(23, yy(68.5, 6.3) + 0.4), col = clay, lwd = 0.75)
    }
  }
  draw <- function() {
    grid.newpage()
    pushViewport(viewport(xscale = c(0, W), yscale = c(-1, H - 1)))
    grid.rect(gp = gpar(fill = "white", col = NA))
    txt("b", 1.7, 37.8, size = 8.2, face = "bold")
    txt("BAG and BAU", 27.5, 37.5, size = 6.9, face = "bold")
    txt("Same BAG, different BAU", 78, 37.5, size = 6.5, face = "bold")
    view(7, 48, FALSE)
    view(56, 100, TRUE)
    popViewport()
  }
  save_diagram(draw, "figure_1b", W, H)
})

local({
  folds <- jsonlite::read_json(file.path(root, "config/center_folds.json"), simplifyVector = TRUE)
  assignment <- folds$outer_center_to_fold
  centre_key <- data.frame(
    outer_fold = as.integer(assignment),
    assessment_centre = as.integer(names(assignment))
  )
  centre_key <- centre_key[order(centre_key$outer_fold, centre_key$assessment_centre), ]
  # Centre 10003 has no NMR measurements.
  centre_key$nmr_available <- centre_key$assessment_centre != 10003L
  centre_key$proteomic_available <- TRUE
  WIDTH_MM <- 183
  HEIGHT_MM <- 82
  FONT_FAMILY <- "Helvetica"

  COL_INK <- "#171717"
  COL_TEXT2 <- "#4E5557"
  COL_RULE <- "#C9CDCE"
  COL_FAINT <- "#ECEEEE"
  COL_SOFT <- "#F6F7F7"
  COL_FIT <- "#565D60"
  COL_CAL <- "#CBD0D1"
  COL_HELD <- "#6D8790"
  COL_NMR <- "#2F7F78"
  COL_PROTEIN <- "#80669D"
  COL_BAG <- "#285F8F"
  COL_BAU <- "#C1644A"


  draw_text <- function(label, x, y, size = 6, face = "plain", hjust = 0.5,
                        vjust = 0.5, colour = COL_INK, lineheight = 1.0,
                        rot = 0) {
    grid.text(
      label,
      x = unit(x, "mm"), y = unit(y, "mm"),
      just = c(hjust, vjust), rot = rot,
      gp = gpar(
        fontfamily = FONT_FAMILY, fontsize = size, fontface = face,
        col = colour, lineheight = lineheight
      )
    )
  }

  draw_rect <- function(x1, x2, y1, y2, fill = "white", colour = COL_RULE,
                        line_width = 0.55) {
    grid.rect(
      x = unit((x1 + x2) / 2, "mm"),
      y = unit((y1 + y2) / 2, "mm"),
      width = unit(x2 - x1, "mm"), height = unit(y2 - y1, "mm"),
      gp = gpar(
        fill = fill, col = colour, lwd = line_width,
        linejoin = "mitre"
      )
    )
  }

  draw_line <- function(x1, x2, y1, y2, colour = COL_TEXT2,
                        line_width = 0.65, line_type = 1) {
    grid.lines(
      x = unit(c(x1, x2), "mm"), y = unit(c(y1, y2), "mm"),
      gp = gpar(
        col = colour, lwd = line_width, lty = line_type,
        lineend = "butt", linejoin = "mitre"
      )
    )
  }

  draw_arrow <- function(x1, x2, y1, y2, colour = COL_TEXT2,
                         line_width = 0.72, arrow_length = 0.72) {
    grid.lines(
      x = unit(c(x1, x2), "mm"), y = unit(c(y1, y2), "mm"),
      gp = gpar(col = colour, lwd = line_width, lineend = "butt"),
      arrow = arrow(
        angle = 25, length = unit(arrow_length, "mm"),
        ends = "last", type = "open"
      )
    )
  }

  draw_dot <- function(x, y, radius = 1.02, fill = COL_INK,
                       colour = fill, line_width = 0.45) {
    grid.circle(
      x = unit(x, "mm"), y = unit(y, "mm"),
      r = unit(radius, "mm"),
      gp = gpar(fill = fill, col = colour, lwd = line_width)
    )
  }

  draw_segmented_rule <- function(x1, x2, y, colours, proportions,
                                  line_width = 1.45) {
    proportions <- proportions / sum(proportions)
    cuts <- x1 + c(0, cumsum(proportions)) * (x2 - x1)
    for (i in seq_along(colours)) {
      draw_line(
        cuts[[i]], cuts[[i + 1L]], y, y,
        colour = colours[[i]], line_width = line_width
      )
    }
  }

  draw_predictive_distribution <- function(cx, baseline, span = 25,
                                           sd_mm = 4.2, peak = 5.5) {
    x <- seq(cx - span / 2, cx + span / 2, length.out = 160)
    y <- baseline + peak * exp(-0.5 * ((x - cx) / sd_mm)^2)
    grid.polygon(
      x = unit(c(x, rev(x)), "mm"),
      y = unit(c(y, rep(baseline, length(y))), "mm"),
      gp = gpar(fill = adjustcolor(COL_HELD, alpha.f = 0.13), col = NA)
    )
    draw_line(cx - span / 2, cx + span / 2, baseline, baseline, COL_RULE, 0.45)
    grid.lines(
      x = unit(x, "mm"), y = unit(y, "mm"),
      gp = gpar(col = COL_HELD, lwd = 0.85, lineend = "round")
    )
    draw_line(cx, cx, baseline, baseline + peak, COL_BAG, 0.8)
    draw_line(
      cx - sd_mm, cx + sd_mm, baseline - 1.25, baseline - 1.25,
      COL_BAU, 0.8
    )
    draw_line(
      cx - sd_mm, cx - sd_mm, baseline - 1.85, baseline - 0.65,
      COL_BAU, 0.8
    )
    draw_line(
      cx + sd_mm, cx + sd_mm, baseline - 1.85, baseline - 0.65,
      COL_BAU, 0.8
    )
  }

  draw_fold_row <- function(modality_key, y, colour) {
    fold_starts <- 48 + (0:4) * 25.65
    box_width <- 23.45

    for (fold in 1:5) {
      x1 <- fold_starts[[fold]]
      x2 <- x1 + box_width
      draw_rect(x1, x2, y - 3.0, y + 3.0,
        fill = "white",
        colour = COL_RULE, line_width = 0.5
      )

      fold_centres <- centre_key[centre_key$outer_fold == fold, ]
      n_slot <- nrow(fold_centres)
      dot_x <- seq(x1 + 3.3, x2 - 3.3, length.out = n_slot)
      available <- if (modality_key == "protein") {
        fold_centres$proteomic_available
      } else {
        fold_centres$nmr_available
      }
      for (j in seq_len(n_slot)) {
        if (available[[j]]) {
          draw_dot(dot_x[[j]], y, fill = colour, colour = colour)
        } else {
          draw_dot(dot_x[[j]], y,
            fill = "white", colour = COL_TEXT2,
            line_width = 0.65
          )
          draw_line(
            dot_x[[j]] - 0.60, dot_x[[j]] + 0.60,
            y - 0.60, y + 0.60, COL_TEXT2, 0.55
          )
        }
      }
    }
  }


  draw_figure <- function() {
    grid.newpage()
    grid.rect(gp = gpar(fill = "white", col = NA))

    draw_text("a", 2.4, 79.1, size = 8.7, face = "bold", hjust = 0)
    draw_text("Centre-held-out outer cross-validation", 8.2, 77.9,
      size = 7.0, face = "bold", hjust = 0
    )
    draw_text(
      "22 UK Biobank assessment centres assigned to five prespecified folds",
      8.2, 73.8,
      size = 5.25, hjust = 0, colour = COL_TEXT2
    )

    fold_starts <- 48 + (0:4) * 25.65
    fold_centres_x <- fold_starts + 23.45 / 2
    for (fold in 1:5) {
      draw_text(paste0("Fold ", fold), fold_centres_x[[fold]], 68.6,
        size = 5.15, face = "bold"
      )
    }
    draw_text("Plasma proteomics", 8.2, 62.7,
      size = 5.6,
      face = "bold", hjust = 0
    )
    draw_text("22 centres", 8.2, 59.8,
      size = 4.8, hjust = 0,
      colour = COL_PROTEIN
    )
    draw_fold_row("protein", 61.6, COL_PROTEIN)

    draw_text("NMR metabolomics", 8.2, 53.5,
      size = 5.6,
      face = "bold", hjust = 0
    )
    draw_text("21 centres", 8.2, 50.6,
      size = 4.8, hjust = 0,
      colour = COL_NMR
    )
    draw_fold_row("nmr", 52.4, COL_NMR)

    draw_dot(48.9, 46.8,
      radius = 0.82, fill = "white",
      colour = COL_TEXT2, line_width = 0.6
    )
    draw_line(48.4, 49.4, 46.3, 47.3, COL_TEXT2, 0.5)
    draw_text("Centre 10003: no available NMR measurements", 51.3, 46.8,
      size = 4.65, hjust = 0, colour = COL_TEXT2
    )
    draw_line(6.0, 180.0, 43.5, 43.5, COL_RULE, 0.6)

    # Independent model-fit and calibration partitions.
    draw_text("Fold-specific model development and application", 8.2, 40.6,
      size = 6.2,
      face = "bold", hjust = 0
    )

    draw_rect(6.0, 83.0, 12.2, 37.8,
      fill = "white",
      colour = COL_RULE, line_width = 0.6
    )
    draw_line(6.0, 83.0, 37.8, 37.8, COL_FIT, 1.5)
    draw_text("Development centres", 8.4, 35.3,
      size = 5.8, face = "bold", hjust = 0
    )
    draw_text("Baseline-healthy reference sample", 8.4, 32.3,
      size = 4.9, hjust = 0, colour = COL_TEXT2
    )

    split_x1 <- 8.4
    split_x2 <- 80.6
    split_cut <- split_x1 + 0.8 * (split_x2 - split_x1)
    draw_rect(split_x1, split_cut, 28.2, 30.3,
      fill = COL_FIT,
      colour = COL_FIT, line_width = 0.4
    )
    draw_rect(split_cut, split_x2, 28.2, 30.3,
      fill = COL_CAL,
      colour = COL_CAL, line_width = 0.4
    )
    draw_text("80%", (split_x1 + split_cut) / 2, 29.25,
      size = 4.5, face = "bold", colour = "white"
    )
    draw_text("20%", (split_cut + split_x2) / 2, 29.25,
      size = 4.5, face = "bold", colour = COL_INK
    )

    draw_line(48.1, 48.1, 14.8, 26.5, COL_FAINT, 0.65)
    draw_text("Model fitting", 8.4, 25.4,
      size = 5.2, face = "bold", hjust = 0
    )
    draw_text("Predictors: omic features and sex", 8.4, 22.1,
      size = 4.75, hjust = 0
    )
    draw_text("Target: chronological age", 8.4, 19.2,
      size = 4.75, hjust = 0
    )
    draw_text("Preprocessing and hyperparameter tuning", 8.4, 16.3,
      size = 4.55, hjust = 0, colour = COL_TEXT2
    )

    draw_text("Calibration reference", 50.4, 25.4,
      size = 5.05, face = "bold", hjust = 0
    )
    draw_text("Mean calibration by sex", 50.4, 22.1,
      size = 4.55, hjust = 0
    )
    draw_text("Variance scaling by sex", 50.4, 19.2,
      size = 4.55, hjust = 0
    )
    draw_text("Sex-specific healthy age reference", 50.4, 16.3,
      size = 4.55, hjust = 0
    )

    draw_arrow(84.0, 88.6, 25.0, 25.0, COL_TEXT2, 0.72, 0.67)
    draw_rect(89.5, 128.5, 12.2, 37.8,
      fill = "white",
      colour = COL_RULE, line_width = 0.6
    )
    draw_line(89.5, 128.5, 37.8, 37.8, COL_HELD, 1.5)
    draw_text("Held-out centres", 92.0, 35.3,
      size = 5.8,
      face = "bold", hjust = 0
    )
    draw_text("All eligible participants", 92.0, 32.3,
      size = 5.0,
      face = "bold", hjust = 0
    )
    draw_text("Frozen fold-specific model and calibrators", 109.0, 27.5,
      size = 4.45, colour = COL_TEXT2
    )
    draw_predictive_distribution(109.0, 18.7,
      span = 25.0,
      sd_mm = 4.1, peak = 5.0
    )
    draw_text("Calibrated predicted-age mean and variance", 109.0, 14.0,
      size = 4.45, colour = COL_TEXT2
    )

    draw_line(128.5, 133.0, 25.0, 25.0, COL_TEXT2, 0.7)
    draw_line(133.0, 133.0, 18.6, 31.4, COL_TEXT2, 0.7)
    draw_arrow(133.0, 137.0, 31.4, 31.4, COL_BAG, 0.78, 0.65)
    draw_arrow(133.0, 137.0, 18.6, 18.6, COL_BAU, 0.78, 0.65)

    draw_rect(138.0, 179.0, 25.2, 37.8,
      fill = "white",
      colour = COL_RULE, line_width = 0.6
    )
    draw_line(138.0, 138.0, 25.2, 37.8, COL_BAG, 1.65)
    draw_text("BAG", 141.0, 34.7,
      size = 5.8, face = "bold",
      hjust = 0, colour = COL_BAG
    )
    draw_text("Calibrated age gap", 141.0, 30.9,
      size = 4.75, hjust = 0
    )
    draw_text("standardized to the healthy reference", 141.0, 28.2,
      size = 4.75, hjust = 0
    )

    draw_rect(138.0, 179.0, 12.2, 24.0,
      fill = "white",
      colour = COL_RULE, line_width = 0.6
    )
    draw_line(138.0, 138.0, 12.2, 24.0, COL_BAU, 1.65)
    draw_text("BAU", 141.0, 21.1,
      size = 5.8, face = "bold",
      hjust = 0, colour = COL_BAU
    )
    draw_text("Log predictive variance", 141.0, 17.5,
      size = 4.75, hjust = 0
    )
    draw_text("standardized to the healthy reference", 141.0, 14.8,
      size = 4.75, hjust = 0
    )

    draw_segmented_rule(6.0, 179.0, 9.4, c(COL_BAG, COL_BAU), c(0.5, 0.5), 1.4)
    draw_text("Five outer-test sets pooled to yield one OOF BAG/BAU pair per participant",
      92.5, 6.5,
      size = 4.85, face = "bold"
    )

    draw_line(179.0, 181.0, 31.4, 31.4, COL_TEXT2, 0.6)
    draw_line(179.0, 181.0, 18.6, 18.6, COL_TEXT2, 0.6)
    draw_line(181.0, 181.0, 6.5, 31.4, COL_TEXT2, 0.6)
    draw_arrow(181.0, 179.4, 6.5, 6.5, COL_TEXT2, 0.66, 0.52)
  }

  save_diagram(draw_figure, "supplementary_figure_1a", WIDTH_MM, HEIGHT_MM)
})
