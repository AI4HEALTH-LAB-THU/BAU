#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(splines)
  library(sandwich)
})
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- Sys.getenv("BAU_CODE_ROOT", normalizePath(file.path(dirname(script), "../..")))
data_root <- Sys.getenv("BAU_DATA_ROOT", "/home/dnanexus/bau-work/data")
results_root <- Sys.getenv("BAU_RESULTS_ROOT", "/home/dnanexus/bau-work/results")
output <- file.path(results_root, "phenotypes/disagreement")
dir.create(output, recursive = TRUE, showWarnings = FALSE)

analyse <- function(data, name, extra = character()) {
  data[, `:=`(sex = factor(Sex_male), centre = factor(AssessmentCenter))]
  variables <- c("H", "BAG", "Age")
  knots <- setNames(lapply(variables, function(variable) {
    list(
      knots = as.numeric(quantile(data[[variable]], c(0.25, 0.5, 0.75))),
      boundary = range(data[[variable]])
    )
  }), variables)
  basis <- function(values, spec) ns(values, knots = spec$knots, Boundary.knots = spec$boundary)
  for (variable in variables) {
    data[, (paste0(variable, "_s", 1:4)) := as.data.table(basis(data[[variable]], knots[[variable]]))]
  }
  exposure_terms <- paste0("H_s", 1:4)
  terms <- c(
    exposure_terms, paste0("BAG_s", 1:4), paste0("Age_s", 1:4),
    "sex", "centre", "input_missing_fraction", extra
  )
  fit <- lm(reformulate(terms, response = "BAU"), data = data, x = TRUE, na.action = na.fail)
  beta <- coef(fit)
  design <- fit$x
  centres <- uniqueN(data$centre)
  covariance <- vcovCL(fit, cluster = data$centre, type = "HC1", cadjust = TRUE)
  stopifnot(fit$rank == ncol(design), all(is.finite(beta)), nobs(fit) == nrow(data))
  critical <- qt(0.975, centres - 1)
  grid <- seq(quantile(data$H, 0.01), quantile(data$H, 0.99), length.out = 201)
  quartiles <- as.numeric(quantile(data$H, c(0.25, 0.75)))
  marginal_design <- function(values) {
    result <- matrix(rep(colMeans(design), each = length(values)),
      nrow = length(values),
      dimnames = list(NULL, names(beta))
    )
    result[, exposure_terms] <- basis(values, knots$H)
    result
  }
  grid_design <- marginal_design(grid)
  estimate <- drop(grid_design %*% beta)
  se <- sqrt(rowSums((grid_design %*% covariance) * grid_design))
  curves <- data.table(
    H = grid, adjusted_mean = estimate,
    ci_low = estimate - critical * se, ci_high = estimate + critical * se,
    n = nrow(data), n_centres = centres
  )
  quartile_design <- marginal_design(quartiles)
  contrast <- quartile_design[2, ] - quartile_design[1, ]
  delta <- sum(contrast * beta)
  delta_se <- sqrt(drop(t(contrast) %*% covariance %*% contrast))
  statistic <- drop(t(beta[exposure_terms]) %*% solve(
    covariance[exposure_terms, exposure_terms],
    beta[exposure_terms]
  )) / 4
  effects <- data.table(
    n = nrow(data), n_centres = centres, H_q25 = quartiles[1], H_q75 = quartiles[2],
    difference_q75_q25 = delta, se = delta_se,
    ci_low = delta - critical * delta_se, ci_high = delta + critical * delta_se,
    omnibus_F = statistic, p = pf(statistic, 4, centres - 1, lower.tail = FALSE)
  )
  fwrite(curves, file.path(output, paste0(name, "_adjusted_curves.tsv")), sep = "\t")
  fwrite(effects, file.path(output, paste0(name, "_adjusted_contrasts.tsv")), sep = "\t")
  distribution <- hist(data$H,
    breaks = seq(min(data$H), max(data$H) + 1e-10, length.out = 101), plot = FALSE
  )
  histogram <- data.table(
    bin_left = head(distribution$breaks, -1),
    bin_right = tail(distribution$breaks, -1), n = distribution$counts
  )
  fwrite(histogram, file.path(output, paste0(name, "_histogram.tsv")), sep = "\t")
  saveRDS(list(
    coefficients = beta, covariance = covariance, knots = knots,
    mean_design = colMeans(design)
  ), file.path(output, paste0(name, "_model.rds")))
}

for (modality in c("protein", "nmr")) {
  data <- as.data.table(read_feather(file.path(results_root, "phenotypes/molecular", modality, "scores.feather")))
  analyse(data, modality)
}
organ <- as.data.table(read_feather(file.path(results_root, "phenotypes/organ_disagreement.feather")))
protein <- as.data.table(read_feather(file.path(results_root, "phenotypes/molecular/protein/scores.feather")))
organ <- merge(organ, protein[, .(eid, input_missing_fraction)], by = "eid")
participants <- as.data.table(read_feather(file.path(data_root, "ukb/support/participants.feather"),
  col_select = c("eid", "Age")
))
stopifnot(!anyDuplicated(participants$eid), all(organ$eid %in% participants$eid))
organ[, Age := participants$Age[match(eid, participants$eid)]]
organ[, H := organ_disagreement]
analyse(organ, "organ", "organ_missing")
