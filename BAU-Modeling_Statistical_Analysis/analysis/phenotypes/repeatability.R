#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
default_root <- normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
root <- Sys.getenv("BAU_CODE_ROOT", default_root)
results_root <- Sys.getenv("BAU_RESULTS_ROOT", "/home/dnanexus/bau-work/results")
pairs_path <- file.path(results_root, "phenotypes/repeat_pairs.feather")
source_dir <- file.path(results_root, "phenotypes/repeatability")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

dat <- as.data.table(read_feather(pairs_path))
setorder(dat, eid)
required <- c(
  "eid", "baseline_bag", "repeat_bag", "baseline_bau", "repeat_bau",
  "baseline_center", "visit_interval_years"
)
stopifnot(all(required %in% names(dat)))
stopifnot(nrow(dat) == 13620L, uniqueN(dat$eid) == nrow(dat))
stopifnot(all(vapply(dat[, ..required], function(x) all(!is.na(x)), logical(1))))

phenotype_statistics <- function(x, phenotype) {
  baseline <- x[[paste0("baseline_", phenotype)]]
  repeat_value <- x[[paste0("repeat_", phenotype)]]
  c(spearman_r = cor(baseline, repeat_value, method = "spearman"))
}

all_statistics <- function(x) {
  bag <- phenotype_statistics(x, "bag")
  bau <- phenotype_statistics(x, "bau")
  c(
    setNames(bag, paste0("bag__", names(bag))),
    setNames(bau, paste0("bau__", names(bau)))
  )
}

observed <- all_statistics(dat)
cluster_indices <- split(seq_len(nrow(dat)), dat$baseline_center)
cluster_names <- names(cluster_indices)
stopifnot(length(cluster_names) == 19L)

set.seed(20260827)
n_boot <- 2000L
boot <- matrix(NA_real_, nrow = n_boot, ncol = length(observed))
colnames(boot) <- names(observed)
for (b in seq_len(n_boot)) {
  sampled_centres <- sample(cluster_names, length(cluster_names), replace = TRUE)
  indices <- unlist(cluster_indices[sampled_centres], use.names = FALSE)
  boot[b, ] <- all_statistics(dat[indices])
}

ci <- apply(boot, 2, quantile, probs = c(0.025, 0.975), type = 8, na.rm = TRUE)
metric_names <- names(observed)
metric_table <- data.table(
  phenotype = toupper(sub("__.*$", "", metric_names)),
  metric = sub("^.*__", "", metric_names),
  estimate = as.numeric(observed),
  ci_low = as.numeric(ci[1, ]),
  ci_high = as.numeric(ci[2, ]),
  n = nrow(dat),
  n_baseline_centres = length(cluster_names),
  bootstrap_replicates = n_boot,
  interval = "95% baseline-assessment-centre cluster-bootstrap percentile CI"
)
fwrite(metric_table, file.path(source_dir, "persistence_metrics.tsv"), sep = "\t")

density_rows <- list()
for (phenotype in c("bag", "bau")) {
  baseline_name <- paste0("baseline_", phenotype)
  repeat_name <- paste0("repeat_", phenotype)
  limit <- ceiling(max(abs(c(dat[[baseline_name]], dat[[repeat_name]]))) * 2) / 2
  breaks <- seq(-limit, limit, length.out = 49L)
  working <- dat[, .(
    baseline_bin = as.integer(cut(get(baseline_name), breaks, include.lowest = TRUE)),
    repeat_bin = as.integer(cut(get(repeat_name), breaks, include.lowest = TRUE))
  )]
  density <- working[!is.na(baseline_bin) & !is.na(repeat_bin), .(count = .N), by = .(baseline_bin, repeat_bin)]
  density[, `:=`(
    phenotype = toupper(phenotype),
    axis_limit = limit,
    baseline_low = breaks[baseline_bin],
    baseline_high = breaks[baseline_bin + 1L],
    repeat_low = breaks[repeat_bin],
    repeat_high = breaks[repeat_bin + 1L]
  )]
  density_rows[[phenotype]] <- density
}
density_table <- rbindlist(density_rows, use.names = TRUE)
setcolorder(density_table, c(
  "phenotype", "axis_limit", "baseline_bin", "repeat_bin",
  "baseline_low", "baseline_high", "repeat_low", "repeat_high", "count"
))
fwrite(density_table, file.path(source_dir, "paired_density.tsv"), sep = "\t")

cohort <- data.table(
  quantity = c(
    "Primary paired participants", "Baseline assessment centres",
    "Median visit interval (years)", "Visit interval IQR lower (years)",
    "Visit interval IQR upper (years)", "Baseline BAG mean", "Repeat BAG mean",
    "Mean BAG change", "Baseline BAU mean", "Repeat BAU mean", "Mean BAU change"
  ),
  value = c(
    nrow(dat), length(cluster_names), median(dat$visit_interval_years),
    quantile(dat$visit_interval_years, 0.25, type = 8),
    quantile(dat$visit_interval_years, 0.75, type = 8),
    mean(dat$baseline_bag), mean(dat$repeat_bag), mean(dat$repeat_bag - dat$baseline_bag),
    mean(dat$baseline_bau), mean(dat$repeat_bau), mean(dat$repeat_bau - dat$baseline_bau)
  )
)
fwrite(cohort, file.path(source_dir, "cohort_summary.tsv"), sep = "\t")

capture.output(sessionInfo(), file = file.path(source_dir, "R_session_info.txt"))
cat("Persistence metrics\n")
print(metric_table)
