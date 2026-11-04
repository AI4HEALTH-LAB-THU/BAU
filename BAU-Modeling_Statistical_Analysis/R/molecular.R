suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(digest)
  library(matrixStats)
  library(splines)
})

collapse_ethnicity <- function(x) {
  fcase(
    x %in% c(1L, 1001L, 1002L, 1003L), "White",
    x %in% c(2L, 2001L, 2002L, 2003L, 2004L), "Mixed",
    x %in% c(3L, 3001L, 3002L, 3003L, 3004L), "Asian",
    x %in% c(4L, 4001L, 4002L, 4003L), "Black",
    x == 5L, "Chinese",
    x == 6L, "Other",
    default = "Unknown"
  )
}

read_feature_order <- function(path) {
  x <- trimws(readLines(path, warn = FALSE))
  x[nzchar(x)]
}

read_molecular_predictions <- function(code_root, modality) {
  results_root <- Sys.getenv("BAU_RESULTS_ROOT", "/home/dnanexus/bau-work/results")
  path <- file.path(
    results_root, "models", modality, "ngboost", "all_participant_predictions.csv"
  )
  columns <- c(
    "eid_int", "Age", "Sex_male", "Baseline_unhealthy",
    "AssessmentCenter", "OuterFold", "BAG", "BAU"
  )
  x <- fread(path, select = columns, showProgress = FALSE)
  setnames(x, "eid_int", "eid")
  x[, eid := as.integer(eid)]
  if (anyDuplicated(x$eid) || !setequal(x$OuterFold, 1:5)) {
    stop(paste("Invalid held-centre predictions for", modality))
  }
  x
}

read_ethnicity <- function(ukb_root) {
  path <- file.path(
    ukb_root,
    "ProcessedData/02-AssessmentCenter/02-Touchscreen/sociodemographics_0.csv"
  )
  x <- fread(path, select = c("eid", "Ethnic background"), showProgress = FALSE)
  setnames(x, "Ethnic background", "ethnicity_code")
  x[, eid := as.integer(eid)]
  if (anyDuplicated(x$eid)) stop("Repeated participant identifier in ethnicity data")
  x
}

add_molecular_metadata <- function(meta, dates, ethnicity) {
  dates <- copy(dates)
  dates[, eid := as.integer(eid)]
  keep <- c("eid", "Date_birth", "Date_attend")
  meta <- merge(meta, dates[, ..keep], by = "eid", all.x = TRUE, sort = FALSE)
  meta <- merge(meta, ethnicity, by = "eid", all.x = TRUE, sort = FALSE)
  attendance <- as.IDate(meta$Date_attend)
  birth <- as.IDate(meta$Date_birth)
  meta[, `:=`(
    Age_exact = as.numeric(attendance - birth) / 365.2425,
    ethnicity_group = factor(
      collapse_ethnicity(as.integer(ethnicity_code)),
      levels = c("White", "Mixed", "Asian", "Black", "Chinese", "Other", "Unknown")
    ),
    AssessmentCenter = factor(as.character(AssessmentCenter)),
    Sex_male = as.integer(Sex_male)
  )]
  if (any(!is.finite(meta$Age_exact)) || max(abs(meta$Age_exact - meta$Age)) > 0.051) {
    stop("Date-derived age does not agree with model age")
  }
  meta
}

filter_target_matrix <- function(meta, y, feature_names, threshold = 0.20) {
  storage.mode(y) <- "double"
  participant_missing <- rowMeans(is.na(y))
  keep_participant <- participant_missing < threshold
  meta <- meta[keep_participant]
  y <- y[keep_participant, , drop = FALSE]
  feature_missing <- colMeans(is.na(y))
  keep_feature <- feature_missing < threshold
  y <- y[, keep_feature, drop = FALSE]
  feature_names <- feature_names[keep_feature]
  colnames(y) <- feature_names
  meta[, target_missing_fraction := rowMeans(is.na(y))]
  list(
    meta = meta,
    y = y,
    features = feature_names,
    feature_missingness = feature_missing[keep_feature]
  )
}

load_protein_targets <- function(code_root, ukb_root, modality, exclude_inputs = FALSE) {
  results_root <- Sys.getenv("BAU_RESULTS_ROOT", "/home/dnanexus/bau-work/results")
  data_root <- Sys.getenv("BAU_DATA_ROOT", "/home/dnanexus/bau-work/data")
  predictions <- read_molecular_predictions(code_root, modality)
  candidate_path <- file.path(data_root, "ukb/modeling", paste0(modality, ".feather"))
  dates <- as.data.table(read_feather(
    candidate_path,
    col_select = c("eid", "Date_birth", "Date_attend")
  ))
  olink_path <- file.path(
    ukb_root, "ProcessedData/03-BiologicalSamples/ukb_protein_full.csv"
  )
  olink <- fread(olink_path, showProgress = TRUE)
  olink[, eid := as.integer(eid)]
  feature_names <- setdiff(names(olink), "eid")
  common <- intersect(predictions$eid, olink$eid)
  meta <- predictions[eid %in% common]
  setorder(meta, eid)
  olink <- olink[match(meta$eid, eid)]
  if (exclude_inputs) {
    clean_path <- file.path(
      results_root, "molecular/cross_omic/nmr_to_protein_transformations.tsv"
    )
    clean_features <- fread(clean_path, select = "feature")$feature
    y <- as.matrix(olink[, ..clean_features])
    storage.mode(y) <- "double"
    keep_participant <- rowMeans(is.na(y)) < 0.20
    meta <- meta[keep_participant]
    y <- y[keep_participant, , drop = FALSE]
    feature_missingness <- colMeans(is.na(y))
    if (any(feature_missingness >= 0.20)) {
      stop("A retained protein has at least 20% missing values")
    }
    meta[, target_missing_fraction := rowMeans(is.na(y))]
    clock_inputs <- read_feature_order(file.path(
      code_root, "config/features/protein.txt"
    ))
    keep <- !clean_features %in% clock_inputs
    dataset <- list(
      meta = meta,
      y = y[, keep, drop = FALSE],
      features = clean_features[keep],
      feature_missingness = feature_missingness[keep]
    )
  } else {
    y <- as.matrix(olink[, ..feature_names])
    dataset <- filter_target_matrix(meta, y, feature_names)
  }
  dataset$meta <- add_molecular_metadata(
    dataset$meta, dates, read_ethnicity(ukb_root)
  )
  dataset$target_layer <- "protein"
  dataset$source_n <- nrow(predictions)
  dataset
}

load_nmr_targets <- function(code_root, ukb_root, modality = "protein") {
  data_root <- Sys.getenv("BAU_DATA_ROOT", "/home/dnanexus/bau-work/data")
  predictions <- read_molecular_predictions(code_root, modality)
  feature_names <- read_feature_order(file.path(
    code_root, "config/features/nmr.txt"
  ))
  candidate <- as.data.table(read_feather(
    file.path(data_root, "ukb/modeling/nmr.feather")
  ))
  candidate[, eid := as.integer(eid)]
  common <- intersect(predictions$eid, candidate$eid)
  meta <- predictions[eid %in% common]
  setorder(meta, eid)
  candidate <- candidate[match(meta$eid, eid)]
  y <- as.matrix(candidate[, ..feature_names])
  dates <- candidate[, .(eid, Date_birth, Date_attend)]
  dataset <- filter_target_matrix(meta, y, feature_names)
  dates <- dates[match(dataset$meta$eid, eid)]
  dataset$meta <- add_molecular_metadata(
    dataset$meta, dates, read_ethnicity(ukb_root)
  )
  dataset$target_layer <- "nmr"
  dataset$source_n <- nrow(predictions)
  dataset
}

transform_molecular_targets <- function(dataset, discovery_rows, tail_fraction = 0.005) {
  y <- dataset$y
  limits <- colQuantiles(
    y[discovery_rows, , drop = FALSE],
    probs = c(tail_fraction, 1 - tail_fraction),
    na.rm = TRUE,
    type = 8
  )
  for (j in seq_len(ncol(y))) {
    observed <- is.finite(y[, j])
    y[observed, j] <- pmax(limits[j, 1], pmin(limits[j, 2], y[observed, j]))
  }
  healthy_rows <- discovery_rows[
    dataset$meta$Baseline_unhealthy[discovery_rows] == 0L
  ]
  reference_mean <- colMeans2(y[healthy_rows, , drop = FALSE], na.rm = TRUE)
  reference_sd <- colSds(y[healthy_rows, , drop = FALSE], na.rm = TRUE)
  if (any(!is.finite(reference_sd)) || any(reference_sd <= 1e-8)) {
    stop("A molecular outcome has no usable reference variance")
  }
  y <- sweep(y, 2, reference_mean, "-")
  y <- sweep(y, 2, reference_sd, "/")
  colnames(y) <- dataset$features
  transform <- data.table(
    feature = dataset$features,
    missing_fraction = as.numeric(dataset$feature_missingness),
    winsor_lower = limits[, 1],
    winsor_upper = limits[, 2],
    healthy_reference_mean = reference_mean,
    healthy_reference_sd = reference_sd
  )
  list(y = y, parameters = transform)
}

make_molecular_design <- function(meta, discovery_rows, bag_form = "linear") {
  attendance <- as.numeric(as.IDate(meta$Date_attend))
  attendance_z <- (
    attendance - mean(attendance[discovery_rows])
  ) / sd(attendance[discovery_rows])
  day <- as.integer(format(as.IDate(meta$Date_attend), "%j"))
  age_fit <- ns(meta$Age_exact[discovery_rows], df = 4, Boundary.knots = c(40, 70))
  age_basis <- predict(age_fit, newx = meta$Age_exact)
  colnames(age_basis) <- paste0("age_spline_", seq_len(ncol(age_basis)))
  design <- data.frame(
    BAG = as.numeric(meta$BAG),
    BAU = as.numeric(meta$BAU),
    BAG_abs = abs(as.numeric(meta$BAG)),
    Sex_male = as.numeric(meta$Sex_male),
    ethnicity_group = meta$ethnicity_group,
    AssessmentCenter = meta$AssessmentCenter,
    attendance_z = attendance_z,
    season_sin = sin(2 * pi * day / 365.2425),
    season_cos = cos(2 * pi * day / 365.2425),
    target_missing_fraction = meta$target_missing_fraction,
    age_basis,
    check.names = FALSE
  )
  age_terms <- paste(colnames(age_basis), collapse = " + ")
  adjustment <- paste(
    "Sex_male + ethnicity_group + AssessmentCenter + attendance_z +",
    "season_sin + season_cos + target_missing_fraction +", age_terms
  )
  formula <- switch(bag_form,
    linear = as.formula(paste("~ BAG + BAU +", adjustment)),
    absolute = as.formula(paste("~ BAG + BAG_abs + BAU +", adjustment)),
    spline = {
      bag_fit <- ns(meta$BAG[discovery_rows], df = 4)
      bag_basis <- predict(bag_fit, newx = meta$BAG)
      colnames(bag_basis) <- paste0("bag_spline_", seq_len(ncol(bag_basis)))
      design <- cbind(design, bag_basis)
      as.formula(paste(
        "~ BAU +", paste(colnames(bag_basis), collapse = " + "), "+", adjustment
      ))
    },
    stop("Unknown BAG form")
  )
  model.matrix(formula, design)
}

fit_molecular_split <- function(y, rows, design, centres) {
  split_y <- y[rows, , drop = FALSE]
  split_x <- design[rows, , drop = FALSE]
  split_centres <- as.character(centres[rows])
  missing_keys <- vapply(seq_len(ncol(split_y)), function(j) {
    digest(which(!is.finite(split_y[, j])), algo = "xxhash64")
  }, character(1))
  groups <- split(seq_len(ncol(split_y)), missing_keys)
  output <- vector("list", length(groups))
  counter <- 0L

  for (columns in groups) {
    counter <- counter + 1L
    observed <- is.finite(split_y[, columns[[1]]])
    x <- split_x[observed, , drop = FALSE]
    outcomes <- split_y[observed, columns, drop = FALSE]
    cluster <- split_centres[observed]
    qr_x <- qr(x, tol = 1e-9, LAPACK = FALSE)
    kept <- sort(qr_x$pivot[seq_len(qr_x$rank)])
    x <- x[, kept, drop = FALSE]
    if (!all(c("BAG", "BAU") %in% colnames(x))) {
      stop("BAG or BAU is aliased in a molecular model")
    }
    fit <- lm.fit(x, outcomes)
    coefficient <- fit$coefficients
    residual <- fit$residuals
    if (is.null(dim(coefficient))) coefficient <- matrix(coefficient, ncol = 1L)
    if (is.null(dim(residual))) residual <- matrix(residual, ncol = 1L)
    rownames(coefficient) <- colnames(x)
    bread <- tryCatch(
      chol2inv(chol(crossprod(x))),
      error = function(e) solve(crossprod(x))
    )
    target <- match(c("BAG", "BAU"), colnames(x))
    weights <- bread[, target, drop = FALSE]
    groups_observed <- unique(cluster)
    scores <- matrix(0, nrow = 2L, ncol = length(groups_observed) * ncol(outcomes))
    dim(scores) <- c(2L, length(groups_observed), ncol(outcomes))
    for (i in seq_along(groups_observed)) {
      take <- cluster == groups_observed[[i]]
      score <- crossprod(x[take, , drop = FALSE], residual[take, , drop = FALSE])
      scores[, i, ] <- crossprod(weights, score)
    }
    n <- nrow(x)
    rank <- ncol(x)
    g <- length(groups_observed)
    correction <- (g / (g - 1)) * ((n - 1) / (n - rank))
    scores_bag <- matrix(scores[1, , ], nrow = g, ncol = ncol(outcomes))
    scores_bau <- matrix(scores[2, , ], nrow = g, ncol = ncol(outcomes))
    variance_bag <- correction * colSums(scores_bag^2)
    variance_bau <- correction * colSums(scores_bau^2)
    covariance <- correction * colSums(scores_bag * scores_bau)
    beta_bag <- coefficient["BAG", ]
    beta_bau <- coefficient["BAU", ]
    beta_difference <- beta_bau - beta_bag
    se_bag <- sqrt(variance_bag)
    se_bau <- sqrt(variance_bau)
    se_difference <- sqrt(pmax(variance_bau + variance_bag - 2 * covariance, 0))
    degrees <- g - 1L
    critical <- qt(0.975, degrees)
    output[[counter]] <- data.table(
      feature = colnames(outcomes),
      n = n,
      centres = g,
      model_rank = rank,
      cluster_df = degrees,
      beta_BAG = as.numeric(beta_bag),
      se_BAG = se_bag,
      ci_lower_BAG = beta_bag - critical * se_bag,
      ci_upper_BAG = beta_bag + critical * se_bag,
      p_BAG = 2 * pt(-abs(beta_bag / se_bag), degrees),
      beta_BAU = as.numeric(beta_bau),
      se_BAU = se_bau,
      ci_lower_BAU = beta_bau - critical * se_bau,
      ci_upper_BAU = beta_bau + critical * se_bau,
      p_BAU = 2 * pt(-abs(beta_bau / se_bau), degrees),
      cov_BAG_BAU = covariance,
      beta_difference = as.numeric(beta_difference),
      se_difference = se_difference,
      ci_lower_difference = beta_difference - critical * se_difference,
      ci_upper_difference = beta_difference + critical * se_difference,
      p_difference = 2 * pt(-abs(beta_difference / se_difference), degrees)
    )
  }
  result <- rbindlist(output)
  result[match(colnames(y), feature)]
}

fit_single_term_split <- function(y, rows, design, centres, term) {
  split_y <- y[rows, , drop = FALSE]
  split_x <- design[rows, , drop = FALSE]
  split_centres <- as.character(centres[rows])
  missing_keys <- vapply(seq_len(ncol(split_y)), function(j) {
    digest(which(!is.finite(split_y[, j])), algo = "xxhash64")
  }, character(1))
  groups <- split(seq_len(ncol(split_y)), missing_keys)
  output <- vector("list", length(groups))
  counter <- 0L
  for (columns in groups) {
    counter <- counter + 1L
    observed <- is.finite(split_y[, columns[[1]]])
    x <- split_x[observed, , drop = FALSE]
    outcomes <- split_y[observed, columns, drop = FALSE]
    cluster <- split_centres[observed]
    qr_x <- qr(x, tol = 1e-9, LAPACK = FALSE)
    kept <- sort(qr_x$pivot[seq_len(qr_x$rank)])
    x <- x[, kept, drop = FALSE]
    if (!term %in% colnames(x)) stop(paste(term, "is aliased"))
    fit <- lm.fit(x, outcomes)
    coefficient <- fit$coefficients
    residual <- fit$residuals
    if (is.null(dim(coefficient))) coefficient <- matrix(coefficient, ncol = 1L)
    if (is.null(dim(residual))) residual <- matrix(residual, ncol = 1L)
    rownames(coefficient) <- colnames(x)
    bread <- tryCatch(
      chol2inv(chol(crossprod(x))),
      error = function(e) solve(crossprod(x))
    )
    weight <- bread[, match(term, colnames(x))]
    groups_observed <- unique(cluster)
    centre_scores <- matrix(0, length(groups_observed), ncol(outcomes))
    for (i in seq_along(groups_observed)) {
      take <- cluster == groups_observed[[i]]
      score <- crossprod(x[take, , drop = FALSE], residual[take, , drop = FALSE])
      centre_scores[i, ] <- as.numeric(crossprod(weight, score))
    }
    n <- nrow(x)
    rank <- ncol(x)
    g <- length(groups_observed)
    correction <- (g / (g - 1)) * ((n - 1) / (n - rank))
    se <- sqrt(correction * colSums(centre_scores^2))
    beta <- coefficient[term, ]
    degrees <- g - 1L
    output[[counter]] <- data.table(
      feature = colnames(outcomes),
      n = n,
      centres = g,
      model_rank = rank,
      cluster_df = degrees,
      beta = as.numeric(beta),
      se = se,
      p_value = 2 * pt(-abs(beta / se), degrees)
    )
  }
  result <- rbindlist(output)
  result[match(colnames(y), feature)]
}

suffix_columns <- function(x, suffix) {
  columns <- setdiff(names(x), "feature")
  setnames(x, columns, paste0(columns, "_", suffix))
  x
}

combine_molecular_splits <- function(discovery, validation, fdr = 0.05) {
  for (term in c("BAG", "BAU", "difference")) {
    discovery[, (paste0("q_", term)) := p.adjust(get(paste0("p_", term)), "BH")]
  }
  result <- merge(
    suffix_columns(copy(discovery), "discovery"),
    suffix_columns(copy(validation), "validation"),
    by = "feature", sort = FALSE
  )
  for (term in c("BAG", "BAU", "difference")) {
    beta_discovery <- result[[paste0("beta_", term, "_discovery")]]
    beta_validation <- result[[paste0("beta_", term, "_validation")]]
    se_validation <- result[[paste0("se_", term, "_validation")]]
    selected <- result[[paste0("q_", term, "_discovery")]] < fdr
    directional_p <- pt(
      beta_validation / se_validation * sign(beta_discovery),
      df = result$cluster_df_validation,
      lower.tail = FALSE
    )
    directional_q <- rep(NA_real_, nrow(result))
    directional_q[selected] <- p.adjust(directional_p[selected], "BH")
    result[, (paste0("discovery_hit_", term)) := selected]
    result[, (paste0("validation_directional_p_", term)) := directional_p]
    result[, (paste0("validation_directional_q_", term)) := directional_q]
    result[, (paste0("replicated_", term)) :=
      selected & sign(beta_discovery) == sign(beta_validation) &
        !is.na(directional_q) & directional_q < fdr]
  }
  result
}
