suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(survival)
})

standardize <- function(x) {
  x <- as.numeric(x)
  (x - mean(x)) / sd(x)
}

prepare_longitudinal_change <- function(x) {
  x <- copy(x)
  x[, `:=`(
    change_bag = standardize(repeat_bag - baseline_bag),
    change_bau = standardize(repeat_bau - baseline_bau),
    interval_z = standardize(visit_interval_years)
  )]
  x
}

read_predictions <- function(code_root, modality, model = "ngboost") {
  results_root <- Sys.getenv("BAU_RESULTS_ROOT", "/home/dnanexus/bau-work/results")
  path <- file.path(
    results_root, "models", modality, model, "all_participant_predictions.csv"
  )
  x <- fread(path)
  setnames(x, "eid_int", "eid")
  x[, eid := as.integer(eid)]
  if (anyDuplicated(x$eid)) stop("Repeated participant identifier in predictions")
  if (!all(is.finite(x$BAG)) || !all(is.finite(x$BAU))) {
    stop("BAG and BAU must be finite")
  }
  x
}

read_ukb_support <- function(code_root) {
  data_root <- Sys.getenv("BAU_DATA_ROOT", "/home/dnanexus/bau-work/data")
  support_dir <- file.path(data_root, "ukb/support")
  covariates <- as.data.table(read_feather(
    file.path(support_dir, "covariates.feather"),
    col_select = c(
      "eid", "Income", "Employment", "Education", "Current smoker",
      "Daily alcohol intake", "Daily moderate to vigorous physical activity",
      "Daily healthy food", "Frequently processed meat", "Frequently red meat",
      "Frequently salt intake", "Family cardiovascular disease history",
      "Family diabetes history", "Hypotensive drugs", "Insulin", "BMI",
      "Systolic blood pressure", "Diabetes history", "Cholesterol",
      "HDL cholesterol"
    )
  ))
  outcomes <- as.data.table(read_feather(file.path(support_dir, "outcomes_primary.feather")))
  history <- as.data.table(read_feather(file.path(support_dir, "disease_history.feather")))
  h7 <- merge(
    as.data.table(read_feather(file.path(support_dir, "h7_events.feather"))),
    as.data.table(read_feather(file.path(support_dir, "h7_baseline.feather"))),
    by = "eid",
    all = TRUE,
    sort = FALSE
  )
  for (x in list(covariates, outcomes, history, h7)) x[, eid := as.integer(eid)]
  list(covariates = covariates, outcomes = outcomes, history = history, h7 = h7)
}

prepare_survival_data <- function(code_root, modality, model = "ngboost", support = NULL) {
  if (is.null(support)) support <- read_ukb_support(code_root)
  x <- read_predictions(code_root, modality, model)
  x <- merge(x, support$covariates, by = "eid", all.x = TRUE, sort = FALSE)
  x <- merge(x, support$outcomes, by = "eid", all.x = TRUE, sort = FALSE)
  history <- copy(support$history)
  if ("Baseline_unhealthy" %in% names(x)) history[, Baseline_unhealthy := NULL]
  x <- merge(x, history, by = "eid", all.x = TRUE, sort = FALSE)
  x <- merge(x, support$h7, by = "eid", all.x = TRUE, sort = FALSE)
  if (anyDuplicated(x$eid)) stop("Duplicate participant IDs after merging supporting tables")

  mvpa <- as.numeric(as.character(x[["Daily moderate to vigorous physical activity"]]))
  x[, `:=`(
    age_z = standardize(Age),
    mvpa_z = standardize(mvpa),
    smoker = as.numeric(as.character(`Current smoker`)),
    daily_alcohol = as.numeric(as.character(`Daily alcohol intake`)),
    education = factor(as.character(Education)),
    income = factor(as.character(Income)),
    employment = factor(as.character(Employment)),
    centre = factor(as.character(AssessmentCenter))
  )]
  required <- c(
    "Age", "Sex_male", "BAG", "BAU", "age_z", "mvpa_z", "smoker",
    "daily_alcohol", "education", "income", "employment", "centre"
  )
  if (any(!complete.cases(x[, ..required]))) stop("Incomplete analysis covariates")
  x
}

assign_bag_bau_state <- function(x) {
  x <- copy(x)
  x[, state := fifelse(
    BAG < 0 & BAU < 0, "LL",
    fifelse(BAG >= 0 & BAU < 0, "HL", fifelse(BAG < 0, "LH", "HH"))
  )]
  x[, state := factor(state, levels = c("LL", "HL", "LH", "HH"))]
  x
}

make_risk_set <- function(x, outcome) {
  if (outcome == "all-cause mortality") {
    duration <- "death duration"
    event <- "death diagnose"
    risk <- copy(x[is.finite(as.numeric(get(duration))) & as.numeric(get(duration)) > 0])
    risk[, `:=`(
      time_days = as.numeric(get(duration)),
      event = as.integer(get(event) == 1L),
      competing_event = 0L
    )]
  } else if (outcome == "loss of healthspan") {
    duration <- "h7_healthspan_duration"
    event <- "h7_healthspan_diagnose"
    risk <- copy(x[is.finite(as.numeric(get(duration))) & as.numeric(get(duration)) > 0])
    risk[, `:=`(
      time_days = as.numeric(get(duration)),
      event = as.integer(get(event) == 1L),
      competing_event = 0L
    )]
  } else if (outcome == "major disease multimorbidity") {
    duration <- "h7_multimorbidity_duration"
    event <- "h7_multimorbidity_diagnose"
    keep <- is.finite(as.numeric(x[[duration]])) & as.numeric(x[[duration]]) > 0 &
      x$baseline_h7_domain_count < 2L
    risk <- copy(x[keep])
    risk[, competing_event := as.integer(
      get(event) == 0L & `death diagnose` == 1L &
        as.numeric(`death duration`) > 0 &
        as.numeric(`death duration`) <= as.numeric(get(duration))
    )]
    risk[, `:=`(
      time_days = pmin(
        as.numeric(get(duration)),
        fifelse(competing_event == 1L, as.numeric(`death duration`), Inf)
      ),
      event = as.integer(get(event) == 1L),
      baseline_h7 = factor(baseline_h7_domain_count, levels = c(0L, 1L))
    )]
  } else {
    duration <- paste0(outcome, " duration")
    event <- paste0(outcome, " diagnose")
    valid <- is.finite(as.numeric(x[[duration]])) & as.numeric(x[[duration]]) > 0 &
      !is.na(x[[event]])
    risk <- copy(x[valid])
    risk[, target_event := as.integer(get(event) == 1L)]
    risk[, death_observed := as.integer(
      `death diagnose` == 1L & as.numeric(`death duration`) > 0
    )]
    risk[, competing_event := as.integer(
      death_observed == 1L &
        ((target_event == 0L & as.numeric(`death duration`) == as.numeric(get(duration))) |
          as.numeric(`death duration`) < as.numeric(get(duration)))
    )]
    risk[, `:=`(
      time_days = pmin(
        as.numeric(get(duration)),
        fifelse(death_observed == 1L, as.numeric(`death duration`), Inf)
      ),
      event = as.integer(target_event == 1L & competing_event == 0L)
    )]
  }
  risk[, time_years := time_days / 365.25]
  risk[, `:=`(
    centre = droplevels(centre),
    education = droplevels(education),
    income = droplevels(income),
    employment = droplevels(employment)
  )]
  if (!all(is.finite(risk$time_years) & risk$time_years > 0)) {
    stop(paste("Invalid follow-up for", outcome))
  }
  risk
}

make_longitudinal_risk_set <- function(base, outcome, type) {
  x <- copy(base)
  x[, landmark_days := visit_interval_years * 365.25]
  x[, death_observed := as.integer(`death diagnose` == 1L & `death duration` > 0)]
  if (type == "mortality") {
    x[, `:=`(
      target_event = as.integer(`death diagnose` == 1L),
      formal_days = as.numeric(`death duration`),
      competing_death = 0L
    )]
    eligible <- is.finite(x$formal_days) & x$formal_days > 0
  } else if (type == "healthspan") {
    x[, `:=`(
      target_event = as.integer(h7_healthspan_diagnose == 1L),
      formal_days = as.numeric(h7_healthspan_duration),
      competing_death = 0L
    )]
    eligible <- is.finite(x$formal_days) & x$formal_days > 0
  } else if (type == "multimorbidity") {
    x[, target_event := as.integer(h7_multimorbidity_diagnose == 1L)]
    x[, competing_death := as.integer(
      target_event == 0L & death_observed == 1L &
        as.numeric(`death duration`) <= as.numeric(h7_multimorbidity_duration)
    )]
    x[, formal_days := pmin(
      as.numeric(h7_multimorbidity_duration),
      fifelse(competing_death == 1L, as.numeric(`death duration`), Inf)
    )]
    eligible <- is.finite(as.numeric(x$h7_multimorbidity_duration)) &
      as.numeric(x$h7_multimorbidity_duration) > 0 &
      x$baseline_h7_domain_count < 2L
  } else {
    duration <- paste0(outcome, " duration")
    diagnosis <- paste0(outcome, " diagnose")
    x[, target_event := as.integer(get(diagnosis) == 1L)]
    x[, competing_death := as.integer(death_observed == 1L & (
      (target_event == 0L & as.numeric(`death duration`) == as.numeric(get(duration))) |
        as.numeric(`death duration`) < as.numeric(get(duration))
    ))]
    x[, formal_days := pmin(
      as.numeric(get(duration)),
      fifelse(competing_death == 1L, as.numeric(`death duration`), Inf)
    )]
    eligible <- is.finite(as.numeric(x[[duration]])) & as.numeric(x[[duration]]) > 0
  }
  x[, event := as.integer(target_event == 1L & competing_death == 0L)]
  risk <- x[eligible & formal_days > landmark_days]
  risk <- risk[
    as.character(education) != "unknown" & as.character(employment) != "unknown"
  ]
  risk[, `:=`(
    time_days = formal_days - landmark_days,
    education = droplevels(education),
    income = droplevels(income),
    employment = droplevels(employment),
    centre = droplevels(centre)
  )]
  risk
}


joint_cox_formula <- function(multimorbidity = FALSE) {
  terms <- c(
    "BAG", "BAU", "age_z", "Sex_male", "smoker", "daily_alcohol",
    "mvpa_z", "education", "income", "employment", "strata(centre)",
    if (multimorbidity) "strata(baseline_h7)" else character(),
    "cluster(centre)"
  )
  as.formula(paste("Surv(time_years, event) ~", paste(terms, collapse = " + ")))
}

fit_joint_cox <- function(risk, outcome, modality, model = "ngboost") {
  fit <- coxph(
    joint_cox_formula(outcome == "major disease multimorbidity"),
    data = risk,
    ties = "efron",
    x = TRUE,
    y = TRUE
  )
  estimates <- summary(fit)$coefficients
  se_name <- if ("robust se" %in% colnames(estimates)) "robust se" else "se(coef)"
  rows <- rbindlist(lapply(c("BAG", "BAU"), function(term) {
    beta <- estimates[term, "coef"]
    se <- estimates[term, se_name]
    data.table(
      modality = modality,
      model = model,
      outcome = outcome,
      exposure = term,
      n_risk = nrow(risk),
      n_events = sum(risk$event),
      follow_up_q1_years = as.numeric(quantile(risk$time_years, 0.25)),
      median_follow_up_years = median(risk$time_years),
      follow_up_q3_years = as.numeric(quantile(risk$time_years, 0.75)),
      beta = beta,
      se = se,
      hr = exp(beta),
      ci_lower = exp(beta - qnorm(0.975) * se),
      ci_upper = exp(beta + qnorm(0.975) * se),
      p_value = 2 * pnorm(-abs(beta / se))
    )
  }))
  zph <- tryCatch(cox.zph(fit, terms = TRUE, global = TRUE), error = function(e) NULL)
  list(estimates = rows, fit = fit, ph = zph)
}

extract_km_curve <- function(risk, outcome) {
  fit <- survfit(Surv(time_years, event) ~ state, data = risk, conf.type = "log-log")
  ends <- cumsum(fit$strata)
  starts <- c(1L, head(ends, -1L) + 1L)
  rbindlist(lapply(seq_along(ends), function(i) {
    index <- starts[[i]]:ends[[i]]
    state <- sub("^state=", "", names(fit$strata)[[i]])
    data.table(
      outcome = outcome,
      state = state,
      time_years = c(0, fit$time[index]),
      estimate = c(0, 1 - fit$surv[index]),
      lower = c(0, 1 - fit$upper[index]),
      upper = c(0, 1 - fit$lower[index])
    )
  }))
}

extract_aj_curve <- function(risk, outcome) {
  x <- copy(risk)
  x[, status := factor(
    fifelse(event == 1L, "target", fifelse(competing_event == 1L, "death", "censor")),
    levels = c("censor", "target", "death")
  )]
  fit <- survfit(Surv(time_years, status) ~ state, data = x, conf.type = "log-log")
  target <- match("target", fit$states)
  ends <- cumsum(fit$strata)
  starts <- c(1L, head(ends, -1L) + 1L)
  rbindlist(lapply(seq_along(ends), function(i) {
    index <- starts[[i]]:ends[[i]]
    state <- sub("^state=", "", names(fit$strata)[[i]])
    data.table(
      outcome = outcome,
      state = state,
      time_years = c(0, fit$time[index]),
      estimate = c(0, fit$pstate[index, target]),
      lower = c(0, fit$lower[index, target]),
      upper = c(0, fit$upper[index, target])
    )
  }))
}

risk_table <- function(risk, outcome, times = c(0, 5, 10, 15)) {
  rbindlist(lapply(levels(risk$state), function(group) {
    x <- risk[state == group]
    data.table(
      outcome = outcome,
      state = group,
      time_years = times,
      n_risk = vapply(times, function(time) sum(x$time_years >= time), integer(1))
    )
  }))
}
