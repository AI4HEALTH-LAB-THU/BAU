#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(haven)
  library(readr)
  library(tidyr)
})

options(dplyr.summarise.inform = FALSE)

args_all <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
default_root <- normalizePath(file.path(dirname(script_file), "../.."), mustWork = TRUE)
CODE_ROOT <- Sys.getenv("BAU_CODE_ROOT", default_root)
EXTERNAL_ROOT <- Sys.getenv("EXTERNAL_DATA_ROOT", "/mnt/project/external-cohorts")
DATA_ROOT <- Sys.getenv("BAU_DATA_ROOT", "/home/dnanexus/bau-work/data")
OUTPUT_ROOT <- file.path(DATA_ROOT, "external")
NHANES_ROOT <- file.path(EXTERNAL_ROOT, "NHANES")
NHANES_RAW <- file.path(NHANES_ROOT, "RawData")
NHANES_MERGE <- file.path(NHANES_ROOT, "ProcessedData", "Merge_Data")

for (subdir in c("NHANES", "CHARLS", "CHNS", "QC", "config")) {
  dir.create(file.path(OUTPUT_ROOT, subdir), recursive = TRUE, showWarnings = FALSE)
}

safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

pick_num <- function(data, candidates) {
  hit <- candidates[candidates %in% names(data)]
  if (!length(hit)) {
    return(rep(NA_real_, nrow(data)))
  }
  safe_num(data[[hit[[1]]]])
}

clean_code <- function(x, missing_codes = c(7, 9, 77, 99, 777, 999)) {
  z <- safe_num(x)
  z[z %in% missing_codes] <- NA_real_
  z
}

status_yes_no <- function(x, borderline = FALSE) {
  z <- safe_num(x)
  case_when(
    z == 1 ~ "yes",
    z == 2 ~ "no",
    borderline & z == 3 ~ "borderline",
    TRUE ~ "unknown"
  )
}

charls_yes_no <- function(x) {
  z <- trimws(as.character(x))
  case_when(
    grepl("^1", z) ~ "yes",
    grepl("^2", z) ~ "no",
    TRUE ~ "unknown"
  )
}

strict_health_status <- function(data, status_columns, positive_values = c("yes", "borderline")) {
  status_data <- data[status_columns]
  positive <- rowSums(as.data.frame(lapply(status_data, function(x) x %in% positive_values)), na.rm = TRUE) > 0
  all_no <- rowSums(as.data.frame(lapply(status_data, function(x) x == "no")), na.rm = TRUE) == length(status_columns)
  ifelse(positive, "disease", ifelse(all_no, "healthy", "unknown"))
}

row_mean_safe <- function(data, columns) {
  available <- columns[columns %in% names(data)]
  if (length(available) == 0) {
    return(rep(NA_real_, nrow(data)))
  }
  values <- as.data.frame(lapply(data[available], safe_num))
  out <- rowMeans(values, na.rm = TRUE)
  out[!is.finite(out)] <- NA_real_
  out
}

date_from_parts <- function(year, month = 7, day = 1, fallback = as.Date(NA)) {
  y <- safe_num(year)
  m <- safe_num(month)
  d <- safe_num(day)
  valid <- !is.na(y) & y >= 1900 & y <= 2100 & !is.na(m) & m >= 1 & m <= 12 & !is.na(d) & d >= 1 & d <= 31
  out <- rep(as.Date(NA), length(y))
  if (length(fallback) == 1) out[] <- fallback else out <- as.Date(fallback)
  candidate <- rep(as.Date(NA), length(y))
  if (any(valid)) {
    candidate[valid] <- suppressWarnings(as.Date(sprintf("%04d-%02d-%02d", y[valid], m[valid], d[valid]), format = "%Y-%m-%d"))
  }
  out[valid & !is.na(candidate)] <- candidate[valid & !is.na(candidate)]
  out
}

date_from_yyyymmdd <- function(x) {
  z <- sprintf("%08.0f", safe_num(x))
  suppressWarnings(as.Date(z, format = "%Y%m%d"))
}

month_start <- function(year, month) date_from_parts(year, month, 1)

month_end <- function(year, month) {
  start <- month_start(year, month)
  next_month <- as.Date(format(start, "%Y-%m-01"))
  next_month <- as.Date(ifelse(format(next_month, "%m") == "12",
    paste0(as.integer(format(next_month, "%Y")) + 1, "-01-01"),
    sprintf("%04d-%02d-01", as.integer(format(next_month, "%Y")), as.integer(format(next_month, "%m")) + 1)
  ))
  next_month - 1
}

mid_date <- function(lower, upper) {
  lower + floor(as.numeric(upper - lower) / 2)
}

years_between <- function(start, end) as.numeric(as.Date(end) - as.Date(start)) / 365.25

missing_profile <- function(data, columns, cohort, stage) {
  data %>%
    summarise(across(all_of(columns), ~ mean(is.na(.x)))) %>%
    pivot_longer(everything(), names_to = "feature", values_to = "missing_fraction") %>%
    mutate(cohort = cohort, stage = stage, n = nrow(data), .before = 1)
}

apply_transform <- function(x, transform) {
  x <- safe_num(x)
  switch(transform,
    identity = x,
    x_times_10 = x * 10,
    x_div_10 = x / 10,
    bun_to_urea = x * 0.3571,
    calcium_mgdl_to_mmoll = x * 0.2495,
    creatinine_mgdl_to_umoll = x * 88.4,
    bilirubin_mgdl_to_umoll = x * 17.104,
    urate_mgdl_to_umoll = x * 59.48,
    cholesterol_mgdl_to_mmoll = x / 38.67,
    triglycerides_mgdl_to_mmoll = x / 88.57,
    glucose_mgdl_to_mmoll = x / 18.018,
    hba1c_ngsp_to_ifcc = (x - 2.152) / 0.09148,
    stop("Unsupported transform: ", transform)
  )
}

align_features <- function(source, specification) {
  out <- source %>% select(
    cohort, participant_id, age, sex, baseline_date,
    blood_missing_n, blood_feature_n, blood_missing_fraction,
    blood_eligible_lt20, blood_complete_case
  )
  for (i in seq_len(nrow(specification))) {
    out[[specification$ukb_feature[i]]] <- apply_transform(source[[specification$source_field[i]]], specification$transform[i])
  }
  out
}

write_bundle <- function(cohort, candidate, survival, analysis) {
  root <- file.path(OUTPUT_ROOT, cohort)
  write_feather(candidate, file.path(root, "candidate_qc.feather"))
  write_feather(survival, file.path(root, "survival_outcomes.feather"))
  write_feather(analysis, file.path(root, "analysis_wide.feather"))
}


feature_specification <- tribble(
  ~cohort, ~source_field, ~source_unit, ~ukb_feature, ~target_unit, ~transform, ~source_variable,
  "NHANES", "albumin_g_dl", "g/dL", "Albumin", "g/L", "x_times_10", "LBXSAL",
  "NHANES", "alt_u_l", "U/L", "Alanine aminotransferase (ALT)", "U/L", "identity", "LBXSATSI",
  "NHANES", "ast_u_l", "U/L", "Aspartate aminotransferase (AST)", "U/L", "identity", "LBXSASSI",
  "NHANES", "alp_u_l", "U/L", "Alkaline phosphatase", "U/L", "identity", "LBXSAPSI/LBDSAPSI",
  "NHANES", "bun_mg_dl", "mg/dL", "Urea", "mmol/L", "bun_to_urea", "LBXSBU",
  "NHANES", "calcium_mg_dl", "mg/dL", "Calcium", "mmol/L", "calcium_mgdl_to_mmoll", "LBXSCA",
  "NHANES", "creatinine_calibrated_mg_dl", "mg/dL", "Creatinine", "umol/L", "creatinine_mgdl_to_umoll", "LBXSCR/LBDSCR calibrated",
  "NHANES", "ggt_u_l", "U/L", "Gamma glutamyltransferase", "U/L", "identity", "LBXSGTSI",
  "NHANES", "bilirubin_mg_dl", "mg/dL", "Total bilirubin", "umol/L", "bilirubin_mgdl_to_umoll", "LBXSTB/LBDSTB",
  "NHANES", "total_protein_g_dl", "g/dL", "Total protein", "g/L", "x_times_10", "LBXSTP",
  "NHANES", "urate_mg_dl", "mg/dL", "Urate", "umol/L", "urate_mgdl_to_umoll", "LBXSUA",
  "NHANES", "wbc_10e9_l", "10^9/L", "White blood cell (leukocyte) count", "10^9/L", "identity", "LBXWBCSI",
  "NHANES", "lymphocyte_pct", "%", "Lymphocyte percentage", "%", "identity", "LBXLYPCT",
  "NHANES", "neutrophil_pct", "%", "Neutrophill percentage", "%", "identity", "LBXNEPCT",
  "NHANES", "monocyte_pct", "%", "Monocyte percentage", "%", "identity", "LBXMOPCT",
  "NHANES", "hemoglobin_g_dl", "g/dL", "Haemoglobin concentration", "g/dL", "identity", "LBXHGB",
  "NHANES", "mcv_fl", "fL", "Mean corpuscular volume", "fL", "identity", "LBXMCVSI",
  "NHANES", "rdw_pct", "%", "Red blood cell (erythrocyte) distribution width", "%", "identity", "LBXRDW",
  "NHANES", "rbc_10e12_l", "10^12/L", "Red blood cell (erythrocyte) count", "10^12/L", "identity", "LBXRBCSI",
  "NHANES", "platelet_10e9_l", "10^9/L", "Platelet count", "10^9/L", "identity", "LBXPLTSI",
  "NHANES", "mpv_fl", "fL", "Mean platelet (thrombocyte) volume", "fL", "identity", "LBXMPSI",
  "NHANES", "cholesterol_mg_dl", "mg/dL", "Cholesterol", "mmol/L", "cholesterol_mgdl_to_mmoll", "LBXTC",
  "NHANES", "hba1c_ngsp_pct", "%", "Glycated haemoglobin (HbA1c)", "mmol/mol", "hba1c_ngsp_to_ifcc", "LBXGH",
  "NHANES", "hdl_mg_dl", "mg/dL", "HDL cholesterol", "mmol/L", "cholesterol_mgdl_to_mmoll", "LBDHDD/LBDHDL/LBXHDD",
  "CHNS", "urea_mmol_l", "mmol/L", "Urea", "mmol/L", "identity", "UREA",
  "CHNS", "urate_umol_l", "umol/L", "Urate", "umol/L", "identity", "UA",
  "CHNS", "creatinine_umol_l", "umol/L", "Creatinine", "umol/L", "identity", "CRE",
  "CHNS", "apob_g_l", "g/L", "Apolipoprotein B", "g/L", "identity", "APO_B",
  "CHNS", "hdl_mmol_l", "mmol/L", "HDL cholesterol", "mmol/L", "identity", "HDL_C",
  "CHNS", "ldl_mmol_l", "mmol/L", "LDL direct", "mmol/L", "identity", "LDL_C",
  "CHNS", "cholesterol_mmol_l", "mmol/L", "Cholesterol", "mmol/L", "identity", "TC",
  "CHNS", "triglycerides_mmol_l", "mmol/L", "Triglycerides", "mmol/L", "identity", "TG",
  "CHNS", "crp_mg_l", "mg/L", "C-reactive protein", "mg/L", "identity", "HS_CRP",
  "CHNS", "glucose_mmol_l", "mmol/L", "Glucose", "mmol/L", "identity", "GLUCOSE",
  "CHNS", "hba1c_ngsp_pct", "%", "Glycated haemoglobin (HbA1c)", "mmol/mol", "hba1c_ngsp_to_ifcc", "HbA1c",
  "CHNS", "alt_u_l", "U/L", "Alanine aminotransferase (ALT)", "U/L", "identity", "ALT",
  "CHNS", "total_protein_g_l", "g/L", "Total protein", "g/L", "identity", "TP",
  "CHNS", "albumin_g_l", "g/L", "Albumin", "g/L", "identity", "ALB",
  "CHNS", "hemoglobin_g_l", "g/L", "Haemoglobin concentration", "g/dL", "x_div_10", "HGB",
  "CHNS", "wbc_10e9_l", "10^9/L", "White blood cell (leukocyte) count", "10^9/L", "identity", "WBC",
  "CHNS", "rbc_10e12_l", "10^12/L", "Red blood cell (erythrocyte) count", "10^12/L", "identity", "RBC",
  "CHNS", "platelet_10e9_l", "10^9/L", "Platelet count", "10^9/L", "identity", "PLT",
  "CHARLS", "wbc_10e9_l", "10^9/L", "White blood cell (leukocyte) count", "10^9/L", "identity", "qc1_vb002",
  "CHARLS", "hemoglobin_g_dl", "g/dL", "Haemoglobin concentration", "g/dL", "identity", "qc1_vb004",
  "CHARLS", "hematocrit_pct", "%", "Haematocrit percentage", "%", "identity", "qc1_vb005",
  "CHARLS", "mcv_fl", "fL", "Mean corpuscular volume", "fL", "identity", "qc1_vb006",
  "CHARLS", "platelet_10e9_l", "10^9/L", "Platelet count", "10^9/L", "identity", "qc1_vb009",
  "CHARLS", "bun_mg_dl", "mg/dL", "Urea", "mmol/L", "bun_to_urea", "newbun",
  "CHARLS", "glucose_mg_dl", "mg/dL", "Glucose", "mmol/L", "glucose_mgdl_to_mmoll", "newglu",
  "CHARLS", "creatinine_mg_dl", "mg/dL", "Creatinine", "umol/L", "creatinine_mgdl_to_umoll", "newcrea",
  "CHARLS", "cholesterol_mg_dl", "mg/dL", "Cholesterol", "mmol/L", "cholesterol_mgdl_to_mmoll", "newcho",
  "CHARLS", "triglycerides_mg_dl", "mg/dL", "Triglycerides", "mmol/L", "triglycerides_mgdl_to_mmoll", "newtg",
  "CHARLS", "hdl_mg_dl", "mg/dL", "HDL cholesterol", "mmol/L", "cholesterol_mgdl_to_mmoll", "newhdl",
  "CHARLS", "ldl_mg_dl", "mg/dL", "LDL direct", "mmol/L", "cholesterol_mgdl_to_mmoll", "newldl",
  "CHARLS", "crp_mg_l", "mg/L", "C-reactive protein", "mg/L", "identity", "newcrp",
  "CHARLS", "hba1c_ngsp_pct", "%", "Glycated haemoglobin (HbA1c)", "mmol/mol", "hba1c_ngsp_to_ifcc", "newhba1c",
  "CHARLS", "urate_mg_dl", "mg/dL", "Urate", "umol/L", "urate_mgdl_to_umoll", "newua",
  "CHARLS", "cystatin_c_mg_l", "mg/L", "Cystatin C", "mg/L", "identity", "cystatinc"
)

write_tsv(feature_specification, file.path(OUTPUT_ROOT, "config", "feature_specification_all_candidates.tsv"))

sample_flow_rows <- list()
missing_rows <- list()
event_rows <- list()

add_flow <- function(cohort, step, n) {
  sample_flow_rows[[length(sample_flow_rows) + 1]] <<- tibble(cohort = cohort, step = step, n = as.integer(n))
}


message("Building NHANES...")

nh_cycles <- c(
  "1999-2000", "2001-2002", "2003-2004", "2005-2006", "2007-2008",
  "2009-2010", "2011-2012", "2013-2014", "2015-2016", "2017-2018"
)
nh_suffix <- c("", "_B", "_C", "_D", "_E", "_F", "_G", "_H", "_I", "_J")
nh_demo_files <- paste0("DEMO", nh_suffix, ".xpt")
nh_bio_files <- c("LAB18.xpt", "L40_B.xpt", "L40_C.xpt", paste0("BIOPRO", nh_suffix[4:10], ".xpt"))
nh_cbc_files <- c("LAB25.xpt", "L25_B.xpt", "L25_C.xpt", paste0("CBC", nh_suffix[4:10], ".xpt"))
nh_tc_files <- c("LAB13.xpt", "L13_B.xpt", "L13_C.xpt", paste0("TCHOL", nh_suffix[4:10], ".xpt"))
nh_a1c_files <- c("LAB10.xpt", "L10_B.xpt", "L10_C.xpt", paste0("GHB", nh_suffix[4:10], ".xpt"))
nh_hdl_files <- c("LAB13.xpt", "L13_B.xpt", "L13_C.xpt", paste0("HDL", nh_suffix[4:10], ".xpt"))

read_nh_demo <- function(cycle, file) {
  x <- read_xpt(file.path(NHANES_RAW, cycle, "Demographic", file))
  age <- pick_num(x, "RIDAGEYR")
  sex_code <- pick_num(x, "RIAGENDR")
  preg <- pick_num(x, "RIDEXPRG")
  sex <- case_when(sex_code == 1 ~ "male", sex_code == 2 ~ "female", TRUE ~ "unknown")
  pregnancy_status <- case_when(
    sex == "male" ~ "not_applicable",
    sex == "female" & (age < 20 | age > 44) ~ "not_applicable",
    preg == 1 ~ "pregnant",
    preg == 2 ~ "not_pregnant",
    TRUE ~ "unknown"
  )
  tibble(
    cohort = "NHANES", participant_id = as.character(safe_num(x$SEQN)), SEQN = safe_num(x$SEQN),
    cycle = cycle, age = age, sex = sex, baseline_date = as.Date(NA),
    pregnancy_raw = preg, pregnancy_status = pregnancy_status,
    race_ethnicity_raw = pick_num(x, c("RIDRETH3", "RIDRETH1")),
    education_raw = pick_num(x, "DMDEDUC2"), income_poverty_ratio = pick_num(x, "INDFMPIR"),
    citizenship_raw = pick_num(x, "DMDCITZN")
  )
}

read_nh_bio <- function(cycle, file) {
  x <- read_xpt(file.path(NHANES_RAW, cycle, "Laboratory", file))
  raw_creatinine <- pick_num(x, c("LBXSCR", "LBDSCR"))
  if (cycle == "1999-2000") {
    calibrated <- 1.013 * raw_creatinine + 0.147
    method <- "CDC: 1.013*x + 0.147 mg/dL"
  } else if (cycle == "2005-2006") {
    calibrated <- -0.016 + 0.978 * raw_creatinine
    method <- "CDC: -0.016 + 0.978*x mg/dL"
  } else {
    calibrated <- raw_creatinine
    method <- "No cycle-specific correction required"
  }
  tibble(
    SEQN = safe_num(x$SEQN), cycle = cycle,
    albumin_g_dl = pick_num(x, "LBXSAL"), alt_u_l = pick_num(x, "LBXSATSI"),
    ast_u_l = pick_num(x, "LBXSASSI"), alp_u_l = pick_num(x, c("LBXSAPSI", "LBDSAPSI")),
    bun_mg_dl = pick_num(x, "LBXSBU"), calcium_mg_dl = pick_num(x, "LBXSCA"),
    creatinine_raw_mg_dl = raw_creatinine, creatinine_calibrated_mg_dl = calibrated,
    creatinine_calibration_method = method,
    ggt_u_l = pick_num(x, "LBXSGTSI"), bilirubin_mg_dl = pick_num(x, c("LBXSTB", "LBDSTB")),
    total_protein_g_dl = pick_num(x, "LBXSTP"), urate_mg_dl = pick_num(x, "LBXSUA")
  )
}

read_nh_cbc <- function(cycle, file) {
  x <- read_xpt(file.path(NHANES_RAW, cycle, "Laboratory", file))
  tibble(
    SEQN = safe_num(x$SEQN), cycle = cycle,
    wbc_10e9_l = pick_num(x, "LBXWBCSI"), lymphocyte_pct = pick_num(x, "LBXLYPCT"),
    neutrophil_pct = pick_num(x, "LBXNEPCT"), monocyte_pct = pick_num(x, "LBXMOPCT"),
    hemoglobin_g_dl = pick_num(x, "LBXHGB"), mcv_fl = pick_num(x, "LBXMCVSI"),
    rdw_pct = pick_num(x, "LBXRDW"), rbc_10e12_l = pick_num(x, "LBXRBCSI"),
    platelet_10e9_l = pick_num(x, "LBXPLTSI"), mpv_fl = pick_num(x, "LBXMPSI")
  )
}

read_nh_single <- function(cycle, file, type) {
  x <- read_xpt(file.path(NHANES_RAW, cycle, "Laboratory", file))
  if (type == "tc") {
    return(tibble(SEQN = safe_num(x$SEQN), cycle = cycle, cholesterol_mg_dl = pick_num(x, "LBXTC")))
  }
  if (type == "a1c") {
    return(tibble(SEQN = safe_num(x$SEQN), cycle = cycle, hba1c_ngsp_pct = pick_num(x, "LBXGH")))
  }
  tibble(SEQN = safe_num(x$SEQN), cycle = cycle, hdl_mg_dl = pick_num(x, c("LBDHDD", "LBDHDL", "LBXHDD")))
}

nh_demo <- bind_rows(Map(read_nh_demo, nh_cycles, nh_demo_files))
nh_lab <- bind_rows(Map(read_nh_bio, nh_cycles, nh_bio_files)) %>%
  full_join(bind_rows(Map(read_nh_cbc, nh_cycles, nh_cbc_files)), by = c("SEQN", "cycle")) %>%
  full_join(bind_rows(Map(read_nh_single, nh_cycles, nh_tc_files, MoreArgs = list(type = "tc"))), by = c("SEQN", "cycle")) %>%
  full_join(bind_rows(Map(read_nh_single, nh_cycles, nh_a1c_files, MoreArgs = list(type = "a1c"))), by = c("SEQN", "cycle")) %>%
  full_join(bind_rows(Map(read_nh_single, nh_cycles, nh_hdl_files, MoreArgs = list(type = "hdl"))), by = c("SEQN", "cycle"))

nh_questionnaire_files <- tibble(
  cycle = nh_cycles,
  bpq = paste0("BPQ", nh_suffix, ".xpt"), diq = paste0("DIQ", nh_suffix, ".xpt"),
  kiq = c("KIQ.xpt", paste0("KIQ_U", nh_suffix[2:10], ".xpt")),
  mcq = paste0("MCQ", nh_suffix, ".xpt"), smq = paste0("SMQ", nh_suffix, ".xpt"),
  alq = paste0("ALQ", nh_suffix, ".xpt"), bmx = paste0("BMX", nh_suffix, ".xpt"),
  bpx = paste0("BPX", nh_suffix, ".xpt")
)

read_nh_questionnaires <- function(cycle, bpq, diq, kiq, mcq, smq, alq, bmx, bpx) {
  qroot <- file.path(NHANES_RAW, cycle, "Questionnaire")
  eroot <- file.path(NHANES_RAW, cycle, "Examination")
  bp <- read_xpt(file.path(qroot, bpq))
  di <- read_xpt(file.path(qroot, diq))
  ki <- read_xpt(file.path(qroot, kiq))
  mc <- read_xpt(file.path(qroot, mcq))
  sm <- read_xpt(file.path(qroot, smq))
  al <- read_xpt(file.path(qroot, alq))
  bm <- read_xpt(file.path(eroot, bmx))
  bx <- read_xpt(file.path(eroot, bpx))
  out <- tibble(SEQN = safe_num(bp$SEQN), cycle = cycle, hypertension_raw = pick_num(bp, "BPQ020")) %>%
    full_join(tibble(SEQN = safe_num(di$SEQN), cycle = cycle, diabetes_raw = pick_num(di, "DIQ010")), by = c("SEQN", "cycle")) %>%
    full_join(tibble(SEQN = safe_num(ki$SEQN), cycle = cycle, kidney_raw = pick_num(ki, c("KIQ020", "KIQ022"))), by = c("SEQN", "cycle")) %>%
    full_join(tibble(
      SEQN = safe_num(mc$SEQN), cycle = cycle,
      heart_failure_raw = pick_num(mc, "MCQ160B"), coronary_heart_disease_raw = pick_num(mc, "MCQ160C"),
      angina_raw = pick_num(mc, "MCQ160D"), heart_attack_raw = pick_num(mc, "MCQ160E"),
      stroke_raw = pick_num(mc, "MCQ160F"), liver_disease_raw = pick_num(mc, "MCQ160L"),
      cancer_raw = pick_num(mc, "MCQ220")
    ), by = c("SEQN", "cycle")) %>%
    full_join(tibble(
      SEQN = safe_num(sm$SEQN), cycle = cycle,
      ever_smoked_raw = pick_num(sm, "SMQ020"), current_smoking_raw = pick_num(sm, "SMQ040")
    ), by = c("SEQN", "cycle")) %>%
    full_join(tibble(SEQN = safe_num(al$SEQN), cycle = cycle, alcohol_drinks_day = clean_code(pick_num(al, "ALQ130"))), by = c("SEQN", "cycle")) %>%
    full_join(tibble(
      SEQN = safe_num(bm$SEQN), cycle = cycle,
      height_cm = pick_num(bm, "BMXHT"), weight_kg = pick_num(bm, "BMXWT"),
      bmi = pick_num(bm, "BMXBMI"), waist_cm = pick_num(bm, "BMXWAIST")
    ), by = c("SEQN", "cycle")) %>%
    full_join(tibble(
      SEQN = safe_num(bx$SEQN), cycle = cycle,
      systolic_bp = row_mean_safe(bx, c("BPXSY1", "BPXSY2", "BPXSY3", "BPXSY4")),
      diastolic_bp = row_mean_safe(bx, c("BPXDI1", "BPXDI2", "BPXDI3", "BPXDI4")),
      pulse = pick_num(bx, "BPXPLS")
    ), by = c("SEQN", "cycle"))
  out
}

nh_cov <- bind_rows(purrr::pmap(nh_questionnaire_files, read_nh_questionnaires))

nh_mort <- readRDS(file.path(NHANES_MERGE, "Death", "Death_1999_2018_merge.rds")) %>%
  transmute(
    SEQN = safe_num(SEQN), mortality_eligible = safe_num(eligstat) == 1,
    event = safe_num(mortstat) == 1, cause_code = safe_num(ucod_leading),
    followup_months = coalesce(safe_num(permth_exm), safe_num(permth_int)),
    followup_years = followup_months / 12,
    cardiovascular_death = safe_num(cardiovascular_death), cancer_death = safe_num(cancer_death)
  )

nh_source_fields <- feature_specification %>%
  filter(cohort == "NHANES") %>%
  pull(source_field)
nh_candidates <- nh_demo %>%
  filter(age >= 20, age < 80, sex %in% c("male", "female")) %>%
  inner_join(nh_lab, by = c("SEQN", "cycle")) %>%
  left_join(nh_cov, by = c("SEQN", "cycle")) %>%
  inner_join(nh_mort, by = "SEQN") %>%
  mutate(
    hypertension_status = status_yes_no(hypertension_raw),
    diabetes_status = status_yes_no(diabetes_raw, borderline = TRUE),
    kidney_disease_status = status_yes_no(kidney_raw),
    heart_failure_status = status_yes_no(heart_failure_raw),
    coronary_heart_disease_status = status_yes_no(coronary_heart_disease_raw),
    angina_status = status_yes_no(angina_raw), heart_attack_status = status_yes_no(heart_attack_raw),
    stroke_status = status_yes_no(stroke_raw), liver_disease_status = status_yes_no(liver_disease_raw),
    cancer_status = status_yes_no(cancer_raw),
    smoking_status = case_when(
      ever_smoked_raw == 2 ~ "never",
      ever_smoked_raw == 1 & current_smoking_raw %in% c(1, 2) ~ "current",
      ever_smoked_raw == 1 & current_smoking_raw == 3 ~ "former",
      TRUE ~ "unknown"
    ),
    blood_missing_n = rowSums(is.na(across(all_of(nh_source_fields)))),
    blood_feature_n = length(nh_source_fields), blood_missing_fraction = blood_missing_n / blood_feature_n,
    blood_eligible_lt20 = blood_missing_fraction < 0.20,
    blood_complete_case = blood_missing_n == 0,
    pregnancy_confirmed_nonpregnant = pregnancy_status %in% c("not_pregnant", "not_applicable"),
    pregnancy_eligible_main = pregnancy_status != "pregnant",
    followup_positive = !is.na(followup_months) & followup_months > 0
  )

nh_health_cols <- c(
  "hypertension_status", "diabetes_status", "kidney_disease_status", "heart_failure_status",
  "coronary_heart_disease_status", "angina_status", "heart_attack_status", "stroke_status",
  "liver_disease_status", "cancer_status"
)
nh_candidates$baseline_health_status <- strict_health_status(nh_candidates, nh_health_cols)

add_flow("NHANES", "Age 20-79, known sex, laboratory record, mortality link", nrow(nh_candidates))
add_flow("NHANES", "Exclude known pregnant; retain unknown with flag", sum(nh_candidates$pregnancy_eligible_main))
add_flow("NHANES", "Blood missing fraction <20%", sum(nh_candidates$pregnancy_eligible_main & nh_candidates$blood_eligible_lt20))
add_flow("NHANES", "Positive mortality follow-up (final)", sum(nh_candidates$pregnancy_eligible_main & nh_candidates$blood_eligible_lt20 & nh_candidates$followup_positive))

missing_rows[[length(missing_rows) + 1]] <- missing_profile(nh_candidates, nh_source_fields, "NHANES", "age-sex-lab-mortality candidates")

nh_final_source <- nh_candidates %>%
  filter(pregnancy_eligible_main, blood_eligible_lt20, followup_positive) %>%
  select(
    cohort, participant_id, SEQN, cycle, age, sex, baseline_date,
    blood_missing_n, blood_feature_n, blood_missing_fraction, blood_eligible_lt20, blood_complete_case,
    all_of(nh_source_fields), creatinine_raw_mg_dl, creatinine_calibration_method,
    race_ethnicity_raw, education_raw, income_poverty_ratio, citizenship_raw,
    height_cm, weight_kg, bmi, waist_cm, systolic_bp, diastolic_bp, pulse,
    smoking_status, alcohol_drinks_day, pregnancy_status, pregnancy_confirmed_nonpregnant,
    baseline_health_status, all_of(nh_health_cols)
  )

nh_aligned <- align_features(nh_final_source, feature_specification %>% filter(cohort == "NHANES"))
nh_health <- nh_final_source %>%
  select(
    cohort, participant_id, baseline_health_status, pregnancy_status, pregnancy_confirmed_nonpregnant,
    all_of(nh_health_cols)
  )
nh_survival <- nh_candidates %>%
  filter(participant_id %in% nh_final_source$participant_id) %>%
  transmute(
    cohort, participant_id, event,
    time_midpoint_years = followup_years,
    time_lower_years = followup_years, time_upper_years = followup_years,
    time_upper_is_infinite = FALSE, time_precision = "exact_month_count",
    last_alive_date = as.Date(NA), death_lower_date = as.Date(NA), death_upper_date = as.Date(NA),
    exact_death_date = as.Date(NA), censor_date = as.Date(NA),
    survival_eligible = followup_positive, outcome_qc_flag = "ok",
    cause_code, cardiovascular_death, cancer_death
  )
nh_wide <- nh_final_source %>%
  select(
    -all_of(nh_source_fields), -creatinine_raw_mg_dl, -creatinine_calibration_method,
    -all_of(nh_health_cols), -baseline_health_status, -pregnancy_status, -pregnancy_confirmed_nonpregnant
  ) %>%
  left_join(nh_aligned, by = c(
    "cohort", "participant_id", "age", "sex", "baseline_date", "blood_missing_n",
    "blood_feature_n", "blood_missing_fraction", "blood_eligible_lt20", "blood_complete_case"
  )) %>%
  left_join(nh_health, by = c("cohort", "participant_id")) %>%
  left_join(nh_survival, by = c("cohort", "participant_id"))
write_bundle("NHANES", nh_candidates, nh_survival, nh_wide)


message("Building CHARLS...")

charls_root <- file.path(EXTERNAL_ROOT, "CHARLS", "RawData")
charls_names <- jsonlite::fromJSON(file.path(CODE_ROOT, "config/charls_names.json"))
read_charls_file <- function(wave, kind) {
  read_dta(file.path(
    charls_root, charls_names$folders[[as.character(wave)]],
    charls_names$files[[kind]]
  ))
}
charls_csv <- read_csv(file.path(EXTERNAL_ROOT, "CHARLS", "Derived", "baseline_2011.csv"),
  col_types = cols(.default = col_character()), show_col_types = FALSE
)
source_columns <- unlist(charls_names$columns)
matched_columns <- match(source_columns, names(charls_csv))
if (anyNA(matched_columns)) stop("Missing CHARLS baseline columns")
names(charls_csv)[matched_columns] <- names(source_columns)

charls_basic <- read_charls_file(2011, "basic")
charls_health_raw <- read_charls_file(2011, "health")
charls_blood_raw <- read_charls_file(2011, "blood")

charls_id_map <- tibble(
  ID_old = as.character(charls_basic$ID), household_id = as.character(charls_basic$householdID)
) %>%
  mutate(last_two = substr(ID_old, nchar(ID_old) - 1, nchar(ID_old)), participant_id = paste0(household_id, "0", last_two)) %>%
  select(ID_old, participant_id)

charls_baseline_dates <- charls_health_raw %>%
  transmute(
    ID_old = as.character(ID),
    baseline_date_exact = date_from_parts(dc001s1, dc001s2, dc001s3)
  )

charls_demo <- charls_csv %>%
  transmute(
    cohort = "CHARLS", ID_old = as.character(participant_id_source), age = safe_num(age_source),
    sex = case_when(grepl("^1", sex_source) ~ "male", grepl("^2", sex_source) ~ "female", TRUE ~ "unknown"),
    survey_year = safe_num(survey_year_source), survey_month = safe_num(survey_month_source), education_raw = education_source,
    height_cm = safe_num(height_cm_source), weight_kg = safe_num(weight_kg_source), bmi = safe_num(BMI),
    waist_cm = safe_num(waist_cm_source),
    systolic_bp = row_mean_safe(charls_csv, c("systolic_bp_1", "systolic_bp_2", "systolic_bp_3")),
    diastolic_bp = row_mean_safe(charls_csv, c("diastolic_bp_1", "diastolic_bp_2", "diastolic_bp_3")),
    pulse = row_mean_safe(charls_csv, c("pulse_1", "pulse_2", "pulse_3")),
    smoking_raw = smoking_source, alcohol_ever_raw = alcohol_source,
    hypertension_raw = hypertension_source, dyslipidemia_raw = dyslipidemia_source,
    diabetes_raw = diabetes_source, cancer_raw = cancer_source,
    lung_disease_raw = lung_disease_source, liver_disease_raw = liver_disease_source,
    heart_problem_raw = heart_problem_source, stroke_raw = stroke_source, kidney_disease_raw = kidney_disease_source,
    asthma_raw = asthma_source
  ) %>%
  inner_join(charls_id_map, by = "ID_old") %>%
  left_join(charls_baseline_dates, by = "ID_old") %>%
  mutate(
    baseline_date = coalesce(baseline_date_exact, date_from_parts(survey_year, survey_month, 15, as.Date("2011-07-01"))),
    participant_id = as.character(participant_id),
    hypertension_status = charls_yes_no(hypertension_raw), dyslipidemia_status = charls_yes_no(dyslipidemia_raw),
    diabetes_status = charls_yes_no(diabetes_raw), cancer_status = charls_yes_no(cancer_raw),
    lung_disease_status = charls_yes_no(lung_disease_raw), liver_disease_status = charls_yes_no(liver_disease_raw),
    heart_problem_status = charls_yes_no(heart_problem_raw), stroke_status = charls_yes_no(stroke_raw),
    kidney_disease_status = charls_yes_no(kidney_disease_raw), asthma_status = charls_yes_no(asthma_raw),
    smoking_status = case_when(grepl("^1", smoking_raw) ~ "ever", grepl("^2", smoking_raw) ~ "never", TRUE ~ "unknown"),
    alcohol_ever_status = case_when(grepl("^1", alcohol_ever_raw) ~ "yes", grepl("^2", alcohol_ever_raw) ~ "no", TRUE ~ "unknown")
  )

charls_blood <- charls_blood_raw %>%
  transmute(
    ID_old = as.character(ID),
    wbc_10e9_l = safe_num(qc1_vb002), hemoglobin_g_dl = safe_num(qc1_vb004),
    hematocrit_pct = safe_num(qc1_vb005), mcv_fl = safe_num(qc1_vb006), platelet_10e9_l = safe_num(qc1_vb009),
    bun_mg_dl = safe_num(newbun), glucose_mg_dl = safe_num(newglu), creatinine_mg_dl = safe_num(newcrea),
    cholesterol_mg_dl = safe_num(newcho), triglycerides_mg_dl = safe_num(newtg),
    hdl_mg_dl = safe_num(newhdl), ldl_mg_dl = safe_num(newldl), crp_mg_l = safe_num(newcrp),
    hba1c_ngsp_pct = safe_num(newhba1c), urate_mg_dl = safe_num(newua), cystatin_c_mg_l = safe_num(cystatinc)
  )

charls_specification_initial <- feature_specification %>% filter(cohort == "CHARLS")
charls_initial_fields <- charls_specification_initial$source_field
charls_candidates <- charls_demo %>%
  filter(age >= 35, age <= 85, sex %in% c("male", "female")) %>%
  inner_join(charls_blood, by = "ID_old") %>%
  mutate(
    blood_missing_n_initial = rowSums(is.na(across(all_of(charls_initial_fields)))),
    blood_feature_n_initial = length(charls_initial_fields),
    blood_missing_fraction_initial = blood_missing_n_initial / blood_feature_n_initial,
    blood_eligible_initial_lt20 = blood_missing_fraction_initial < 0.20
  )

charls_specification <- charls_specification_initial
charls_fields <- charls_specification$source_field
charls_candidates <- charls_candidates %>%
  mutate(
    blood_missing_n = rowSums(is.na(across(all_of(charls_fields)))),
    blood_feature_n = length(charls_fields), blood_missing_fraction = blood_missing_n / blood_feature_n,
    blood_eligible_lt20 = blood_missing_fraction < 0.20,
    blood_complete_case = blood_missing_n == 0
  )

write_tsv(charls_specification, file.path(OUTPUT_ROOT, "config", "charls_feature_specification_final.tsv"))

charls_health_cols <- c(
  "hypertension_status", "dyslipidemia_status", "diabetes_status", "cancer_status",
  "lung_disease_status", "liver_disease_status", "heart_problem_status", "stroke_status",
  "kidney_disease_status", "asthma_status"
)
charls_candidates$baseline_health_status <- strict_health_status(charls_candidates, charls_health_cols)

read_charls_status <- function(wave) {
  x <- read_charls_file(wave, "status")
  tibble(
    participant_id = as.character(x$ID), wave = as.numeric(wave), died = safe_num(x$died),
    report_date = date_from_parts(x$iyear, x$imonth, 15, as.Date(paste0(wave, "-07-01")))
  )
}
charls_status_long <- bind_rows(lapply(c(2015, 2018, 2020), read_charls_status))

charls_2013_health <- read_charls_file(2013, "health") %>%
  transmute(participant_id = as.character(ID), alive_date = date_from_parts(dc001s1, dc001s2, dc001s3, as.Date("2013-07-01")))
charls_alive_long <- bind_rows(
  charls_candidates %>% transmute(participant_id, alive_date = baseline_date, alive_wave = 2011),
  charls_2013_health %>% mutate(alive_wave = 2013),
  charls_status_long %>% filter(died == 0) %>% transmute(participant_id, alive_date = report_date, alive_wave = wave)
)

charls_exit_2013 <- read_charls_file(2013, "exit_2013") %>%
  transmute(
    participant_id = as.character(ID), death_year = safe_num(exb001_1), death_month = safe_num(exb001_2),
    death_day = NA_real_, death_source = "2013_exit", report_wave = 2013,
    report_date = as.Date("2013-12-31")
  )
charls_exit_2020 <- read_charls_file(2020, "exit_2020") %>%
  transmute(
    participant_id = as.character(ID), death_year = safe_num(exb001_1), death_month = safe_num(exb001_2),
    death_day = safe_num(exb001_3), death_source = "2020_exit", report_wave = 2020,
    report_date = as.Date("2020-12-31")
  )
charls_status_deaths <- charls_status_long %>%
  filter(died == 1) %>%
  transmute(participant_id,
    death_year = NA_real_, death_month = NA_real_, death_day = NA_real_,
    death_source = paste0("status_", wave), report_wave = wave, report_date
  )

charls_death_long <- bind_rows(charls_exit_2013, charls_exit_2020, charls_status_deaths) %>%
  arrange(participant_id, report_date) %>%
  group_by(participant_id) %>%
  summarise(
    first_report_date = min(report_date, na.rm = TRUE), first_report_wave = min(report_wave, na.rm = TRUE),
    exact_year = {
      z <- death_year[!is.na(death_year) & death_year >= 2011 & death_year <= 2020]
      if (length(z)) z[1] else NA_real_
    },
    exact_month = {
      i <- which(!is.na(death_year) & death_year >= 2011 & death_year <= 2020 & !is.na(death_month) & death_month >= 1 & death_month <= 12)
      if (length(i)) death_month[i[1]] else NA_real_
    },
    exact_day = {
      i <- which(!is.na(death_year) & death_year >= 2011 & death_year <= 2020 & !is.na(death_month) & death_month >= 1 & death_month <= 12 & !is.na(death_day) & death_day >= 1 & death_day <= 31)
      if (length(i)) death_day[i[1]] else NA_real_
    },
    exact_source = {
      i <- which(!is.na(death_year) & death_year >= 2011 & death_year <= 2020)
      if (length(i)) death_source[i[1]] else NA_character_
    },
    .groups = "drop"
  )

build_charls_survival_row <- function(participant_id, baseline_date) {
  d <- charls_death_long %>% filter(.data$participant_id == !!participant_id)
  alive <- charls_alive_long %>% filter(.data$participant_id == !!participant_id, alive_date >= baseline_date)
  if (nrow(d) == 0) {
    censor <- if (nrow(alive)) max(alive$alive_date, na.rm = TRUE) else baseline_date
    time <- years_between(baseline_date, censor)
    return(tibble(
      participant_id,
      event = FALSE, last_alive_date = censor, death_lower_date = as.Date(NA),
      death_upper_date = as.Date(NA), exact_death_date = as.Date(NA), censor_date = censor,
      time_lower_years = time, time_upper_years = NA_real_, time_upper_is_infinite = TRUE,
      time_midpoint_years = time, time_upper_bound_years = time, time_precision = "right_censored",
      outcome_qc_flag = ifelse(time > 0, "ok", "no_postbaseline_followup"), survival_eligible = time > 0,
      first_death_report_wave = NA_real_
    ))
  }
  report_date <- d$first_report_date[1]
  alive_before <- alive %>% filter(alive_date < report_date)
  last_alive <- if (nrow(alive_before)) max(alive_before$alive_date, na.rm = TRUE) else baseline_date
  if (!is.na(d$exact_year[1]) && !is.na(d$exact_month[1])) {
    if (!is.na(d$exact_day[1])) {
      lower <- upper <- date_from_parts(d$exact_year[1], d$exact_month[1], d$exact_day[1])
      precision <- "exact_day"
    } else {
      lower <- month_start(d$exact_year[1], d$exact_month[1])
      upper <- month_end(d$exact_year[1], d$exact_month[1])
      precision <- "month_interval"
    }
    lower <- max(lower, last_alive)
  } else {
    lower <- last_alive
    upper <- report_date
    precision <- "wave_interval"
  }
  inconsistent <- is.na(lower) || is.na(upper) || upper <= baseline_date || lower > upper
  midpoint <- if (!inconsistent) mid_date(lower, upper) else as.Date(NA)
  tibble(
    participant_id,
    event = TRUE, last_alive_date = last_alive, death_lower_date = lower,
    death_upper_date = upper, exact_death_date = ifelse(precision == "exact_day", as.character(lower), NA_character_) %>% as.Date(),
    censor_date = as.Date(NA),
    time_lower_years = ifelse(inconsistent, NA_real_, years_between(baseline_date, lower)),
    time_upper_years = ifelse(inconsistent, NA_real_, years_between(baseline_date, upper)),
    time_upper_is_infinite = FALSE,
    time_midpoint_years = ifelse(inconsistent, NA_real_, years_between(baseline_date, midpoint)),
    time_upper_bound_years = ifelse(inconsistent, NA_real_, years_between(baseline_date, upper)),
    time_precision = precision, outcome_qc_flag = ifelse(inconsistent, "inconsistent_event_interval", "ok"),
    survival_eligible = !inconsistent && years_between(baseline_date, midpoint) > 0,
    first_death_report_wave = d$first_report_wave[1]
  )
}

charls_final_source <- charls_candidates %>%
  filter(blood_eligible_lt20) %>%
  select(
    cohort, participant_id, ID_old, age, sex, baseline_date,
    blood_missing_n, blood_feature_n, blood_missing_fraction, blood_eligible_lt20, blood_complete_case,
    all_of(charls_fields), education_raw, height_cm, weight_kg, bmi, waist_cm,
    systolic_bp, diastolic_bp, pulse, smoking_status, alcohol_ever_status,
    baseline_health_status, all_of(charls_health_cols)
  )
charls_survival <- bind_rows(Map(build_charls_survival_row, charls_final_source$participant_id, charls_final_source$baseline_date)) %>%
  mutate(cohort = "CHARLS", .before = 1)

charls_aligned <- align_features(charls_final_source, charls_specification)
charls_health <- charls_final_source %>% select(cohort, participant_id, baseline_health_status, all_of(charls_health_cols))
charls_wide <- charls_final_source %>%
  select(-all_of(charls_fields), -all_of(charls_health_cols), -baseline_health_status) %>%
  left_join(charls_aligned, by = c(
    "cohort", "participant_id", "age", "sex", "baseline_date", "blood_missing_n",
    "blood_feature_n", "blood_missing_fraction", "blood_eligible_lt20", "blood_complete_case"
  )) %>%
  left_join(charls_health, by = c("cohort", "participant_id")) %>%
  left_join(charls_survival, by = c("cohort", "participant_id"))
add_flow("CHARLS", "Age 35-85, known sex, raw blood DTA", nrow(charls_candidates))
add_flow("CHARLS", "Initial 16-feature missing fraction <20%", sum(charls_candidates$blood_eligible_initial_lt20))
add_flow("CHARLS", paste0("Final ", length(charls_fields), "-feature missing fraction <20%"), nrow(charls_final_source))
add_flow("CHARLS", "Complete cases in final retained panel", sum(charls_wide$blood_complete_case, na.rm = TRUE))
add_flow("CHARLS", "Cox-eligible after outcome reconstruction", sum(charls_wide$survival_eligible, na.rm = TRUE))
missing_rows[[length(missing_rows) + 1]] <- missing_profile(charls_candidates, charls_initial_fields, "CHARLS", "initial 16-feature candidates")
missing_rows[[length(missing_rows) + 1]] <- missing_profile(charls_final_source, charls_fields, "CHARLS", "final retained panel")
event_rows[[length(event_rows) + 1]] <- charls_survival %>%
  count(event, time_precision, outcome_qc_flag, name = "n") %>%
  mutate(cohort = "CHARLS", .before = 1)

write_bundle("CHARLS", charls_candidates, charls_survival, charls_wide)


message("Building CHNS...")

chns_root <- file.path(EXTERNAL_ROOT, "CHNS", "RawData")
chns_bio_raw <- read_sas(file.path(chns_root, "Biomarker", "biomarker_09.sas7bdat"))
chns_master_raw <- read_sas(file.path(chns_root, "Individual_Data", "mast_pub_12.sas7bdat"))
chns_exam_all <- read_sas(file.path(chns_root, "Individual_Data", "pexam_00.sas7bdat"))
chns_edu_raw <- read_sas(file.path(chns_root, "Individual_Data", "educ_12.sas7bdat"))
chns_job_raw <- read_sas(file.path(chns_root, "Individual_Data", "jobs_12.sas7bdat"))
chns_health_raw <- read_sas(file.path(chns_root, "Individual_Data", "hlth_12.sas7bdat"))
chns_ins_raw <- read_sas(file.path(chns_root, "Individual_Data", "ins_12.sas7bdat"))
chns_surveys_raw <- read_sas(file.path(chns_root, "Individual_Data", "surveys_pub_12.sas7bdat"))

chns_blood <- chns_bio_raw %>%
  transmute(
    participant_id = as.character(IDind), urea_mmol_l = safe_num(UREA), urate_umol_l = safe_num(UA),
    creatinine_umol_l = safe_num(CRE), apob_g_l = safe_num(APO_B), hdl_mmol_l = safe_num(HDL_C),
    ldl_mmol_l = safe_num(LDL_C), cholesterol_mmol_l = safe_num(TC), triglycerides_mmol_l = safe_num(TG),
    crp_mg_l = safe_num(HS_CRP), glucose_mmol_l = safe_num(GLUCOSE), hba1c_ngsp_pct = safe_num(HbA1c),
    alt_u_l = safe_num(ALT), total_protein_g_l = safe_num(TP), albumin_g_l = safe_num(ALB),
    hemoglobin_g_l = safe_num(HGB), wbc_10e9_l = safe_num(WBC), rbc_10e12_l = safe_num(RBC),
    platelet_10e9_l = safe_num(PLT)
  )

chns_master <- chns_master_raw %>%
  transmute(
    participant_id = as.character(Idind), sex = case_when(GENDER == 1 ~ "male", GENDER == 2 ~ "female", TRUE ~ "unknown"),
    birth_year = safe_num(WEST_DOB_Y), ethnicity_raw = safe_num(NATIONALITY),
    death_year_raw = safe_num(DOD_Y), death_report_year_raw = safe_num(DOD_RPT)
  ) %>%
  mutate(age = 2009 - birth_year)

chns_exam_2009 <- chns_exam_all %>%
  filter(WAVE == 2009) %>%
  transmute(
    participant_id = as.character(IDind), exam_date = date_from_yyyymmdd(T7_PE), urban_rural_raw = safe_num(T2),
    height_cm = safe_num(HEIGHT), weight_kg = safe_num(WEIGHT),
    bmi = ifelse(height_cm > 0 & weight_kg > 0, weight_kg / (height_cm / 100)^2, NA_real_),
    waist_cm = safe_num(U10), hip_cm = safe_num(U9),
    systolic_bp = row_mean_safe(pick(everything()), c("SYSTOL1", "SYSTOL2", "SYSTOL3")),
    diastolic_bp = row_mean_safe(pick(everything()), c("DIASTOL1", "DIASTOL2", "DIASTOL3")),
    hypertension_raw = safe_num(U22), diabetes_raw = safe_num(U24A), mi_raw = safe_num(U24J),
    stroke_raw = safe_num(U24L), asthma_raw = safe_num(U24Q), cancer_raw_2009 = safe_num(U24W),
    antihypertensive_raw = safe_num(U24), pregnancy_raw = safe_num(U56),
    ever_smoked_raw = safe_num(U25), alcohol_frequency_raw = safe_num(U41)
  )

chns_edu <- chns_edu_raw %>%
  filter(WAVE == 2009) %>%
  transmute(participant_id = as.character(IDind), education_raw = safe_num(A12))
chns_job <- chns_job_raw %>%
  filter(wave == 2009) %>%
  transmute(participant_id = as.character(IDind), working_raw = safe_num(B2), occupation_raw = safe_num(B4))
chns_ins <- chns_ins_raw %>%
  filter(wave == 2009) %>%
  transmute(participant_id = as.character(IDind), medical_insurance_raw = safe_num(M1))
chns_health <- chns_health_raw %>%
  filter(wave == 2009) %>%
  transmute(
    participant_id = as.character(IDind), sick_4w_raw = safe_num(M23), illness_severity_raw = safe_num(M25),
    fever_cough_4w_raw = safe_num(M24B_1), diarrhea_stomachache_4w_raw = safe_num(M24B_2)
  )

chns_cancer_2011 <- chns_exam_all %>%
  filter(WAVE == 2011) %>%
  transmute(participant_id = as.character(IDind), cancer_2011_raw = safe_num(U24W), cancer_age_first_2011 = safe_num(U24X))

chns_surveys <- chns_surveys_raw %>%
  transmute(participant_id = as.character(Idind), wave = safe_num(wave), survey_date = date_from_yyyymmdd(t7)) %>%
  filter(wave >= 2009)
chns_baseline_survey <- chns_surveys %>%
  filter(wave == 2009) %>%
  group_by(participant_id) %>%
  summarise(baseline_survey_date = min(survey_date, na.rm = TRUE), .groups = "drop") %>%
  mutate(baseline_survey_date = ifelse(is.infinite(as.numeric(baseline_survey_date)), NA, as.numeric(baseline_survey_date)) %>% as.Date(origin = "1970-01-01"))

chns_fields <- feature_specification %>%
  filter(cohort == "CHNS") %>%
  pull(source_field)
chns_candidates <- chns_master %>%
  filter(age >= 20, age <= 80, sex %in% c("male", "female")) %>%
  inner_join(chns_blood, by = "participant_id") %>%
  left_join(chns_exam_2009, by = "participant_id") %>%
  left_join(chns_edu, by = "participant_id") %>%
  left_join(chns_job, by = "participant_id") %>%
  left_join(chns_ins, by = "participant_id") %>%
  left_join(chns_health, by = "participant_id") %>%
  left_join(chns_cancer_2011, by = "participant_id") %>%
  left_join(chns_baseline_survey, by = "participant_id") %>%
  mutate(
    cohort = "CHNS", baseline_date = coalesce(baseline_survey_date, exam_date, as.Date("2009-07-01")),
    blood_missing_n = rowSums(is.na(across(all_of(chns_fields)))), blood_feature_n = length(chns_fields),
    blood_missing_fraction = blood_missing_n / blood_feature_n, blood_eligible_lt20 = blood_missing_fraction < 0.20,
    blood_complete_case = blood_missing_n == 0,
    hypertension_status = status_yes_no(hypertension_raw), diabetes_status = status_yes_no(diabetes_raw),
    mi_status = status_yes_no(mi_raw), stroke_status = status_yes_no(stroke_raw), asthma_status = status_yes_no(asthma_raw),
    cancer_status = case_when(
      cancer_2011_raw == 1 & !is.na(cancer_age_first_2011) & cancer_age_first_2011 >= 0 & cancer_age_first_2011 < age ~ "yes_retrospective_prebaseline",
      cancer_2011_raw == 1 & !is.na(cancer_age_first_2011) & cancer_age_first_2011 >= 0 & cancer_age_first_2011 == age ~ "yes_retrospective_baseline_timing_ambiguous",
      cancer_2011_raw == 1 & !is.na(cancer_age_first_2011) & cancer_age_first_2011 > age ~ "yes_retrospective_postbaseline",
      cancer_2011_raw == 1 ~ "yes_retrospective_diagnosis_time_unknown",
      TRUE ~ "unknown_not_collected_in_2009"
    ),
    pregnancy_status = case_when(
      sex == "male" ~ "not_applicable", pregnancy_raw == 1 ~ "pregnant", pregnancy_raw == 0 ~ "not_pregnant", TRUE ~ "unknown"
    ),
    smoking_status = case_when(ever_smoked_raw == 1 ~ "ever", ever_smoked_raw == 0 ~ "never", TRUE ~ "unknown"),
    alcohol_frequency = case_when(
      alcohol_frequency_raw %in% c(1, 2) ~ "frequent", alcohol_frequency_raw == 3 ~ "monthly",
      alcohol_frequency_raw %in% c(4, 5) ~ "rare_or_never", TRUE ~ "unknown"
    )
  )

chns_measured_health_cols <- c("hypertension_status", "diabetes_status", "mi_status", "stroke_status", "asthma_status")
chns_candidates$health_status_excluding_cancer <- strict_health_status(chns_candidates, chns_measured_health_cols)
chns_candidates$baseline_health_status <- ifelse(
  chns_candidates$health_status_excluding_cancer == "disease" | chns_candidates$cancer_status == "yes_retrospective_prebaseline",
  "disease", "unknown_due_unmeasured_baseline_cancer"
)

build_chns_survival_row <- function(participant_id, baseline_date, death_year_raw, death_report_year_raw) {
  alive <- chns_surveys %>% filter(.data$participant_id == !!participant_id, survey_date >= baseline_date)
  report_year <- safe_num(death_report_year_raw)
  event <- !is.na(report_year) && report_year >= 2009 && report_year <= 2015
  prebaseline_death_record <- (!is.na(report_year) && report_year < 2009) ||
    (!is.na(death_year_raw) && death_year_raw != 20 && death_year_raw < 2009)
  if (!event) {
    censor <- if (nrow(alive)) max(alive$survey_date, na.rm = TRUE) else baseline_date
    time <- years_between(baseline_date, censor)
    return(tibble(
      participant_id,
      event = ifelse(prebaseline_death_record, NA, FALSE), last_alive_date = censor, death_lower_date = as.Date(NA), death_upper_date = as.Date(NA),
      exact_death_date = as.Date(NA), censor_date = censor, time_lower_years = time, time_upper_years = NA_real_,
      time_upper_is_infinite = TRUE, time_midpoint_years = time, time_upper_bound_years = time,
      time_precision = ifelse(prebaseline_death_record, "invalid_prebaseline_record", "right_censored"),
      outcome_qc_flag = ifelse(prebaseline_death_record, "death_record_prebaseline", ifelse(time > 0, "ok", "no_postbaseline_followup")),
      survival_eligible = !prebaseline_death_record && time > 0, death_year_usable = FALSE, death_report_year = report_year
    ))
  }
  report_rows <- alive %>% filter(wave == report_year)
  report_date <- if (nrow(report_rows) && any(!is.na(report_rows$survey_date))) min(report_rows$survey_date, na.rm = TRUE) else as.Date(paste0(report_year, "-12-31"))
  alive_before <- alive %>% filter(wave < report_year, survey_date < report_date)
  last_alive <- if (nrow(alive_before)) max(alive_before$survey_date, na.rm = TRUE) else baseline_date
  valid_year <- !is.na(death_year_raw) && death_year_raw >= 2009 && death_year_raw <= report_year
  prebaseline_record <- !is.na(death_year_raw) && death_year_raw != 20 && death_year_raw < 2009
  if (valid_year) {
    lower <- max(as.Date(paste0(death_year_raw, "-01-01")), last_alive)
    upper <- min(as.Date(paste0(death_year_raw, "-12-31")), report_date)
    precision <- "year_interval"
  } else {
    lower <- last_alive
    upper <- report_date
    precision <- "report_wave_interval"
  }
  inconsistent <- prebaseline_record || is.na(lower) || is.na(upper) || upper <= baseline_date || lower > upper
  midpoint <- if (!inconsistent) mid_date(lower, upper) else as.Date(NA)
  tibble(
    participant_id,
    event = TRUE, last_alive_date = last_alive, death_lower_date = lower, death_upper_date = upper,
    exact_death_date = as.Date(NA), censor_date = as.Date(NA),
    time_lower_years = ifelse(inconsistent, NA_real_, years_between(baseline_date, lower)),
    time_upper_years = ifelse(inconsistent, NA_real_, years_between(baseline_date, upper)),
    time_upper_is_infinite = FALSE,
    time_midpoint_years = ifelse(inconsistent, NA_real_, years_between(baseline_date, midpoint)),
    time_upper_bound_years = ifelse(inconsistent, NA_real_, years_between(baseline_date, upper)),
    time_precision = precision,
    outcome_qc_flag = ifelse(prebaseline_record, "death_year_prebaseline", ifelse(inconsistent, "inconsistent_event_interval", "ok")),
    survival_eligible = !inconsistent && years_between(baseline_date, midpoint) > 0,
    death_year_usable = valid_year, death_report_year = report_year
  )
}

chns_final_source <- chns_candidates %>%
  filter(blood_eligible_lt20) %>%
  select(
    cohort, participant_id, age, sex, baseline_date, birth_year,
    blood_missing_n, blood_feature_n, blood_missing_fraction, blood_eligible_lt20, blood_complete_case,
    all_of(chns_fields), ethnicity_raw, urban_rural_raw, working_raw, occupation_raw, education_raw,
    medical_insurance_raw, height_cm, weight_kg, bmi, waist_cm, hip_cm, systolic_bp, diastolic_bp,
    smoking_status, alcohol_frequency, pregnancy_status, baseline_health_status, health_status_excluding_cancer,
    all_of(chns_measured_health_cols), cancer_status, cancer_2011_raw, cancer_age_first_2011,
    death_year_raw, death_report_year_raw
  )

chns_survival <- bind_rows(Map(
  build_chns_survival_row, chns_final_source$participant_id, chns_final_source$baseline_date,
  chns_final_source$death_year_raw, chns_final_source$death_report_year_raw
)) %>% mutate(cohort = "CHNS", .before = 1)

chns_aligned <- align_features(chns_final_source, feature_specification %>% filter(cohort == "CHNS"))
chns_health_labels <- chns_final_source %>%
  select(
    cohort, participant_id, baseline_health_status, health_status_excluding_cancer,
    all_of(chns_measured_health_cols), cancer_status, pregnancy_status
  )
chns_wide <- chns_final_source %>%
  select(
    -all_of(chns_fields), -all_of(chns_measured_health_cols), -baseline_health_status,
    -health_status_excluding_cancer, -cancer_status, -pregnancy_status,
    -death_year_raw, -death_report_year_raw
  ) %>%
  left_join(chns_aligned, by = c(
    "cohort", "participant_id", "age", "sex", "baseline_date", "blood_missing_n",
    "blood_feature_n", "blood_missing_fraction", "blood_eligible_lt20", "blood_complete_case"
  )) %>%
  left_join(chns_health_labels, by = c("cohort", "participant_id")) %>%
  left_join(chns_survival, by = c("cohort", "participant_id"))
add_flow("CHNS", "Age 20-80, known sex, 2009 biomarker record", nrow(chns_candidates))
add_flow("CHNS", "18-feature missing fraction <20%", nrow(chns_final_source))
add_flow("CHNS", "Complete cases in retained panel", sum(chns_wide$blood_complete_case, na.rm = TRUE))
add_flow("CHNS", "Cox-eligible after outcome reconstruction", sum(chns_wide$survival_eligible, na.rm = TRUE))
missing_rows[[length(missing_rows) + 1]] <- missing_profile(chns_candidates, chns_fields, "CHNS", "18-feature candidates")
event_rows[[length(event_rows) + 1]] <- chns_survival %>%
  count(event, time_precision, outcome_qc_flag, death_year_usable, name = "n") %>%
  mutate(cohort = "CHNS", .before = 1)

write_bundle("CHNS", chns_candidates, chns_survival, chns_wide)


sample_flow <- bind_rows(sample_flow_rows)
feature_missingness <- bind_rows(missing_rows)
event_qc <- bind_rows(event_rows)

write_tsv(sample_flow, file.path(OUTPUT_ROOT, "QC", "sample_flow.tsv"))
write_tsv(feature_missingness, file.path(OUTPUT_ROOT, "QC", "feature_missingness.tsv"))
write_tsv(event_qc, file.path(OUTPUT_ROOT, "QC", "survival_event_qc.tsv"))

message("External cohorts written to ", OUTPUT_ROOT)
print(sample_flow, n = Inf)
