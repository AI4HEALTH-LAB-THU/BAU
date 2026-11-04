suppressPackageStartupMessages({
  library(data.table)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
default_root <- normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
project_root <- Sys.getenv("BAU_CODE_ROOT", default_root)
results_root <- Sys.getenv("BAU_RESULTS_ROOT", "/home/dnanexus/bau-work/results")
model_root <- file.path(results_root, "models")
output_dir <- file.path(results_root, "phenotypes/algorithm_robustness")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

modalities <- c("nmr", "protein")
algorithms <- c("ngboost", "lightgbmlss", "xgboostlss")
algorithm_labels <- c(
  ngboost = "NGBoost",
  lightgbmlss = "LightGBMLSS",
  xgboostlss = "XGBoostLSS"
)
expected_n <- c(nmr = 272672L, protein = 43867L)

load_prediction <- function(modality, algorithm) {
  path <- file.path(
    model_root, modality, algorithm, "all_participant_predictions.csv"
  )
  x <- fread(
    path,
    select = c(
      "eid_int", "Age", "Sex_male", "AssessmentCenter", "OuterFold",
      "BAG", "BAU"
    )
  )
  setnames(x, "eid_int", "eid")
  x[, eid := as.integer(eid)]
  setorder(x, eid)
  stopifnot(
    nrow(x) == expected_n[[modality]],
    !anyDuplicated(x$eid),
    all(is.finite(x$BAG)), all(is.finite(x$BAU)),
    identical(sort(unique(x$OuterFold)), 1:5)
  )
  x
}

predictions <- list()
for (modality in modalities) {
  for (algorithm in algorithms) {
    predictions[[paste(modality, algorithm, sep = "/")]] <-
      load_prediction(modality, algorithm)
  }
}

for (modality in modalities) {
  reference <- predictions[[paste(modality, "ngboost", sep = "/")]]
  for (algorithm in algorithms[-1L]) {
    candidate <- predictions[[paste(modality, algorithm, sep = "/")]]
    stopifnot(
      identical(reference$eid, candidate$eid),
      isTRUE(all.equal(reference$Age, candidate$Age, tolerance = 0)),
      identical(reference$Sex_male, candidate$Sex_male),
      identical(as.character(reference$AssessmentCenter), as.character(candidate$AssessmentCenter)),
      identical(reference$OuterFold, candidate$OuterFold)
    )
  }
}

rank_agreement <- rbindlist(lapply(modalities, function(modality) {
  rbindlist(lapply(algorithms, function(algorithm_1) {
    x <- predictions[[paste(modality, algorithm_1, sep = "/")]]
    rbindlist(lapply(algorithms, function(algorithm_2) {
      y <- predictions[[paste(modality, algorithm_2, sep = "/")]]
      data.table(
        modality = modality,
        algorithm_1 = algorithm_1,
        algorithm_1_label = algorithm_labels[[algorithm_1]],
        algorithm_2 = algorithm_2,
        algorithm_2_label = algorithm_labels[[algorithm_2]],
        n = nrow(x),
        spearman_r = suppressWarnings(cor(x$BAU, y$BAU, method = "spearman"))
      )
    }))
  }))
}))

setorder(rank_agreement, modality, algorithm_1, algorithm_2)
fwrite(rank_agreement, file.path(output_dir, "bau_rank_agreement.tsv"), sep = "\t")
