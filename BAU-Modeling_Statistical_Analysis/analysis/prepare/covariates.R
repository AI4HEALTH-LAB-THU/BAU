#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(arrow))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: export_covariates.R participant_index source_rds output_feather")
}

index <- as.data.frame(read_feather(args[[1L]], col_select = "eid"))
source <- as.data.frame(readRDS(args[[2L]]))
to_eid <- function(x) as.integer(round(as.numeric(as.character(x))))

index$eid <- to_eid(index$eid)
source$eid <- to_eid(source$eid)
if (anyNA(index$eid) || anyDuplicated(index$eid)) stop("Invalid participant index")
if (anyNA(source$eid) || anyDuplicated(source$eid)) stop("Invalid covariate identifiers")

row_index <- match(index$eid, source$eid)
if (anyNA(row_index)) stop("Covariates do not cover the participant index")
out <- source[row_index, , drop = FALSE]
out$eid <- index$eid
write_feather(out, args[[3L]], compression = "zstd")
