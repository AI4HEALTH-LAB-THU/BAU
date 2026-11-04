#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
population_arg <- grep("^--population=", args, value = TRUE)
population <- if (length(population_arg)) sub("^--population=", "", population_arg[[1]]) else "primary"
if (!population %in% c("primary", "healthy")) stop("Population must be primary or healthy")
args_all <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
default_root <- normalizePath(file.path(dirname(script_file), "../.."), mustWork = TRUE)
code_root <- Sys.getenv("BAU_CODE_ROOT", default_root)
data_root <- Sys.getenv("BAU_DATA_ROOT", "/home/dnanexus/bau-work/data")
results_root <- Sys.getenv("BAU_RESULTS_ROOT", "/home/dnanexus/bau-work/results")
reference_dir <- Sys.getenv("GENE_SET_ROOT", file.path(data_root, "references/msigdb"))
source(file.path(code_root, "R/enrichment.R"))
reference <- read_gene_sets(reference_dir)
display_config <- fromJSON(file.path(
  code_root, "config/displayed_pathways.json"
))
result_dir <- file.path(results_root, "molecular/enrichment", population)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
estimands <- c("BAG", "BAU", "difference")
stages <- c("discovery", "validation")
collections <- c("Hallmark", "Reactome")
gsea_reference <- rbind(reference$hallmark, reference$reactome, fill = TRUE)

for (modality_index in 1:2) {
  modality <- c("nmr", "protein")[[modality_index]]
  analysis_id <- if (modality == "nmr") "nmr_to_protein" else "protein_to_noninput_protein"
  association_dir <- if (population == "primary") "cross_omic" else "healthy"
  association <- fread(file.path(
    results_root, "molecular", association_dir,
    paste0(analysis_id, "_associations.tsv")
  ))
  rank_table <- rbindlist(lapply(estimands, function(estimand) {
    data.table(
      feature = association$feature,
      estimand = estimand,
      discovery_statistic = association[[paste0("beta_", estimand, "_discovery")]] /
        association[[paste0("se_", estimand, "_discovery")]],
      validation_statistic = association[[paste0("beta_", estimand, "_validation")]] /
        association[[paste0("se_", estimand, "_validation")]]
    )
  }))

  run_one <- function(estimand_index, stage_index, offset = 0L) {
    term <- estimands[[estimand_index]]
    stage_name <- stages[[stage_index]]
    ranks <- rank_table[
      estimand == term,
      setNames(get(paste0(stage_name, "_statistic")), feature)
    ]
    if (population == "primary" && modality == "nmr") {
      seed <- c(20261000L, 20262000L)[[stage_index]] + estimand_index
      result <- run_gsea_table(gsea_reference, ranks, seed, 10000L, 101L)
    } else {
      result <- rbindlist(lapply(1:2, function(collection_index) {
        collection_name <- collections[[collection_index]]
        seed <- if (population == "primary") {
          c(20262500L, 20263500L)[[stage_index]] + estimand_index +
            collection_index * 100000L + offset
        } else {
          20268000L + modality_index * 1000L + estimand_index * 100L +
            stage_index * 10L + collection_index + offset
        }
        run_gsea_table(
          gsea_reference[collection_label == collection_name], ranks, seed,
          if (term == "BAG") 150000L else 50000L,
          501L,
          retry = population == "healthy" && term == "BAG"
        )
      }))
    }
    result[, `:=`(estimand = term, stage = stage_name)]
    result
  }

  message("Running ", population, " enrichment for ", modality)
  results <- rbindlist(lapply(seq_along(estimands), function(i) {
    rbind(run_one(i, 1L), run_one(i, 2L))
  }), fill = TRUE)
  stability <- NULL
  if (population == "healthy" || modality == "protein") {
    production <- copy(results[estimand == "BAG" & stage == "discovery"])
    production[, seed_offset := 0L]
    seed_results <- rbindlist(c(list(production), lapply(c(1009L, 2027L), function(offset) {
      run_one(1L, 1L, offset)[, seed_offset := offset]
    })), fill = TRUE)
    stability <- seed_results[, .(
      stable_discovery_significant = .N == 3L &&
        all(is.finite(padj) & padj < 0.05) &&
        uniqueN(sign(NES[is.finite(NES)])) == 1L
    ), by = .(collection_label, gene_set)]
    fwrite(seed_results, file.path(result_dir, paste0(modality, "_bag_seed_results.tsv")), sep = "\t")
  }
  paired <- gsea_replication(results, stability)
  selected_terms <- as.data.table(display_config[[paste0(modality, "_gsea")]])
  selected_terms[, display_order := seq_len(.N)]
  displayed <- merge(
    paired, selected_terms,
    by = "gene_set", all = FALSE, sort = FALSE
  )
  fwrite(rank_table, file.path(result_dir, paste0(modality, "_rank_statistics.tsv")), sep = "\t")
  fwrite(results, file.path(result_dir, paste0(modality, "_gsea_results.tsv")), sep = "\t")
  fwrite(paired, file.path(result_dir, paste0(modality, "_gsea_replication.tsv")), sep = "\t")
  fwrite(displayed, file.path(result_dir, paste0(modality, "_gsea_displayed.tsv")), sep = "\t")

  if (population == "primary") {
    direction_sets <- list(
      "BAG positive" = association[replicated_BAG & beta_BAG_discovery > 0, feature],
      "BAG negative" = association[replicated_BAG & beta_BAG_discovery < 0, feature],
      "BAU positive" = association[replicated_BAU & beta_BAU_discovery > 0, feature],
      "BAU negative" = association[replicated_BAU & beta_BAU_discovery < 0, feature],
      "BAU > BAG" = association[replicated_difference & beta_difference_discovery > 0, feature],
      "BAU < BAG" = association[replicated_difference & beta_difference_discovery < 0, feature]
    )
    go_tables <- list(BP = reference$go_bp, CC = reference$go_cc, MF = reference$go_mf)
    go_results <- rbindlist(lapply(names(direction_sets), function(direction) {
      rbindlist(lapply(names(go_tables), function(ontology) {
        go_overrepresentation(
          go_tables[[ontology]], association$feature,
          direction_sets[[direction]], direction, ontology
        )
      }))
    }))
    go_terms <- as.data.table(display_config[[paste0(modality, "_go")]])
    go_terms[, display_order := seq_len(.N)]
    go_displayed <- merge(
      go_results[fdr < 0.05 & overlap >= 3L & fold_enrichment > 1],
      go_terms,
      by = "gene_set", all = FALSE, sort = FALSE
    )
    fwrite(go_results, file.path(result_dir, paste0(modality, "_go_results.tsv")), sep = "\t")
    fwrite(go_displayed, file.path(result_dir, paste0(modality, "_go_displayed.tsv")), sep = "\t")
  } else {
    primary <- fread(file.path(
      results_root, "molecular/enrichment/primary",
      paste0(modality, "_gsea_replication.tsv")
    ))
    keep <- c("estimand", "collection_label", "gene_set", "NES_discovery", "replicated")
    comparison <- merge(
      primary[, ..keep], paired[, ..keep],
      by = c("estimand", "collection_label", "gene_set"),
      suffixes = c("_primary", "_healthy")
    )
    fwrite(comparison, file.path(result_dir, paste0(modality, "_primary_comparison.tsv")), sep = "\t")
  }
}
