suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
})

read_gene_sets <- function(reference_dir) {
  names <- c("hallmark", "reactome", "go_bp", "go_cc", "go_mf")
  reference <- setNames(lapply(names, function(name) {
    fread(file.path(reference_dir, paste0("MSigDB_2026.1_Hs_", name, ".tsv.gz")))
  }), names)
  versions <- unique(unlist(lapply(reference, function(x) x$db_version)))
  if (!identical(versions, "2026.1.Hs")) stop("Gene-set reference version differs")
  reference$hallmark[, collection_label := "Hallmark"]
  reference$reactome[, collection_label := "Reactome"]
  reference
}

run_gsea_table <- function(annotation, ranks, seed, permutations, sample_size, retry = FALSE) {
  annotation <- unique(annotation[gene_symbol %in% names(ranks)])
  membership <- unique(annotation[, .(gene_set, gene_symbol)])
  pathways <- split(membership$gene_symbol, membership$gene_set)
  ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)
  set.seed(seed)
  result <- as.data.table(fgseaMultilevel(
    pathways = pathways,
    stats = ranks,
    sampleSize = sample_size,
    minSize = 10,
    maxSize = 500,
    eps = 1e-20,
    nPermSimple = permutations,
    scoreType = "std",
    nproc = 1
  ))
  unresolved <- result[!is.finite(pval), pathway]
  if (retry && length(unresolved)) {
    set.seed(seed + 500000L)
    repeated <- as.data.table(fgseaMultilevel(
      pathways = pathways[unresolved],
      stats = ranks,
      sampleSize = sample_size,
      minSize = 10,
      maxSize = 500,
      eps = 1e-20,
      nPermSimple = 1500000L,
      scoreType = "std",
      nproc = 1
    ))
    result <- rbindlist(list(result[!pathway %in% unresolved], repeated))
  }
  setnames(result, "pathway", "gene_set")
  result[, leading_edge := vapply(leadingEdge, function(x) paste(x, collapse = ";"), character(1))]
  result[, leadingEdge := NULL]
  result <- merge(
    result, unique(annotation[, .(gene_set, collection_label, description)]),
    by = "gene_set", all.x = TRUE
  )
  result[, padj := p.adjust(pval, "BH"), by = collection_label]
  result[, seed := seed]
  result
}

gsea_replication <- function(results, stability = NULL) {
  paired <- merge(
    results[stage == "discovery"], results[stage == "validation"],
    by = c("estimand", "collection_label", "gene_set", "description"),
    suffixes = c("_discovery", "_validation")
  )
  paired[, numerically_stable := TRUE]
  if (!is.null(stability)) {
    paired <- merge(paired, stability, by = c("collection_label", "gene_set"), all.x = TRUE)
    paired[estimand == "BAG", numerically_stable :=
      !is.na(stable_discovery_significant) & stable_discovery_significant]
  }
  paired[, discovery_selected :=
    is.finite(padj_discovery) & padj_discovery < 0.05 & numerically_stable]
  paired[, validation_directional_p := fifelse(
    is.finite(NES_discovery) & is.finite(NES_validation) &
      sign(NES_discovery) == sign(NES_validation),
    pval_validation, 1
  )]
  paired[, validation_directional_q := NA_real_]
  paired[discovery_selected == TRUE,
    validation_directional_q := p.adjust(validation_directional_p, "BH"),
    by = .(estimand, collection_label)
  ]
  paired[, replicated := discovery_selected &
    sign(NES_discovery) == sign(NES_validation) &
    is.finite(validation_directional_q) & validation_directional_q < 0.05]
  paired[is.na(replicated), replicated := FALSE]
  paired
}

go_overrepresentation <- function(annotation, measured, selected, direction, ontology) {
  universe <- intersect(measured, unique(annotation$gene_symbol))
  foreground <- intersect(selected, universe)
  annotation <- unique(annotation[
    gene_symbol %in% universe, .(gene_set, description, gene_symbol)
  ])
  result <- annotation[, .(
    term_size = .N,
    overlap = sum(gene_symbol %in% foreground),
    overlapping_proteins = paste(sort(gene_symbol[gene_symbol %in% foreground]), collapse = ";")
  ), by = .(gene_set, description)][term_size >= 10L & term_size <= 500L]
  result[, p_value := phyper(
    overlap - 1L, term_size, length(universe) - term_size,
    length(foreground),
    lower.tail = FALSE
  )]
  result[, `:=`(
    fdr = p.adjust(p_value, "BH"),
    fold_enrichment = (overlap / length(foreground)) / (term_size / length(universe)),
    foreground_n = length(foreground),
    background_n = length(universe),
    direction = direction,
    ontology = ontology
  )]
  result[, log2_enrichment := log2(fold_enrichment)]
  result
}
