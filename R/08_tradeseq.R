source("R/00_utils.R")

# ── GAM fitting ────────────────────────────────────────────────────────────────

prepare_tradeseq_counts <- function(obj, sds_obj, assay = "RNA") {
  variable_genes <- Seurat::VariableFeatures(obj)
  if (length(variable_genes) == 0) {
    warning("No variable features found — using all genes for tradeSeq")
    variable_genes <- rownames(obj)
  }
  counts <- Seurat::GetAssayData(obj, assay = assay, layer = "counts")
  counts[variable_genes, , drop = FALSE]
}

run_fitgam <- function(filt_counts, sds_obj, params, audit_dir = NULL, label = "") {
  # sds_obj must be the SlingshotDataSet from which getCurves was already run.
  # Original code referenced 'curves' (undefined) and 'sce2' (undefined).
  # Fixed: sds_obj is the explicit argument; sce_post is the return value.
  bp <- BiocParallel::MulticoreParam(workers = params$tradeseq$cores)
  BiocParallel::register(bp)

  message(label, ": Fitting GAMs with tradeSeq (nknots=",
          params$tradeseq$nknots, ", cores=", params$tradeseq$cores, ") ...")
  sce <- tradeSeq::fitGAM(
    counts   = as.matrix(filt_counts),
    sds      = sds_obj,
    nknots   = params$tradeseq$nknots,
    BPPARAM  = bp
  )

  # Store pseudotime per lineage in colData
  pt_gam <- slingshot::slingPseudotime(sds_obj)
  common  <- intersect(colnames(sce), rownames(pt_gam))
  pt_sub  <- pt_gam[common, , drop = FALSE]
  for (j in seq_len(ncol(pt_sub))) {
    SummarizedExperiment::colData(sce)[common, paste0("slingshot_pseudotime_", j)] <-
      pt_sub[, j]
  }

  if (!is.null(audit_dir)) {
    audit <- list(
      module    = paste0("tradeseq_fitgam_", label),
      timestamp = as.character(Sys.time()),
      n_genes   = nrow(sce),
      n_cells   = ncol(sce),
      nknots    = params$tradeseq$nknots,
      cores     = params$tradeseq$cores,
      n_lineages = ncol(pt_gam),
      warnings  = "sds_obj must come from getCurves on the correct object (see fix log)"
    )
    write_audit(audit, name = paste0(label, "_tradeseq_fitgam"), out_dir = audit_dir)
  }

  message(label, ": GAM fitting complete")
  sce
}

# ── Statistical tests ──────────────────────────────────────────────────────────

run_association_test <- function(sce) {
  message("Running associationTest ...")
  tradeSeq::associationTest(sce)
}

run_startvsend_test <- function(sce, lineages = TRUE) {
  # Original code referenced 'sce2' — corrected to use the passed sce object
  message("Running startVsEndTest ...")
  res <- tradeSeq::startVsEndTest(sce, lineages = lineages)
  res %>%
    tibble::rownames_to_column(var = "gene") %>%
    dplyr::filter(!grepl("^MT-", gene))
}

run_pattern_test <- function(sce) {
  message("Running patternTest ...")
  res <- tradeSeq::patternTest(sce, pairwise = TRUE)
  res %>%
    tibble::rownames_to_column(var = "gene") %>%
    dplyr::filter(!grepl("^MT-", gene))
}

run_earlyDE_test <- function(sce) {
  message("Running earlyDETest ...")
  tradeSeq::earlyDETest(sce, knots = c(1, 2))
}

run_diffend_test <- function(sce) {
  message("Running diffEndTest ...")
  tradeSeq::diffEndTest(sce)
}

# ── Results summaries ──────────────────────────────────────────────────────────

top_genes_per_lineage <- function(startRes, n_lineages = 4, top_n = 100) {
  lapply(seq_len(n_lineages), function(k) {
    colk <- paste0("waldStat_lineage", k)
    if (!colk %in% colnames(startRes)) return(NULL)
    startRes %>%
      dplyr::arrange(dplyr::desc(.data[[colk]])) %>%
      dplyr::slice_head(n = top_n) %>%
      dplyr::select(gene = gene, waldStat = .data[[colk]])
  }) |> setNames(paste0("Lineage", seq_len(n_lineages)))
}

build_transient_score <- function(patternRes, endRes) {
  compare <- merge(
    patternRes %>% dplyr::select(Gene = gene, pattern = waldStat),
    data.frame(Gene = rownames(endRes), end = endRes$waldStat),
    by = "Gene", all = FALSE
  )
  compare$transientScore <-
    rank(-compare$end, ties.method = "min")^2 +
    rank(compare$pattern, ties.method = "random")^2
  compare
}

# ── Heatmaps ───────────────────────────────────────────────────────────────────

scale_rows <- function(x) {
  m <- rowMeans(x, na.rm = TRUE)
  s <- apply(x, 1, stats::sd, na.rm = TRUE)
  s[s == 0 | is.na(s)] <- 1
  sweep(sweep(x, 1, m, "-"), 1, s, "/")
}

plot_tradeseq_heatmaps <- function(sce, avg_exp, top_by_lin, annotations,
                                    out_dir, prefix) {
  col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("darkgreen", "white", "red"))
  avg_norm <- setNames(colnames(avg_exp), normalize_labels(colnames(avg_exp)))

  lineage_names <- names(top_by_lin)
  ht_list <- lapply(seq_along(lineage_names), function(k) {
    lin_name    <- lineage_names[k]
    genes_k     <- intersect(top_by_lin[[lin_name]], rownames(avg_exp))
    wanted_disp <- annotations$tradeseq_lineage_celltypes[[lin_name]]
    if (is.null(wanted_disp) || length(genes_k) == 0) return(NULL)

    wanted_norm  <- normalize_labels(wanted_disp)
    mapped       <- unname(avg_norm[wanted_norm])
    missing_mask <- is.na(mapped)

    avg_ext <- avg_exp
    if (any(missing_mask)) {
      miss_names <- wanted_disp[missing_mask]
      filler     <- matrix(
        NA_real_, nrow = nrow(avg_exp), ncol = sum(missing_mask),
        dimnames = list(rownames(avg_exp), miss_names)
      )
      avg_ext    <- cbind(avg_exp, filler)
      mapped[missing_mask] <- miss_names
    }

    mat_k      <- as.matrix(avg_ext[genes_k, mapped, drop = FALSE])
    colnames(mat_k) <- wanted_disp
    mat_kz     <- scale_rows(mat_k)

    ComplexHeatmap::Heatmap(
      mat_kz, name = paste0("Z (", lin_name, ")"), col = col_fun,
      na_col = "grey95", cluster_rows = FALSE, cluster_columns = FALSE,
      show_row_names = TRUE, show_column_names = TRUE,
      column_title = paste("Cell types:", lin_name),
      row_title    = lin_name,
      row_title_gp = grid::gpar(fontface = "bold"),
      heatmap_legend_param = list(at = c(-2, 0, 2), labels = c("Low", "Mid", "High"))
    )
  })
  ht_list <- Filter(Negate(is.null), ht_list)
  ht_list
}

save_tradeseq_results <- function(results_list, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  for (nm in names(results_list)) {
    path <- file.path(out_dir, paste0(nm, ".csv"))
    utils::write.csv(results_list[[nm]], path, row.names = TRUE)
    message("Saved: ", path)
  }
}
