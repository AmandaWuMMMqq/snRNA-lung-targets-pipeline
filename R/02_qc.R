source("R/00_utils.R")

run_qc_check <- function(obj, params, label = "", audit_dir = NULL) {
  qc <- params$qc
  n_before <- ncol(obj)
  cells_keep <- colnames(obj)

  if (!is.null(qc$min_features)) {
    nf <- Matrix::colSums(obj[["RNA"]]$counts > 0)
    cells_keep <- intersect(cells_keep, names(nf)[nf >= qc$min_features])
  }
  if (!is.null(qc$max_percent_mt) && "percent.mt" %in% colnames(obj@meta.data)) {
    cells_keep <- intersect(
      cells_keep,
      rownames(obj@meta.data)[obj@meta.data$percent.mt <= qc$max_percent_mt]
    )
  }

  n_after <- length(cells_keep)
  n_removed <- n_before - n_after

  if (n_removed > 0) {
    message(label, ": QC removed ", n_removed, " cells (",
            round(100 * n_removed / n_before, 1), "%)")
    obj <- subset(obj, cells = cells_keep)
  } else {
    message(label, ": QC — no cells removed (thresholds not binding)")
  }

  if (!is.null(audit_dir)) {
    audit <- list(
      module     = "qc",
      label      = label,
      timestamp  = as.character(Sys.time()),
      n_before   = n_before,
      n_after    = n_after,
      n_removed  = n_removed,
      qc_params  = qc,
      warnings   = if (n_removed == 0) "No cells removed — check QC thresholds" else NULL
    )
    write_audit(audit, name = paste0(label, "_qc"), out_dir = audit_dir)
  }

  obj
}

compute_qc_metrics_table <- function(obj) {
  meta <- obj@meta.data
  data.frame(
    n_cells      = ncol(obj),
    n_genes      = nrow(obj),
    median_nFeature = if ("nFeature_RNA" %in% colnames(meta)) median(meta$nFeature_RNA, na.rm = TRUE) else NA,
    median_pct_mt   = if ("percent.mt"   %in% colnames(meta)) median(meta$percent.mt,   na.rm = TRUE) else NA,
    n_batches    = if ("batch" %in% colnames(meta)) length(unique(meta$batch)) else NA
  )
}
