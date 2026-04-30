source("R/00_utils.R")

# ── Contaminant filtering ──────────────────────────────────────────────────────

build_contaminant_filter_expr <- function(contaminant_genes) {
  parts <- vapply(contaminant_genes, function(g) {
    sprintf("%s < %s", g$gene, g$threshold)
  }, character(1))
  paste(parts, collapse = " & ")
}

filter_contaminants <- function(obj, contaminant_genes, lineage_col, lineage_label) {
  # Build explicit subset call using Seurat::subset with a filter expression.
  # The original code had bare expressions without assignment — this version
  # explicitly assigns the result. (Fix for flagged issue #1.)
  n_before <- ncol(obj)

  filter_str <- sprintf(
    "%s == '%s' & %s",
    lineage_col,
    lineage_label,
    build_contaminant_filter_expr(contaminant_genes)
  )

  obj_filtered <- Seurat::subset(
    obj,
    subset = eval(str2expression(filter_str))
  )

  n_after <- ncol(obj_filtered)
  message("Contaminant filter: ", n_before, " → ", n_after, " cells removed: ",
          n_before - n_after)
  obj_filtered
}

# ── Pre-merge epithelial subsetting ───────────────────────────────────────────

subset_epithelial_pre <- function(fetal_obj, pediatric_obj, params, audit_dir = NULL) {
  contaminants <- params$subsetting$contaminant_genes
  lineage_col  <- "lineage_lvl1"
  lineage_label <- "Epithelium"

  check_required_metadata(fetal_obj,     lineage_col)
  check_required_metadata(pediatric_obj, lineage_col)

  message("Pre-merge: subsetting epithelium from fetal object ...")
  fetal_epi <- filter_contaminants(
    fetal_obj, contaminants, lineage_col, lineage_label
  )

  message("Pre-merge: subsetting epithelium from pediatric object ...")
  ped_epi <- filter_contaminants(
    pediatric_obj, contaminants, lineage_col, lineage_label
  )

  message("Fetal epithelial cells: ",    ncol(fetal_epi))
  message("Pediatric epithelial cells: ", ncol(ped_epi))

  merged <- merge(fetal_epi, ped_epi)
  Seurat::DefaultAssay(merged) <- "RNA"

  if (!is.null(audit_dir)) {
    audit <- list(
      module        = "subsetting_pre",
      timestamp     = as.character(Sys.time()),
      fetal_n_in    = ncol(fetal_obj),
      fetal_n_out   = ncol(fetal_epi),
      ped_n_in      = ncol(pediatric_obj),
      ped_n_out     = ncol(ped_epi),
      merged_n      = ncol(merged),
      contaminants  = contaminants,
      lineage_col   = lineage_col,
      lineage_label = lineage_label
    )
    write_audit(audit, name = "pre_epi_subsetting", out_dir = audit_dir)
  }

  merged
}

# ── Post-merge epithelial subsetting ──────────────────────────────────────────

subset_epithelial_post <- function(joint_obj, params, lineage_col = "lineage_lvl2",
                                    audit_dir = NULL) {
  contaminants  <- params$subsetting$contaminant_genes
  lineage_label <- "Epithelium"

  check_required_metadata(joint_obj, lineage_col)

  n_before <- ncol(joint_obj)

  message("Post-merge: subsetting epithelium from joint object ...")
  # Explicit assignment — fixes flagged issue #1 (original had bare expressions)
  obj_filtered <- filter_contaminants(
    joint_obj, contaminants, lineage_col, lineage_label
  )

  message("Post-merge epithelial cells: ", ncol(obj_filtered))

  if (!is.null(audit_dir)) {
    audit <- list(
      module        = "subsetting_post",
      timestamp     = as.character(Sys.time()),
      n_before      = n_before,
      n_after       = ncol(obj_filtered),
      lineage_col   = lineage_col,
      lineage_label = lineage_label,
      contaminants  = contaminants
    )
    write_audit(audit, name = "post_epi_subsetting", out_dir = audit_dir)
  }

  obj_filtered
}

# ── Alveolar subsetting ────────────────────────────────────────────────────────

subset_alveolar <- function(epi_obj, cell_type_col, alveolar_labels,
                             audit_dir = NULL, label = "") {
  check_required_metadata(epi_obj, cell_type_col)

  n_before <- ncol(epi_obj)
  alve_obj <- subset(epi_obj, subset = .data[[cell_type_col]] %in% alveolar_labels)
  n_after  <- ncol(alve_obj)

  message(label, " alveolar subset: ", n_before, " → ", n_after, " cells")

  if (!is.null(audit_dir)) {
    audit <- list(
      module         = paste0("alveolar_subsetting_", label),
      timestamp      = as.character(Sys.time()),
      n_before       = n_before,
      n_after        = n_after,
      cell_type_col  = cell_type_col,
      alveolar_labels = alveolar_labels
    )
    write_audit(audit, name = paste0(label, "_alveolar_subsetting"), out_dir = audit_dir)
  }

  alve_obj
}
