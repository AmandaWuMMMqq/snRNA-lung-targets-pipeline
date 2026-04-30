source("R/00_utils.R")

integrate_seurat <- function(obj, params, split_by = "batch", audit_dir = NULL, label = "") {
  dims       <- seq_len(params$seurat$dims)
  resolution <- params$seurat$resolution
  n_before   <- ncol(obj)

  message(label, ": Running SCTransform ...")
  Seurat::DefaultAssay(obj) <- "RNA"
  obj <- Seurat::JoinLayers(obj)
  obj[["RNA"]] <- Seurat::split(obj[["RNA"]], f = obj@meta.data[[split_by]])

  obj <- Seurat::SCTransform(
    obj,
    vars.to.regress = params$seurat$vars_to_regress,
    verbose = FALSE
  )

  message(label, ": Running PCA ...")
  obj <- Seurat::RunPCA(obj, dims = dims, verbose = FALSE)

  message(label, ": Integrating layers (RPCA) ...")
  obj <- Seurat::IntegrateLayers(
    object               = obj,
    method               = Seurat::RPCAIntegration,
    normalization.method = params$seurat$normalization_method,
    verbose              = FALSE
  )

  obj[["RNA"]] <- Seurat::JoinLayers(obj[["RNA"]])

  message(label, ": FindNeighbors + FindClusters + RunUMAP ...")
  has_integrated_dr <- "integrated.dr" %in% names(obj@reductions)
  reduction_use <- if (has_integrated_dr) "integrated.dr" else "pca"

  obj <- Seurat::FindNeighbors(obj, dims = dims, reduction = reduction_use, verbose = FALSE)
  obj <- Seurat::FindClusters(obj, resolution = resolution)
  obj <- Seurat::RunUMAP(obj, dims = dims, reduction = reduction_use)

  if (!is.null(audit_dir)) {
    audit <- list(
      module        = "integration",
      label         = label,
      timestamp     = as.character(Sys.time()),
      n_cells       = ncol(obj),
      n_genes       = nrow(obj),
      assays        = names(obj@assays),
      reductions    = names(obj@reductions),
      meta_cols     = colnames(obj@meta.data),
      dims_used     = params$seurat$dims,
      resolution    = resolution,
      split_by      = split_by,
      reduction_used = reduction_use
    )
    write_audit(audit, name = paste0(label, "_integration"), out_dir = audit_dir)
  }

  message(label, ": Integration complete — ", ncol(obj), " cells")
  obj
}

downsample_seurat <- function(obj, n, seed = 123) {
  if (ncol(obj) <= n) {
    message("Object has ", ncol(obj), " cells — no downsampling needed (target: ", n, ")")
    return(obj)
  }
  set.seed(seed)
  cells <- sample(colnames(obj), n)
  message("Downsampled from ", ncol(obj), " to ", n, " cells")
  subset(obj, cells = cells)
}

merge_objects <- function(obj1, obj2) {
  merged <- merge(obj1, obj2)
  Seurat::DefaultAssay(merged) <- "RNA"
  merged
}
