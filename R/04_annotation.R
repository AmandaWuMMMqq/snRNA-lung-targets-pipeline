source("R/00_utils.R")

apply_cluster_annotation <- function(obj, mapping, new_col, cluster_col = "seurat_clusters") {
  if (!cluster_col %in% colnames(obj@meta.data)) {
    if (cluster_col == "seurat_clusters") {
      cluster_vals <- as.character(Seurat::Idents(obj))
    } else {
      stop("Column not found: ", cluster_col)
    }
  } else {
    cluster_vals <- as.character(obj@meta.data[[cluster_col]])
  }

  annotated <- mapping[cluster_vals]
  n_unmapped <- sum(is.na(annotated))
  if (n_unmapped > 0) {
    warning(n_unmapped, " cells could not be mapped in '", new_col,
            "' — check annotation map for missing cluster IDs")
  }

  obj@meta.data[[new_col]] <- annotated
  message("Applied annotation '", new_col, "': ",
          length(unique(na.omit(annotated))), " unique labels")
  obj
}

compute_annotation_crosswalk <- function(obj, cluster_col, ref_col) {
  if (!all(c(cluster_col, ref_col) %in% colnames(obj@meta.data))) {
    stop("Required columns missing: ", paste(c(cluster_col, ref_col), collapse = ", "))
  }
  obj@meta.data %>%
    dplyr::filter(!is.na(.data[[ref_col]])) %>%
    dplyr::group_by(.data[[cluster_col]], .data[[ref_col]]) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(.data[[cluster_col]]) %>%
    dplyr::mutate(percentage = count / sum(count) * 100) %>%
    dplyr::slice_max(order_by = percentage, n = 1, with_ties = FALSE) %>%
    dplyr::arrange(as.numeric(as.character(.data[[cluster_col]])))
}

find_top_markers <- function(obj, assay = "RNA", n = 10) {
  if (!requireNamespace("presto", quietly = TRUE)) stop("Package 'presto' required")
  res <- presto::wilcoxauc(obj, seurat_assay = assay)
  presto::top_markers(res, n = n)
}

store_cluster_ids <- function(obj) {
  obj@meta.data$cluster <- Seurat::Idents(obj)
  obj
}
