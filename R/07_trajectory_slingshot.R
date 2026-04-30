source("R/00_utils.R")

prepare_slingshot_inputs <- function(obj, reduction = "umap", assay = "RNA",
                                      cluster_col = "cell_type2") {
  check_required_metadata(obj, cluster_col)
  if (!reduction %in% names(obj@reductions)) {
    stop("Reduction '", reduction, "' not found in object. Available: ",
         paste(names(obj@reductions), collapse = ", "))
  }
  list(
    dimred     = obj@reductions[[reduction]]@cell.embeddings,
    clustering = obj@meta.data[[cluster_col]],
    counts     = Seurat::GetAssayData(obj, assay = assay, layer = "counts")
  )
}

extract_seurat_colors <- function(obj, group_col) {
  Seurat::Idents(obj) <- group_col
  p          <- Seurat::DimPlot(obj, group.by = group_col, label = TRUE)
  pb         <- ggplot2::ggplot_build(p)
  umap_cols  <- unique(pb$data[[1]]$colour)
  col_levels <- levels(Seurat::Idents(obj))
  setNames(umap_cols, col_levels)
}

run_slingshot_lineages <- function(dimred, clustering, start_clus, dist_method,
                                    params) {
  tp <- params$trajectory
  message("Running getLineages: start='", start_clus, "', dist='", dist_method, "' ...")
  lin <- slingshot::getLineages(
    data          = dimred,
    clusterLabels = clustering,
    start.clus    = start_clus,
    dist.method   = dist_method
  )
  message("Running getCurves ...")
  sds <- slingshot::getCurves(
    lin,
    approx_points = tp$approx_points,
    thresh        = tp$thresh,
    stretch       = tp$stretch,
    allow.breaks  = tp$allow_breaks,
    shrink        = tp$shrink
  )
  slingshot::SlingshotDataSet(sds)
}

subset_slingshot <- function(sds, lineage_indices) {
  sub         <- sds
  sub@lineages <- sds@lineages[lineage_indices]
  sub@curves   <- sds@curves[lineage_indices]
  names(sub@lineages) <- names(sub@curves)
  sub
}

compute_pseudotime <- function(sds, cells) {
  combined_pt_from_sds(sds, cells)
}

plot_slingshot_lineages <- function(dimred, clustering, sds_list, colors_map,
                                     out_path, line_colors = NULL, line_lty = NULL,
                                     legend_labels = NULL) {
  grDevices::pdf(out_path, width = 12, height = 10, useDingbats = FALSE)
  graphics::plot(
    dimred[, 1:2],
    col  = colors_map[as.character(clustering)],
    cex  = 0.5, pch = 16,
    xlab = "", ylab = "", xaxt = "n", yaxt = "n", axes = FALSE
  )
  if (is.null(line_colors)) line_colors <- rep("steelblue", length(sds_list))
  if (is.null(line_lty))    line_lty    <- seq_along(sds_list)

  for (i in seq_along(sds_list)) {
    slingshot::lines(
      sds_list[[i]],
      lwd = 5,
      col = grDevices::adjustcolor(line_colors[i], alpha.f = 0.8),
      lty = line_lty[i]
    )
  }

  unique_clusters <- sort(unique(as.character(clustering)))
  for (cl in unique_clusters) {
    idx <- which(as.character(clustering) == cl)
    if (length(idx) > 0)
      graphics::text(mean(dimred[idx, 1]), mean(dimred[idx, 2]),
                     labels = cl, font = 2, cex = 1)
  }

  if (!is.null(legend_labels)) {
    graphics::legend("topleft", legend = legend_labels,
                     col = line_colors, lwd = 5, lty = line_lty, bty = "n")
  }
  grDevices::dev.off()
  message("Saved slingshot plot: ", out_path)
}

plot_pseudotime_umap <- function(dimred, pt_combined, clustering, sds_list,
                                  palette_fn, out_path) {
  n_colors <- 200
  pal <- palette_fn(n_colors)
  idx <- 1 + floor(pt_combined * (n_colors - 1))
  idx[idx < 1 | idx > n_colors | is.na(idx)] <- 1
  pt_cols <- pal[idx]

  grDevices::pdf(out_path, width = 12, height = 9, useDingbats = FALSE)
  graphics::plot(
    dimred[, 1:2], col = pt_cols, pch = 16, cex = 0.5,
    xlab = "", ylab = "", xaxt = "n", yaxt = "n", axes = FALSE
  )
  unique_clusters <- sort(unique(as.character(clustering)))
  for (cl in unique_clusters) {
    idx2 <- which(clustering == cl)
    if (length(idx2) > 0)
      graphics::text(mean(dimred[idx2, 1]), mean(dimred[idx2, 2]),
                     labels = cl, font = 2, cex = 1.2, col = "darkred")
  }
  for (s in sds_list) {
    slingshot::lines(s, lwd = 2, col = grDevices::adjustcolor("mediumblue", 0.5))
  }
  grDevices::dev.off()
  message("Saved pseudotime plot: ", out_path)
}

plot_pseudotime_ridges <- function(pt_all, clustering, annotations, palette_fn,
                                    out_dir, prefix, plotting_cfg) {
  pt_long <- as.data.frame(pt_all) %>%
    dplyr::mutate(cell_type = clustering) %>%
    tidyr::pivot_longer(
      cols      = colnames(pt_all),
      names_to  = "lineage",
      values_to = "pt"
    ) %>%
    dplyr::filter(!is.na(pt))

  keep_list <- annotations$pseudotime_ridge_keep

  pt_filt <- pt_long %>%
    dplyr::filter(mapply(function(ln, ct) ct %in% keep_list[[ln]], lineage, cell_type))

  meds <- pt_filt %>%
    dplyr::group_by(lineage, cell_type) %>%
    dplyr::summarize(med = median(pt), .groups = "drop")

  pt_filt <- pt_filt %>%
    dplyr::left_join(meds, by = c("lineage", "cell_type")) %>%
    dplyr::mutate(cell_type_ord = tidytext::reorder_within(cell_type, med, within = lineage))

  p <- ggplot2::ggplot(pt_filt, ggplot2::aes(x = pt, y = cell_type_ord, fill = ..x..)) +
    ggridges::geom_density_ridges_gradient(
      scale = 2.2, rel_min_height = 0.01, size = 0.2, color = "white"
    ) +
    ggplot2::scale_fill_gradientn(colors = palette_fn(200), name = "Pseudotime") +
    tidytext::scale_y_reordered() +
    ggplot2::facet_wrap(~ lineage, ncol = 1, scales = "free_y") +
    ggplot2::labs(x = "Pseudotime", y = NULL,
                  title = "Pseudotime by specified cell types per lineage") +
    ggridges::theme_ridges() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   panel.grid = ggplot2::element_blank())

  save_plot(p, paste0(prefix, "_pseudotime_ridges.png"), out_dir,
            width = 8, height = 12, dpi = plotting_cfg$dimensions$dpi)
  p
}
