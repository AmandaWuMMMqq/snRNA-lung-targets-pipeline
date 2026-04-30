source("R/00_utils.R")

plot_umap_panel <- function(obj, group_vars, out_dir, prefix, plotting_cfg,
                              label_col = NULL) {
  plots <- lapply(group_vars, function(gv) {
    if (!gv %in% colnames(obj@meta.data) && gv != "seurat_clusters") {
      message("Skipping missing column: ", gv)
      return(NULL)
    }
    p <- Seurat::DimPlot(
      obj, group.by = gv, reduction = "umap",
      label = plotting_cfg$umap$label,
      raster = plotting_cfg$umap$raster,
      repel = plotting_cfg$umap$repel
    ) +
      ggplot2::theme_void() +
      ggplot2::ggtitle(paste0(prefix, " – ", gv)) +
      ggplot2::theme(
        legend.text  = ggplot2::element_text(size = plotting_cfg$umap$legend_text_size),
        legend.title = ggplot2::element_text(
          size = plotting_cfg$umap$legend_title_size, face = "bold"
        )
      )
    if (!is.null(label_col) && gv == label_col) {
      p <- Seurat::LabelClusters(plot = p, id = gv,
                                  size = plotting_cfg$umap$label_size)
    }
    p
  })
  plots <- Filter(Negate(is.null), plots)

  if (length(plots) == 0) return(invisible(NULL))
  combined <- patchwork::wrap_plots(plots, ncol = 2)
  save_plot(
    combined,
    filename = paste0(prefix, "_umap_panel.png"),
    out_dir  = out_dir,
    width    = plotting_cfg$dimensions$umap_width,
    height   = ceiling(length(plots) / 2) * 5,
    dpi      = plotting_cfg$dimensions$dpi
  )
  invisible(plots)
}

plot_deg_heatmap <- function(obj, cell_type_col, top_n_per_type = 10,
                               out_dir, prefix, plotting_cfg) {
  Seurat::Idents(obj) <- cell_type_col

  deg_list <- lapply(levels(Seurat::Idents(obj)), function(ct) {
    tryCatch(
      Seurat::FindMarkers(obj, ident.1 = ct, only.pos = TRUE,
                           min.pct = 0.25, logfc.threshold = 0.25),
      error = function(e) { message("FindMarkers failed for ", ct); NULL }
    )
  })
  names(deg_list) <- levels(Seurat::Idents(obj))

  top_genes <- unique(unlist(lapply(deg_list, function(deg) {
    if (is.null(deg) || nrow(deg) == 0) return(NULL)
    head(rownames(deg[order(deg$avg_log2FC, decreasing = TRUE), ]), top_n_per_type)
  })))

  if (length(top_genes) == 0) {
    warning("No genes found for heatmap")
    return(invisible(NULL))
  }

  expr_data <- Seurat::FetchData(obj, vars = c(top_genes, cell_type_col)) %>%
    tidyr::pivot_longer(cols = -dplyr::all_of(cell_type_col),
                        names_to = "gene", values_to = "expression")

  heatmap_data <- expr_data %>%
    dplyr::group_by(.data[[cell_type_col]], gene) %>%
    dplyr::summarise(avg_expression = mean(expression, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = cell_type_col, values_from = avg_expression) %>%
    tibble::column_to_rownames(var = "gene")

  scaled_data <- t(scale(t(as.matrix(heatmap_data))))

  ht <- ComplexHeatmap::Heatmap(
    as.matrix(scaled_data), name = "Expression",
    cluster_rows = TRUE, cluster_columns = TRUE,
    row_title = "Genes", column_title = "Cell Types",
    col = viridis::viridis(100),
    show_row_names = TRUE, show_column_names = TRUE,
    row_names_gp    = grid::gpar(fontsize = plotting_cfg$heatmap$row_names_fontsize),
    column_names_gp = grid::gpar(fontsize = plotting_cfg$heatmap$col_names_fontsize, rot = 45),
    heatmap_legend_param = list(
      title    = "Expression Level",
      title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 10)
    )
  )

  save_pdf(
    ComplexHeatmap::draw(ht),
    filename = paste0(prefix, "_deg_heatmap.pdf"),
    out_dir  = out_dir,
    width    = 10, height = 20
  )
  ht
}

plot_dimplot_labeled <- function(obj, group_col, out_dir, prefix, plotting_cfg) {
  p <- Seurat::DimPlot(
    obj, group.by = group_col,
    label = FALSE, repel = TRUE, raster = FALSE
  ) + ggplot2::ggtitle(NULL) + ggplot2::theme_void() +
    ggplot2::theme(
      legend.text  = ggplot2::element_text(size = plotting_cfg$umap$legend_text_size),
      legend.title = ggplot2::element_text(
        size = plotting_cfg$umap$legend_title_size, face = "bold"
      )
    )
  p_lab <- Seurat::LabelClusters(plot = p, id = group_col,
                                   size = plotting_cfg$umap$label_size)
  save_plot(p_lab, paste0(prefix, "_", group_col, "_labeled.png"), out_dir,
            width = plotting_cfg$dimensions$umap_width,
            height = plotting_cfg$dimensions$umap_height,
            dpi = plotting_cfg$dimensions$dpi)
  p_lab
}
