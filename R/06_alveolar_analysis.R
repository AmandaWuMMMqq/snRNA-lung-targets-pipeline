source("R/00_utils.R")

plot_alveolar_umap_panel <- function(obj, group_vars, out_dir, prefix,
                                      plotting_cfg) {
  plots <- lapply(group_vars, function(gv) {
    if (!gv %in% colnames(obj@meta.data)) {
      message("Skipping UMAP for missing column: ", gv)
      return(NULL)
    }
    Seurat::DimPlot(obj, group.by = gv, reduction = "umap",
                    label = FALSE, raster = FALSE) +
      ggplot2::theme_void() +
      ggplot2::ggtitle(paste0(prefix, " – ", gv))
  })
  plots <- Filter(Negate(is.null), plots)

  combined <- patchwork::wrap_plots(plots, ncol = 2)
  save_plot(
    combined,
    filename = paste0(prefix, "_umap_panel.png"),
    out_dir  = out_dir,
    width    = plotting_cfg$dimensions$umap_width,
    height   = plotting_cfg$dimensions$umap_height,
    dpi      = plotting_cfg$dimensions$dpi
  )
  plots
}

plot_age_gradient_umap <- function(obj, out_dir, prefix, plotting_cfg) {
  meta    <- obj@meta.data %>% tibble::rownames_to_column("cell")
  umap_df <- Seurat::Embeddings(obj, "umap") %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell")
  plot_df <- dplyr::left_join(umap_df, meta, by = "cell") %>%
    dplyr::filter(group %in% c("fetal", "pediatric"), !is.na(age_days))

  fetal_pal <- RColorBrewer::brewer.pal(9, plotting_cfg$color_palettes$age_gradient$fetal_palette)
  ped_pal   <- RColorBrewer::brewer.pal(9, plotting_cfg$color_palettes$age_gradient$pediatric_palette)

  p <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = dplyr::filter(plot_df, group == "fetal"),
      ggplot2::aes(umap_1, umap_2, color = age_days), size = 0.5
    ) +
    ggplot2::scale_color_gradientn(
      colors = fetal_pal, name = "Fetal Age (days)",
      limits = range(dplyr::filter(plot_df, group == "fetal")$age_days, na.rm = TRUE),
      na.value = "transparent"
    ) +
    ggnewscale::new_scale_color() +
    ggplot2::geom_point(
      data = dplyr::filter(plot_df, group == "pediatric"),
      ggplot2::aes(umap_1, umap_2, color = age_days), size = 0.5
    ) +
    ggplot2::scale_color_gradientn(
      colors = ped_pal, name = "Pediatric Age (days)",
      limits = range(dplyr::filter(plot_df, group == "pediatric")$age_days, na.rm = TRUE),
      na.value = "transparent"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.text  = ggplot2::element_text(size = 16),
      legend.title = ggplot2::element_text(size = 18, face = "bold")
    )

  save_plot(p, paste0(prefix, "_age_gradient_umap.png"), out_dir,
            width = 8, height = 4, dpi = plotting_cfg$dimensions$dpi)
  p
}

plot_feature_panel <- function(obj, genes, group_name, out_dir, prefix,
                                my_colors_scale, plotting_cfg, ncol = 2) {
  Seurat::DefaultAssay(obj) <- "SCT"
  plots <- Seurat::FeaturePlot(obj, features = genes,
                                pt.size = plotting_cfg$feature_plot$pt_size,
                                raster = FALSE, combine = FALSE)
  plots <- lapply(seq_along(plots), function(i) {
    plots[[i]] + my_colors_scale + ggplot2::theme_void() +
      ggplot2::ggtitle(genes[i]) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 12))
  })
  combined <- patchwork::wrap_plots(plots, ncol = ncol) +
    patchwork::plot_annotation(
      title = paste(group_name, "Markers"),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 20, face = "bold"))
    )
  save_plot(combined, paste0(prefix, "_feature_panel_", group_name, ".png"),
            out_dir, width = 8, height = 12, dpi = plotting_cfg$dimensions$dpi)
  combined
}

compute_dotplot_data <- function(obj, genes_of_interest, cell_type_col) {
  if (!cell_type_col %in% colnames(obj@meta.data)) stop("Column not found: ", cell_type_col)
  Seurat::FetchData(obj, vars = c(cell_type_col, genes_of_interest)) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(genes_of_interest),
                        names_to = "gene", values_to = "expression") %>%
    dplyr::group_by(.data[[cell_type_col]], gene) %>%
    dplyr::summarise(
      avg_expression = mean(expression, na.rm = TRUE),
      pct_cells      = mean(expression > 0) * 100,
      .groups        = "drop"
    ) %>%
    dplyr::rename(cell_type = 1)
}

plot_dotplot <- function(dotplot_data, genes_of_interest, out_dir, prefix,
                          plotting_cfg) {
  dotplot_data$gene <- factor(dotplot_data$gene, levels = genes_of_interest)

  dp <- ggplot2::ggplot(
    dotplot_data,
    ggplot2::aes(x = gene, y = cell_type, size = pct_cells, color = avg_expression)
  ) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradientn(
      colors = c("lightgrey", "#EDF8B1", "#C7E9B4", "#7FCDBB",
                  "#41B6C4", "#1D91C0", "#225EA8", "#253494", "#081D58")
    ) +
    ggplot2::scale_size_area(max_size = 9) +
    ggplot2::labs(x = "Genes", y = "Cell Type",
                  size = "% Expressing Cells", color = "Avg Expression",
                  title = "Genes of Interest") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text     = ggplot2::element_text(size = 14),
      axis.title    = ggplot2::element_text(size = 16),
      axis.text.x   = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid    = ggplot2::element_blank()
    )

  save_plot(dp, paste0(prefix, "_dotplot.png"), out_dir,
            width = 12, height = 6, dpi = plotting_cfg$dimensions$dpi)
  dp
}

find_deg_markers <- function(obj, ident_col, ident1, ident2 = NULL) {
  Seurat::Idents(obj) <- obj@meta.data[[ident_col]]
  if (is.null(ident2)) {
    Seurat::FindMarkers(obj, ident.1 = ident1,
                         only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  } else {
    Seurat::FindMarkers(obj, ident.1 = ident1, ident.2 = ident2)
  }
}
