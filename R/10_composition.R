source("R/00_utils.R")

compute_composition_tables <- function(obj, cell_type_col = "cell_type3",
                                        group_col = "group",
                                        batch_col = "batch",
                                        age_col = "age_days") {
  check_required_metadata(obj, c(cell_type_col, group_col, batch_col))

  meta <- obj@meta.data
  total_df <- meta %>%
    dplyr::group_by(.data[[batch_col]]) %>%
    dplyr::summarise(total_cells = dplyr::n(), .groups = "drop")

  pct_table <- meta %>%
    dplyr::filter(.data[[group_col]] %in% c("fetal", "pediatric")) %>%
    dplyr::group_by(.data[[batch_col]], .data[[group_col]],
                    .data[[age_col]], .data[[cell_type_col]]) %>%
    dplyr::summarise(celltype_cells = dplyr::n(), .groups = "drop") %>%
    dplyr::left_join(total_df, by = batch_col) %>%
    dplyr::mutate(pct = 100 * celltype_cells / total_cells) %>%
    dplyr::filter(!is.na(.data[[age_col]]), !is.na(pct))

  list(
    total_per_batch = total_df,
    pct_table       = pct_table
  )
}

plot_age_distributions <- function(obj, celltypes, cell_type_col = "cell_type3",
                                    group_col = "group", age_col = "age_days",
                                    title, out_dir, prefix, plotting_cfg) {
  meta <- obj@meta.data

  plot_df <- meta[
    meta[[cell_type_col]] %in% celltypes & meta[[group_col]] %in% c("fetal", "pediatric"),
  ]
  plot_df[[cell_type_col]] <- factor(plot_df[[cell_type_col]], levels = celltypes)
  plot_df[[group_col]]     <- factor(plot_df[[group_col]], levels = c("fetal", "pediatric"))

  p <- ggplot2::ggplot(plot_df,
                        ggplot2::aes(x = .data[[age_col]],
                                     fill = .data[[group_col]])) +
    ggplot2::geom_histogram(bins = 25, color = "black") +
    ggplot2::facet_grid(
      stats::as.formula(paste(cell_type_col, "~", group_col)),
      scales = "free_x"
    ) +
    ggplot2::scale_fill_manual(
      values = c("fetal" = plotting_cfg$color_palettes$group_fetal,
                  "pediatric" = plotting_cfg$color_palettes$group_pediatric)
    ) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::labs(title = title, x = "Age (days)", y = "Number of cells") +
    ggplot2::guides(fill = "none")

  save_plot(p, paste0(prefix, "_age_distribution.png"), out_dir,
            width = plotting_cfg$dimensions$default_width,
            height = plotting_cfg$dimensions$tall_height,
            dpi = plotting_cfg$dimensions$dpi)
  p
}

plot_composition_bars <- function(pct_table, cell_type_col = "cell_type3",
                                   group_col = "group", batch_col = "batch",
                                   age_col = "age_days",
                                   title, out_dir, prefix, plotting_cfg,
                                   filter_celltype = NULL) {
  plot_df <- pct_table
  if (!is.null(filter_celltype)) {
    plot_df <- dplyr::filter(plot_df, .data[[cell_type_col]] == filter_celltype)
  }

  plot_df <- plot_df %>%
    dplyr::arrange(.data[[cell_type_col]], .data[[group_col]], .data[[age_col]]) %>%
    dplyr::mutate(sample_age = paste0(.data[[batch_col]], " (", .data[[age_col]], ")")) %>%
    dplyr::group_by(.data[[cell_type_col]], .data[[group_col]]) %>%
    dplyr::mutate(sample_age = factor(sample_age, levels = unique(sample_age))) %>%
    dplyr::ungroup()

  p <- ggplot2::ggplot(plot_df,
                        ggplot2::aes(x = sample_age, y = pct,
                                     fill = .data[[group_col]])) +
    ggplot2::geom_col(color = "black", width = 0.85) +
    ggplot2::facet_grid(stats::as.formula(paste(".", "~", group_col)),
                        scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_manual(
      values = c("fetal" = plotting_cfg$color_palettes$group_fetal,
                  "pediatric" = plotting_cfg$color_palettes$group_pediatric)
    ) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                        vjust = 0.5, hjust = 1)) +
    ggplot2::labs(title = title, x = "Sample (age_days)", y = "Cells (%)") +
    ggplot2::guides(fill = "none")

  save_plot(p, paste0(prefix, "_composition_bars.png"), out_dir,
            width = plotting_cfg$dimensions$default_width,
            height = plotting_cfg$dimensions$default_height,
            dpi = plotting_cfg$dimensions$dpi)
  p
}

plot_lineage_composition <- function(obj, sds_obj, lineage_name, celltypes,
                                      cell_type_col = "cell_type3",
                                      group_col = "group",
                                      batch_col = "batch",
                                      age_col = "age_days",
                                      out_dir, prefix, plotting_cfg) {
  w_sling <- SummarizedExperiment::assay(sds_obj, "weights")
  if (!lineage_name %in% colnames(w_sling)) {
    warning("Lineage '", lineage_name, "' not found in slingshot weights — skipping")
    return(NULL)
  }
  cells_lin <- rownames(w_sling)[w_sling[, lineage_name] > 0]

  meta_df <- obj@meta.data %>%
    dplyr::mutate(cell = rownames(obj@meta.data)) %>%
    dplyr::filter(cell %in% cells_lin,
                  .data[[cell_type_col]] %in% celltypes,
                  .data[[group_col]] %in% c("fetal", "pediatric")) %>%
    dplyr::group_by(.data[[cell_type_col]], .data[[group_col]],
                    .data[[batch_col]], .data[[age_col]]) %>%
    dplyr::summarise(n_cells = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(.data[[cell_type_col]], .data[[group_col]], .data[[age_col]]) %>%
    dplyr::group_by(.data[[cell_type_col]], .data[[group_col]]) %>%
    dplyr::mutate(
      sample_age = factor(
        paste0(.data[[batch_col]], " (", .data[[age_col]], ")"),
        levels = unique(paste0(.data[[batch_col]], " (", .data[[age_col]], ")"))
      )
    ) %>%
    dplyr::ungroup()

  meta_df[[cell_type_col]] <- factor(meta_df[[cell_type_col]], levels = celltypes)
  meta_df[[group_col]]     <- factor(meta_df[[group_col]], levels = c("fetal", "pediatric"))

  p <- ggplot2::ggplot(meta_df,
                        ggplot2::aes(x = sample_age, y = n_cells,
                                     fill = .data[[group_col]])) +
    ggplot2::geom_col(color = "black", width = 0.85) +
    ggplot2::facet_grid(
      stats::as.formula(paste(cell_type_col, "~", group_col)),
      scales = "free_x", space = "free_x"
    ) +
    ggplot2::scale_fill_manual(
      values = c("fetal" = plotting_cfg$color_palettes$group_fetal,
                  "pediatric" = plotting_cfg$color_palettes$group_pediatric)
    ) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                        vjust = 0.5, hjust = 1)) +
    ggplot2::labs(
      title = paste0("Cell counts per sample (", lineage_name, " by weights)"),
      x = "Sample (age_days)", y = "Number of cells"
    ) +
    ggplot2::guides(fill = "none")

  save_plot(p, paste0(prefix, "_", lineage_name, "_lineage_composition.png"), out_dir,
            width = plotting_cfg$dimensions$default_width,
            height = plotting_cfg$dimensions$default_height,
            dpi = plotting_cfg$dimensions$dpi)
  p
}
