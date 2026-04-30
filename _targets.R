library(targets)
library(tarchetypes)

tar_option_set(
  packages = c(
    "dplyr", "tidyr", "tibble", "ggplot2", "patchwork",
    "Seurat", "SeuratWrappers", "sctransform",
    "slingshot", "SingleCellExperiment", "tradeSeq", "BiocParallel",
    "RColorBrewer", "viridis", "ComplexHeatmap", "pheatmap",
    "ggridges", "ggrepel", "ggnewscale", "cowplot", "circlize",
    "presto", "writexl",
    "CellChat", "NMF", "ggalluvial",
    "yaml", "R.utils", "tidytext"
  )
)

# Source all R modules
lapply(list.files("R", pattern = "\\.R$", full.names = TRUE), source)

list(
  # ── Config ─────────────────────────────────────────────────────────────────
  tar_target(paths,       load_config("config/paths.yaml")),
  tar_target(params,      load_config("config/parameters.yaml")),
  tar_target(markers,     load_config("config/markers.yaml")),
  tar_target(annotations, load_config("config/annotations.yaml")),
  tar_target(plotting,    load_config("config/plotting.yaml")),

  # Ensure output directories exist
  tar_target(output_dirs, check_output_dirs(paths), deployment = "main"),

  # Write the code fix log once at the start
  tar_target(fix_log, write_fix_log(paths$output$audit), deployment = "main"),

  # ── Input data ─────────────────────────────────────────────────────────────
  tar_target(
    fetal_raw,
    load_seurat_object(paths$input$fetal_rds, label = "fetal", validate = TRUE)
  ),
  tar_target(
    pediatric_raw,
    load_seurat_object(paths$input$pediatric_rds, label = "pediatric", validate = TRUE)
  ),
  tar_target(
    input_audit, {
      audit <- list(
        module    = "input",
        timestamp = as.character(Sys.time()),
        fetal_n_cells = ncol(fetal_raw),
        fetal_n_genes = nrow(fetal_raw),
        fetal_batches = paste(unique(fetal_raw$batch), collapse = ", "),
        ped_n_cells   = ncol(pediatric_raw),
        ped_n_genes   = nrow(pediatric_raw),
        ped_batches   = paste(unique(pediatric_raw$batch), collapse = ", ")
      )
      write_audit(audit, "input", paths$output$audit)
    }
  ),

  # ── Individual integration ─────────────────────────────────────────────────
  tar_target(
    fetal_integrated,
    integrate_seurat(fetal_raw, params, split_by = "batch",
                     audit_dir = paths$output$audit, label = "fetal")
  ),
  tar_target(
    pediatric_downsampled,
    downsample_seurat(pediatric_raw, n = params$pediatric$downsample_n,
                      seed = params$random_seed)
  ),
  tar_target(
    pediatric_integrated,
    integrate_seurat(pediatric_downsampled, params, split_by = "batch",
                     audit_dir = paths$output$audit, label = "pediatric")
  ),
  tar_target(
    individual_integration_audit, {
      audit <- list(
        module    = "individual_integration",
        timestamp = as.character(Sys.time()),
        fetal_n_cells     = ncol(fetal_integrated),
        fetal_reductions  = paste(names(fetal_integrated@reductions), collapse = ", "),
        ped_n_cells       = ncol(pediatric_integrated),
        ped_reductions    = paste(names(pediatric_integrated@reductions), collapse = ", "),
        dims              = params$seurat$dims,
        resolution        = params$seurat$resolution
      )
      write_audit(audit, "individual_integration", paths$output$audit)
    }
  ),

  # ── Lineage annotation ─────────────────────────────────────────────────────
  tar_target(
    fetal_annotated, {
      obj <- apply_cluster_annotation(
        fetal_integrated,
        mapping  = unlist(annotations$fetal_lineage),
        new_col  = "lineage_lvl1"
      )
      obj <- store_cluster_ids(obj)
      obj
    }
  ),
  tar_target(
    pediatric_annotated, {
      obj <- apply_cluster_annotation(
        pediatric_integrated,
        mapping  = unlist(annotations$pediatric_lineage),
        new_col  = "lineage_lvl1"
      )
      obj <- store_cluster_ids(obj)
      obj
    }
  ),
  tar_target(
    annotation_audit, {
      audit <- list(
        module         = "annotation",
        timestamp      = as.character(Sys.time()),
        fetal_lineages = paste(unique(na.omit(fetal_annotated$lineage_lvl1)), collapse = ", "),
        ped_lineages   = paste(unique(na.omit(pediatric_annotated$lineage_lvl1)), collapse = ", ")
      )
      write_audit(audit, "annotation", paths$output$audit)
    }
  ),

  # ── Joint fetal–pediatric integration ─────────────────────────────────────
  tar_target(
    joint_integrated, {
      merged <- merge_objects(fetal_annotated, pediatric_annotated)
      integrate_seurat(merged, params, split_by = "batch",
                       audit_dir = paths$output$audit, label = "joint")
    }
  ),
  tar_target(
    joint_annotated, {
      obj <- apply_cluster_annotation(
        joint_integrated,
        mapping = unlist(annotations$joint_celltype),
        new_col = "cell_type"
      )
      obj <- apply_cluster_annotation(
        obj,
        mapping = unlist(annotations$joint_lineage),
        new_col = "lineage_lvl2"
      )
      obj <- store_cluster_ids(obj)
      obj
    }
  ),
  tar_target(
    joint_integration_audit, {
      audit <- list(
        module     = "joint_integration",
        timestamp  = as.character(Sys.time()),
        n_cells    = ncol(joint_annotated),
        cell_types = paste(unique(na.omit(joint_annotated$cell_type)), collapse = "; "),
        lineages   = paste(unique(na.omit(joint_annotated$lineage_lvl2)), collapse = ", ")
      )
      write_audit(audit, "joint_integration", paths$output$audit)
    }
  ),

  # ── Epithelial branches ────────────────────────────────────────────────────

  # Pre-merge: subset from individual objects, then merge + integrate
  tar_target(
    pre_epi_merged,
    subset_epithelial_pre(fetal_annotated, pediatric_annotated, params,
                           audit_dir = paths$output$audit)
  ),
  tar_target(
    pre_epi_integrated,
    integrate_seurat(pre_epi_merged, params, split_by = "batch",
                     audit_dir = paths$output$audit, label = "pre_epi")
  ),
  tar_target(
    pre_epi, {
      obj <- apply_cluster_annotation(
        pre_epi_integrated,
        mapping = unlist(annotations$epi_post_celltype),
        new_col = "cell_type2"
      )
      obj <- store_cluster_ids(obj)
      obj
    }
  ),

  # Post-merge: subset from joint integrated object
  tar_target(
    post_epi_raw,
    subset_epithelial_post(joint_annotated, params,
                            lineage_col = "lineage_lvl2",
                            audit_dir   = paths$output$audit)
  ),
  tar_target(
    post_epi_integrated,
    integrate_seurat(post_epi_raw, params, split_by = "batch",
                     audit_dir = paths$output$audit, label = "post_epi")
  ),
  tar_target(
    post_epi, {
      obj <- apply_cluster_annotation(
        post_epi_integrated,
        mapping = unlist(annotations$epi_post_celltype),
        new_col = "cell_type2"
      )
      obj <- store_cluster_ids(obj)
      obj
    }
  ),
  tar_target(
    epithelial_branch_audit, {
      audit <- list(
        module        = "epithelial_branch",
        timestamp     = as.character(Sys.time()),
        pre_n_cells   = ncol(pre_epi),
        post_n_cells  = ncol(post_epi),
        pre_cell_types  = paste(unique(na.omit(pre_epi$cell_type2)),  collapse = "; "),
        post_cell_types = paste(unique(na.omit(post_epi$cell_type2)), collapse = "; ")
      )
      write_audit(audit, "epithelial_branch", paths$output$audit)
    }
  ),

  # ── Alveolar branches ─────────────────────────────────────────────────────
  tar_target(
    pre_alveolar, {
      obj <- apply_cluster_annotation(
        pre_epi,
        mapping = unlist(annotations$alveolar_pre_celltype2),
        new_col = "cell_type2",
        cluster_col = "seurat_clusters"
      )
      subset_alveolar(obj, "cell_type2",
                       alveolar_labels = params$subsetting$alveolar_labels,
                       audit_dir = paths$output$audit, label = "pre")
    }
  ),
  tar_target(
    post_alveolar, {
      obj <- apply_cluster_annotation(
        post_epi,
        mapping = unlist(annotations$alveolar_post_alve_celltype2),
        new_col = "cell_type2",
        cluster_col = "seurat_clusters"
      )
      subset_alveolar(obj, "cell_type2",
                       alveolar_labels = params$subsetting$alveolar_labels,
                       audit_dir = paths$output$audit, label = "post")
    }
  ),
  tar_target(
    alveolar_audit, {
      audit <- list(
        module        = "alveolar",
        timestamp     = as.character(Sys.time()),
        pre_n_cells   = ncol(pre_alveolar),
        post_n_cells  = ncol(post_alveolar),
        pre_cell_types  = paste(unique(na.omit(pre_alveolar$cell_type2)),  collapse = "; "),
        post_cell_types = paste(unique(na.omit(post_alveolar$cell_type2)), collapse = "; ")
      )
      write_audit(audit, "alveolar", paths$output$audit)
    }
  ),

  # ── Slingshot trajectory ───────────────────────────────────────────────────
  tar_target(
    pre_slingshot_inputs,
    prepare_slingshot_inputs(pre_alveolar, reduction = "umap",
                              cluster_col = "cell_type2")
  ),
  tar_target(
    pre_slingshot, {
      sds <- run_slingshot_lineages(
        pre_slingshot_inputs$dimred,
        pre_slingshot_inputs$clustering,
        start_clus  = params$trajectory$start_cluster,
        dist_method = params$trajectory$dist_method_primary,
        params       = params
      )
      sds
    }
  ),
  tar_target(
    post_slingshot_inputs,
    prepare_slingshot_inputs(post_alveolar, reduction = "umap",
                              cluster_col = "cell_type2")
  ),
  tar_target(
    post_slingshot, {
      sds <- run_slingshot_lineages(
        post_slingshot_inputs$dimred,
        post_slingshot_inputs$clustering,
        start_clus  = params$trajectory$start_cluster,
        dist_method = params$trajectory$dist_method_primary,
        params       = params
      )
      sds
    }
  ),
  tar_target(
    trajectory_audit, {
      audit <- list(
        module          = "trajectory",
        timestamp       = as.character(Sys.time()),
        start_cluster   = params$trajectory$start_cluster,
        dist_method     = params$trajectory$dist_method_primary,
        pre_n_lineages  = length(pre_slingshot@lineages),
        post_n_lineages = length(post_slingshot@lineages),
        pre_lineage_names  = paste(names(pre_slingshot@lineages),  collapse = ", "),
        post_lineage_names = paste(names(post_slingshot@lineages), collapse = ", ")
      )
      write_audit(audit, "trajectory", paths$output$audit)
    }
  ),

  # ── tradeSeq GAM fitting ───────────────────────────────────────────────────
  tar_target(
    pre_tradeseq_counts,
    prepare_tradeseq_counts(pre_alveolar, pre_slingshot)
  ),
  tar_target(
    pre_tradeseq,
    run_fitgam(pre_tradeseq_counts, pre_slingshot, params,
               audit_dir = paths$output$audit, label = "pre")
  ),
  tar_target(
    post_tradeseq_counts,
    prepare_tradeseq_counts(post_alveolar, post_slingshot)
  ),
  tar_target(
    post_tradeseq,
    run_fitgam(post_tradeseq_counts, post_slingshot, params,
               audit_dir = paths$output$audit, label = "post")
  ),
  tar_target(
    pre_tradeseq_results, {
      list(
        association = run_association_test(pre_tradeseq),
        startvsend  = run_startvsend_test(pre_tradeseq),
        pattern     = run_pattern_test(pre_tradeseq),
        earlyDE     = run_earlyDE_test(pre_tradeseq),
        diffend     = run_diffend_test(pre_tradeseq)
      )
    }
  ),
  tar_target(
    post_tradeseq_results, {
      list(
        association = run_association_test(post_tradeseq),
        startvsend  = run_startvsend_test(post_tradeseq),
        pattern     = run_pattern_test(post_tradeseq),
        earlyDE     = run_earlyDE_test(post_tradeseq),
        diffend     = run_diffend_test(post_tradeseq)
      )
    }
  ),
  tar_target(
    tradeseq_audit, {
      audit <- list(
        module    = "tradeseq",
        timestamp = as.character(Sys.time()),
        pre_n_genes  = nrow(pre_tradeseq),
        post_n_genes = nrow(post_tradeseq),
        nknots    = params$tradeseq$nknots,
        cores     = params$tradeseq$cores,
        fix_note  = "sce2 variable fixed; sds_obj passed explicitly; see code_fix_log.yaml"
      )
      write_audit(audit, "tradeseq", paths$output$audit)
    }
  ),

  # ── CellChat ───────────────────────────────────────────────────────────────
  tar_target(
    pre_cellchat, {
      cc <- create_cellchat_object(pre_alveolar, label_col = "cell_type2",
                                    sample_col = "batch", params = params)
      run_cellchat_pipeline(cc, params, audit_dir = paths$output$audit)
    }
  ),
  tar_target(
    post_cellchat, {
      cc <- create_cellchat_object(post_alveolar, label_col = "cell_type2",
                                    sample_col = "batch", params = params)
      run_cellchat_pipeline(cc, params, audit_dir = paths$output$audit)
    }
  ),
  tar_target(
    cellchat_audit, {
      audit <- list(
        module    = "cellchat",
        timestamp = as.character(Sys.time()),
        pre_pathways  = paste(pre_cellchat@netP$pathways,  collapse = ", "),
        post_pathways = paste(post_cellchat@netP$pathways, collapse = ", "),
        database  = params$cellchat$database,
        min_cells = params$cellchat$min_cells
      )
      write_audit(audit, "cellchat", paths$output$audit)
    }
  ),

  # ── Cell composition ───────────────────────────────────────────────────────
  tar_target(
    composition_tables, {
      list(
        pre  = compute_composition_tables(pre_alveolar),
        post = compute_composition_tables(post_alveolar)
      )
    }
  ),
  tar_target(
    composition_plots, {
      out_dir <- file.path(paths$output$plots, "composition")
      for (branch in c("pre", "post")) {
        obj <- if (branch == "pre") pre_alveolar else post_alveolar
        tbl <- composition_tables[[branch]]

        plot_composition_bars(
          pct_table = tbl$pct_table,
          title     = paste0(branch, " – all cell types"),
          out_dir   = out_dir,
          prefix    = branch,
          plotting_cfg = plotting
        )

        for (ct in unique(na.omit(obj$cell_type3 %||% obj$cell_type2))) {
          tryCatch(
            plot_composition_bars(
              pct_table     = tbl$pct_table,
              filter_celltype = ct,
              title         = paste0(ct, " as % total (", branch, ")"),
              out_dir       = out_dir,
              prefix        = paste0(branch, "_", gsub("[^a-zA-Z0-9]", "_", ct)),
              plotting_cfg  = plotting
            ),
            error = function(e) message("Composition plot failed for ", ct, ": ", e$message)
          )
        }
      }
      out_dir
    }
  ),
  tar_target(
    composition_audit, {
      audit <- list(
        module    = "composition",
        timestamp = as.character(Sys.time()),
        pre_n_cells  = ncol(pre_alveolar),
        post_n_cells = ncol(post_alveolar)
      )
      write_audit(audit, "composition", paths$output$audit)
    }
  ),

  # ── Final plots and tables ─────────────────────────────────────────────────
  tar_target(
    final_plots, {
      # UMAP panels for both branches
      plot_umap_panel(pre_alveolar,  c("cell_type2", "group", "batch"),
                      file.path(paths$output$plots, "alveolar"), "pre_alveolar", plotting)
      plot_umap_panel(post_alveolar, c("cell_type2", "group", "batch"),
                      file.path(paths$output$plots, "alveolar"), "post_alveolar", plotting)
      "done"
    }
  ),
  tar_target(
    final_tables, {
      out_dir <- paths$output$tables
      # Save tradeSeq top genes per lineage
      for (branch in c("pre", "post")) {
        res <- if (branch == "pre") pre_tradeseq_results else post_tradeseq_results
        top <- top_genes_per_lineage(res$startvsend)
        save_tradeseq_results(top, file.path(out_dir, "tradeseq",
                                              paste0(branch, "_top_genes")))
        utils::write.csv(res$startvsend,
                          file.path(out_dir, "tradeseq",
                                    paste0(branch, "_tradeseq_startvsend_results.csv")))
        utils::write.csv(res$association,
                          file.path(out_dir, "tradeseq",
                                    paste0(branch, "_tradeseq_association_results.csv")))
        utils::write.csv(res$pattern,
                          file.path(out_dir, "tradeseq",
                                    paste0(branch, "_tradeseq_pattern_results.csv")))
      }
      "done"
    }
  ),

  # ── Reports ────────────────────────────────────────────────────────────────
  tar_render(
    module_audit_report,
    "reports/module_audit_report.qmd",
    output_dir = paths$output$reports
  ),
  tar_render(
    final_report,
    "reports/final_report.qmd",
    output_dir = paths$output$reports
  )
)
