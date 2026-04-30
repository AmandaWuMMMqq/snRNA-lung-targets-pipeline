source("R/00_utils.R")

create_cellchat_object <- function(obj, label_col, sample_col = "batch",
                                    params) {
  check_required_metadata(obj, c(label_col, sample_col))
  Seurat::DefaultAssay(obj) <- "SCT"
  data.input <- Seurat::GetAssayData(obj, assay = "SCT", layer = "data")

  labels <- obj@meta.data[[label_col]]
  batch  <- as.factor(obj@meta.data[[sample_col]])
  meta   <- data.frame(labels = labels, samples = batch,
                        row.names = colnames(obj))

  message("Creating CellChat object ...")
  cc <- CellChat::createCellChat(object = data.input, meta = meta)

  db <- if (params$cellchat$database == "human") {
    CellChat::CellChatDB.human
  } else {
    CellChat::CellChatDB.mouse
  }
  cc@DB <- db
  cc
}

run_cellchat_pipeline <- function(cc, params, audit_dir = NULL) {
  message("CellChat: subsetData ...")
  cc <- CellChat::subsetData(cc)

  message("CellChat: identifyOverExpressedGenes ...")
  future::plan("multisession", workers = params$cellchat$workers %||% 4)
  cc <- CellChat::identifyOverExpressedGenes(cc)
  cc <- CellChat::identifyOverExpressedInteractions(cc)
  cc <- CellChat::smoothData(cc, adj = CellChat::PPI.human)

  message("CellChat: computeCommunProb ...")
  future::plan("sequential")
  options(future.globals.maxSize = 8 * 1024^3)
  cc <- CellChat::computeCommunProb(
    cc, type = params$cellchat$compute_prob_type,
    trim = NULL, raw.use = FALSE
  )

  message("CellChat: filterCommunication + aggregateNet ...")
  cc <- CellChat::filterCommunication(cc, min.cells = params$cellchat$min_cells)
  cc <- CellChat::computeCommunProbPathway(cc)
  cc <- CellChat::aggregateNet(cc)

  message("CellChat: netAnalysis_computeCentrality ...")
  cc <- CellChat::netAnalysis_computeCentrality(cc, slot.name = "netP")

  message("CellChat: computeNetSimilarity (functional + structural) ...")
  cc <- CellChat::computeNetSimilarity(cc, type = "functional")
  cc <- CellChat::netEmbedding(cc,   type = "functional")
  cc <- CellChat::netClustering(cc,  type = "functional")

  cc <- CellChat::computeNetSimilarity(cc, type = "structural")
  cc <- CellChat::netEmbedding(cc,  type = "structural")
  cc <- CellChat::netClustering(cc, type = "structural")

  if (!is.null(audit_dir)) {
    audit <- list(
      module    = "cellchat",
      timestamp = as.character(Sys.time()),
      n_pathways = length(cc@netP$pathways),
      pathways  = cc@netP$pathways,
      params    = params$cellchat
    )
    write_audit(audit, name = "cellchat", out_dir = audit_dir)
  }

  cc
}

reorder_cellchat_clusters <- function(cc, cluster_order) {
  CellChat::updateClusterLabels(
    object           = cc,
    old.cluster.name = cluster_order,
    new.cluster.name = cluster_order,
    new.order        = cluster_order,
    new.cluster.metaname = "celltype_ordered"
  )
}

add_seurat_umap_to_cellchat <- function(cc, seurat_obj) {
  dimred_cc <- seurat_obj@reductions$umap@cell.embeddings
  cc        <- CellChat::addReduction(cc, dr = dimred_cc, dr.name = "umap", force.add = TRUE)
  cc@meta$samples <- factor("all")
  cc
}

identify_outgoing_patterns <- function(cc, n_patterns, params) {
  CellChat::selectK(cc, pattern = "outgoing")
  CellChat::identifyCommunicationPatterns(cc, pattern = "outgoing", k = n_patterns)
}

identify_incoming_patterns <- function(cc, n_patterns, params) {
  CellChat::selectK(cc, pattern = "incoming")
  CellChat::identifyCommunicationPatterns(cc, pattern = "incoming", k = n_patterns)
}

plot_cellchat_pathway <- function(cc, pathway, out_dir, prefix) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Circle plot
  p_path <- file.path(out_dir, paste0(prefix, "_", pathway, "_circle.png"))
  grDevices::png(p_path, width = 800, height = 800)
  CellChat::netVisual_aggregate(cc, signaling = pathway, layout = "circle")
  graphics::title(paste(pathway, "Signaling – Circle"), line = 1, cex.main = 1.2)
  grDevices::dev.off()
  message("Saved: ", p_path)

  # Heatmap
  p_heat <- CellChat::netVisual_heatmap(cc, signaling = pathway, color.heatmap = "Reds")
  save_plot(p_heat, paste0(prefix, "_", pathway, "_heatmap.png"), out_dir)
}

plot_cellchat_bubble_panel <- function(cc, sources, targets, pathways,
                                        out_dir, prefix, plotting_cfg) {
  for (pw in pathways) {
    p <- CellChat::netVisual_bubble(
      cc,
      sources.use    = sources,
      targets.use    = targets,
      signaling      = pw,
      remove.isolate = FALSE
    )
    save_plot(p, paste0(prefix, "_bubble_", pw, ".png"), out_dir,
              width = plotting_cfg$dimensions$umap_width,
              height = plotting_cfg$dimensions$default_height,
              dpi = plotting_cfg$dimensions$dpi)
  }
}

`%||%` <- function(x, y) if (is.null(x)) y else x
