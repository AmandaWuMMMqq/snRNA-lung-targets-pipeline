library(yaml)

# ── Color helpers ──────────────────────────────────────────────────────────────

make_sequential_scale <- function(plotting_cfg) {
  ggplot2::scale_color_gradientn(colors = plotting_cfg$color_palettes$sequential)
}

make_sequential_palette_fn <- function(plotting_cfg) {
  grDevices::colorRampPalette(rev(plotting_cfg$color_palettes$sequential[-1]))
}

# ── Numeric helpers ────────────────────────────────────────────────────────────

scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[2] == rng[1])
    return(rep(NA_real_, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

combined_pt_from_sds <- function(sds, cells) {
  pt <- slingshot::slingPseudotime(sds)
  w  <- slingshot::slingCurveWeights(sds)
  pt <- pt[cells, , drop = FALSE]
  w  <- w[cells,  , drop = FALSE]
  pt_scaled <- apply(pt, 2, scale01)
  w_sum     <- rowSums(w, na.rm = TRUE)
  pt_one    <- rowSums(pt_scaled * w, na.rm = TRUE) / w_sum
  pt_one[!is.finite(pt_one)] <- NA_real_
  pt_one
}

normalize_labels <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

# ── Config loaders ─────────────────────────────────────────────────────────────

load_config <- function(file) {
  if (!file.exists(file)) stop("Config file not found: ", file)
  yaml::read_yaml(file)
}

load_all_configs <- function(config_dir = "config") {
  list(
    paths      = load_config(file.path(config_dir, "paths.yaml")),
    params     = load_config(file.path(config_dir, "parameters.yaml")),
    markers    = load_config(file.path(config_dir, "markers.yaml")),
    annotations = load_config(file.path(config_dir, "annotations.yaml")),
    plotting   = load_config(file.path(config_dir, "plotting.yaml"))
  )
}

# ── Validation helpers ─────────────────────────────────────────────────────────

check_required_metadata <- function(obj, required_cols) {
  missing <- setdiff(required_cols, colnames(obj@meta.data))
  if (length(missing) > 0) {
    stop("Missing required metadata columns: ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

check_required_assays <- function(obj, required_assays) {
  present <- names(obj@assays)
  missing <- setdiff(required_assays, present)
  if (length(missing) > 0) {
    stop("Missing required assays: ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

check_output_dirs <- function(paths) {
  dirs <- unlist(paths$output)
  for (d in dirs) {
    if (!dir.exists(d)) {
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
      message("Created output directory: ", d)
    }
  }
  invisible(TRUE)
}

check_config_keys <- function(cfg, required_keys, cfg_name = "config") {
  missing <- setdiff(required_keys, names(cfg))
  if (length(missing) > 0) {
    stop(cfg_name, " missing required keys: ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

# ── Plot helpers ───────────────────────────────────────────────────────────────

save_plot <- function(plot, filename, out_dir, width = 8, height = 6, dpi = 300) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(out_dir, filename)
  ggplot2::ggsave(path, plot = plot, width = width, height = height, dpi = dpi)
  message("Saved: ", path)
  invisible(path)
}

save_pdf <- function(expr, filename, out_dir, width = 12, height = 10) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(out_dir, filename)
  grDevices::pdf(path, width = width, height = height, useDingbats = FALSE)
  force(expr)
  grDevices::dev.off()
  message("Saved PDF: ", path)
  invisible(path)
}

# ── Object summary ─────────────────────────────────────────────────────────────

summarize_seurat <- function(obj, label = "") {
  assays_present    <- names(obj@assays)
  reductions_present <- names(obj@reductions)
  meta_cols         <- colnames(obj@meta.data)
  list(
    label      = label,
    n_cells    = ncol(obj),
    n_genes    = nrow(obj),
    assays     = assays_present,
    reductions = reductions_present,
    meta_cols  = meta_cols
  )
}
