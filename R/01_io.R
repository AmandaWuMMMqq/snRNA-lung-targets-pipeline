source("R/00_utils.R")

REQUIRED_META_COLS <- c(
  "batch", "sex", "age", "age_days", "group",
  "TF_annotation", "TF_annotation_lvl1", "TF_annotation_lvl2",
  "TF_annotation_lvl3", "percent.mt"
)

load_seurat_object <- function(path, label = "", validate = TRUE) {
  if (!file.exists(path)) stop("RDS file not found: ", path)
  message("Loading ", label, " from: ", path)
  obj <- readRDS(path)
  if (validate) {
    missing_cols <- setdiff(REQUIRED_META_COLS, colnames(obj@meta.data))
    if (length(missing_cols) > 0) {
      warning("Missing metadata columns in ", label, ": ",
              paste(missing_cols, collapse = ", "))
    }
  }
  message("  Loaded: ", ncol(obj), " cells, ", nrow(obj), " genes")
  obj
}

save_seurat_object <- function(obj, filename, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(out_dir, filename)
  message("Saving to: ", path)
  saveRDS(obj, file = path)
  invisible(path)
}

inspect_objects <- function(...) {
  objs <- list(...)
  nms  <- as.character(match.call())[-1]
  for (i in seq_along(objs)) {
    cat("\n--- ", nms[i], " ---\n")
    print(objs[[i]])
    cat("Metadata columns:\n")
    print(colnames(objs[[i]]@meta.data))
    cat("Cells per batch:\n")
    if ("batch" %in% colnames(objs[[i]]@meta.data))
      print(table(objs[[i]]$batch))
  }
  invisible(NULL)
}
