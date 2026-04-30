library(yaml)

write_audit <- function(audit_list, name, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # YAML
  yaml_path <- file.path(out_dir, paste0(name, "_audit.yaml"))
  yaml::write_yaml(audit_list, yaml_path)

  # CSV (flatten one level)
  flat <- lapply(audit_list, function(v) {
    if (length(v) > 1) paste(v, collapse = "; ") else as.character(v)
  })
  csv_path <- file.path(out_dir, paste0(name, "_audit.csv"))
  utils::write.csv(as.data.frame(flat), csv_path, row.names = FALSE)

  # Markdown
  md_lines <- c(
    paste0("# Audit: ", name),
    paste0("Generated: ", Sys.time()),
    "",
    "| Key | Value |",
    "|-----|-------|"
  )
  for (nm in names(audit_list)) {
    val <- audit_list[[nm]]
    if (is.null(val)) val <- "null"
    if (length(val) > 1) val <- paste(val, collapse = ", ")
    md_lines <- c(md_lines, paste0("| ", nm, " | ", val, " |"))
  }
  md_path <- file.path(out_dir, paste0(name, "_audit.md"))
  writeLines(md_lines, md_path)

  message("Audit written: ", name, " → ", out_dir)
  invisible(list(yaml = yaml_path, csv = csv_path, md = md_path))
}

read_audit <- function(name, out_dir) {
  yaml_path <- file.path(out_dir, paste0(name, "_audit.yaml"))
  if (!file.exists(yaml_path)) stop("Audit file not found: ", yaml_path)
  yaml::read_yaml(yaml_path)
}

list_audits <- function(out_dir) {
  if (!dir.exists(out_dir)) return(character(0))
  files <- list.files(out_dir, pattern = "_audit\\.yaml$", full.names = FALSE)
  gsub("_audit\\.yaml$", "", files)
}

consolidate_audits <- function(audit_dir) {
  names <- list_audits(audit_dir)
  if (length(names) == 0) {
    message("No audit files found in ", audit_dir)
    return(invisible(NULL))
  }
  all_audits <- lapply(names, function(nm) {
    a <- read_audit(nm, audit_dir)
    a$audit_name <- nm
    a
  })
  names(all_audits) <- names
  all_audits
}

write_fix_log <- function(out_dir = "outputs/audit") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  fixes <- list(
    list(
      issue = "Contaminant gene filter — post-merge epithelial subsetting",
      location = "R/05_subsetting.R: subset_epithelial_post()",
      original = "Original code had bare expressions like `CLDN5 < 0.1 & ...` without assignment — filters were not applied",
      fix = "Wrapped in filter_contaminants() which explicitly returns and assigns the subsetted object"
    ),
    list(
      issue = "Undefined variable 'curves' in tradeSeq fitGAM",
      location = "R/08_tradeseq.R: run_fitgam()",
      original = "Original code referenced 'curves' which was not defined; should be getCurves(sds_BT)",
      fix = "run_fitgam() takes sds_obj as an explicit argument — the caller passes getCurves result"
    ),
    list(
      issue = "Undefined variable 'sce2' in startVsEndTest",
      location = "R/08_tradeseq.R: run_startvsend_test()",
      original = "Original code referenced 'sce2' which should have been 'sce_post'",
      fix = "run_startvsend_test() takes sce as an explicit argument — no reliance on global variable"
    )
  )
  yaml_path <- file.path(out_dir, "code_fix_log.yaml")
  yaml::write_yaml(fixes, yaml_path)
  md_lines <- c("# Code Fix Log", "", "Issues from original Rmd fixed during modularization:", "")
  for (f in fixes) {
    md_lines <- c(
      md_lines,
      paste0("## ", f$issue),
      paste0("**Location:** `", f$location, "`"),
      paste0("**Original issue:** ", f$original),
      paste0("**Fix applied:** ", f$fix),
      ""
    )
  }
  md_path <- file.path(out_dir, "code_fix_log.md")
  writeLines(md_lines, md_path)
  message("Fix log written: ", yaml_path)
  invisible(list(yaml = yaml_path, md = md_path))
}
