#!/usr/bin/env Rscript
# Safely clean selected outputs.
# Does NOT delete raw data, archive/, or config/.
# Usage: Rscript scripts/clean_outputs.R [--targets] [--plots] [--tables] [--audit] [--reports] [--all]

args <- commandArgs(trailingOnly = TRUE)

safe_delete_dir <- function(path, label) {
  if (dir.exists(path)) {
    unlink(path, recursive = TRUE)
    cat("Removed:", path, "\n")
  } else {
    cat("Not found (skipped):", path, "\n")
  }
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  cat("Recreated empty directory:", path, "\n")
}

if ("--targets" %in% args || "--all" %in% args) {
  if (dir.exists("_targets")) {
    cat("WARNING: Deleting _targets/ will force a full rerun of the pipeline.\n")
    cat("Proceed? (yes/no): ")
    ans <- readLines("stdin", n = 1)
    if (tolower(trimws(ans)) == "yes") {
      unlink("_targets", recursive = TRUE)
      cat("Removed _targets/\n")
    } else {
      cat("Skipped _targets/\n")
    }
  }
}

if ("--plots" %in% args || "--all" %in% args) {
  safe_delete_dir("outputs/plots", "plots")
}
if ("--tables" %in% args || "--all" %in% args) {
  safe_delete_dir("outputs/tables", "tables")
}
if ("--audit" %in% args || "--all" %in% args) {
  safe_delete_dir("outputs/audit", "audit")
}
if ("--reports" %in% args || "--all" %in% args) {
  safe_delete_dir("outputs/reports", "reports")
}

if (length(args) == 0) {
  cat("No flags provided. Usage:\n")
  cat("  Rscript scripts/clean_outputs.R [--targets] [--plots] [--tables] [--audit] [--reports] [--all]\n")
  cat("\nNOTE: Raw data, archive/, and config/ are never deleted by this script.\n")
}
