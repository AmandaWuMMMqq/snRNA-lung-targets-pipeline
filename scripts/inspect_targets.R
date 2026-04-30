#!/usr/bin/env Rscript
# Inspect target status, outdated targets, and dependency summary.
# Usage: Rscript scripts/inspect_targets.R

library(targets)

cat("=== Target manifest ===\n")
manifest <- tar_manifest()
print(manifest)

cat("\n=== Outdated targets ===\n")
outdated <- tar_outdated()
if (length(outdated) == 0) {
  cat("All targets are up to date.\n")
} else {
  cat(paste(outdated, collapse = "\n"), "\n")
}

cat("\n=== Target status ===\n")
status <- tar_progress()
print(status)

cat("\n=== Validation ===\n")
tryCatch(
  tar_validate(),
  error = function(e) cat("Validation error:", e$message, "\n")
)
cat("Validation complete.\n")
