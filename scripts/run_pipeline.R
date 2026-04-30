#!/usr/bin/env Rscript
# Run the full targets pipeline.
# Usage: Rscript scripts/run_pipeline.R

library(targets)

options(future.globals.maxSize = 40 * 1024^3)
set.seed(123)

cat("Starting full pipeline...\n")
tar_make()
cat("Pipeline complete.\n")
