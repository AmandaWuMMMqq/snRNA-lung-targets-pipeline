library(testthat)
library(yaml)

# Project root is two levels up from tests/testthat/
.proj <- normalizePath(file.path(test_path(), "..", ".."), mustWork = FALSE)
.cfg  <- function(f) file.path(.proj, "config", f)

test_that("config files can be loaded", {
  expect_no_error(yaml::read_yaml(.cfg("paths.yaml")))
  expect_no_error(yaml::read_yaml(.cfg("parameters.yaml")))
  expect_no_error(yaml::read_yaml(.cfg("markers.yaml")))
  expect_no_error(yaml::read_yaml(.cfg("annotations.yaml")))
  expect_no_error(yaml::read_yaml(.cfg("plotting.yaml")))
})

test_that("paths.yaml has required keys", {
  cfg <- yaml::read_yaml(.cfg("paths.yaml"))
  expect_true("input"  %in% names(cfg))
  expect_true("output" %in% names(cfg))
  expect_true("fetal_rds"     %in% names(cfg$input))
  expect_true("pediatric_rds" %in% names(cfg$input))
  expect_true("root"  %in% names(cfg$output))
  expect_true("rds"   %in% names(cfg$output))
  expect_true("plots" %in% names(cfg$output))
  expect_true("audit" %in% names(cfg$output))
})

test_that("parameters.yaml has required keys", {
  cfg <- yaml::read_yaml(.cfg("parameters.yaml"))
  expect_true("seurat"     %in% names(cfg))
  expect_true("trajectory" %in% names(cfg))
  expect_true("tradeseq"   %in% names(cfg))
  expect_true("cellchat"   %in% names(cfg))
  expect_true("dims"       %in% names(cfg$seurat))
  expect_true("resolution" %in% names(cfg$seurat))
  expect_gt(cfg$seurat$dims, 0)
  expect_gt(cfg$seurat$resolution, 0)
})

test_that("markers.yaml has expected lineage markers", {
  cfg <- yaml::read_yaml(.cfg("markers.yaml"))
  expect_true("lineage_markers" %in% names(cfg))
  expect_true("epithelium"  %in% names(cfg$lineage_markers))
  expect_true("mesenchyme"  %in% names(cfg$lineage_markers))
  expect_true("immune"      %in% names(cfg$lineage_markers))
  expect_true("endothelium" %in% names(cfg$lineage_markers))
})

test_that("annotations.yaml has required mapping tables", {
  cfg <- yaml::read_yaml(.cfg("annotations.yaml"))
  expect_true("fetal_lineage"     %in% names(cfg))
  expect_true("pediatric_lineage" %in% names(cfg))
  expect_true("joint_celltype"    %in% names(cfg))
  expect_true("joint_lineage"     %in% names(cfg))
})

test_that("plotting.yaml has required keys", {
  cfg <- yaml::read_yaml(.cfg("plotting.yaml"))
  expect_true("dimensions"     %in% names(cfg))
  expect_true("color_palettes" %in% names(cfg))
  expect_gt(cfg$dimensions$dpi, 0)
})
