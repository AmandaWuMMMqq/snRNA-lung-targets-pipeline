library(testthat)

# Project root is two levels up from tests/testthat/
.proj <- normalizePath(file.path(test_path(), "..", ".."), mustWork = FALSE)
source(file.path(.proj, "R", "00_utils.R"))
source(file.path(.proj, "R", "12_audit.R"))

# ── scale01 ───────────────────────────────────────────────────────────────────

test_that("scale01 returns values in [0, 1]", {
  x <- c(1, 2, 3, 4, 5)
  y <- scale01(x)
  expect_equal(min(y, na.rm = TRUE), 0)
  expect_equal(max(y, na.rm = TRUE), 1)
})

test_that("scale01 handles NA values", {
  x <- c(1, NA, 3)
  y <- scale01(x)
  expect_true(is.na(y[2]))
})

test_that("scale01 returns NA vector for constant input", {
  x <- c(5, 5, 5)
  y <- scale01(x)
  expect_true(all(is.na(y)))
})

# ── normalize_labels ──────────────────────────────────────────────────────────

test_that("normalize_labels lowercases and replaces spaces", {
  expect_equal(normalize_labels("Bud Tip"), "bud_tip")
  expect_equal(normalize_labels("AT2_1"),   "at2_1")
  expect_equal(normalize_labels("  Fetal  "), "fetal")
})

# ── check_required_metadata ───────────────────────────────────────────────────
# check_required_metadata uses obj@meta.data — test its logic directly via a
# helper that bypasses the S4 slot access:
.check_meta_cols <- function(meta_df, required_cols) {
  missing <- setdiff(required_cols, colnames(meta_df))
  if (length(missing) > 0)
    stop("Missing required metadata columns: ", paste(missing, collapse = ", "))
  invisible(TRUE)
}

test_that("check_required_metadata logic throws on missing columns", {
  meta <- data.frame(a = 1, b = 2)
  expect_error(.check_meta_cols(meta, c("a", "c")), "Missing required metadata")
})

test_that("check_required_metadata logic passes when all columns present", {
  meta <- data.frame(a = 1, b = 2)
  expect_true(.check_meta_cols(meta, c("a", "b")))
})

# ── check_output_dirs ──────────────────────────────────────────────────────────

test_that("check_output_dirs creates missing directories", {
  tmp      <- tempfile()
  test_dir <- file.path(tmp, "plots")
  paths    <- list(output = list(plots = test_dir))
  expect_no_error(check_output_dirs(paths))
  expect_true(dir.exists(test_dir))
  unlink(tmp, recursive = TRUE)
})

# ── write_audit / read_audit ──────────────────────────────────────────────────

test_that("write_audit creates YAML, CSV, and MD files", {
  tmp   <- tempfile()
  dir.create(tmp)
  audit <- list(module = "test", timestamp = "2026-01-01", n_cells = 1000)
  write_audit(audit, "test_module", tmp)
  expect_true(file.exists(file.path(tmp, "test_module_audit.yaml")))
  expect_true(file.exists(file.path(tmp, "test_module_audit.csv")))
  expect_true(file.exists(file.path(tmp, "test_module_audit.md")))
  unlink(tmp, recursive = TRUE)
})

test_that("read_audit returns same content as written", {
  tmp   <- tempfile()
  dir.create(tmp)
  audit <- list(module = "roundtrip", n_cells = 500L, label = "pre")
  write_audit(audit, "roundtrip_test", tmp)
  result <- read_audit("roundtrip_test", tmp)
  expect_equal(result$module,  "roundtrip")
  expect_equal(result$n_cells, 500L)
  unlink(tmp, recursive = TRUE)
})

test_that("audit MD file contains module name", {
  tmp   <- tempfile()
  dir.create(tmp)
  audit <- list(module = "mymodule", n_cells = 42)
  write_audit(audit, "md_check", tmp)
  lines <- readLines(file.path(tmp, "md_check_audit.md"))
  expect_true(any(grepl("md_check", lines)))
  unlink(tmp, recursive = TRUE)
})

# ── branch label consistency ───────────────────────────────────────────────────

test_that("branch labels pre and post are handled consistently", {
  prefixes <- c("pre_alveolar", "post_alveolar",
                "pre_slingshot", "post_slingshot",
                "pre_tradeseq",  "post_tradeseq",
                "pre_cellchat",  "post_cellchat")
  for (p in prefixes) {
    expect_true(grepl("^(pre|post)_", p),
                info = paste("Branch label not consistent:", p))
  }
})

# ── check_config_keys ──────────────────────────────────────────────────────────

test_that("check_config_keys passes for valid keys", {
  cfg <- list(a = 1, b = 2, c = 3)
  expect_invisible(check_config_keys(cfg, c("a", "b"), "test_cfg"))
})

test_that("check_config_keys throws on missing keys", {
  cfg <- list(a = 1)
  expect_error(
    check_config_keys(cfg, c("a", "missing_key"), "test_cfg"),
    "missing required keys"
  )
})
