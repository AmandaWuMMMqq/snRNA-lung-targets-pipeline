# Set the working directory to the project root for all tests.
# testthat changes cwd to the test file directory; this restores it.
proj_root <- normalizePath(file.path(dirname(dirname(getwd())), ".."), mustWork = FALSE)

# Walk up from tests/testthat until we find _targets.R
find_proj_root <- function(start) {
  path <- normalizePath(start, mustWork = FALSE)
  for (i in seq_len(10)) {
    if (file.exists(file.path(path, "_targets.R"))) return(path)
    parent <- dirname(path)
    if (parent == path) break
    path <- parent
  }
  stop("Could not find project root (_targets.R) from: ", start)
}

.proj_root <- tryCatch(
  find_proj_root(getwd()),
  error = function(e) NULL
)

if (!is.null(.proj_root)) {
  setwd(.proj_root)
}
