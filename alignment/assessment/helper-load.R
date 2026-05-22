#' helper-load — testthat helper: source alignment-data-sources.R for tests.
#'
#' Auto-loaded by testthat before each test file in `alignment/assessment/`.
#' `testthat::test_dir()` walks up from the test file to find `helper-*.R`
#' files, so this works whether you run tests from the repo root or
#' anywhere else.
#'
#' @section Reads:
#'   - `alignment/alignment-data-sources.R` (sourced relative to this file)

# Auto-loaded by testthat before each test file is run.
# testthat::test_dir() walks up from the test file to find helper-*.R files,
# so this works whether you run tests from the repo root or anywhere else.

# Resolve the alignment dir relative to this helper, which testthat sources
# with a known path.
.assessment_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile, mustWork = FALSE)),
  error = function(e) getwd()
)
.alignment_dir <- normalizePath(file.path(.assessment_dir, ".."),
                                mustWork = FALSE)

source(file.path(.alignment_dir, "alignment-data-sources.R"))
