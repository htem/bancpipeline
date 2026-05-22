#' test-data-sources — testthat suite for alignment-data-sources.R.
#'
#' Covers `parse_alignment_data_args` argument parsing (all CLI forms +
#' env-var fallback + validation errors). `resolve_alignment_paths` tests
#' are TODO — would need a `withr::with_tempdir` filesystem mock because
#' the function pokes real O2 paths.
#'
#' @section CLI:
#'   Run with: `Rscript -e 'testthat::test_dir("alignment/assessment")'`

###############################################################################
# testthat tests for alignment/alignment-data-sources.R
#
# Run with:
#   Rscript -e 'testthat::test_dir("alignment/assessment")'
#
# Covers parse_alignment_data_args argument parsing (all CLI forms + env var
# fallback + validation errors). resolve_alignment_paths tests are left as
# TODO — need a withr::with_tempdir filesystem mock because the function
# pokes real O2 paths.
###############################################################################
if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("Install testthat to run this suite: install.packages('testthat')")
}
library(testthat)
# helper-load.R (auto-sourced by testthat) has already loaded
# alignment-data-sources.R from the parent directory.

# -------------------------------------------------------------------------
# parse_alignment_data_args
# -------------------------------------------------------------------------
context("parse_alignment_data_args")

test_that("space-separated --source arg is parsed", {
  out <- parse_alignment_data_args(c("--source", "local"))
  expect_equal(out$source, "local")
  expect_equal(out$syn_source, "v3")        # default
})

test_that("equals-form --source=gcs is parsed", {
  out <- parse_alignment_data_args(c("--source=gcs"))
  expect_equal(out$source, "gcs")
})

test_that("--local-root is captured", {
  out <- parse_alignment_data_args(c("--source=local", "--local-root=/tmp/foo"))
  expect_equal(out$local_root, "/tmp/foo")
})

test_that("--syn-source=v2 round-trips", {
  out <- parse_alignment_data_args(c("--source=local", "--syn-source=v2"))
  expect_equal(out$syn_source, "v2")
})

test_that("positional args pass through", {
  out <- parse_alignment_data_args(c("some_target", "--source=local", "other"))
  expect_equal(out$positional, c("some_target", "other"))
})

test_that("BANC_ALIGNMENT_SOURCE env var supplies missing --source", {
  withr::with_envvar(c(BANC_ALIGNMENT_SOURCE = "gcs"), {
    out <- parse_alignment_data_args(character(0))
    expect_equal(out$source, "gcs")
  })
})

test_that("BANC_SYN_SOURCE env var supplies missing --syn-source", {
  withr::with_envvar(c(BANC_SYN_SOURCE = "v2", BANC_ALIGNMENT_SOURCE = "local"), {
    out <- parse_alignment_data_args(character(0))
    expect_equal(out$syn_source, "v2")
  })
})

test_that("unset source triggers explicit error", {
  withr::with_envvar(c(BANC_ALIGNMENT_SOURCE = ""), {
    expect_error(parse_alignment_data_args(character(0)),
                 regexp = "Data source must be specified")
  })
})

test_that("invalid --source value is rejected", {
  expect_error(parse_alignment_data_args(c("--source=remote")),
               regexp = "must be 'local' or 'gcs'")
})

test_that("invalid --syn-source value is rejected", {
  expect_error(parse_alignment_data_args(c("--source=local", "--syn-source=v1")),
               regexp = "must be 'v2' or 'v3'")
})

test_that("--source without value errors", {
  expect_error(parse_alignment_data_args(c("--source")),
               regexp = "--source requires a value")
})

# -------------------------------------------------------------------------
# resolve_alignment_paths  —  TODO
# -------------------------------------------------------------------------
# Left deferred: requires mocking banc.versioned.save.path + banc.meta.save.path
# via withr::with_tempdir + an options()/assign stub. Safe to write once the
# alignment sweep lands; until then we only need the parse-arg safety net.
