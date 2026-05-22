#' banc-fix-major-cell-type-errors — Integrate external tracing decisions (ORN + others).
#'
#' Manual spot-fixer; not part of any automated run.
#'
#' @section Reads:
#'   - `banc_edgelist()`, SeaTable `banc_meta`
#'
#' @section Writes:
#'   - SeaTable `banc_meta` (per-block targeted columns)
#'
#' @section Notes:
#'   - Review each block before uncommenting the push.

#############################################
### BANC INTEGRATE EXTERNAL TRACING WORK ###
#############################################
source("banc/banc-startup.R")
redo <- FALSE
version <- banc.nblast.version

# Get edgelist!
banc.el <- banc_edgelist() 

# Get meta!
banc.meta <- banctable_query()

#################
### ORN FIXES ###
#################