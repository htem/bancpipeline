#' banc-updateids — Refresh root_ids in SeaTable via bancr::banctable_updateids.
#'
#' Thin wrapper around `bancr:::banctable_updateids()` that resolves stale
#' root_ids in `banc_meta` and pushes the refresh back to SeaTable.
#'
#' @section Reads:
#'   - SeaTable `banc_meta`
#'   - CAVE (root_id resolution)
#'
#' @section Writes:
#'   - SeaTable `banc_meta`: col `root_id`
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_updateids.sh`

###################################
### Update BANC IDS in seatable ###
###################################
source("banc/banc-startup.R")
bancr:::banctable_updateids()