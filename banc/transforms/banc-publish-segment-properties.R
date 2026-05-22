#!/usr/bin/env Rscript
#' banc-publish-segment-properties — Refresh per-dataset segment_properties JSON only.
#'
#' Sources just the helper functions from `banc-ngl-upload.R`, reads each
#' compiled_data feather from GCS (FANC uses local CSV), builds the
#' tag-rich `segment_properties` JSON, and pushes to the per-dataset GCS
#' precomputed layer. Cheap alternative to `banc-ngl-upload.R` when you
#' don't need to redo mesh uploads — skips the SeaTable + franken_meta
#' pulls that cost ~5 min on login nodes.
#'
#' @section Reads:
#'   - GCS `compiled_data/banc_<ver>/*.feather`
#'   - FANC: local `<banc.meta.save.path>/fanc_meta.csv`
#'
#' @section Writes:
#'   - GCS `gs://lee-lab_..._/precomputed/<dataset>/segment_properties/info`
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_publish_segment_properties.sh`

###############################################################################
### Standalone segment_properties publisher                                 ###
###                                                                         ###
### Sources just the helper functions from banc-ngl-upload.R, reads each    ###
### compiled_data feather from GCS (FANC uses local CSV), builds tag-rich   ###
### segment_properties JSON, and pushes to the per-dataset GCS layer.       ###
###                                                                         ###
### Use this when you only need to refresh annotations without doing mesh   ###
### uploads. The full banc-ngl-upload.R also does this, but pulls in        ###
### banctable + franken_meta which costs ~5 min on login nodes.             ###
###                                                                         ###
### Run:  Rscript banc/transforms/banc-publish-segment-properties.R         ###
###############################################################################
suppressPackageStartupMessages({
  library(arrow)
  library(jsonlite)
  library(dplyr)
  library(readr)
})

# Minimal startup: only what's needed for path resolution + helpers.
banc.meta.save.path <- "/n/data1/hms/neurobio/wilson/banc/meta"
version <- "elastix_tpsreg_240721"        # banc.nblast.version
malecns.version <- "navis_tpsreg_250206"  # banc.nblast.malecns.version

# Pull helpers from banc-ngl-upload.R (single source of truth).
.script <- "banc/transforms/banc-ngl-upload.R"
.src    <- readLines(.script)
.start  <- grep("^build_segment_properties <- function", .src)[1]
.end    <- grep("^OTHER_TAG_COLS <- ", .src)[1]
eval(parse(text = paste(.src[.start:.end], collapse = "\n")))

# Per-dataset config (compact table; one row per neuroglancer layer).
DATASETS <- list(
  list(name = "BANC v888 (named)",
       meta_key = "banc_888", id_col = "root_888",
       layer = "imported_meshes/banc_meshes",
       sp_subdir = "segment_properties_v888",
       tag_cols = BANC_TAG_COLS),
  list(name = "BANC v888 (default)",
       meta_key = "banc_888", id_col = "root_888",
       layer = "imported_meshes/banc_meshes",
       sp_subdir = "segment_properties",
       tag_cols = BANC_TAG_COLS),
  list(name = "MaleCNS",
       meta_key = "malecns_09", id_col = "malecns_09_id",
       layer = sprintf("imported_meshes/malecns_v0.9_meshes_%s", malecns.version),
       sp_subdir = "segment_properties",
       tag_cols = OTHER_TAG_COLS),
  list(name = "FAFB",
       meta_key = "fafb_783", id_col = "fafb_783_id",
       layer = sprintf("imported_meshes/fafb_783_meshes_%s", version),
       sp_subdir = "segment_properties",
       tag_cols = OTHER_TAG_COLS),
  list(name = "MANC",
       meta_key = "manc_121", id_col = "manc_121_id",
       layer = sprintf("imported_meshes/manc_v1.2.1_meshes_%s", version),
       sp_subdir = "segment_properties",
       tag_cols = OTHER_TAG_COLS),
  list(name = "Hemibrain",
       meta_key = "hemibrain_121", id_col = "hemibrain_121_id",
       layer = sprintf("imported_meshes/hemibrain_v1.2.1_meshes_%s", version),
       sp_subdir = "segment_properties",
       tag_cols = OTHER_TAG_COLS),
  list(name = "FANC (local meta)",
       meta_key = "fanc", id_col = "cell_id",
       layer = sprintf("imported_meshes/fanc_1116_meshes_%s", version),
       sp_subdir = "segment_properties",
       tag_cols = OTHER_TAG_COLS)
)

GCS_ROOT <- "gs://lee-lab_brain-and-nerve-cord-fly-connectome"

# Cache compiled meta reads so the two BANC entries don't double-fetch.
meta_cache <- new.env(parent = emptyenv())
get_meta <- function(key) {
  if (!exists(key, envir = meta_cache, inherits = FALSE)) {
    message("  Reading meta: ", key)
    assign(key, read_compiled_meta(key), envir = meta_cache)
  }
  get(key, envir = meta_cache)
}

for (d in DATASETS) {
  message("\n### ", d$name, " ###")
  layer_gs <- file.path(GCS_ROOT, d$layer)
  ensure_layer_info(layer_gs)
  meta <- get_meta(d$meta_key)
  json <- build_segment_properties(meta, d$id_col, "cell_type", d$tag_cols)
  gs_info  <- file.path(layer_gs, d$sp_subdir, "info")
  loc_info <- file.path("setup", d$layer, d$sp_subdir, "info")
  ok <- push_segment_properties(json, gs_info, loc_info)
  if (!isTRUE(ok)) {
    message("  WARN: gsutil cp returned non-zero for ", gs_info)
  } else {
    message("  Pushed: ", gs_info)
  }
}

cat("\n=== Done ===\n")
