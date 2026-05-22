#' banc-export-skeletons — Package per-neuron SWCs into the versioned release tree.
#'
#' Prefers the detailed skeleton; falls back to L2 if no detailed SWC exists.
#'
#' @section Reads:
#'   - `banc_<ver>_meta.feather`, detailed + L2 SWCs
#'
#' @section Writes:
#'   - `banc_banc_space_swc/<root_id>_{skeleton,l2}.swc` (+ GCS mirror)
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_export_skeletons.sh`.
#'
#' @section Schema:
#'   banc_swc_skeletons.md.

###########################################################
### Export per-neuron SWCs for compiled_data/banc_<ver>/
###
### For each valid banc_<ver> root_id (from banc_<ver>_meta.feather):
###   - If a detailed skeleton exists in banc.swc.save.path,
###     copy as <root_id>_skeleton.swc
###   - Else if an L2 skeleton exists in banc.l2swc.save.path,
###     copy as <root_id>_l2.swc
###   - Else skip
###
### Output dir: <banc_<ver>>/banc_banc_space_swc/
### GCS sync:   gs://lee-lab_brain-and-nerve-cord-fly-connectome/
###             compiled_data/banc_<ver>/banc_banc_space_swc/
###
### Runs standalone after banc-data.R has refreshed
### banc_<ver>_meta.feather.
###
### Usage: Rscript banc/share/banc-export-skeletons.R
###########################################################
source("banc/banc-startup.R")

local({

message("### banc: exporting per-neuron SWCs ###")
t_start <- Sys.time()

save.path <- file.path("/n/data1/hms/neurobio/wilson/connectomes/banc",
                       paste0("banc_", banc.version))
out.dir <- file.path(save.path, "banc_banc_space_swc")
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

meta.file <- file.path(save.path, paste0("banc_", banc.version, "_meta.feather"))
if (!file.exists(meta.file)) {
  stop(sprintf("Meta feather not found: %s — run banc-data.R first", meta.file))
}

meta <- arrow::read_feather(meta.file)
version_id_col <- paste0("banc_", banc.version, "_id")
if (!version_id_col %in% names(meta)) {
  stop(sprintf("Meta missing column %s", version_id_col))
}

valid.ids <- as.character(meta[[version_id_col]])
valid.ids <- valid.ids[!is.na(valid.ids) & nzchar(valid.ids) & valid.ids != "0"]
valid.ids <- unique(valid.ids)
message(sprintf("Valid banc_%s ids: %d", banc.version, length(valid.ids)))

skel.files <- list.files(banc.swc.save.path,   pattern = "\\.swc$", full.names = FALSE)
l2.files   <- list.files(banc.l2swc.save.path, pattern = "\\.swc$", full.names = FALSE)
skel.ids   <- gsub("\\.swc$", "", skel.files)
l2.ids     <- gsub("\\.swc$", "", l2.files)
message(sprintf("Available: %d detailed skeletons, %d L2 skeletons",
                length(skel.ids), length(l2.ids)))

have.skel <- intersect(valid.ids, skel.ids)
need.l2   <- setdiff(valid.ids, skel.ids)
have.l2   <- intersect(need.l2, l2.ids)
missing   <- setdiff(need.l2, l2.ids)

message(sprintf("Plan: %d detailed + %d L2 fallback (%d valid ids have no skeleton)",
                length(have.skel), length(have.l2), length(missing)))

# Build source -> destination manifest
src <- c(file.path(banc.swc.save.path,   paste0(have.skel, ".swc")),
         file.path(banc.l2swc.save.path, paste0(have.l2,   ".swc")))
dst <- c(file.path(out.dir, paste0(have.skel, "_skeleton.swc")),
         file.path(out.dir, paste0(have.l2,   "_l2.swc")))

# Skip already-exported (idempotent for resume)
existing.dst <- list.files(out.dir, pattern = "\\.swc$", full.names = TRUE)
todo <- !(dst %in% existing.dst)
message(sprintf("Already in %s: %d files. To copy: %d.",
                basename(out.dir), sum(!todo), sum(todo)))

if (any(todo)) {
  src.todo <- src[todo]
  dst.todo <- dst[todo]
  # Drop sources that disappeared between listing and copy
  exists.src <- file.exists(src.todo)
  if (any(!exists.src)) {
    message(sprintf("  WARNING: %d source files vanished — skipping", sum(!exists.src)))
    src.todo <- src.todo[exists.src]
    dst.todo <- dst.todo[exists.src]
  }
  message(sprintf("Copying %d files...", length(src.todo)))
  ok <- file.copy(src.todo, dst.todo, overwrite = FALSE, copy.date = TRUE)
  message(sprintf("  Copied: %d/%d", sum(ok), length(ok)))
  if (any(!ok)) {
    failed <- src.todo[!ok][seq_len(min(5, sum(!ok)))]
    message("  First few failures:")
    for (f in failed) message("    ", f)
  }
}

# Sanity-check: no stale files (root_ids that are no longer valid)
final.files <- list.files(out.dir, pattern = "\\.swc$", full.names = FALSE)
final.ids   <- sub("_(skeleton|l2)\\.swc$", "", final.files)
stale       <- setdiff(unique(final.ids), valid.ids)
if (length(stale)) {
  stale.files <- final.files[final.ids %in% stale]
  message(sprintf("  Removing %d stale SWCs (root_ids not in current meta)",
                  length(stale.files)))
  file.remove(file.path(out.dir, stale.files))
}

# Final tally
final.files <- list.files(out.dir, pattern = "\\.swc$", full.names = FALSE)
n.skel <- sum(grepl("_skeleton\\.swc$", final.files))
n.l2   <- sum(grepl("_l2\\.swc$", final.files))
message(sprintf("Final: %d SWCs in %s (%d detailed + %d L2)",
                length(final.files), basename(out.dir), n.skel, n.l2))

# GCS sync
message("\n=== GCS sync ===")
gcs.dst <- sprintf("gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/banc_%s/banc_banc_space_swc/",
                   banc.version)
sync.cmd <- sprintf("gsutil -m rsync -r %s %s", out.dir, gcs.dst)
message("  ", sync.cmd)
rc <- system(sync.cmd)
if (rc != 0) message(sprintf("  WARNING: gsutil rsync exited with code %d", rc))

message(sprintf("\n### Done [%s] ###",
                format(round(difftime(Sys.time(), t_start, units = "mins"), 1))))

})
