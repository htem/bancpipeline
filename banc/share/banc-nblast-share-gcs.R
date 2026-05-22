#' banc-nblast-share-gcs — Publish per-dataset NBLAST feathers + reviewed CSVs to GCS.
#'
#' @section Reads:
#'   - `banc_<target>_<ver>_nblast.feather` (FAFB, MANC, FANC, hemibrain, maleCNS, lr, native)
#'   - reviewed match CSVs
#'
#' @section Writes:
#'   - GCS `nblast/`
#'
#' @section Invoked by:
#'   production v888 rebuild chain.
#'
#' @section Used by:
#'   BANC-project/R/figures/panels_proofread_matching.R (loads
#'   `banc_{fafb_783,manc,hemibrain,fanc,malecns}_nblast.feather` via `read_feather_gcs`).
#'
#' @section Paper:
#'   Methods §"NBLAST cross-dataset matching".

###########################################################
### Push NBLAST feather + reviewed match CSV files to
### Google Cloud Storage for public sharing.
###
### Master bucket: gs://lee-lab_brain-and-nerve-cord-fly-connectome/
### NBLAST files live at <bucket>/nblast/.
###########################################################
source("banc/banc-startup.R")

local({

message("### banc: pushing NBLAST data to GCS ###")
t_start <- Sys.time()

gcs_nblast <- "gs://lee-lab_brain-and-nerve-cord-fly-connectome/nblast"

# Files to push
nblast_files <- c(
  "banc_fafb_783_nblast.feather",
  "banc_manc_v1.2.1_nblast.feather",
  "banc_hemibrain_v1.2.1_nblast.feather",
  "banc_fanc_1116_nblast.feather",
  "banc_mirror_nblast.feather",
  "banc_native_nblast.feather",
  "banc_malecns_v0.9_nblast.feather",
  "banc_manc_reviewed_matches.csv",
  "banc_fafb_reviewed_matches.csv",
  "banc_hemibrain_reviewed_matches.csv",
  "banc_mirror_reviewed_matches.csv",
  "banc_fanc_reviewed_matches.csv",
  "banc_malecns_reviewed_matches.csv"
)

message(sprintf("Pushing %d files to %s ...", length(nblast_files), gcs_nblast))
for (f in nblast_files) {
  src <- file.path(banc.meta.save.path, f)
  if (file.exists(src)) {
    ret <- system(sprintf("gsutil cp %s %s/", src, gcs_nblast))
    if (ret != 0) message(sprintf("  WARNING: gsutil failed for %s (exit %d)", f, ret))
  } else {
    message(sprintf("  Skipping %s (not found)", f))
  }
}

message(sprintf("### banc: GCS push complete [%s] ###",
                format(round(difftime(Sys.time(), t_start, units = "mins"), 1))))

})
