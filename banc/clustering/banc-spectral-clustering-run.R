#' banc-spectral-clustering-run — Dispatcher: invoke `banc-spectral-clustering.py` from R.
#'
#' Helper that auto-locates the BANC data directory (O2 versioned tree or
#' local GCS mirror) and shells out to the Python spectral-clustering
#' script with the right `--data-dir` / `--banc-version` flags.
#'
#' @section Reads:
#'   - same inputs as `banc-spectral-clustering.py` (data dir auto-detected)
#'
#' @section Writes:
#'   - same outputs as `banc-spectral-clustering.py`
#'
#' @section CLI:
#'   --gcs   use the GCS-mirror data dir rather than the O2 versioned tree

###########################################################
### Run banc-spectral-clustering.py from R
###
### Helper to invoke the Python spectral clustering script
### from an R session, using the same data-dir
### auto-detection as the rest of the pipeline.
###
### Usage:
###   Rscript banc/clustering/banc-spectral-clustering-run.R
###   Rscript banc/clustering/banc-spectral-clustering-run.R --gcs
###########################################################
source("banc/banc-startup.R")

local({
  message("### banc: running spectral clustering (Python) ###")

  cli_args <- commandArgs(trailingOnly = TRUE)
  use_gcs <- "--gcs" %in% cli_args

  # banc.version set in banc-startup.R

  # Locate data directory (same candidates as Python script)
  data_candidates <- c(
    file.path("/n/data1/hms/neurobio/wilson/connectomes/banc",
              paste0("banc_", banc.version)),
    file.path("lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data",
              paste0("banc_", banc.version))
  )
  data_dir <- NULL
  for (d in data_candidates) {
    if (dir.exists(d)) { data_dir <- d; break }
  }
  if (is.null(data_dir)) stop("Could not find BANC ", banc.version, " data directory")

  # Build command
  script <- "banc/clustering/banc-spectral-clustering.py"
  stopifnot(file.exists(script))

  cmd <- sprintf("python3 %s --data-dir %s --banc-version %s",
                 shQuote(script), shQuote(data_dir), banc.version)
  if (use_gcs) cmd <- paste(cmd, "--gcs")

  message(sprintf("Running: %s", cmd))
  ret <- system(cmd)

  if (ret != 0) {
    stop("banc-spectral-clustering.py exited with code ", ret)
  }

  message("### banc: spectral clustering (Python) complete ###")
})
