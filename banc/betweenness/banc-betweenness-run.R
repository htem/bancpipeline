#' banc-betweenness-run — Dispatcher: run `banc-betweenness.py` from R with data-dir auto-detection.
#'
#' Helper that invokes the Python betweenness centrality script using the
#' same data-dir auto-detection as the rest of the pipeline. Parses
#' `--source v2|v3` (or `BANC_SYN_SOURCE`) and a positional task name
#' (`all_to_all` | `afferent_efferent` | `both`).
#'
#' @section Reads:
#'   - `<banc.versioned.save.path>/banc_<ver>_edgelist_simple_<src>.feather`
#'     (resolved by the Python child)
#'   - env var `BANC_SYN_SOURCE`
#'
#' @section Writes:
#'   - `<banc.versioned.save.path>/banc_<ver>_betweenness_<mode>_<src>.csv`
#'     via the Python child
#'
#' @section CLI:
#'   [task]              all_to_all | afferent_efferent | both (default both)
#'   --source {v2,v3}    default banc.synapse.source.default
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_betweenness.sh`,
#'   `o2/production/o2_banc_v888_rebuild.sh`,
#'   `o2/oneshots/o2_banc_v890_rebuild.sh`
#'
#' @section Schema:
#'   BANC-project/manuscript/print/dataverse/documentation/banc_888_betweenness.md
#'
#' @section Paper:
#'   Methods §"Betweenness centrality".

###########################################################
### Run banc-betweenness.py from R
###
### Helper to invoke the Python betweenness centrality
### script from an R session, using the same data-dir
### auto-detection as the rest of the pipeline.
###
### Usage:
###   Rscript banc/betweenness/banc-betweenness-run.R
###   Rscript banc/betweenness/banc-betweenness-run.R afferent_efferent
###   Rscript banc/betweenness/banc-betweenness-run.R all_to_all
###   Rscript banc/betweenness/banc-betweenness-run.R both --source v2
###########################################################
source("banc/banc-startup.R")

local({
  message("### banc: running betweenness centrality (Python) ###")

  cli_args <- commandArgs(trailingOnly = TRUE)

  # Parse --source flag (strip it out before positional-arg lookup)
  i <- which(cli_args == "--source")
  if (length(i) == 1 && length(cli_args) >= i + 1) {
    syn_source <- tolower(cli_args[i + 1])
    cli_args <- cli_args[-c(i, i + 1)]
  } else {
    env <- Sys.getenv("BANC_SYN_SOURCE", unset = NA_character_)
    syn_source <- if (!is.na(env) && nzchar(env)) tolower(env)
                  else if (exists("banc.synapse.source.default"))
                        tolower(banc.synapse.source.default)
                  else "v3"
  }
  stopifnot(syn_source %in% c("v2", "v3"))

  # Parse optional mode argument (afferent_efferent | all_to_all | both)
  mode <- if (length(cli_args) >= 1) cli_args[1] else "both"
  stopifnot(mode %in% c("afferent_efferent", "all_to_all", "both"))

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
  script <- "banc/betweenness/banc-betweenness.py"
  stopifnot(file.exists(script))

  cmd <- sprintf("python3 %s --data-dir %s --banc-version %s --mode %s --source %s",
                 shQuote(script), shQuote(data_dir), banc.version, mode, syn_source)

  message(sprintf("Running: %s", cmd))
  ret <- system(cmd)

  if (ret != 0) {
    stop("banc-betweenness.py exited with code ", ret)
  }

  message("### banc: betweenness centrality complete ###")
})
