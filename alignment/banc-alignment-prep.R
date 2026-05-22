#' banc-alignment-prep — Prepare data for the cross-dataset alignment algorithm.
#'
#' Dispatcher that loads a region preset and runs the preset's prep logic.
#' Each preset bundles the region filter, the reference-dataset vocabulary
#' rules, the seed-labelling policy, and the per-type capacity formula —
#' those differ enough between regions (e.g. optic lobe vs whole brain) that
#' a single parameterised script would obscure the logic. Presets live at
#' `alignment/presets/<region>/prep.R`.
#'
#' @section CLI:
#'   --region {optic-lobe,whole-brain}  required; selects the preset
#'   [side]                              positional: right | left | both
#'   --source {local,gcs}                where to load BANC + FAFB inputs from
#'   --local-root PATH                   override local data root
#'   --syn-source {v2,v3}                synapse-version suffix for edgelists
#'
#' @section Adding a new region:
#'   1. Create `alignment/presets/<your-region>/prep.R` (copy an existing
#'      preset and adapt filters / output dir / seed rules).
#'   2. Run `Rscript alignment/banc-alignment-prep.R --region <your-region> ...`.
#'
#' @section Paper:
#'   Methods §"Automated typing by morphology and connectivity". The paper
#'   uses the optic-lobe preset. The whole-brain preset extends the same
#'   algorithm to brain-wide cell-type transfer.

args <- commandArgs(trailingOnly = TRUE)
region_idx <- which(args == "--region")
if (length(region_idx) == 0 || region_idx == length(args))
  stop("--region {optic-lobe,whole-brain} required; see `Rscript alignment/banc-alignment-prep.R --help` (or the dispatcher header)")
region <- args[region_idx + 1]

preset_path <- sprintf("alignment/presets/%s/prep.R", region)
if (!file.exists(preset_path))
  stop(sprintf("unknown region '%s' — no preset at %s. Available: %s",
               region, preset_path,
               paste(list.files("alignment/presets/", full.names = FALSE), collapse = ", ")))

# Strip the --region flag + its value before delegating, so preset scripts
# see the same argv they used to see when invoked directly.
commandArgs <- function(trailingOnly = FALSE) {
  base <- base::commandArgs(trailingOnly = trailingOnly)
  drop_idx <- which(base == "--region")
  if (length(drop_idx) == 0) return(base)
  base[-c(drop_idx, drop_idx + 1)]
}

message(sprintf("[banc-alignment-prep] dispatching to preset: %s", region))
source(preset_path, local = FALSE)
