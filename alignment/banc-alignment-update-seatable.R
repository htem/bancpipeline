#' banc-alignment-update-seatable — Push alignment results to SeaTable.
#'
#' Dispatcher that loads a region preset and runs the preset's seatable-push
#' logic. Each preset's update-seatable script differs in which SeaTable
#' columns it touches and how it reconciles conflicts with existing curation.
#' Presets live at `alignment/presets/<region>/update-seatable.R`.
#'
#' @section CLI:
#'   --region {optic-lobe,whole-brain}  required; selects the preset
#'   ...                                 preset-specific args follow

args <- commandArgs(trailingOnly = TRUE)
region_idx <- which(args == "--region")
if (length(region_idx) == 0 || region_idx == length(args))
  stop("--region {optic-lobe,whole-brain} required")
region <- args[region_idx + 1]

preset_path <- sprintf("alignment/presets/%s/update-seatable.R", region)
if (!file.exists(preset_path))
  stop(sprintf("unknown region '%s' — no preset at %s. Available: %s",
               region, preset_path,
               paste(list.files("alignment/presets/", full.names = FALSE), collapse = ", ")))

commandArgs <- function(trailingOnly = FALSE) {
  base <- base::commandArgs(trailingOnly = trailingOnly)
  drop_idx <- which(base == "--region")
  if (length(drop_idx) == 0) return(base)
  base[-c(drop_idx, drop_idx + 1)]
}

message(sprintf("[banc-alignment-update-seatable] dispatching to preset: %s", region))
source(preset_path, local = FALSE)
