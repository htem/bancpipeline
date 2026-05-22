#' banc-make-proofread-ids — Emit a txt file of proofread BANC root_ids.
#'
#' Used by the proofread-redo workflow: each downstream NBLAST script
#' honours `BANC_TEST_IDS_FILE` and processes only the listed IDs.
#'
#' @section Reads:
#'   - `banc_meta.csv` (override path via `BANC_META_PATH`)
#'
#' @section Writes:
#'   - `banc_proofread_ids.txt` (override path via `BANC_PROOFREAD_IDS_FILE`)
#'
#' @section Invoked by:
#'   `o2/oneshots/o2_banc_proofread_redo.sh`.

###########################################################
### Emit a text file of root_ids whose proofread == TRUE
### in banc_meta.csv. Used by the proofread-redo workflow
### (set BANC_TEST_IDS_FILE to the output path before each
### NBLAST script, with BANC_NBLAST_REDO=TRUE).
###########################################################
suppressPackageStartupMessages({
  library(readr); library(dplyr)
})

meta_path  <- Sys.getenv("BANC_META_PATH",
                         unset = "/n/data1/hms/neurobio/wilson/banc/meta/banc_meta.csv")
out_path   <- Sys.getenv("BANC_PROOFREAD_IDS_FILE",
                         unset = "/n/data1/hms/neurobio/wilson/banc/meta/banc_proofread_ids.txt")

stopifnot(file.exists(meta_path))

m <- readr::read_csv(meta_path, col_types = cols(.default = "c"),
                     show_col_types = FALSE)
ids <- m %>%
  dplyr::filter(proofread == "TRUE", !is.na(root_id), root_id != "") %>%
  dplyr::pull(root_id) %>% unique()

writeLines(ids, out_path)
message(sprintf("Wrote %d proofread root_ids to %s",
                length(ids), out_path))
