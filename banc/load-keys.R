#' load-keys — Load private Google / Drive / Sheet identifiers into `banc.keys`.
#'
#' Mirrors `BANC-project/R/startup/load-keys.R`. Reads `data/private/keys.csv`
#' (gitignored) into the named list `banc.keys` exposed to the global
#' environment. Each row is `name,value,description`. Downstream scripts
#' reference `banc.keys$<name>`; missing keys cause the consuming script
#' to skip its Drive-bound step rather than fail hard.
#'
#' @section Reads:
#'   - `data/private/keys.csv` (NOT committed)
#'
#' @section Notes:
#'   - To populate locally: create `data/private/keys.csv` with the rows
#'     documented in CLAUDE.md.

###########################################################################
### Load private identifiers (Google Sheet / Doc / Drive IDs, tokens etc.)
###
### Mirrors BANC-project's R/startup/load-keys.R. Reads
### `data/private/keys.csv` (gitignored) into a named list `banc.keys`
### exposed to the global environment. Each row is `name,value,description`.
###
### Downstream scripts (e.g. banc/annotations/banc-tracing-*.R) reference
### `banc.keys$<name>`. If a key is missing, the script should skip the
### Drive-bound step rather than fail hard.
###
### To populate locally: create `data/private/keys.csv` with the rows
### documented in CLAUDE.md / the project README. The file is never
### committed.
###########################################################################

banc.keys <- list()

.keys_path <- file.path("data", "private", "keys.csv")
if (file.exists(.keys_path)) {
  .keys_df <- tryCatch(
    readr::read_csv(.keys_path, show_col_types = FALSE,
                    col_types = readr::cols(.default = readr::col_character())),
    error = function(e) {
      warning("Could not parse ", .keys_path, ": ", conditionMessage(e))
      NULL
    }
  )
  if (!is.null(.keys_df) && all(c("name", "value") %in% colnames(.keys_df))) {
    banc.keys <- setNames(as.list(.keys_df$value), .keys_df$name)
    message(sprintf("Loaded %d private keys from %s",
                    length(banc.keys), .keys_path))
  } else {
    warning(.keys_path, " missing required columns `name` and `value`; ",
            "banc.keys left empty.")
  }
  rm(.keys_df)
} else {
  message(sprintf(
    "No %s found — banc.keys is empty. Scripts that need Drive/Sheet IDs (e.g. banc/annotations/banc-tracing-*.R) will skip those steps. To enable them, create the file with rows: name,value,description (see CLAUDE.md for the canonical entries).",
    .keys_path))
}
rm(.keys_path)
