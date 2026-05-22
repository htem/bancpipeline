#' alignment-data-sources — Shared CLI / path resolution for alignment prep.
#'
#' Sourced by `alignment/presets/<region>/prep.R` to decide where the BANC
#' edgelist, target (FAFB / MANC / hemibrain / maleCNS) edgelist and
#' BANC-target NBLAST feather come from: the current BANC version's data
#' tree (on O2 or a local root) or the public GCS bucket. Exposes
#' `parse_alignment_data_args`, `resolve_alignment_paths`, `gcs_cache`,
#' and `load_target_meta_for_alignment`.
#'
#' @section Reads:
#'   - env vars `BANC_ALIGNMENT_SOURCE`, `BANC_ALIGNMENT_LOCAL_ROOT`,
#'     `BANC_SYN_SOURCE` (gated CLI overrides)
#'   - `banc.versioned.save.path` / `banc.connectivity.save.path` /
#'     `banc.meta.save.path` (from `banc/banc-startup.R`)
#'
#' @section Notes:
#'   - Target-agnostic via `load_target_meta_for_alignment()`: dataset-
#'     specific ID columns are renamed to a neutral `target_*` schema so the
#'     preset scripts don't need branching on dataset.

###########################################################
### Alignment data-source resolution
###
### Sourced by banc-alignment-prep.R --region optic-lobe and banc-alignment-prep.R --region whole-brain
### to decide where the BANC edgelist, FAFB edgelist and BANC-FAFB
### NBLAST feather come from: the current BANC version's data tree
### (on O2 or a local root), or the public GCS bucket.
###
### Expects `banc.versioned.save.path`, `banc.connectivity.save.path`
### and `banc.meta.save.path` to be defined (set by banc-startup.R).
###
### Public API:
###   parse_alignment_data_args(args)
###     Strips --source / --local-root / --syn-source from a CLI arg
###     vector and returns list(positional, source, local_root, syn_source).
###     Honours BANC_ALIGNMENT_SOURCE / BANC_ALIGNMENT_LOCAL_ROOT /
###     BANC_SYN_SOURCE env vars. Errors if --source not specified.
###     syn_source defaults to banc.synapse.source.default ('v3').
###
###   resolve_alignment_paths(source, local_root, banc_version,
###                           nblast_version = "783",
###                           syn_source = "v3")
###     Returns list(banc_edgelist, fafb_edgelist, nblast). BANC edgelist
###     carries the syn_source suffix (v2 or v3). Both edgelists are
###     hard requirements; nblast is soft — NA if a local file is
###     missing or a GCS fetch fails.
###
###   gcs_cache(gcs_path)
###     Downloads `gcs_path` to a per-session tempdir cache and
###     returns the local path. Errors if `gsutil cp` fails.
###
###   alignment_path(stage, query, target, ..., date = NULL, ext = "csv",
###                  dir = ".", scheme = "auto")
###     D.6.2 self-describing filename builder. Returns a path of the
###     form
###       {stage}__{query}_{target}[__tags...]__{YYYYMMDD}.{ext}
###     where optional tags (region, side, vq, vt, syn, extra) are
###     appended in canonical order ONLY when provided. See
###     alignment/RENAME_NOTES.md §D.6.2 for the full grammar.
###
###     scheme = "auto" (default) emits the new "v2" name. scheme =
###     "legacy" calls the legacy-name builder; scripts can opt into
###     dual-write by checking Sys.getenv("BANC_ALIGNMENT_NAMING").
###########################################################

#' Load the target dataset's metadata table with normalised column names.
#'
#' Returns the target connectome's neuron metadata tibble with the dataset-
#' specific identifier columns renamed to a neutral `target_*` schema. This
#' keeps the alignment algorithm code target-agnostic — preset scripts only
#' touch `target.meta$target_id`, `$target_cell_type`, etc. rather than
#' `$fafb_id`, so adding MANC / hemibrain / maleCNS support means slotting
#' a new branch into this function, not touching the algorithm.
#'
#' Currently only `target = "fafb"` is wired up (the paper run). MANC,
#' hemibrain, and maleCNS branches are stubbed for future work.
#'
#' @param target  string; one of "fafb", "manc", "hemibrain", "malecns".
#'                Defaults to "fafb" — the paper-run target.
#' @return        tibble with columns: target_id, target_region, target_side,
#'                target_super_class, target_cell_class, target_cell_sub_class,
#'                target_cell_type, target_top_nt. `target_id` is character.
load_target_meta <- function(target = "fafb") {
  target <- tolower(target)
  if (target == "fafb") {
    franken_meta(
      "SELECT fafb_id, region, side, super_class, cell_class, cell_sub_class, cell_type, top_nt FROM franken_meta",
      base = "cns_meta") %>%
      dplyr::filter(!is.na(fafb_id)) %>%
      dplyr::rename(
        target_id            = fafb_id,
        target_region        = region,
        target_side          = side,
        target_super_class   = super_class,
        target_cell_class    = cell_class,
        target_cell_sub_class = cell_sub_class,
        target_cell_type     = cell_type,
        target_top_nt        = top_nt
      ) %>%
      dplyr::mutate(target_id = as.character(target_id))
  } else if (target %in% c("manc", "hemibrain", "malecns")) {
    stop(sprintf(
      "target='%s' not yet wired up — add a branch to load_target_meta() ",
      target),
      "in alignment/alignment-data-sources.R. The schema returned must have ",
      "columns target_id, target_region, target_side, target_super_class, ",
      "target_cell_class, target_cell_sub_class, target_cell_type, target_top_nt.")
  } else {
    stop("Unknown target dataset: '", target,
         "'. Supported: fafb (others: see load_target_meta source).")
  }
}

#' Build a self-describing alignment-output filename (D.6.2).
#'
#' Implements the schema documented in alignment/RENAME_NOTES.md §D.6.2:
#'
#'   {stage}__{query}_{target}[__region-…][__side-…][__vq-…][__vt-…]
#'         [__syn-…][__<extra>]__{YYYYMMDD}.{ext}
#'
#' Required: stage, query, target. The trailing date defaults to today
#' (UTC, YYYYMMDD). All other tags are optional and only emitted when
#' non-NULL / non-empty.
#'
#' Naming-scheme switch:
#'   BANC_ALIGNMENT_NAMING=legacy  →  call legacy_alignment_path() instead
#'   BANC_ALIGNMENT_NAMING=v2 (or unset, default) → new scheme
#'
#' @param stage     character(1); pipeline stage tag, e.g. "prep-banc-meta",
#'                  "align-cosine-tier", "validate-holdout-accuracy".
#'                  Dash-internal, no double-underscores.
#' @param query     query dataset name (lowercase-with-dashes), e.g. "banc".
#' @param target    target dataset name, e.g. "fafb", "manc", "hemibrain".
#' @param region    optional region tag (e.g. "whole-brain", "optic-lobe").
#'                  Emitted as `region-<value>`.
#' @param side      optional `right`/`left`/`both`. Emitted as `side-<value>`.
#' @param vq        optional query-dataset version pin (e.g. "888").
#'                  Emitted as `vq-<value>`.
#' @param vt        optional target-dataset version pin (e.g. "783").
#'                  Emitted as `vt-<value>`.
#' @param syn       optional synapse-table version (e.g. "v2", "v3").
#'                  Emitted as `syn-<value>`.
#' @param extra     optional character vector of preset-defined extra tags,
#'                  each lowercase-with-dashes, emitted in order.
#' @param date      ISO date or YYYYMMDD. Defaults to today UTC.
#' @param ext       file extension, default "csv". No leading dot.
#' @param dir       output directory; default "." (current).
#' @param scheme    "auto" | "v2" | "legacy". Defaults to honour
#'                  BANC_ALIGNMENT_NAMING env var (legacy / v2 / unset → v2).
#' @return          character(1) path.
alignment_path <- function(stage, query, target,
                           region = NULL, side = NULL,
                           vq = NULL, vt = NULL, syn = NULL,
                           extra = NULL,
                           date = NULL, ext = "csv",
                           dir = ".", scheme = "auto") {
  if (identical(scheme, "auto")) {
    scheme <- tolower(Sys.getenv("BANC_ALIGNMENT_NAMING", unset = "v2"))
  }
  if (identical(scheme, "legacy")) {
    return(legacy_alignment_path(stage = stage, query = query, target = target,
                                  region = region, side = side,
                                  vq = vq, vt = vt, syn = syn, extra = extra,
                                  date = date, ext = ext, dir = dir))
  }
  if (!identical(scheme, "v2")) {
    stop("alignment_path(): scheme must be 'v2' or 'legacy'; got '", scheme, "'.")
  }
  # Stage must be a single dash-internal token (no double-underscores).
  if (grepl("__", stage)) {
    stop("alignment_path(): stage may not contain '__' (field separator).")
  }
  # Compose tag fields in canonical order.
  parts <- c(stage, paste0(query, "_", target))
  tag <- function(prefix, value) {
    if (is.null(value)) return(NULL)
    if (is.na(value) || !nzchar(as.character(value))) return(NULL)
    paste0(prefix, "-", as.character(value))
  }
  for (p in list(c("region", region), c("side", side),
                 c("vq", vq), c("vt", vt), c("syn", syn))) {
    t <- tag(p[1], p[2])
    if (!is.null(t)) parts <- c(parts, t)
  }
  if (!is.null(extra) && length(extra) > 0) {
    extra <- extra[!is.na(extra) & nzchar(as.character(extra))]
    if (any(grepl("__", extra))) {
      stop("alignment_path(): extra tag may not contain '__'.")
    }
    parts <- c(parts, extra)
  }
  # Date: default today UTC, accept Date / POSIXct / YYYY-MM-DD / YYYYMMDD.
  if (is.null(date)) {
    date_str <- format(Sys.time(), "%Y%m%d", tz = "UTC")
  } else if (inherits(date, c("Date", "POSIXt"))) {
    date_str <- format(date, "%Y%m%d", tz = "UTC")
  } else {
    s <- gsub("-", "", as.character(date))
    if (!grepl("^[0-9]{8}$", s)) {
      stop("alignment_path(): date must be YYYY-MM-DD or YYYYMMDD; got '",
           date, "'.")
    }
    date_str <- s
  }
  parts <- c(parts, date_str)
  fname <- paste0(paste(parts, collapse = "__"), ".", ext)
  if (identical(dir, ".") || is.null(dir)) fname else file.path(dir, fname)
}

#' Legacy alignment-output path (pre-D.6.2 naming).
#'
#' Maps the new (stage, query, target, region, side, ...) parameters back
#' to the legacy filename convention so producers can dual-write during the
#' migration. Only the stage→legacy-name mappings used by the paper-run set
#' are implemented; unrecognised stages raise an error to flag the gap.
#'
#' @inheritParams alignment_path
#' @return character(1) path matching the pre-D.6.2 naming.
legacy_alignment_path <- function(stage, query, target,
                                  region = NULL, side = NULL,
                                  vq = NULL, vt = NULL, syn = NULL,
                                  extra = NULL,
                                  date = NULL, ext = "csv", dir = ".") {
  # Legacy infix for region. Optic-lobe used "optic", whole-brain used "brain".
  region_legacy_infix <- function(r) {
    if (is.null(r)) return(NULL)
    switch(tolower(r),
      "optic-lobe" = "optic",
      "whole-brain" = "brain",
      stop("legacy_alignment_path(): no legacy infix for region='", r, "'."))
  }
  infix <- region_legacy_infix(region)
  side_l <- if (is.null(side)) "" else paste0(side)
  # Build the {prefix}_{side}_{stage_legacy} bones per legacy convention.
  q_t_infix_side <- function() {
    paste(c(query, infix, side_l), collapse = "_")
  }
  q_target_infix_side <- function() {
    paste(c(query, target, infix, side_l), collapse = "_")
  }
  # Map known stages to their legacy filename.
  fname <- switch(stage,
    "prep-banc-meta"               = sprintf("%s_meta.csv", q_t_infix_side()),
    "prep-target-meta"             = sprintf("%s_%s_meta.csv", target, paste(c(infix, side_l), collapse = "_")),
    "prep-banc-edges"              = sprintf("%s_edgelist.feather", q_t_infix_side()),
    "prep-target-edges"            = sprintf("%s_%s_edgelist.feather", target, paste(c(infix, side_l), collapse = "_")),
    "prep-seeds"                   = sprintf("%s_seeds.csv", q_t_infix_side()),
    "prep-seeds-production"        = sprintf("%s_seeds_production.csv", q_t_infix_side()),
    "prep-nblast"                  = sprintf("%s_nblast.csv", q_target_infix_side()),
    "prep-capacity"                = sprintf("%s_type_capacity.csv", target),
    "prep-forbidden"               = "forbidden-matches.csv",
    "validate-holdout-accuracy"    = sprintf("%s_holdout_accuracy.csv", q_t_infix_side()),
    "validate-holdout-confusions"  = sprintf("%s_holdout_confusions.csv", q_t_infix_side()),
    "validate-nt-mismatches"       = sprintf("%s_nt_mismatches.csv", q_t_infix_side()),
    "validate-mismatches"          = sprintf("%s_mismatches.csv", q_t_infix_side()),
    "validate-discrepancies"       = {
      extra_tag <- if (length(extra) > 0) paste0("_", paste(extra, collapse = "_")) else ""
      sprintf("%s_discrepancies%s.csv", q_t_infix_side(), extra_tag)
    },
    {
      # Bare "align" stage: <q>_<infix>_<side>_alignment.<ext> (no variant rest).
      if (identical(stage, "align")) {
        sprintf("%s_alignment.%s", q_t_infix_side(), ext)
      } else if (grepl("^align-", stage)) {
        # align-* legacy names: <q>_<infix>_<side>_alignment_<rest>.csv
        rest <- sub("^align-", "", stage)
        sprintf("%s_alignment_%s.%s", q_t_infix_side(), gsub("-", "_", rest), ext)
      } else if (grepl("^ntac", stage)) {
        rest <- sub("^ntac", "ntac", stage)
        sprintf("%s_%s.%s", q_t_infix_side(), gsub("-", "_", rest), ext)
      } else {
        stop("legacy_alignment_path(): no legacy mapping for stage='", stage, "'.")
      }
    }
  )
  if (identical(dir, ".") || is.null(dir)) fname else file.path(dir, fname)
}

gcs_cache <- function(gcs_path,
                      cache_dir = file.path(tempdir(), "gcs_cache")) {
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  local_file <- file.path(cache_dir, basename(gcs_path))
  if (!file.exists(local_file)) {
    message("  Downloading ", basename(gcs_path), " from GCS...")
    system2("gsutil", c("cp", gcs_path, local_file),
            stdout = FALSE, stderr = FALSE)
    if (!file.exists(local_file)) stop("gsutil cp failed for ", gcs_path)
  } else {
    message("  Using cached ", basename(local_file))
  }
  local_file
}

parse_alignment_data_args <- function(args) {
  source_val <- NULL
  local_root <- NULL
  syn_source <- NULL
  positional <- character(0)
  i <- 1L
  while (i <= length(args)) {
    a <- args[i]
    if (identical(a, "--source")) {
      if (i == length(args)) stop("--source requires a value (local|gcs)")
      source_val <- args[i + 1L]; i <- i + 2L
    } else if (startsWith(a, "--source=")) {
      source_val <- sub("^--source=", "", a); i <- i + 1L
    } else if (identical(a, "--local-root")) {
      if (i == length(args)) stop("--local-root requires a path")
      local_root <- args[i + 1L]; i <- i + 2L
    } else if (startsWith(a, "--local-root=")) {
      local_root <- sub("^--local-root=", "", a); i <- i + 1L
    } else if (identical(a, "--syn-source")) {
      if (i == length(args)) stop("--syn-source requires a value (v2|v3)")
      syn_source <- args[i + 1L]; i <- i + 2L
    } else if (startsWith(a, "--syn-source=")) {
      syn_source <- sub("^--syn-source=", "", a); i <- i + 1L
    } else {
      positional <- c(positional, a); i <- i + 1L
    }
  }
  if (is.null(source_val)) {
    env_src <- Sys.getenv("BANC_ALIGNMENT_SOURCE", unset = "")
    if (nzchar(env_src)) source_val <- env_src
  }
  if (is.null(local_root)) {
    env_root <- Sys.getenv("BANC_ALIGNMENT_LOCAL_ROOT", unset = "")
    if (nzchar(env_root)) local_root <- env_root
  }
  if (is.null(syn_source)) {
    env_syn <- Sys.getenv("BANC_SYN_SOURCE", unset = "")
    if (nzchar(env_syn)) syn_source <- env_syn
    else if (exists("banc.synapse.source.default")) syn_source <- banc.synapse.source.default
    else syn_source <- "v3"
  }
  if (is.null(source_val) || !nzchar(source_val)) {
    stop("Data source must be specified. Pass --source {local|gcs} or set ",
         "BANC_ALIGNMENT_SOURCE={local|gcs}.")
  }
  source_val <- tolower(source_val)
  if (!source_val %in% c("local", "gcs")) {
    stop("--source must be 'local' or 'gcs'; got '", source_val, "'")
  }
  syn_source <- tolower(syn_source)
  if (!syn_source %in% c("v2", "v3")) {
    stop("--syn-source must be 'v2' or 'v3'; got '", syn_source, "'")
  }
  list(positional = positional, source = source_val, local_root = local_root,
       syn_source = syn_source)
}

resolve_alignment_paths <- function(source, local_root = NULL,
                                    banc_version,
                                    nblast_version = "783",
                                    syn_source = NULL) {
  if (is.null(syn_source)) {
    env_syn <- Sys.getenv("BANC_SYN_SOURCE", unset = "")
    syn_source <- if (nzchar(env_syn)) env_syn
                  else if (exists("banc.synapse.source.default"))
                    banc.synapse.source.default
                  else "v3"
  }
  syn_source <- tolower(syn_source)
  if (!syn_source %in% c("v2", "v3")) {
    stop("syn_source must be 'v2' or 'v3'; got '", syn_source, "'")
  }

  gcs_base        <- "gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data"
  gcs_nblast_base <- "gs://lee-lab_brain-and-nerve-cord-fly-connectome/nblast"

  banc_el_name    <- sprintf("banc_%s_edgelist_simple_%s.feather",
                             banc_version, syn_source)
  fafb_el_name    <- sprintf("fafb_%s_simple_edgelist.feather", nblast_version)
  nblast_name     <- sprintf("banc_fafb_%s_nblast.feather", nblast_version)

  # Local FAFB edgelist lives under wilson/connectomes/fafb on O2, NOT under
  # banc.connectivity.save.path. Hard-coded until we add a proper variable
  # to banc-startup.R.
  fafb_local_dir <- "/n/data1/hms/neurobio/wilson/connectomes/fafb"

  if (source == "local") {
    if (is.null(local_root) || !nzchar(local_root)) {
      banc_el  <- file.path(banc.versioned.save.path, banc_el_name)
      fafb_el  <- file.path(fafb_local_dir, fafb_el_name)
      nblast_f <- file.path(banc.meta.save.path, nblast_name)
      src_desc <- "default O2 paths"
    } else {
      # Convention under --local-root:
      #   <root>/banc_<ver>/banc_<ver>_edgelist_simple.feather
      #   <root>/connectivity/fafb_<ver>_simple_edgelist.feather
      #   <root>/meta/banc_fafb_<ver>_nblast.feather
      banc_el  <- file.path(local_root, sprintf("banc_%s", banc_version), banc_el_name)
      fafb_el  <- file.path(local_root, "connectivity", fafb_el_name)
      nblast_f <- file.path(local_root, "meta", nblast_name)
      src_desc <- sprintf("--local-root %s", local_root)
    }
    message(sprintf("  Data source: local (%s)", src_desc))

    missing_el <- !file.exists(c(banc_el, fafb_el))
    if (any(missing_el)) {
      stop("Local edgelist files missing:\n  ",
           paste(c(banc_el, fafb_el)[missing_el], collapse = "\n  "),
           "\nUse --source gcs to fetch from GCS, or fix --local-root.")
    }
    if (!file.exists(nblast_f)) {
      message("  NBLAST file not found locally: ", nblast_f)
      message("  Proceeding without NBLAST (Tier 2 seeds will be empty).")
      nblast_f <- NA_character_
    }
    return(list(banc_edgelist = banc_el,
                fafb_edgelist = fafb_el,
                nblast        = nblast_f))
  }

  # --source gcs
  message("  Data source: gcs")
  # Layout under compiled_data is `<dataset>_<version>/<file>`.
  banc_el <- gcs_cache(file.path(gcs_base,
    sprintf("banc_%s", banc_version), banc_el_name))
  fafb_el <- gcs_cache(file.path(gcs_base,
    sprintf("fafb_%s", nblast_version), fafb_el_name))
  nblast_f <- tryCatch(
    gcs_cache(file.path(gcs_nblast_base, nblast_name)),
    error = function(e) {
      message("  GCS NBLAST fetch failed: ", conditionMessage(e))
      message("  Proceeding without NBLAST (Tier 2 seeds will be empty).")
      NA_character_
    }
  )
  list(banc_edgelist = banc_el,
       fafb_edgelist = fafb_el,
       nblast        = nblast_f)
}
