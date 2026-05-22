#' banc-calculate-completion — Synapse-capture rates across v1, v2, v3 detection rounds.
#'
#' Classifies each synapse's pre/post as proofread / identified / fragment
#' and emits per-source capture-rate CSVs + a 3-way comparison.
#'
#' @section Reads:
#'   - GCS v1/v2/v3 synapse parquets, `banc_<ver>_synapses_v3.parquet`
#'   - SeaTable `banc_meta`, `banc_backbone_proofread()` (CAVE)
#'
#' @section Writes (under `capture_rates/`):
#'   - `banc_<ver>_<source>_<gross|inout|region|neuropil>_capture_rates.csv`
#'   - `banc_<ver>_v1_v2_v3_{summary,region_summary}.csv`
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_completion.sh`, production v888 rebuild chain.
#'
#' @section Notes:
#'   - The v3 sample block hit the SLURM-only arrow C++ URI bug; see
#'     `banc-v3-synapse-sample.R` for the standalone replacement.

###############################################
### Synapse completion metrics (v1 vs v2 vs v3)
###
### For a target BANC version (default banc.version, override with
### BANC_V3_TARGET_VERSION env var), pulls the v888-rooted synapse parquets
### for v1, v2, v3 from GCS, attaches spatial neuropil/region/side, classifies
### each synapse's pre/post root_id as proofread / identified / fragment, and
### writes per-source capture-rate CSVs plus a 3-way comparison summary so we
### can judge which synapse-detection version captures the most identified
### connectivity.
###
### Sources at v<ver> on gs://lee-lab_brain-and-nerve-cord-fly-connectome/v<ver>/
###   v1: synapses_v1_human_readable_id_size_prerootid_postrootid_prex_prey_prez_neuropil.parquet
###       (already carries `neuropil`; region/side parsed from string)
###   v2: synapses_v2_human_readable_id_size_prerootid_postrootid_prex_prey_prez.parquet
###       joined on `id` with banc.connectivity.save.path/banc_synapses_to_neuropils_v2.csv
###   v3: read locally from banc.synapses.v3.save.path/banc_<ver>_synapses_v3.parquet
###       (built by banc/metrics/banc-synapses-v3-optimised.R)
###
### Status classification (matches legacy/banc-assess-synapses.R semantics):
###   proofread  = pre/post_root_id in banc_backbone_proofread() (CAVE)
###   identified = pre/post_root_id in banctable_query()$root_<ver>
###                (or fall back to local banc_<ver>_meta.feather if SeaTable unreachable)
###   fragment   = neither
###
### Outputs (all under banc.synapses.v3.save.path/capture_rates/):
###   banc_<ver>_<source>_<gross|inout|region|neuropil>_capture_rates.csv  (12 files)
###   banc_<ver>_v1_v2_v3_summary.csv          (3-way gross + inout summary)
###   banc_<ver>_v1_v2_v3_region_summary.csv   (3-way per-region+side capture %)
###############################################
source("banc/banc-startup.R")
library(knitr)

local({

  # -------------------------------------------------------
  # Config
  # -------------------------------------------------------
  target.version <- {
    v <- Sys.getenv("BANC_V3_TARGET_VERSION", unset = "")
    if (!nzchar(v)) banc.version else v
  }

  # When TRUE, skip GCS downloads if the cached parquet already exists.
  USE_CACHE <- TRUE

  # Size thresholds: capture rates are computed once per threshold, per source.
  # 0L = no filter (output filenames have no suffix); non-zero values append
  # `_size_thresh_<N>` to every CSV name produced for that pass.
  # Project convention historically is >= 5 (matches banc-calculate-connectivity.R).
  # 10L is added because Zetta v3 source has min size = 10, so this gives the
  # cleanest apples-to-apples slice across all three sources (v1+v2 also
  # filtered down to v3's natural floor).
  SIZE_THRESHOLDS <- c(0L, 5L, 10L)

  # v3 stratified neuroglancer sample (post-loop, NO size threshold applied).
  # 100 syns per (np_primary × side), drawn evenly across N_BINS mean_score bins.
  V3_SAMPLE_PER_GROUP <- 100L
  V3_SAMPLE_BINS      <- 10L

  # Public release layout: published synapse tables live under
  # neuron_connectivity/v<ver>/ (mirrors the v626 preprint layout). Was a
  # stray flat path pre-2026-05-13; now cleaned up.
  GCS_BASE <- sprintf(
    "gs://lee-lab_brain-and-nerve-cord-fly-connectome/neuron_connectivity/v%s",
    target.version
  )

  banc.synapses.v3.save.path <- file.path(banc.save.path, "synapses_v3")

  # rates_dir: BANC-project's data/completion/ so outputs land where the
  # paper-figure scripts expect them and can be committed/pushed straight
  # from BANC-project. Override via BANC_PROJECT_DATA_DIR env var (e.g.
  # to redirect to a checkout on another machine).
  banc.project.data.dir <- Sys.getenv(
    "BANC_PROJECT_DATA_DIR",
    unset = file.path(banc.save.path, "BANC-project", "data"))
  rates_dir <- file.path(banc.project.data.dir, "completion")

  cache_dir <- file.path(banc.synapses.v3.save.path, "cache",
                         sprintf("v%s_inputs", target.version))
  for (d in c(rates_dir, cache_dir)) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  }

  message("### banc: synapse completion metrics (v1 / v2 / v3) ###")
  message(sprintf("  target_version=%s (banc.version=%s)",
                  target.version, banc.version))
  message(sprintf("  rates_dir=%s", rates_dir))
  message(sprintf("  cache_dir=%s", cache_dir))

  # -------------------------------------------------------
  # GCS pull (idempotent)
  # -------------------------------------------------------
  pull_from_gcs <- function(remote, local) {
    if (USE_CACHE && file.exists(local) && file.size(local) > 0L) {
      message(sprintf("  cached: %s", basename(local)))
      return(invisible(local))
    }
    message(sprintf("  gsutil cp %s -> %s", basename(remote), local))
    st <- system2("gsutil", c("-q", "cp", remote, local),
                  stdout = "", stderr = "")
    if (st != 0L) stop(sprintf("gsutil cp failed (%d): %s", st, remote))
    invisible(local)
  }

  v1_remote <- file.path(
    GCS_BASE,
    "synapses_v1_human_readable_id_size_prerootid_postrootid_prex_prey_prez_neuropil.parquet"
  )
  v2_remote <- file.path(
    GCS_BASE,
    "synapses_v2_human_readable_id_size_prerootid_postrootid_prex_prey_prez.parquet"
  )
  # v3 root_ids come from the GCS export of the ingested CAVE synapses_v3
  # table — those are CAVE-curated and trusted. We DO NOT use any locally
  # produced root_ids for v3 (the local ones came from a contaminated
  # svid_cache and were ~30% spuriously "0"; see tasks.md).
  #
  # As of v890 the curated parquet is no longer staged at GCS; we read the
  # CSV.gz instead (same pattern as v2/v3 in banc-calculate-connectivity.R).
  # Pre-v890 versions can still use the parquet path — try parquet first,
  # fall back to CSV.gz.
  v3_remote_parquet <- file.path(
    GCS_BASE,
    "synapses_v3_human_readable_id_size_prerootid_postrootid_prex_prey_prez.parquet"
  )
  v3_remote_csv <- file.path(GCS_BASE, "synapses_v3_human_readable.csv.gz")
  v1_local <- file.path(cache_dir, basename(v1_remote))
  v2_local <- file.path(cache_dir, basename(v2_remote))
  v3_local_parquet <- file.path(cache_dir, basename(v3_remote_parquet))
  v3_local_csv_gz <- file.path(cache_dir, basename(v3_remote_csv))
  v3_local_csv <- sub("\\.gz$", "", v3_local_csv_gz)

  message("\n--- pulling v1, v2 from GCS ---")
  pull_from_gcs(v1_remote, v1_local)
  pull_from_gcs(v2_remote, v2_local)

  # v3 source resolution: parquet if it exists upstream (legacy versions),
  # otherwise CSV.gz. Records which we got into v3_source for load_v3().
  v3_source <- NA_character_
  v3_local <- NA_character_
  if (file.exists(v3_local_parquet) && file.size(v3_local_parquet) > 0L) {
    message("Using cached v3 parquet: ", v3_local_parquet)
    v3_source  <- "parquet"
    v3_local   <- v3_local_parquet
  } else if (file.exists(v3_local_csv) && file.size(v3_local_csv) > 0L) {
    message("Using cached v3 CSV: ", v3_local_csv)
    v3_source  <- "csv"
    v3_local   <- v3_local_csv
  } else {
    # Try parquet upstream first (cheap probe via gsutil stat)
    has_parquet <- system2("gsutil", c("-q", "stat", v3_remote_parquet),
                           stdout = FALSE, stderr = FALSE) == 0L
    if (has_parquet) {
      message("\n--- pulling v3 parquet from GCS ---")
      pull_from_gcs(v3_remote_parquet, v3_local_parquet)
      v3_source <- "parquet"
      v3_local  <- v3_local_parquet
    } else {
      message("\n--- v3 parquet not at GCS; pulling CSV.gz ---")
      pull_from_gcs(v3_remote_csv, v3_local_csv_gz)
      message("Decompressing v3 CSV...")
      system(sprintf("gunzip %s", v3_local_csv_gz))
      stopifnot(file.exists(v3_local_csv))
      v3_source <- "csv"
      v3_local  <- v3_local_csv
    }
  }

  # Spatial neuropil lookup for v3: produced by Stage 2 of
  # banc/metrics/banc-synapses-v3-optimised.R. We use ONLY the spatial cols
  # (syn_id, region, side, neuropil) — its root_ids cols are corrupt and
  # ignored. syn_id is stable across BANC root-id versions (see tasks.md
  # Section 2), so older versioned spatial parquets are valid for newer
  # rebuilds. Falls back: target.version → 888 → 850.
  .pick_spatial <- function() {
    for (v in unique(c(target.version, "888", "850"))) {
      p <- file.path(banc.synapses.v3.save.path,
                     sprintf("banc_%s_synapses_v3.parquet", v))
      if (file.exists(p)) return(p)
    }
    stop("No v3 spatial parquet found for any of target/888/850.")
  }
  v3_spatial_local <- .pick_spatial()
  message(sprintf("v3 spatial parquet: %s", basename(v3_spatial_local)))

  # -------------------------------------------------------
  # Status classification (proofread / identified / fragment)
  # -------------------------------------------------------
  message("\n--- loading proofread + identified IDs ---")

  proof.ids <- tryCatch({
    bb <- bancr::banc_backbone_proofread()
    unique(as.character(bb$pt_root_id))
  }, error = function(e) {
    warning(sprintf("banc_backbone_proofread() failed (%s); proofread split disabled",
                    e$message))
    character(0)
  })
  proof.ids <- proof.ids[!is.na(proof.ids) & nzchar(proof.ids)]

  # SeaTable banc_meta is the source of truth for root_<ver> mapping. Fall back
  # to the local meta feather only if SeaTable is unreachable.
  identified.ids <- tryCatch({
    bm <- banctable_query()
    col <- paste0("root_", target.version)
    if (!col %in% names(bm)) col <- "root_id"
    unique(as.character(bm[[col]]))
  }, error = function(e) {
    warning(sprintf("banctable_query() failed (%s); falling back to local meta feather",
                    e$message))
    f1 <- file.path(banc.connectivity.save.path,
                    sprintf("banc_%s_meta.feather", target.version))
    f2 <- file.path("/n/data1/hms/neurobio/wilson/connectomes/banc",
                    paste0("banc_", target.version),
                    sprintf("banc_%s_meta.feather", target.version))
    f <- if (file.exists(f1)) f1 else if (file.exists(f2)) f2 else NA_character_
    if (is.na(f)) stop("No identified-IDs source available (SeaTable + meta both missing)")
    bm <- arrow::read_feather(f)
    col <- paste0("root_", target.version)
    if (!col %in% names(bm)) col <- "root_id"
    unique(as.character(bm[[col]]))
  })
  identified.ids <- identified.ids[!is.na(identified.ids) & nzchar(identified.ids)]

  message(sprintf("  identified=%s  proofread=%s (overlap: %s)",
                  format(length(identified.ids), big.mark = ","),
                  format(length(proof.ids),       big.mark = ","),
                  format(length(intersect(identified.ids, proof.ids)),
                         big.mark = ",")))

  # Sanity guard: SeaTable's root_<ver> column is sometimes empty when a new
  # BANC version is in the middle of being staged (the rebuild's banc-ids.R
  # populates it). Without identified.ids, every non-proofread cell falls into
  # "fragment" and the capture-rate CSVs become silently misleading. Fail loud.
  if (length(identified.ids) < 10000L) {
    stop(sprintf(paste0("identified.ids has only %d entries -- SeaTable banc_meta ",
                        "may not have root_%s populated yet. Run banc/update/banc-ids.R ",
                        "(or the v%s rebuild orchestrator) first to populate it. ",
                        "Outputs would be misleading without this; refusing to continue."),
                 length(identified.ids), target.version, target.version))
  }

  classify_status <- function(rid) {
    rid <- as.character(rid)
    out <- rep("fragment", length(rid))
    out[rid %in% identified.ids] <- "identified"
    out[rid %in% proof.ids]      <- "proofread"
    out
  }

  # -------------------------------------------------------
  # v1 region/side parsing
  # -------------------------------------------------------
  # v1's `neuropil` column already encodes side as "_R"/"_L" suffix, and
  # region as a prefix (OL_ optic, CB_ central, VNC_, NK_/neck_, MB_, etc.).
  # Map prefixes to the same {central_brain, optic_lobes, suboesophageal_zone,
  # vnc, neck, outside} buckets used by v2/v3 so the 3-way region comparison
  # is meaningful.
  parse_v1_neuropil <- function(np) {
    np <- ifelse(is.na(np) | np == "", "outside", np)
    side <- dplyr::case_when(
      grepl("_L$", np) ~ "left",
      grepl("_R$", np) ~ "right",
      TRUE             ~ "centre"
    )
    # Strip _L/_R suffix for the bare neuropil name
    np_bare <- sub("_[LR]$", "", np)

    # Pull the {OL, CB, VNC, NK, MB, ...} prefix
    prefix <- sub("_.*$", "", np_bare)
    # SEZ neuropils inside CB: GNG, CAN, FLA, AMMC, SAD, PRW
    sez_mask <- grepl("(^|_)(GNG|CAN|FLA|AMMC|SAD|PRW)(_|$)", np_bare)
    region <- dplyr::case_when(
      np_bare == "outside"            ~ "outside",
      grepl("^outside", np_bare)      ~ "outside",
      sez_mask                        ~ "suboesophageal_zone",
      prefix == "OL"                  ~ "optic_lobes",
      prefix == "VNC"                 ~ "vnc",
      prefix %in% c("NK", "Neck",
                    "neck")           ~ "neck",
      prefix == "CB"                  ~ "central_brain",
      prefix == "MB"                  ~ "central_brain",   # mushroom body lives in CB
      grepl("_unassigned",  np_bare)  ~ "outside",
      TRUE                            ~ "central_brain"
    )
    list(neuropil = np_bare, region = region, side = side)
  }

  # -------------------------------------------------------
  # Loaders — each returns a tibble with the canonical columns:
  #   pre_root_id, post_root_id, region, side, neuropil
  # -------------------------------------------------------
  load_v1 <- function() {
    message(sprintf("\n--- loading v1: %s ---", basename(v1_local)))
    df <- arrow::read_parquet(v1_local,
                              col_select = c("pre_root_id", "post_root_id",
                                             "neuropil", "size"))
    parsed <- parse_v1_neuropil(df$neuropil)
    tibble::tibble(
      pre_root_id  = as.character(df$pre_root_id),
      post_root_id = as.character(df$post_root_id),
      region       = parsed$region,
      side         = parsed$side,
      neuropil     = parsed$neuropil,
      size         = df$size
    )
  }

  load_v2 <- function() {
    message(sprintf("\n--- loading v2: %s + neuropil CSV ---", basename(v2_local)))
    df <- arrow::read_parquet(v2_local,
                              col_select = c("id", "pre_root_id",
                                             "post_root_id", "size",
                                             "pre_x", "pre_y", "pre_z")) %>%
      dplyr::transmute(
        id           = as.integer(id),
        pre_root_id  = as.character(pre_root_id),
        post_root_id = as.character(post_root_id),
        size         = size,
        # GCS v2 export stores pre coords in nm (verified 2026-05-05 from
        # written sample magnitudes: pre_x ~ 4 * Coord1_X). Converted to
        # BANC raw voxel via banc_nm2raw at sample-write time.
        pre_x = pre_x, pre_y = pre_y, pre_z = pre_z
      )
    csv_path <- file.path(banc.connectivity.save.path,
                          "banc_synapses_to_neuropils_v2.csv")
    stopifnot(file.exists(csv_path))
    # X, Y, Z brought in for the v2 stratified sample (post-loop).
    # compute_completion ignores them; they ride along for free.
    np <- readr::read_csv(csv_path,
                          col_types = readr::cols_only(
                            id       = readr::col_integer(),
                            X        = readr::col_double(),
                            Y        = readr::col_double(),
                            Z        = readr::col_double(),
                            side     = readr::col_character(),
                            region   = readr::col_character(),
                            neuropil = readr::col_character()
                          ))
    out <- df %>% dplyr::left_join(np, by = "id")
    n_unmatched <- sum(is.na(out$neuropil))
    if (n_unmatched > 0) {
      message(sprintf("  %s/%s v2 synapses had no neuropil match (%.2f%%); marked outside",
                      format(n_unmatched,  big.mark = ","),
                      format(nrow(out),    big.mark = ","),
                      100 * n_unmatched / nrow(out)))
    }
    out %>%
      dplyr::mutate(
        region   = ifelse(is.na(region),   "outside", region),
        side     = ifelse(is.na(side),     "centre",  side),
        neuropil = ifelse(is.na(neuropil), "outside", neuropil)
      ) %>%
      dplyr::transmute(
        syn_id       = id,
        pre_root_id  = pre_root_id,
        post_root_id = post_root_id,
        region       = region,
        side         = side,
        neuropil     = neuropil,
        size         = size,
        X = X, Y = Y, Z = Z,            # nm centroid (cleft midpoint)
        pre_x = pre_x, pre_y = pre_y, pre_z = pre_z   # nm (converted at sample write)
      )
  }

  load_v3 <- function() {
    message(sprintf("\n--- loading v3 from GCS [%s]: %s ---",
                    v3_source, basename(v3_local)))
    # v3 source: either the curated parquet (legacy versions, schema:
    # id, size, pre_root_id, post_root_id, pre_x/y/z) or the canonical
    # CSV.gz (v890+, schema: 15 cols matching v2 — see banc-calculate-
    # connectivity.R for the column list). Either way, CAVE-ingested
    # trusted root_ids.
    if (v3_source == "parquet") {
      gcs <- arrow::read_parquet(v3_local,
                                 col_select = c("id", "size",
                                                "pre_root_id", "post_root_id"))
    } else {                                          # csv
      column_names <- c('id', 'pre_x', 'pre_y', 'pre_z', 'post_x', 'post_y',
                        'post_z', 'ctr_x', 'ctr_y', 'ctr_z', 'size',
                        'pre_supervoxel_id', 'pre_root_id',
                        'post_supervoxel_id', 'post_root_id')
      col_types <- vroom::cols(
        id           = vroom::col_character(),
        size         = vroom::col_double(),
        pre_root_id  = vroom::col_character(),
        post_root_id = vroom::col_character(),
        .default     = vroom::col_double()
      )
      gcs <- vroom::vroom(v3_local,
                          col_names   = column_names,
                          col_select  = dplyr::all_of(c("id", "size",
                                                        "pre_root_id",
                                                        "post_root_id")),
                          col_types   = col_types,
                          skip        = 0)            # no header in v3 CSV
    }
    n_gcs <- nrow(gcs)
    message(sprintf("  GCS v3: %s rows", format(n_gcs, big.mark = ",")))

    # Spatial lookup from our processed v3 batches (Stage 2 output). Only use
    # syn_id + region/side/neuropil — IGNORE its root_id columns (corrupt).
    if (!file.exists(v3_spatial_local)) {
      stop(sprintf("v3 spatial lookup not found: %s\nRun stage 2 of banc/metrics/banc-synapses-v3-optimised.R first.",
                   v3_spatial_local))
    }
    message(sprintf("  loading spatial lookup: %s",
                    basename(v3_spatial_local)))
    # mean_score + median_score + X/Y/Z brought in for the post-loop
    # neuroglancer sample + score-distribution plots.
    # compute_completion() ignores them; they ride along for free.
    #
    # Use open_dataset + lazy select instead of read_parquet(col_select=...) —
    # job 38760712 SIGFPE'd inside read_parquet on the 15-column read of this
    # 6 GB / 259M-row parquet (column-set widening from 5 → 15 tipped it over
    # an int64-inflation bug, same family as Fix C in banc-synapses-v3-optimised.R
    # Stage 4). open_dataset takes a different code path (Arrow C++ scanner)
    # that doesn't trip the bug.
    spatial <- arrow::open_dataset(v3_spatial_local) %>%
      dplyr::select(syn_id, region, side, neuropil,
                    mean_score, median_score,
                    X, Y, Z,
                    presyn_x, presyn_y, presyn_z,
                    postsyn_x, postsyn_y, postsyn_z) %>%
      dplyr::collect() %>%
      dplyr::distinct(syn_id, .keep_all = TRUE)   # belt-and-braces dedup
    message(sprintf("  spatial lookup: %s rows", format(nrow(spatial),
                                                          big.mark = ",")))

    # Inner join: keep only GCS rows that have a spatial match. Hard-fail
    # if coverage drops below 99% — we'd rather notice than ship silently.
    # GCS `id` is int64 (bit64); our spatial `syn_id` was stored as double
    # — both exactly represent the same integer values (largest seen ~5.5e12,
    # well below 2^53 = 9e15), so coercing GCS id to double is lossless and
    # makes dplyr happy.
    gcs <- gcs %>% dplyr::mutate(id = as.numeric(id))
    joined <- gcs %>%
      dplyr::inner_join(spatial, by = c("id" = "syn_id"))
    coverage <- nrow(joined) / n_gcs
    message(sprintf("  spatial join coverage: %s / %s = %.4f%%",
                    format(nrow(joined), big.mark = ","),
                    format(n_gcs,        big.mark = ","),
                    100 * coverage))
    if (coverage < 0.99) {
      stop(sprintf(paste0("v3 spatial join coverage too low: %.2f%% < 99%%. ",
                          "Likely cause: spatial lookup was built with ",
                          "MIN_SIZE > 0 and GCS contains smaller detections, ",
                          "or syn_id mismatches. Re-run Stage 2 of ",
                          "banc-synapses-v3-optimised.R with MIN_SIZE = 0L."),
                   100 * coverage))
    }

    joined %>% dplyr::transmute(
      syn_id       = id,           # GCS join key, retained for v3 sample export
      pre_root_id  = as.character(pre_root_id),
      post_root_id = as.character(post_root_id),
      region       = ifelse(is.na(region),   "outside", region),
      side         = ifelse(is.na(side),     "centre",  side),
      neuropil     = ifelse(is.na(neuropil) | neuropil == "", "outside", neuropil),
      size         = size,
      mean_score   = mean_score,     # for stratified sample + distribution plots
      median_score = median_score,   # for distribution plots
      X = X, Y = Y, Z = Z,           # nm coords for neuroglancer (centroid)
      # presyn / postsyn in V3 voxel units (16x16x45 nm/voxel) — converted to
      # BANC raw voxel space at sample-write time.
      presyn_x = presyn_x, presyn_y = presyn_y, presyn_z = presyn_z,
      postsyn_x = postsyn_x, postsyn_y = postsyn_y, postsyn_z = postsyn_z
    )
  }

  # -------------------------------------------------------
  # Per-source completion metrics, parameterised by size threshold.
  # threshold > 0 appends `_size_thresh_<N>` to every output filename so
  # unfiltered (threshold = 0) and filtered passes co-exist on disk.
  # -------------------------------------------------------
  thresh_suffix <- function(threshold) {
    if (threshold <= 0L) "" else sprintf("_size_thresh_%d", threshold)
  }

  compute_completion <- function(syns_full, prefix, threshold) {
    if (threshold > 0L) {
      n_pre <- nrow(syns_full)
      syns <- syns_full %>% dplyr::filter(size >= threshold)
      message(sprintf("  [%s] size>=%d filter: kept %s/%s rows (%.2f%%)",
                      prefix, threshold,
                      format(nrow(syns), big.mark = ","),
                      format(n_pre,      big.mark = ","),
                      100 * nrow(syns) / n_pre))
    } else {
      syns <- syns_full
    }

    syns <- syns %>%
      dplyr::mutate(
        pre_status  = classify_status(pre_root_id),
        post_status = classify_status(post_root_id)
      )

    suffix <- thresh_suffix(threshold)

    write_rate <- function(df, name) {
      f <- file.path(rates_dir,
                     sprintf("banc_%s_%s_%s_capture_rates%s.csv",
                             target.version, prefix, name, suffix))
      readr::write_csv(df, f)
      message(sprintf("  [%s] saved %s (%s rows)",
                      prefix, basename(f),
                      format(nrow(df), big.mark = ",")))
      f
    }

    gross <- syns %>%
      dplyr::count(pre_status, post_status) %>%
      dplyr::mutate(prop = round(n / sum(n), 4))
    write_rate(gross, "gross")

    inout <- syns %>%
      dplyr::mutate(in_mesh = ifelse(region == "outside", "outside", "inside")) %>%
      dplyr::group_by(in_mesh) %>%
      dplyr::count(pre_status, post_status, in_mesh) %>%
      dplyr::mutate(prop = round(n / sum(n), 4)) %>%
      dplyr::ungroup()
    write_rate(inout, "inout")

    region_rates <- syns %>%
      dplyr::group_by(region, side) %>%
      dplyr::count(pre_status, post_status, side, region) %>%
      dplyr::mutate(prop = round(n / sum(n), 4)) %>%
      dplyr::ungroup()
    write_rate(region_rates, "region")

    neuropil_rates <- syns %>%
      dplyr::group_by(region, side, neuropil) %>%
      dplyr::count(pre_status, post_status, side, region, neuropil) %>%
      dplyr::mutate(prop = round(n / sum(n), 4)) %>%
      dplyr::ungroup()
    write_rate(neuropil_rates, "neuropil")

    captured <- function(s) s %in% c("proofread", "identified")
    n_total      <- nrow(syns)
    summary_row <- tibble::tibble(
      source              = prefix,
      size_threshold      = threshold,
      n_total             = n_total,
      n_outside           = sum(syns$region == "outside"),
      pct_outside         = round(100 * sum(syns$region == "outside") / n_total, 2),
      n_pre_captured      = sum(captured(syns$pre_status)),
      pct_pre_captured    = round(100 * sum(captured(syns$pre_status)) / n_total, 2),
      n_post_captured     = sum(captured(syns$post_status)),
      pct_post_captured   = round(100 * sum(captured(syns$post_status)) / n_total, 2),
      n_both_captured     = sum(captured(syns$pre_status) & captured(syns$post_status)),
      pct_both_captured   = round(100 * sum(captured(syns$pre_status) &
                                              captured(syns$post_status)) / n_total, 2),
      n_pre_proofread     = sum(syns$pre_status  == "proofread"),
      pct_pre_proofread   = round(100 * sum(syns$pre_status  == "proofread") / n_total, 2),
      n_post_proofread    = sum(syns$post_status == "proofread"),
      pct_post_proofread  = round(100 * sum(syns$post_status == "proofread") / n_total, 2),
      n_both_proofread    = sum(syns$pre_status == "proofread" &
                                  syns$post_status == "proofread"),
      pct_both_proofread  = round(100 * sum(syns$pre_status == "proofread" &
                                              syns$post_status == "proofread") / n_total, 2)
    )

    message(sprintf("\n  [%s thresh>=%d] %s synapses; outside=%.1f%%; captured pre=%.1f%% post=%.1f%% both=%.1f%%; proofread pre=%.1f%% post=%.1f%% both=%.1f%%",
                    prefix, threshold, format(n_total, big.mark = ","),
                    summary_row$pct_outside,
                    summary_row$pct_pre_captured,
                    summary_row$pct_post_captured,
                    summary_row$pct_both_captured,
                    summary_row$pct_pre_proofread,
                    summary_row$pct_post_proofread,
                    summary_row$pct_both_proofread))

    list(summary = summary_row, region = region_rates)
  }

  # -------------------------------------------------------
  # Run all three sources × all thresholds.
  # Each source loaded ONCE; compute_completion applies the threshold to a
  # cheap subset copy and emits suffix-tagged outputs.
  # -------------------------------------------------------
  results <- list()  # results[[src]][[as.character(threshold)]] = list(summary, region)
  for (src in c("v1", "v2", "v3")) {
    syns_full <- switch(src, v1 = load_v1(), v2 = load_v2(), v3 = load_v3())
    message(sprintf("  [%s] loaded %s rows", src,
                    format(nrow(syns_full), big.mark = ",")))
    results[[src]] <- list()
    for (threshold in SIZE_THRESHOLDS) {
      results[[src]][[as.character(threshold)]] <-
        compute_completion(syns_full, src, threshold)
    }

    # -----------------------------------------------------
    # v2 stratified neuroglancer sample. Same recipe as v3 but binned on
    # `size` (v2 has no mean_score; size is the natural quality proxy and
    # spans a wide range, all the way down to size=1 in the v2 detector).
    # No size threshold applied — matches v3 sample's all-detections scope.
    #
    # Coordinate convention (BANC raw voxel space, 4x4x45 nm/voxel):
    #   - Coordinate 1: cleft centroid (best post-synapse proxy available
    #     for v2 — the GCS export only carries presyn coords; centroid sits
    #     ~half-a-cleft from true postsyn, fine for visual review).
    #   - Coordinate 2: presyn position (in nm from GCS export; converted
    #                   to BANC raw voxel via banc_nm2raw to match Coord 1).
    #   Type = "Line" so the pair renders as a pre→post vector.
    #
    # Two files are written for v2:
    #   - banc_<ver>_v2_synapse_sample_<date>.csv                          all sizes
    #   - banc_<ver>_v2_synapse_sample_lessthanorequalto_5_<date>.csv      size <= 5 only
    # -----------------------------------------------------
    if (src == "v2") {
      build_and_save_v2_sample <- function(pool, label) {
        message(sprintf("\n  [%s/%s] drawing stratified sample (%d per (np_primary x side), %d size bins)...",
                        src, label, V3_SAMPLE_PER_GROUP, V3_SAMPLE_BINS))

        # Per-pool size quantile bins. Empty quantile boundaries collapse via
        # unique(); cut() handles a single bin gracefully.
        size_breaks <- unique(stats::quantile(
          pool$size,
          probs = seq(0, 1, length.out = V3_SAMPLE_BINS + 1),
          na.rm = TRUE
        ))
        n_size_bins <- length(size_breaks) - 1L
        message(sprintf("    global size bins (quantile, n=%d): [%g, %g]",
                        n_size_bins, size_breaks[1],
                        size_breaks[length(size_breaks)]))
        message(sprintf("    bin edges: %s",
                        paste(size_breaks, collapse = ", ")))

        set.seed(42)
        if (n_size_bins >= 1L) {
          sampled <- pool %>%
            dplyr::mutate(
              size_bin = cut(size, breaks = size_breaks,
                             include.lowest = TRUE, labels = FALSE)
            ) %>%
            dplyr::group_by(np_primary, side, size_bin) %>%
            dplyr::slice_sample(n = ceiling(V3_SAMPLE_PER_GROUP / V3_SAMPLE_BINS),
                                replace = FALSE) %>%
            dplyr::ungroup() %>%
            as.data.frame(stringsAsFactors = FALSE)
        } else {
          sampled <- pool %>%
            dplyr::mutate(size_bin = 1L) %>%
            dplyr::group_by(np_primary, side) %>%
            dplyr::slice_sample(n = V3_SAMPLE_PER_GROUP, replace = FALSE) %>%
            dplyr::ungroup() %>%
            as.data.frame(stringsAsFactors = FALSE)
        }
        n_groups <- length(unique(paste(sampled$np_primary, sampled$side)))
        message(sprintf("    sampled %s syns from %d (np_primary x side) groups",
                        format(nrow(sampled), big.mark = ","), n_groups))

        # Coordinate 1: cleft centroid in BANC raw voxel (4x4x45 nm/voxel).
        # X/Y/Z are in nm — convert via banc_nm2raw().
        coord1 <- bancr::banc_nm2raw(cbind(sampled$X, sampled$Y, sampled$Z))
        sampled$`Coordinate 1` <- sprintf("(%.0f, %.0f, %.0f)",
                                           coord1[, 1], coord1[, 2], coord1[, 3])
        # Coordinate 2: presyn position. The v2 GCS export stores pre_x/y/z
        # in nm (not raw voxel as previously assumed); convert via banc_nm2raw
        # so Coord 2 is in the same BANC raw voxel space as Coord 1.
        coord2 <- bancr::banc_nm2raw(cbind(sampled$pre_x, sampled$pre_y, sampled$pre_z))
        sampled$`Coordinate 2` <- sprintf("(%.0f, %.0f, %.0f)",
                                           coord2[, 1], coord2[, 2], coord2[, 3])

        ngl_ann <- data.frame(
          `Coordinate 1`         = sampled$`Coordinate 1`,
          `Coordinate 2`         = sampled$`Coordinate 2`,
          `Ellipsoid Dimensions` = "",
          Tags                   = "",
          Description            = sampled$np_primary,
          `Segment IDs`          = "",
          `Parent ID`            = "",
          Type                   = "Line",
          ID                     = sprintf("%.0f", sampled$syn_id),
          syn_id                 = sprintf("%.0f", sampled$syn_id),
          size                   = sampled$size,
          size_bin               = sampled$size_bin,
          region                 = sampled$region,
          side                   = sampled$side,
          neuropil               = sampled$np_primary,
          pre_root_id            = sampled$pre_root_id,
          post_root_id           = sampled$post_root_id,
          pre_status             = sampled$pre_status,
          post_status            = sampled$post_status,
          stringsAsFactors       = FALSE
        )
        colnames(ngl_ann) <- gsub("\\.", " ", colnames(ngl_ann))

        sample_date <- format(Sys.time(), "%Y-%m-%d")
        suffix <- if (label == "all") "" else sprintf("_%s", label)
        ngl_file <- file.path(rates_dir,
                              sprintf("banc_%s_v2_synapse_sample%s_%s.csv",
                                      target.version, suffix, sample_date))
        readr::write_excel_csv(ngl_ann, ngl_file)
        message(sprintf("    saved %s (%s rows)", basename(ngl_file),
                        format(nrow(ngl_ann), big.mark = ",")))
      }

      v2_pool <- syns_full %>%
        dplyr::mutate(
          np_primary = ifelse(neuropil == "" | is.na(neuropil), "outside",
                              sub(",.*", "", neuropil)),
          pre_status  = classify_status(pre_root_id),
          post_status = classify_status(post_root_id)
        )

      # File 1: all sizes
      build_and_save_v2_sample(v2_pool, "all")

      # File 2: size <= 5 only
      v2_pool_small <- v2_pool %>% dplyr::filter(size <= 5)
      message(sprintf("  [%s/lessthanorequalto_5] %s rows out of %s with size <= 5",
                      src,
                      format(nrow(v2_pool_small), big.mark = ","),
                      format(nrow(v2_pool),       big.mark = ",")))
      if (nrow(v2_pool_small) > 0L) {
        build_and_save_v2_sample(v2_pool_small, "lessthanorequalto_5")
      } else {
        message("    no v2 synapses with size <= 5 — skipping size-thresh file")
      }

      rm(v2_pool, v2_pool_small); gc(verbose = FALSE)
    }

    # -----------------------------------------------------
    # v3-only: stratified neuroglancer sample (no size threshold).
    # Outside the threshold loop — runs once per completion job.
    # -----------------------------------------------------
    if (src == "v3") {
      message(sprintf("\n  [%s] drawing stratified sample (%d per (np_primary x side), %d score bins)...",
                      src, V3_SAMPLE_PER_GROUP, V3_SAMPLE_BINS))

      v3_pool <- syns_full %>%
        dplyr::mutate(
          np_primary = ifelse(neuropil == "" | is.na(neuropil), "outside",
                              sub(",.*", "", neuropil)),
          # classify pre/post status for the saved sample (so reviewers see
          # which v3 detections sit on identified vs fragment cells)
          pre_status  = classify_status(pre_root_id),
          post_status = classify_status(post_root_id)
        )

      # GLOBAL score bins — defined once across the whole v3 pool, then
      # applied uniformly to every neuropil x side group. So a synapse with
      # mean_score=X falls in the same bin in OL_ME_R as in CB_AVLP_L.
      #
      # Quantile breaks (was equal-width until 2026-04-20). Reason: v3 mean_score
      # is heavily right-skewed — equal-width bins put 4-5 of 10 bins in
      # near-empty regions of the score range, leaving the high-quality tail
      # (and the threshold-curve high-x end) statistically meaningless. Quantile
      # breaks give ~equal N per bin (~26M syns each) so sampling and the
      # threshold curve are stable across the full distribution.
      # `unique()` collapses any duplicate quantile boundaries from heavy ties;
      # if that drops bin count below V3_SAMPLE_BINS, downstream code adapts.
      score_breaks <- unique(stats::quantile(
        v3_pool$mean_score,
        probs = seq(0, 1, length.out = V3_SAMPLE_BINS + 1),
        na.rm = TRUE
      ))
      n_bins_effective <- length(score_breaks) - 1L
      message(sprintf("    global score bins (quantile, n=%d): [%.4f, %.4f]",
                      n_bins_effective,
                      score_breaks[1], score_breaks[length(score_breaks)]))
      message(sprintf("    bin edges: %s",
                      paste(sprintf("%.4f", score_breaks), collapse = ", ")))

      # Equal-width breaks for the histogram plot (kept separately — equal-width
      # is the right choice for visualizing the distribution shape; using
      # quantile breaks would hide the skew).
      score_breaks_equal <- seq(min(v3_pool$mean_score, na.rm = TRUE),
                                 max(v3_pool$mean_score, na.rm = TRUE),
                                 length.out = V3_SAMPLE_BINS + 1)

      set.seed(42)
      v3_sampled <- v3_pool %>%
        dplyr::mutate(
          score_bin = cut(mean_score, breaks = score_breaks,
                          include.lowest = TRUE, labels = FALSE)
        ) %>%
        dplyr::group_by(np_primary, side, score_bin) %>%
        dplyr::slice_sample(n = ceiling(V3_SAMPLE_PER_GROUP / V3_SAMPLE_BINS),
                            replace = FALSE) %>%
        dplyr::ungroup() %>%
        as.data.frame(stringsAsFactors = FALSE)

      n_groups <- length(unique(paste(v3_sampled$np_primary, v3_sampled$side)))
      message(sprintf("    sampled %s syns from %d (np_primary x side) groups",
                      format(nrow(v3_sampled), big.mark = ","), n_groups))

      # Neuroglancer-formatted CSV (Coordinate columns + metadata).
      # Coordinate convention (BANC raw voxel space, 4x4x45 nm/voxel):
      #   Coordinate 1: TRUE postsynaptic point (postsyn_x/y/z, V3 voxel
      #                  16x16x45 nm — multiply by V3_VOXEL_NM = (16,16,45)
      #                  to get nm, then banc_nm2raw to BANC voxel).
      #   Coordinate 2: TRUE presynaptic point (presyn_x/y/z, same chain).
      # Type = "Line" so the pair renders as a pre→post vector.
      V3_VOXEL_NM <- c(16, 16, 45)
      post_nm <- cbind(v3_sampled$postsyn_x * V3_VOXEL_NM[1],
                       v3_sampled$postsyn_y * V3_VOXEL_NM[2],
                       v3_sampled$postsyn_z * V3_VOXEL_NM[3])
      pre_nm  <- cbind(v3_sampled$presyn_x  * V3_VOXEL_NM[1],
                       v3_sampled$presyn_y  * V3_VOXEL_NM[2],
                       v3_sampled$presyn_z  * V3_VOXEL_NM[3])
      coord1 <- bancr::banc_nm2raw(post_nm)
      coord2 <- bancr::banc_nm2raw(pre_nm)
      v3_sampled$`Coordinate 1` <- sprintf("(%.0f, %.0f, %.0f)",
                                            coord1[, 1], coord1[, 2], coord1[, 3])
      v3_sampled$`Coordinate 2` <- sprintf("(%.0f, %.0f, %.0f)",
                                            coord2[, 1], coord2[, 2], coord2[, 3])
      ngl_ann <- data.frame(
        `Coordinate 1`         = v3_sampled$`Coordinate 1`,
        `Coordinate 2`         = v3_sampled$`Coordinate 2`,
        `Ellipsoid Dimensions` = "",
        Tags                   = "",
        Description            = v3_sampled$np_primary,
        `Segment IDs`          = "",
        `Parent ID`            = "",
        Type                   = "Line",
        ID                     = sprintf("%.0f", v3_sampled$syn_id),
        syn_id                 = sprintf("%.0f", v3_sampled$syn_id),
        size                   = v3_sampled$size,
        mean_score             = round(v3_sampled$mean_score, 4),
        score_bin              = v3_sampled$score_bin,
        region                 = v3_sampled$region,
        side                   = v3_sampled$side,
        neuropil               = v3_sampled$np_primary,
        pre_root_id            = v3_sampled$pre_root_id,
        post_root_id           = v3_sampled$post_root_id,
        pre_status             = v3_sampled$pre_status,
        post_status            = v3_sampled$post_status,
        stringsAsFactors       = FALSE
      )
      colnames(ngl_ann) <- gsub("\\.", " ", colnames(ngl_ann))

      sample_date <- format(Sys.time(), "%Y-%m-%d")
      ngl_file <- file.path(rates_dir,
                            sprintf("banc_%s_v3_synapse_sample_%s.csv",
                                    target.version, sample_date))
      readr::write_excel_csv(ngl_ann, ngl_file)
      message(sprintf("    saved %s (%s rows)", basename(ngl_file),
                      format(nrow(ngl_ann), big.mark = ",")))

      # ---------------------------------------------------
      # mean_score threshold curve: id+proof rate vs threshold,
      # evaluated at each of the global bin edges. Uses cumulative
      # tallies so it's O(N) once + O(B) for the curve.
      # ---------------------------------------------------
      message("\n    building mean_score threshold curve...")
      v3_for_curve <- v3_pool %>%
        dplyr::mutate(
          score_bin    = cut(mean_score, breaks = score_breaks,
                             include.lowest = TRUE, labels = FALSE),
          pre_idproof  = pre_status  %in% c("identified", "proofread"),
          post_idproof = post_status %in% c("identified", "proofread"),
          both_idproof = pre_idproof & post_idproof
        )

      bin_tally <- v3_for_curve %>%
        dplyr::group_by(score_bin) %>%
        dplyr::summarise(
          n      = dplyr::n(),
          n_pre  = sum(pre_idproof),
          n_post = sum(post_idproof),
          n_both = sum(both_idproof),
          .groups = "drop"
        ) %>%
        dplyr::arrange(score_bin)

      # Cumulative from the high-score end down: row b reflects "all syns
      # with score_bin >= b" (i.e. mean_score >= score_breaks[b]).
      curve <- bin_tally %>%
        dplyr::arrange(dplyr::desc(score_bin)) %>%
        dplyr::mutate(
          cum_n    = cumsum(n),
          cum_pre  = cumsum(n_pre),
          cum_post = cumsum(n_post),
          cum_both = cumsum(n_both)
        ) %>%
        dplyr::arrange(score_bin) %>%
        dplyr::mutate(
          mean_score_threshold = score_breaks[score_bin],
          pct_pre_idproof      = round(100 * cum_pre  / cum_n, 3),
          pct_post_idproof     = round(100 * cum_post / cum_n, 3),
          pct_both_idproof     = round(100 * cum_both / cum_n, 3)
        ) %>%
        dplyr::select(score_bin, mean_score_threshold,
                      cum_n, cum_pre, cum_post, cum_both,
                      pct_pre_idproof, pct_post_idproof, pct_both_idproof)

      curve_csv <- file.path(rates_dir,
                              sprintf("banc_%s_v3_score_threshold_curve.csv",
                                      target.version))
      readr::write_csv(curve, curve_csv)
      message(sprintf("    saved %s", basename(curve_csv)))
      print(knitr::kable(curve, format = "pipe"))

      # Plot — three lines, x = mean_score threshold (lower bin edge), y = pct id+proof
      curve_long <- curve %>%
        tidyr::pivot_longer(cols = c(pct_pre_idproof, pct_post_idproof, pct_both_idproof),
                            names_to = "side", values_to = "pct") %>%
        dplyr::mutate(side = factor(
          dplyr::recode(side,
                        pct_pre_idproof  = "pre",
                        pct_post_idproof = "post",
                        pct_both_idproof = "both"),
          levels = c("pre", "post", "both")
        ))

      curve_plot <- ggplot2::ggplot(
        curve_long,
        ggplot2::aes(x = mean_score_threshold, y = pct, colour = side, group = side)
      ) +
        ggplot2::geom_line(linewidth = 0.9) +
        ggplot2::geom_point(size = 2) +
        ggplot2::scale_colour_manual(values = c(pre = "#1f77b4",
                                                  post = "#d62728",
                                                  both = "#2ca02c")) +
        ggplot2::scale_y_continuous(limits = c(0, 100),
                                    breaks = seq(0, 100, 10)) +
        ggplot2::labs(
          title    = sprintf("v3 (banc_%s): id+proof capture vs mean_score threshold",
                              target.version),
          subtitle = sprintf("Cumulative over %d global score bins; each x is the lower edge of a bin",
                              V3_SAMPLE_BINS),
          x        = "mean_score threshold (synapses with mean_score >= x)",
          y        = "id+proof capture rate (%)",
          colour   = "side"
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(legend.position = "bottom",
                       plot.title       = ggplot2::element_text(face = "bold"))

      curve_png <- file.path(rates_dir,
                              sprintf("banc_%s_v3_score_threshold_curve.png",
                                      target.version))
      ggplot2::ggsave(curve_png, curve_plot,
                      width = 7, height = 5, dpi = 150)
      message(sprintf("    saved %s", basename(curve_png)))

      # ---------------------------------------------------
      # Score distributions: histograms + density plots.
      # - mean_score histogram: uses `score_breaks` (same 10 bins as the
      #   stratified sample, for visual consistency).
      # - median_score histogram: its own 10 equal-width bins (different
      #   variable, different range — tying it to mean_score's breaks would
      #   dump most synapses into a few bins or leave them unbinned).
      # - densities: geom_density on ~260M syns is prohibitive (kernel
      #   estimation is O(N)); downsample to 1M for smooth visualization.
      # ---------------------------------------------------
      message("\n    plotting score distributions (histograms + densities)...")

      # Derive median_score breaks for its own histogram
      median_breaks <- seq(min(syns_full$median_score, na.rm = TRUE),
                           max(syns_full$median_score, na.rm = TRUE),
                           length.out = V3_SAMPLE_BINS + 1)

      n_total_v3 <- nrow(syns_full)
      set.seed(17)
      density_n <- min(n_total_v3, 1e6L)
      v3_density_sample <- syns_full %>%
        dplyr::slice_sample(n = density_n) %>%
        dplyr::select(mean_score, median_score)

      save_score_plot <- function(plot, name) {
        f <- file.path(rates_dir,
                       sprintf("banc_%s_v3_%s.png", target.version, name))
        ggplot2::ggsave(f, plot, width = 7, height = 4.5, dpi = 150)
        message(sprintf("    saved %s", basename(f)))
      }

      # Histogram: mean_score with equal-width bins (separate from sample's
      # quantile bins — equal-width is correct for showing the skew).
      p_hist_mean <- ggplot2::ggplot(syns_full,
                                     ggplot2::aes(x = mean_score)) +
        ggplot2::geom_histogram(breaks = score_breaks_equal,
                                fill = "#1f77b4", colour = "white") +
        ggplot2::labs(
          title    = sprintf("v3 (banc_%s): mean_score histogram",
                              target.version),
          subtitle = sprintf("n = %s synapses, %d equal-width bins (sampling uses quantile bins)",
                              format(n_total_v3, big.mark = ","),
                              V3_SAMPLE_BINS),
          x        = "mean_score",
          y        = "count"
        ) +
        ggplot2::theme_minimal(base_size = 12)
      save_score_plot(p_hist_mean, "mean_score_histogram")

      # Histogram: median_score (own bins)
      p_hist_median <- ggplot2::ggplot(syns_full,
                                       ggplot2::aes(x = median_score)) +
        ggplot2::geom_histogram(breaks = median_breaks,
                                fill = "#d62728", colour = "white") +
        ggplot2::labs(
          title    = sprintf("v3 (banc_%s): median_score histogram",
                              target.version),
          subtitle = sprintf("n = %s synapses, %d bins over median_score range",
                              format(n_total_v3, big.mark = ","),
                              V3_SAMPLE_BINS),
          x        = "median_score",
          y        = "count"
        ) +
        ggplot2::theme_minimal(base_size = 12)
      save_score_plot(p_hist_median, "median_score_histogram")

      # Density: mean_score (downsampled)
      p_dens_mean <- ggplot2::ggplot(v3_density_sample,
                                     ggplot2::aes(x = mean_score)) +
        ggplot2::geom_density(fill = "#1f77b4", alpha = 0.5,
                               colour = "#1f77b4") +
        ggplot2::labs(
          title    = sprintf("v3 (banc_%s): mean_score density",
                              target.version),
          subtitle = sprintf("density estimated on %s-row random sample (from %s total)",
                              format(density_n, big.mark = ","),
                              format(n_total_v3, big.mark = ",")),
          x        = "mean_score",
          y        = "density"
        ) +
        ggplot2::theme_minimal(base_size = 12)
      save_score_plot(p_dens_mean, "mean_score_density")

      # Density: median_score (same sample)
      p_dens_median <- ggplot2::ggplot(v3_density_sample,
                                       ggplot2::aes(x = median_score)) +
        ggplot2::geom_density(fill = "#d62728", alpha = 0.5,
                               colour = "#d62728") +
        ggplot2::labs(
          title    = sprintf("v3 (banc_%s): median_score density",
                              target.version),
          subtitle = sprintf("density estimated on %s-row random sample (from %s total)",
                              format(density_n, big.mark = ","),
                              format(n_total_v3, big.mark = ",")),
          x        = "median_score",
          y        = "density"
        ) +
        ggplot2::theme_minimal(base_size = 12)
      save_score_plot(p_dens_median, "median_score_density")

      rm(v3_pool, v3_sampled, ngl_ann, v3_for_curve, bin_tally, curve,
         curve_long, curve_plot, v3_density_sample,
         p_hist_mean, p_hist_median, p_dens_mean, p_dens_median);
      gc(verbose = FALSE)
    }

    rm(syns_full); gc(verbose = FALSE)
  }

  # -------------------------------------------------------
  # 3-way summaries — one per threshold
  # -------------------------------------------------------
  message("\n--- 3-way comparison summaries ---")

  for (threshold in SIZE_THRESHOLDS) {
    suffix <- thresh_suffix(threshold)
    summary_3way <- dplyr::bind_rows(lapply(results, function(src_results) {
      src_results[[as.character(threshold)]]$summary
    }))
    summary_file <- file.path(rates_dir,
                              sprintf("banc_%s_v1_v2_v3_summary%s.csv",
                                      target.version, suffix))
    readr::write_csv(summary_3way, summary_file)
    message(sprintf("  saved %s", basename(summary_file)))
    print(knitr::kable(summary_3way, format = "pipe"))
  }

  # Region+side capture %, joined across sources for the columns that exist in all.
  # Aggregate counts within each (source, region, side) and compute a captured %.
  # One wide-pivoted CSV per threshold.
  for (threshold in SIZE_THRESHOLDS) {
    suffix <- thresh_suffix(threshold)
    region_capt <- dplyr::bind_rows(lapply(names(results), function(src) {
      r <- results[[src]][[as.character(threshold)]]$region
      r %>%
        dplyr::group_by(region, side) %>%
        dplyr::summarise(
          n_total           = sum(n),
          n_pre_captured    = sum(n[pre_status  %in% c("proofread", "identified")]),
          n_post_captured   = sum(n[post_status %in% c("proofread", "identified")]),
          n_both_captured   = sum(n[pre_status  %in% c("proofread", "identified") &
                                     post_status %in% c("proofread", "identified")]),
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          source              = src,
          pct_pre_captured    = round(100 * n_pre_captured  / n_total, 2),
          pct_post_captured   = round(100 * n_post_captured / n_total, 2),
          pct_both_captured   = round(100 * n_both_captured / n_total, 2)
        )
    }))

    region_wide <- region_capt %>%
      dplyr::select(source, region, side,
                    n_total, pct_pre_captured, pct_post_captured, pct_both_captured) %>%
      tidyr::pivot_wider(
        names_from  = source,
        values_from = c(n_total, pct_pre_captured, pct_post_captured, pct_both_captured),
        names_glue  = "{.value}_{source}"
      ) %>%
      dplyr::arrange(region, side)

    region_file <- file.path(rates_dir,
                             sprintf("banc_%s_v1_v2_v3_region_summary%s.csv",
                                     target.version, suffix))
    readr::write_csv(region_wide, region_file)
    message(sprintf("  saved %s", basename(region_file)))
  }

  # -------------------------------------------------------
  # Proportion-on-cell ECDF (Fig 1 supplement)
  # -------------------------------------------------------
  # Cumulative share of synapses on identified (banc_meta) root_ids vs
  # fragments, as a function of per-root_id synapse count. Answers
  # "above what synapse-count threshold is most of the network actually
  # neurons?" Replaces banc/legacy/banc-synapse-proportion-plot.R
  # (hardcoded v821, raw 30 GB CSV, 64 G OOM-prone).
  #
  # Source: banc_<ver>_synapses_v2_enriched.parquet (banc/meta/banc-data.R).
  # That table includes synapses where at least one endpoint is in root.ids
  # — exactly what we need to see the fragment tail.
  message("\n--- Proportion-on-cell ECDF ---")

  # banc-data.R writes the enriched parquet under the versioned save path;
  # an older location (banc.save.path/) is checked first for backward compat.
  ecdf_candidates <- c(
    file.path(banc.save.path,
              sprintf("banc_%s_synapses_v2_enriched.parquet", target.version)),
    file.path("/n/data1/hms/neurobio/wilson/connectomes/banc",
              paste0("banc_", target.version),
              sprintf("banc_%s_synapses_v2_enriched.parquet", target.version))
  )
  ecdf_parquet <- ecdf_candidates[file.exists(ecdf_candidates) &
                                    file.info(ecdf_candidates)$size > 0L][1]
  if (is.na(ecdf_parquet)) {
    message(sprintf("  skipped: enriched v2 parquet not found at any of:\n    - %s\n  (run banc/meta/banc-data.R --source v2 first)",
                    paste(ecdf_candidates, collapse = "\n    - ")))
  } else {
    t0 <- Sys.time()
    message(sprintf("  reading %s ...", basename(ecdf_parquet)))
    syns_v2 <- arrow::read_parquet(ecdf_parquet,
                                    col_select = c("pre_root_id", "post_root_id"))
    syns_v2$pre_root_id  <- as.character(syns_v2$pre_root_id)
    syns_v2$post_root_id <- as.character(syns_v2$post_root_id)
    message(sprintf("  read %s rows in %s",
                    format(nrow(syns_v2), big.mark = ","),
                    format(round(difftime(Sys.time(), t0, units = "mins"), 1))))

    process_root_id <- function(df, root_col) {
      counts <- df %>%
        dplyr::count(root_id = .data[[root_col]], name = "n_syn")
      counts$status <- ifelse(counts$root_id %in% identified.ids,
                              "neuron", "fragment")
      # Each root_id with n_syn = k contributes k synapses to bin k. Sum to
      # get total synapses per (status, n_syn), then cumulate.
      counts %>%
        dplyr::group_by(status, n_syn) %>%
        dplyr::summarise(n = sum(n_syn), .groups = "drop") %>%
        dplyr::arrange(n_syn) %>%
        tidyr::pivot_wider(names_from = status, values_from = n, values_fill = 0) %>%
        dplyr::mutate(
          cum_neuron   = if ("neuron"   %in% names(.)) cumsum(neuron)   else rep(0, dplyr::n()),
          cum_fragment = if ("fragment" %in% names(.)) cumsum(fragment) else rep(0, dplyr::n()),
          cum_total    = cum_neuron + cum_fragment,
          pct_neuron   = dplyr::if_else(cum_total > 0,
                                         cum_neuron / cum_total,
                                         NA_real_),
          root_id_col  = root_col
        )
    }

    dat_pre  <- process_root_id(syns_v2, "pre_root_id")
    dat_post <- process_root_id(syns_v2, "post_root_id")
    syn_comp <- dplyr::bind_rows(dat_pre, dat_post)
    rm(syns_v2); gc(verbose = FALSE)

    g.ecdf <- ggplot2::ggplot(syn_comp,
                              ggplot2::aes(x = n_syn, y = pct_neuron,
                                            color = root_id_col)) +
      ggplot2::geom_line(linewidth = 1.2) +
      ggplot2::scale_x_log10() +
      ggplot2::annotation_logticks(sides = "b") +
      ggplot2::labs(
        x     = "number of synapses per root_id (threshold)",
        y     = "proportion of synapses in neurons (≤ threshold)",
        color = "",
        title = sprintf("cumulative share of synapses on identified neurons (banc_%s, v2)",
                        target.version)
      ) +
      ggplot2::scale_color_manual(values = c(pre_root_id  = paper.cols[["pre"]],
                                              post_root_id = paper.cols[["post"]])) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(legend.position = "none")

    # Persist the long-format table too — small, useful for downstream
    # paper-figure tweaks without re-reading the 10 GB parquet.
    readr::write_csv(syn_comp,
                     file.path(rates_dir,
                               sprintf("banc_%s_synapse_proportion_on_cell.csv",
                                       target.version)))
    for (.ext in c("png", "pdf")) {
      out <- file.path(rates_dir,
                       sprintf("banc_%s_synapse_proportion_on_cell.%s",
                               target.version, .ext))
      ggplot2::ggsave(plot = g.ecdf, filename = out,
                      width = 6, height = 4, dpi = 300, bg = "transparent")
      message(sprintf("  saved %s", basename(out)))
    }
  }

  message("\n### banc: completion metrics done ###")
})
