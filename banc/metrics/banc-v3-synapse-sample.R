#' banc-v3-synapse-sample — Regenerate the v3 synapse spatial sample CSV.
#'
#' Standalone replacement for the v3 sample block in
#' `banc-calculate-completion.R`, which hit the SLURM-only arrow C++ URI bug
#' ("Invalid: Unrecognized filesystem type in URI: file:///_") on the 15 GB
#' spatial parquet. Uses narrow parquet reads; falls back to a 2-pass join
#' if a single 12-col read FPEs.
#'
#' @section Reads:
#'   - local v3 spatial parquet (override with `BANC_V3_SPATIAL_PARQUET`)
#'
#' @section Writes:
#'   - `$BANC_PROJECT_DATA_DIR/completion/banc_<ver>_v3_synapse_sample_<DATE>.csv`
#'
#' @section Invoked by:
#'   `o2/oneshots/o2_banc_v3_synapse_sample.sh`.

###############################################################################
### Standalone v3 synapse sample producer.
###
### Standalone replacement for the v3 sample block inside
### banc/metrics/banc-calculate-completion.R. That block has bombed at v888
### with the SLURM-only arrow C++ URI bug ("Invalid: Unrecognized filesystem
### type in URI: file:///_") when reading the 15-column 15 GB local spatial
### parquet. The on-disk Apr 24 v3 sample is the LAST successful write, and
### its format is from before the May 5 fixes:
###   - Coordinate 1 is in nm, not BANC raw voxel.
###   - Coordinate 2 is empty.
###   - Type is "Point" (should be "Line").
###   - ID is empty (should be syn_id).
###
### This script regenerates the file with the correct format using narrow
### parquet reads that the SLURM-arrow bug doesn't trip on. Falls back to a
### 2-pass-and-join read if a single 12-col read still FPEs.
###
### The spatial info (region/side/neuropil/coordinates/scores) is keyed by
### syn_id, which is stable across BANC root_id versions — so the v888
### spatial parquet's sample IS valid for v890 reads of the same syn_ids.
###
### Usage:
###   Rscript banc/metrics/banc-v3-synapse-sample.R
###
### Output:
###   $BANC_PROJECT_DATA_DIR/completion/banc_<ver>_v3_synapse_sample_<DATE>.csv
###
### Optional env vars:
###   BANC_V3_SPATIAL_PARQUET — override path to local v3 spatial parquet
###     (default: <banc.synapses.v3.save.path>/banc_<ver>_synapses_v3.parquet,
###     falling back to v888 / v850 if banc.version's file doesn't exist).
###   BANC_V3_SAMPLE_VERSION — version label used in the output filename
###     (default: banc.version).
###############################################################################
source("banc/banc-startup.R")

local({

message("### banc v3 synapse sample (standalone) ###")
t_start <- Sys.time()

###########################
### Resolve spatial parquet path
###########################

# banc.synapses.v3.save.path is defined locally inside completion.R; not
# exported by banc-startup.R. Recompute here from the canonical save.path.
banc.synapses.v3.save.path <- file.path(banc.save.path, "synapses_v3")

.pick_spatial <- function() {
  override <- Sys.getenv("BANC_V3_SPATIAL_PARQUET", unset = "")
  if (nzchar(override) && file.exists(override)) return(override)
  for (v in unique(c(banc.version, "888", "850"))) {
    p <- file.path(banc.synapses.v3.save.path,
                   sprintf("banc_%s_synapses_v3.parquet", v))
    if (file.exists(p)) return(p)
  }
  stop("No v3 spatial parquet found for any of banc.version/888/850.",
       call. = FALSE)
}
v3_spatial_local <- .pick_spatial()
sample_version   <- Sys.getenv("BANC_V3_SAMPLE_VERSION",
                               unset = banc.version)
message(sprintf("v3 spatial parquet: %s", v3_spatial_local))
message(sprintf("sample version label: %s", sample_version))

###########################
### Narrow read(s)      ###
###########################

# Goal columns from the spatial parquet. Drop median_score + X/Y/Z (centroid)
# vs the failed 15-col read in completion.R — we only need score for binning
# and presyn/postsyn for Coord 1/2.
cols_meta    <- c("syn_id", "size", "region", "side", "neuropil",
                  "mean_score")
cols_spatial <- c("syn_id", "presyn_x", "presyn_y", "presyn_z",
                  "postsyn_x", "postsyn_y", "postsyn_z")

.do_two_pass <- function() {
  # The 12-col one-pass read OOMs at 48G on 39969154 (SIGKILL bypasses
  # tryCatch); go straight to 2-pass narrow reads + join on syn_id. Each
  # pass peaks ~5-8GB; total well under any reasonable mem-per-cpu budget.
  message("  reading meta cols (6)...")
  meta <- arrow::read_parquet(v3_spatial_local,
                              col_select = dplyr::all_of(cols_meta))
  message(sprintf("  meta: %s rows", format(nrow(meta), big.mark = ",")))
  message("  reading spatial cols (7)...")
  spat <- arrow::read_parquet(v3_spatial_local,
                              col_select = dplyr::all_of(cols_spatial))
  message(sprintf("  spatial: %s rows", format(nrow(spat), big.mark = ",")))
  message("  joining on syn_id...")
  meta %>% dplyr::inner_join(spat, by = "syn_id")
}

spatial <- .do_two_pass()
spatial <- spatial %>% dplyr::distinct(syn_id, .keep_all = TRUE)
message(sprintf("  spatial rows: %s",
                format(nrow(spatial), big.mark = ",")))

###########################
### Filter + sample     ###
###########################

# Drop rows lacking the bucket key or score; collapse comma-joined neuropil
# strings (a single v3 synapse can intersect multiple neuropils) to a single
# "primary" label.
#
# VNC sampling uses ONLY COURT_vnc_* as bucket keys — MANC_vnc_* duplicates
# the same anatomy under a less-canonical naming scheme and was inflating the
# sample (~800 extra rows). Synapses in region=="vnc" that lack any
# COURT_vnc_* label in their neuropil string are dropped from sampling.
# Non-VNC regions fall back to first-comma-entry, which already guarantees
# one neuropil per synapse.
#
# `sez` region is folded into `brain` (only ~2 rows total in v888).
spatial <- spatial %>%
  dplyr::filter(!is.na(side), nzchar(side),
                !is.na(neuropil), nzchar(neuropil),
                !is.na(mean_score)) %>%
  dplyr::mutate(
    region = ifelse(region == "sez", "brain", region),
    has_court_vnc = grepl("COURT_vnc_", neuropil, fixed = TRUE),
    np_primary = ifelse(
      has_court_vnc,
      # extract the COURT_vnc_* segment regardless of its position in the
      # comma-list; lazy `.*?` consumes any leading "MANC_xxx," etc.
      sub(".*?(COURT_vnc_[^,]*).*", "\\1", neuropil),
      # no COURT_vnc_ in this string — take first comma-separated entry
      sub(",.*", "", neuropil)
    ),
    # Some neuropil names are slash-joined (e.g. "ITO_optic_ME/LO_L/R" or
    # similar multi-region labels). Collapse those to the first segment too.
    np_primary = sub("/.*", "", np_primary)
  )

# VNC: require COURT_vnc_* label. Drops MANC_vnc_*-only synapses, which
# inflate the sample with redundant anatomy.
n_pre_vnc_filter <- nrow(spatial)
spatial <- spatial %>%
  dplyr::filter(region != "vnc" | has_court_vnc) %>%
  dplyr::select(-has_court_vnc)
message(sprintf("  dropped %s vnc rows without COURT_vnc_* label",
                format(n_pre_vnc_filter - nrow(spatial), big.mark = ",")))
message(sprintf("  unique np_primary values: %d (before collapse it was %d)",
                length(unique(spatial$np_primary)),
                length(unique(spatial$neuropil))))

# Target: max 5-6k sampled syns total. With ~80-100 (np_primary x side)
# groups after the COURT_vnc_* collapse + MANC_vnc_* drop, 20 per group ×
# 10 bins / 10 = 2 per bin = 20 per group → ~2-3.5k base sample. The
# per-region quotas (optic_lobes=1000, outside=200) top up to ~4-5k total,
# under the V3_SAMPLE_TOTAL_CAP downsample step below.
V3_SAMPLE_BINS      <- 10L
V3_SAMPLE_PER_GROUP <- 20L
V3_SAMPLE_TOTAL_CAP <- 6000L

# Quantile bin edges (mean_score is right-skewed; equal-width wastes bins).
qbreaks <- stats::quantile(spatial$mean_score,
                           probs = seq(0, 1, length.out = V3_SAMPLE_BINS + 1),
                           na.rm = TRUE,
                           names = FALSE)
qbreaks <- unique(qbreaks)
if (length(qbreaks) < 3L) {
  message("  too few unique mean_score values for quantile bins; using one bin")
  qbreaks <- range(spatial$mean_score, na.rm = TRUE)
}

set.seed(42)
sampled <- spatial %>%
  dplyr::mutate(score_bin = cut(mean_score, breaks = qbreaks,
                                include.lowest = TRUE, labels = FALSE)) %>%
  dplyr::group_by(np_primary, side, score_bin) %>%
  dplyr::slice_sample(n = ceiling(V3_SAMPLE_PER_GROUP / V3_SAMPLE_BINS),
                      replace = FALSE) %>%
  dplyr::ungroup() %>%
  as.data.frame(stringsAsFactors = FALSE)

n_groups <- length(unique(paste(sampled$np_primary, sampled$side)))
message(sprintf("  stratified sample: %s syns from %d (np_primary x side) groups",
                format(nrow(sampled), big.mark = ","), n_groups))

###########################
### Per-region minimum quotas
###########################
# Stratifying by (np_primary x side) under-represents regions with few unique
# np_primary values: optic_lobes has only 4-6 large neuropils so the natural
# cap is ~250 rows; outside has even fewer. Top them up with random draws
# from the unsampled remainder of each region. Score-bin stratification is
# NOT preserved in the boost (the goal is coverage, not balanced binning).
spatial$score_bin <- cut(spatial$mean_score, breaks = qbreaks,
                         include.lowest = TRUE, labels = FALSE)
REGION_QUOTA <- list(optic_lobes = 1000L, outside = 200L)
for (rgn in names(REGION_QUOTA)) {
  current <- sum(sampled$region == rgn, na.rm = TRUE)
  needed  <- REGION_QUOTA[[rgn]] - current
  if (needed > 0L) {
    pool <- spatial %>%
      dplyr::filter(region == rgn,
                    !syn_id %in% sampled$syn_id)
    extra <- pool %>%
      dplyr::slice_sample(n = min(needed, nrow(pool)), replace = FALSE) %>%
      as.data.frame(stringsAsFactors = FALSE)
    sampled <- dplyr::bind_rows(sampled, extra)
    message(sprintf("  region %s: topped up +%d rows (was %d, now %d / quota %d)",
                    rgn, nrow(extra), current, current + nrow(extra),
                    REGION_QUOTA[[rgn]]))
  } else {
    message(sprintf("  region %s: already at %d rows (quota %d) — no top-up",
                    rgn, current, REGION_QUOTA[[rgn]]))
  }
}
message(sprintf("  TOTAL after quotas: %s syns",
                format(nrow(sampled), big.mark = ",")))

# Hard cap: if total exceeds V3_SAMPLE_TOTAL_CAP, uniformly downsample while
# preserving the per-region quota floors. Picks the over-allocated regions
# (above their quota / unquoted) and trims them proportionally.
if (nrow(sampled) > V3_SAMPLE_TOTAL_CAP) {
  excess <- nrow(sampled) - V3_SAMPLE_TOTAL_CAP
  message(sprintf("  capping: %s rows over %d cap, downsampling...",
                  format(excess, big.mark = ","), V3_SAMPLE_TOTAL_CAP))
  # Per-region quota floor (0 for regions with no quota — they can be
  # trimmed down to 0 if needed, but in practice each non-quoted region
  # has at most a few hundred rows).
  region_floor <- c(unlist(REGION_QUOTA))
  region_floor <- setNames(rep(0L, length(unique(sampled$region))),
                           unique(sampled$region))
  for (rgn in names(REGION_QUOTA))
    region_floor[rgn] <- REGION_QUOTA[[rgn]]
  # Trimmable count per region = current - floor. Drop `excess` rows in
  # proportion to each region's trimmable surplus.
  region_counts <- table(sampled$region)
  surplus <- pmax(0L, as.integer(region_counts) - region_floor[names(region_counts)])
  names(surplus) <- names(region_counts)
  total_surplus <- sum(surplus)
  if (total_surplus < excess) {
    message(sprintf("  WARNING: trimmable surplus (%d) < excess (%d) — quotas will eat into floor",
                    total_surplus, excess))
  }
  set.seed(43)
  for (rgn in names(surplus)) {
    if (surplus[[rgn]] <= 0L) next
    to_drop <- ceiling(excess * surplus[[rgn]] / max(1L, total_surplus))
    to_drop <- min(to_drop, surplus[[rgn]])
    if (to_drop <= 0L) next
    rgn_idx <- which(sampled$region == rgn)
    drop_idx <- sample(rgn_idx, size = to_drop, replace = FALSE)
    sampled <- sampled[-drop_idx, , drop = FALSE]
  }
  message(sprintf("  after cap: %s syns",
                  format(nrow(sampled), big.mark = ",")))
}

###########################
### Coordinate transform — V3 voxel -> nm -> BANC raw voxel
###########################
# v3 detections live in V3 voxel space (16 x 16 x 45 nm). To get them into
# the BANC raw-voxel coordinate space used by neuroglancer:
#   (i)  multiply V3 voxel coords by V3_VOXEL_NM to get nm
#   (ii) bancr::banc_nm2raw to convert nm to BANC raw voxel (4 x 4 x 45 nm)
V3_VOXEL_NM <- c(16, 16, 45)
post_nm <- cbind(sampled$postsyn_x * V3_VOXEL_NM[1],
                 sampled$postsyn_y * V3_VOXEL_NM[2],
                 sampled$postsyn_z * V3_VOXEL_NM[3])
pre_nm  <- cbind(sampled$presyn_x  * V3_VOXEL_NM[1],
                 sampled$presyn_y  * V3_VOXEL_NM[2],
                 sampled$presyn_z  * V3_VOXEL_NM[3])
coord1 <- bancr::banc_nm2raw(post_nm)   # Coord 1 = postsynaptic
coord2 <- bancr::banc_nm2raw(pre_nm)    # Coord 2 = presynaptic

sampled$`Coordinate 1` <- sprintf("(%.0f, %.0f, %.0f)",
                                  coord1[, 1], coord1[, 2], coord1[, 3])
sampled$`Coordinate 2` <- sprintf("(%.0f, %.0f, %.0f)",
                                  coord2[, 1], coord2[, 2], coord2[, 3])

###########################
### Neuroglancer-format CSV
###########################
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
  mean_score             = round(sampled$mean_score, 4),
  score_bin              = sampled$score_bin,
  region                 = sampled$region,
  side                   = sampled$side,
  neuropil               = sampled$np_primary,
  stringsAsFactors       = FALSE
)
colnames(ngl_ann) <- gsub("\\.", " ", colnames(ngl_ann))

###########################
### Write               ###
###########################
out_dir <- file.path(
  Sys.getenv("BANC_PROJECT_DATA_DIR",
             unset = "/n/data1/hms/neurobio/wilson/banc/BANC-project/data"),
  "completion")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
sample_date <- format(Sys.time(), "%Y-%m-%d")
ngl_file <- file.path(out_dir,
                      sprintf("banc_%s_v3_synapse_sample_%s.csv",
                              sample_version, sample_date))
readr::write_excel_csv(ngl_ann, ngl_file)
message(sprintf("saved %s (%s rows)", basename(ngl_file),
                format(nrow(ngl_ann), big.mark = ",")))

message(sprintf("### banc v3 synapse sample done [%s] ###",
                format(round(difftime(Sys.time(), t_start, units = "mins"), 1))))

})
