#' banc-test — Push 10 neurons through the full pipeline as a smoke test.
#'
#' Selects 10 test neurons (preferring proofread rows with super_class that
#' have not yet been processed) and runs them through phases 1-5: metrics,
#' morphology, cross-dataset NBLAST, NBLAST images, SeaTable push. Does
#' NOT run `banc-ids.R`, `banc-delete.R`, `banc-nblast-compile.R`.
#'
#' @section Reads:
#'   - SeaTable `banc_meta` (to select candidates)
#'
#' @section Writes:
#'   - Per-neuron outputs under `banc.l2split.save.path/`, `banc.swc.save.path/`,
#'     `banc.metrics.save.path/`, `banc.nblast_*` trees (small subset)
#'   - SeaTable `banc_meta` (NT predictions, metrics for the 10 test neurons)
#'
#' @section Notes:
#'   - Mutates SeaTable — do not run against production unless you intend
#'     the side effects.

###########################################################
### Test harness: push 10 neurons through the pipeline
###
### Selects 10 test neurons (preferring proofread, with
### super_class, not yet processed) and runs them through:
###
### Phase 1 (metrics):
###   L2 download → root positions → regions → synapses →
###   L2 metrics → volumes → NT predictions
###
### Phase 2 (morphology):
###   Detailed skeletons → flow centrality splits
###
### Phase 3 (NBLAST — all datasets):
###   FAFB, FANC, MANC, hemibrain, maleCNS, left-right
###
### Phase 4 (NBLAST images — all datasets):
###   FAFB, FANC, MANC, hemibrain, maleCNS, left-right
###
### Phase 5 (SeaTable push):
###   Join per-metric feather files → push to seatable
###   Push NT predictions to seatable
###
### Does NOT run: banc-ids.R, banc-delete.R,
### banc-nblast-compile.R.
###
### Usage: Rscript banc/banc-test.R
###########################################################
source("banc/banc-startup.R")

bancr::choose_banc()

message("========================================")
message("=== BANC PIPELINE TEST MODE          ===")
message("========================================")
t_start <- Sys.time()

###########################
### Select test neurons ###
###########################

bc <- banctable_query("SELECT _id, status, proofread, super_class, root_id, supervoxel_id, position from banc_meta") %>%
  banc_filter_neurons() %>%
  dplyr::filter(!is.na(root_id))

# Read existing metrics to find unprocessed neurons
synapses_file <- file.path(banc.save.path, "banc_detections.feather")
if (file.exists(synapses_file)) {
  processed <- arrow::read_feather(synapses_file) %>%
    dplyr::filter(!is.na(input_connections)) %>%
    dplyr::pull(root_id)
} else {
  processed <- character(0)
}

# Score neurons: proofread=+2, has super_class=+1, unprocessed=+1
candidates <- bc %>%
  dplyr::mutate(
    score = as.integer(!is.na(proofread) & proofread == "t") * 2 +
            as.integer(!is.na(super_class) & super_class != "") +
            as.integer(!root_id %in% processed)
  ) %>%
  dplyr::arrange(dplyr::desc(score))

# Take 10 from top-scoring candidates
set.seed(42)
banc.test.ids <<- candidates %>%
  dplyr::slice_head(n = min(50, nrow(candidates))) %>%
  dplyr::slice_sample(n = min(10, nrow(.))) %>%
  dplyr::pull(root_id) %>%
  as.character()

message(sprintf("Selected %d test neurons:", length(banc.test.ids)))
message(paste(" ", banc.test.ids, collapse = "\n"))

# Helper: source with tryCatch
safe_source <- function(script, label) {
  tryCatch(
    source(script),
    error = function(e) message(sprintf("  %s failed: %s", label, e$message))
  )
}

###########################
### Phase 1: Metrics    ###
###########################

n_step <- 0
n_total <- 23

message("\n========================================")
message("=== Phase 1: Core metrics            ===")
message("========================================")

n_step <- n_step + 1
message(sprintf("\n--- Step %d/%d: L2 skeleton download ---", n_step, n_total))
safe_source("banc/metrics/banc-l2.R", "L2 download")

n_step <- n_step + 1
message(sprintf("\n--- Step %d/%d: Root positions ---", n_step, n_total))
safe_source("banc/metrics/banc-calculate-root-positions.R", "Root positions")

n_step <- n_step + 1
message(sprintf("\n--- Step %d/%d: Regions ---", n_step, n_total))
safe_source("banc/metrics/banc-calculate-regions.R", "Regions")

n_step <- n_step + 1
message(sprintf("\n--- Step %d/%d: Synapses ---", n_step, n_total))
safe_source("banc/metrics/banc-calculate-synapses.R", "Synapses")

n_step <- n_step + 1
message(sprintf("\n--- Step %d/%d: L2 metrics ---", n_step, n_total))
safe_source("banc/metrics/banc-calculate-l2-metrics.R", "L2 metrics")

n_step <- n_step + 1
message(sprintf("\n--- Step %d/%d: Volumes ---", n_step, n_total))
safe_source("banc/metrics/banc-calculate-volumes.R", "Volumes")

n_step <- n_step + 1
#message(sprintf("\n--- Step %d/%d: NT predictions ---", n_step, n_total))
#safe_source("banc/metrics/banc-calculate-ntpred.R", "NT predictions")

###########################
### Phase 2: Morphology ###
###########################

message("\n========================================")
message("=== Phase 2: Morphology              ===")
message("========================================")

n_step <- n_step + 1
message(sprintf("\n--- Step %d/%d: Detailed skeletons ---", n_step, n_total))
safe_source("banc/metrics/banc-calculate-skeletons.R", "Skeletons")

n_step <- n_step + 1
message(sprintf("\n--- Step %d/%d: Flow centrality splits ---", n_step, n_total))
safe_source("banc/metrics/banc-calculate-split.R", "Splits")

###########################
### Phase 3: NBLAST     ###
###########################

message("\n========================================")
message("=== Phase 3: NBLAST (all datasets)   ===")
message("========================================")

nblast_scripts <- list(
  list(script = "banc/nblast/banc-fafb-nblast.R",     label = "FAFB NBLAST"),
  list(script = "banc/nblast/banc-fanc-nblast.R",      label = "FANC NBLAST"),
  list(script = "banc/nblast/banc-manc-nblast.R",      label = "MANC NBLAST"),
  list(script = "banc/nblast/banc-hemibrain-nblast.R",  label = "Hemibrain NBLAST"),
  list(script = "banc/nblast/banc-malecns-nblast.R",   label = "maleCNS NBLAST")
  #list(script = "banc/nblast/banc-nblast-lr.R",        label = "Left-right NBLAST")
)

for (ns in nblast_scripts) {
  n_step <- n_step + 1
  message(sprintf("\n--- Step %d/%d: %s ---", n_step, n_total, ns$label))
  safe_source(ns$script, ns$label)
}

###########################
### Phase 4: NBLAST img ###
###########################

message("\n========================================")
message("=== Phase 4: NBLAST images           ===")
message("========================================")

nblast_image_scripts <- list(
  list(script = "banc/nblast/banc-fafb-nblast-images.R",    label = "FAFB images"),
  list(script = "banc/nblast/banc-fanc-nblast-images.R",     label = "FANC images"),
  list(script = "banc/nblast/banc-manc-nblast-images.R",     label = "MANC images"),
  list(script = "banc/nblast/banc-hemibrain-nblast-images.R", label = "Hemibrain images"),
  list(script = "banc/nblast/banc-malecns-nblast-images.R",  label = "maleCNS images")
  #list(script = "banc/nblast/banc-lr-nblast-images.R",       label = "Left-right images")
)

for (ns in nblast_image_scripts) {
  n_step <- n_step + 1
  message(sprintf("\n--- Step %d/%d: %s ---", n_step, n_total, ns$label))
  safe_source(ns$script, ns$label)
}

###########################
### Phase 5: SeaTable   ###
###########################

message("\n========================================")
message("=== Phase 5: SeaTable push           ===")
message("========================================")

n_step <- n_step + 1
message(sprintf("\n--- Step %d/%d: Push metrics to seatable ---", n_step, n_total))
safe_source("banc/update/banc-update-metrics.R", "SeaTable push")

n_step <- n_step + 1
#message(sprintf("\n--- Step %d/%d: Push NT predictions to seatable ---", n_step, n_total))
#safe_source("banc/update/banc-update-ntpred.R", "NT prediction push")

###########################
### Report results      ###
###########################

message("\n========================================")
message("=== TEST RESULTS                     ===")
message("========================================")

# Phase 1: Per-metric feather files
combined_file <- file.path(banc.save.path, "banc_metrics.feather")
if (file.exists(combined_file)) {
  final <- arrow::read_feather(combined_file) %>%
    dplyr::filter(root_id %in% banc.test.ids) %>%
    dplyr::select(dplyr::any_of(c("root_id", "root_position_nm", "side", "region",
                  "l2_nodes", "l2_cable_length_um",
                  "input_connections", "output_connections",
                  "volume_nm3")))
  message(sprintf("\n[Phase 1] Metrics populated for %d/%d test neurons", nrow(final), length(banc.test.ids)))
  print(as.data.frame(final))

  # Check completeness
  missing <- banc.test.ids[!banc.test.ids %in% final$root_id]
  if (length(missing)) {
    message(sprintf("\nMissing from combined feather: %s", paste(missing, collapse = ", ")))
  }
  na_cols <- colSums(is.na(final[, -1]))
  if (any(na_cols > 0)) {
    message("\nColumns with NA values:")
    for (col in names(na_cols[na_cols > 0])) {
      message(sprintf("  %s: %d/%d missing", col, na_cols[col], nrow(final)))
    }
  }
} else {
  message("WARNING: banc_metrics.feather not found — no combined metrics to report")
  message("  (Individual feather files may still exist — check update-metrics.R)")
}

# Phase 2: Skeletons and splits
message("\n[Phase 2] Morphology outputs:")
swc_files <- list.files(banc.swc.save.path, pattern = "\\.swc$")
swc_ids <- gsub("\\.swc$", "", swc_files)
swc_ids <- gsub("datahmneurobiwilsobanob", "", swc_ids)
n_skels <- sum(banc.test.ids %in% swc_ids)
message(sprintf("  Detailed skeletons: %d/%d test neurons", n_skels, length(banc.test.ids)))

split_swcs <- list.files(file.path(banc.l2split.save.path, "swc"), pattern = "\\.swc$")
split_ids <- gsub("\\.swc$", "", split_swcs)
n_splits <- sum(banc.test.ids %in% split_ids)
message(sprintf("  Flow centrality splits: %d/%d test neurons", n_splits, length(banc.test.ids)))

split_metrics <- list.files(file.path(banc.l2split.save.path, "metrics"), pattern = "\\.csv$")
message(sprintf("  Split metrics CSVs: %d files", length(split_metrics)))

split_images <- list.files(file.path(banc.l2split.save.path, "images"), pattern = "\\.png$")
split_img_ids <- gsub("\\.png$", "", split_images)
n_split_imgs <- sum(banc.test.ids %in% split_img_ids)
message(sprintf("  Split comparison images: %d/%d test neurons", n_split_imgs, length(banc.test.ids)))

# Phase 3 & 4: NBLAST results and images
message("\n[Phase 3 & 4] NBLAST outputs:")

nblast_datasets <- list(
  list(name = "FAFB",      save_path = banc.nblast.fafb.save.path,      version = banc.nblast.version),
  list(name = "FANC",      save_path = banc.nblast.fanc.save.path,      version = banc.nblast.version),
  list(name = "MANC",      save_path = banc.nblast.manc.save.path,      version = banc.nblast.version),
  list(name = "Hemibrain", save_path = banc.nblast.hemibrain.save.path,  version = banc.nblast.version),
  list(name = "maleCNS",   save_path = banc.nblast.malecns.save.path,   version = banc.nblast.malecns.version),
  list(name = "Mirror",    save_path = banc.nblast.mirror.save.path,    version = NULL)
)

for (ds in nblast_datasets) {
  # Results
  if (!is.null(ds$version)) {
    nblast_dir <- file.path(ds$save_path, "results", ds$version)
  } else {
    nblast_dir <- file.path(ds$save_path, "results")
  }
  if (dir.exists(nblast_dir)) {
    nblast_files <- list.files(nblast_dir, pattern = "\\.csv$")
    nblast_rids <- gsub(".*_root_id_|\\.csv", "", nblast_files)
    n_nblast <- sum(banc.test.ids %in% nblast_rids)
    message(sprintf("  %s NBLAST results: %d/%d", ds$name, n_nblast, length(banc.test.ids)))
  } else {
    message(sprintf("  %s NBLAST results dir not found", ds$name))
  }

  # Images
  if (!is.null(ds$version)) {
    img_dir <- file.path(ds$save_path, "images", ds$version, "todo")
  } else {
    img_dir <- file.path(ds$save_path, "images", "todo")
  }
  if (dir.exists(img_dir)) {
    imgs <- list.files(img_dir, pattern = "\\.png$", recursive = TRUE)
    if (length(imgs)) {
      img_rids <- unique(regmatches(imgs, regexpr("(?<=root_id_)\\d+", imgs, perl = TRUE)))
      n_imgs <- sum(banc.test.ids %in% img_rids)
      message(sprintf("  %s NBLAST images: %d/%d", ds$name, n_imgs, length(banc.test.ids)))
    } else {
      message(sprintf("  %s NBLAST images: 0", ds$name))
    }
  }
}

message(sprintf("\n=== Test complete [%s] ===",
                format(round(difftime(Sys.time(), t_start, units = "mins"), 1))))

# Clean up
rm(banc.test.ids, envir = .GlobalEnv)
