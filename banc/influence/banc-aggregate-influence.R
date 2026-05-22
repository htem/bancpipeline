#' banc-aggregate-influence — Aggregate the all-to-all influence by sensory and effector sub-class.
#'
#' Reads the sharded parquet chunks from `banc-build-influence.R` and
#' produces two aggregated parquets: (1) sensory sub_class → all neurons,
#' (2) all neurons → effector sub_class. Then pushes both to GCS.
#'
#' @section Reads:
#'   - `<banc.versioned.save.path>/influence/all_to_all/chunk_*.parquet`
#'   - `<banc.versioned.save.path>/banc_<ver>_meta.feather`
#'
#' @section Writes:
#'   - `<banc.versioned.save.path>/influence_sensory_subclass_to_all.parquet`
#'   - `<banc.versioned.save.path>/influence_all_to_effector_subclass.parquet`
#'   - GCS mirror of both
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_aggregate_influence.sh`,
#'   `o2/production/o2_banc_aggregate_influence_priority.sh`
#'
#' @section Schema:
#'   BANC-project/manuscript/print/dataverse/documentation/influence_all_to_effector_subclass.md
#'   BANC-project/manuscript/print/dataverse/documentation/influence_sensory_subclass_to_all.md
#'
#' @section Paper:
#'   Methods §"Influence".

###############################################
### Aggregate all-to-all influence by sensory
### and effector cell_sub_class.
###
### Reads the 309 sharded parquet chunks from
### banc-build-influence.R and produces two
### parquet files:
###   1. sensory sub_class -> all neurons
###   2. all neurons -> effector sub_class
###
### Pushes outputs to GCS.
###############################################
source("banc/banc-startup.R")
library(data.table)

local({

  message("### banc: aggregate influence by subclass ###")
  t_start <- Sys.time()

  save.path  <- file.path("/n/data1/hms/neurobio/wilson/connectomes/banc",
                          paste0("banc_", banc.version))
  chunk_dir  <- file.path(save.path, "influence", "all_to_all")
  output_dir <- file.path(save.path, "influence")

  meta <- arrow::read_feather(file.path(save.path,
                                        sprintf("banc_%s_meta.feather", banc.version)))

  # --- Sensory lookup: super_class contains "sensory" OR body_part_sensory defined ---
  is_sensory <- grepl("sensory", meta$super_class, ignore.case = TRUE) |
                (!is.na(meta$body_part_sensory) &
                 meta$body_part_sensory != "" &
                 meta$body_part_sensory != "unknown")
  is_sensory[is.na(is_sensory)] <- FALSE
  sensory <- meta[is_sensory & !is.na(meta$cell_sub_class), ]
  sensory_lookup <- setNames(sensory$cell_sub_class,
                             sensory[[paste0("banc_", banc.version, "_id")]])

  # --- Effector lookup: flow == "efferent" ---
  is_effector <- !is.na(meta$flow) & meta$flow == "efferent"
  effector <- meta[is_effector & !is.na(meta$cell_sub_class), ]
  effector_lookup <- setNames(effector$cell_sub_class,
                              effector[[paste0("banc_", banc.version, "_id")]])

  message(sprintf("  Sensory:  %d neurons in %d sub_classes (dropped %d with NA sub_class)",
                  length(sensory_lookup), length(unique(sensory_lookup)),
                  sum(is_sensory) - length(sensory_lookup)))
  message(sprintf("  Effector: %d neurons in %d sub_classes (dropped %d with NA sub_class)",
                  length(effector_lookup), length(unique(effector_lookup)),
                  sum(is_effector) - length(effector_lookup)))

  # --- Chunk files ---
  chunk_files <- sort(list.files(chunk_dir,
                                 pattern = "^chunk_[0-9]+\\.parquet$",
                                 full.names = TRUE))
  message(sprintf("  %d chunk files in %s", length(chunk_files), chunk_dir))

  # --- Accumulators ---
  sensory_acc  <- data.table(source = character(0),
                             target = character(0),
                             influence = numeric(0))
  effector_acc <- data.table(source = character(0),
                             target = character(0),
                             influence = numeric(0))

  # --- Main loop ---
  for (i in seq_along(chunk_files)) {
    chunk <- arrow::read_parquet(chunk_files[i])
    setDT(chunk)

    chunk <- chunk[grepl("^[0-9]+$", upstream_id)]

    # Sensory direction
    s_idx <- chunk$upstream_id %in% names(sensory_lookup)
    if (any(s_idx)) {
      s_chunk <- chunk[s_idx]
      s_chunk[, source := sensory_lookup[upstream_id]]
      partial_s <- s_chunk[, .(influence = sum(raw_influence, na.rm = TRUE)),
                           by = .(source, downstream_id)]
      setnames(partial_s, "downstream_id", "target")
      sensory_acc <- rbindlist(list(sensory_acc, partial_s))[
        , .(influence = sum(influence)), by = .(source, target)]
    }

    # Effector direction
    e_idx <- chunk$downstream_id %in% names(effector_lookup)
    if (any(e_idx)) {
      e_chunk <- chunk[e_idx]
      e_chunk[, target := effector_lookup[downstream_id]]
      partial_e <- e_chunk[, .(influence = sum(raw_influence, na.rm = TRUE)),
                           by = .(upstream_id, target)]
      setnames(partial_e, "upstream_id", "source")
      effector_acc <- rbindlist(list(effector_acc, partial_e))[
        , .(influence = sum(influence)), by = .(source, target)]
    }

    rm(chunk)
    if (i %% 50 == 0) gc(verbose = FALSE)

    if (i %% 10 == 0 || i == length(chunk_files)) {
      elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
      eta <- elapsed / i * (length(chunk_files) - i)
      message(sprintf("  Chunk %d/%d (%.0f%%) | sensory: %s rows | effector: %s rows | %.1fm elapsed | ETA %.1fm",
                      i, length(chunk_files), 100 * i / length(chunk_files),
                      format(nrow(sensory_acc), big.mark = ","),
                      format(nrow(effector_acc), big.mark = ","),
                      elapsed, eta))
    }
  }

  # --- Filter zeros ---
  sensory_acc  <- sensory_acc[influence > 0]
  effector_acc <- effector_acc[influence > 0]

  message(sprintf("  Final sensory:  %s rows (%d sub_classes)",
                  format(nrow(sensory_acc), big.mark = ","),
                  length(unique(sensory_acc$source))))
  message(sprintf("  Final effector: %s rows (%d sub_classes)",
                  format(nrow(effector_acc), big.mark = ","),
                  length(unique(effector_acc$target))))

  # --- Write ---
  sensory_file  <- file.path(output_dir, "influence_sensory_subclass_to_all.parquet")
  effector_file <- file.path(output_dir, "influence_all_to_effector_subclass.parquet")

  write_connectome_data(as.data.frame(sensory_acc), sensory_file, format = "parquet")
  write_connectome_data(as.data.frame(effector_acc), effector_file, format = "parquet")

  message(sprintf("  Wrote %s", sensory_file))
  message(sprintf("  Wrote %s", effector_file))

  # --- GCS sync ---
  gs_base <- sprintf("gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/banc_%s/influence/",
                     banc.version)
  for (f in c(sensory_file, effector_file)) {
    cmd <- sprintf("gsutil cp %s %s", f, gs_base)
    status <- system(cmd)
    if (status == 0) {
      message(sprintf("  Uploaded %s -> %s", basename(f), gs_base))
    } else {
      warning(sprintf("  Failed to upload %s (exit %d)", basename(f), status))
    }
  }

  message(sprintf("### banc: influence aggregation complete (%.1f min) ###",
                  as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
})
