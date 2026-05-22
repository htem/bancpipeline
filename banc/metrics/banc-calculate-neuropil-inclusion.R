#' banc-calculate-neuropil-inclusion — Assign neuropil / region / side to each synapse.
#'
#' Join by stable synapse `id` against `synapse_neuropil_lookup.parquet`
#' first; classify residuals via alpha shapes (pre-computed, fork-parallel).
#'
#' @section Reads:
#'   - `banc_<ver>_synapses_<src>.parquet`, `synapse_neuropil_lookup.parquet`
#'
#' @section Writes:
#'   - `banc_<ver>_synapses_<src>_neuropils.parquet`
#'   - `synapse_neuropil_lookup.parquet` (extended)
#'
#' @section Invoked by:
#'   production v888 rebuild chain.
#'
#' @section Notes:
#'   - Resumable: chunked into temp parquet files.

###############################################
### Calculate neuropil inclusion for BANC   ###
### synapses                                ###
###############################################
### Reads versioned synapse parquet from
### banc-calculate-connectivity.R, classifies
### each synapse into neuropil/region/side.
###
### Optimization 1: loads a version-agnostic
### synapse_neuropil_lookup.parquet and joins
### by synapse `id` first. Only synapses with
### no lookup match are classified via alpha
### shapes. This avoids recomputing stable
### neuropil assignments across version rebuilds.
###
### Optimization 2: pre-computes all alpha shapes
### ONCE before chunking (~5-10 min), then
### uses fork-based parallel::mclapply so
### workers inherit via copy-on-write (zero
### serialization overhead).
###
### Output:
###   banc_{version}_synapses_v2_neuropils.parquet
###   synapse_neuropil_lookup.parquet (updated)
###
### Resumable: processed chunks saved as temp
### parquet files.
###############################################
source("banc/banc-startup.R")

local({

message("### banc: calculating neuropil inclusion ###")
t_total_start <- Sys.time()

#######################
### CONFIGURATION   ###
#######################

# banc.version set in banc-startup.R

# Processing parameters
chunk_size <- 500000   # rows per chunk (500K balances memory vs overhead)
n_cores <- max(1, parallel::detectCores() - 1)
alpha_value <- 50000   # alpha shape parameter

# Input/output paths
# banc-calculate-connectivity.R --source v2 writes `banc_<ver>_synapses_v2.parquet`
# (the `_v2` suffix was added during the v3-source migration). This was
# previously `banc_<ver>_synapses.parquet`; updated 2026-04-21 after the v888
# rebuild crashed on the missing unsuffixed file.
input_parquet <- file.path(banc.connectivity.save.path,
                           sprintf("banc_%s_synapses_v2.parquet", banc.version))
output_parquet <- file.path(banc.connectivity.save.path,
                            sprintf("banc_%s_synapses_v2_neuropils.parquet", banc.version))
temp_dir <- file.path(banc.connectivity.save.path,
                      sprintf("neuropil_chunks_%s", banc.version))
lookup_file <- file.path(banc.connectivity.save.path, "synapse_neuropil_lookup.parquet")

stopifnot(file.exists(input_parquet))
dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

message(sprintf("Input: %s", input_parquet))
message(sprintf("Output: %s", output_parquet))
message(sprintf("Cores: %d, Chunk size: %s", n_cores, format(chunk_size, big.mark = ",")))

##########################################
### READ SYNAPSES & JOIN LOOKUP        ###
##########################################

message("Reading synapse parquet...")
synapses <- arrow::read_parquet(input_parquet)
synapses$id <- as.character(synapses$id)
n_total <- nrow(synapses)
message(sprintf("  Total synapses: %s", format(n_total, big.mark = ",")))

# Try to load existing neuropil lookup
matched <- NULL
to_classify <- synapses

if (file.exists(lookup_file)) {
  message("Loading synapse neuropil lookup...")
  lookup <- arrow::read_parquet(lookup_file)
  lookup$id <- as.character(lookup$id)
  message(sprintf("  Lookup contains %s entries", format(nrow(lookup), big.mark = ",")))

  # Join lookup to synapses by id
  synapses_joined <- dplyr::left_join(synapses, lookup, by = "id")

  # Treat legacy / stale lookup entries as "needs reclassification":
  #   - "outside_*" suffixed labels (e.g. "outside_ME_L") from a deprecated
  #     scheme; we now want plain "outside" for synapses outside every alpha
  #     shape.
  #   - empty strings (would otherwise pass the !is.na test).
  # Anything else with a non-NA neuropil is trusted.
  legacy <- !is.na(synapses_joined$neuropil) &
              (grepl("^outside_", synapses_joined$neuropil) |
                 synapses_joined$neuropil == "")
  if (any(legacy)) {
    n_legacy <- sum(legacy)
    message(sprintf("  Invalidated %s legacy/empty lookup entries — will reclassify",
                    format(n_legacy, big.mark = ",")))
    synapses_joined$neuropil[legacy] <- NA_character_
    synapses_joined$region[legacy]   <- NA_character_
    synapses_joined$side[legacy]     <- NA_character_
  }

  matched <- synapses_joined[!is.na(synapses_joined$neuropil), ]
  to_classify <- synapses_joined[is.na(synapses_joined$neuropil), ]

  # Drop the empty neuropil/region/side columns from unmatched rows
  to_classify$neuropil <- NULL
  to_classify$region <- NULL
  to_classify$side <- NULL

  n_matched <- nrow(matched)
  message(sprintf("  Matched %s/%s from lookup (%.1f%%), %s to classify",
                  format(n_matched, big.mark = ","),
                  format(n_total, big.mark = ","),
                  100 * n_matched / n_total,
                  format(nrow(to_classify), big.mark = ",")))

  rm(synapses_joined, lookup)
  gc()
}

rm(synapses)
gc()

##########################################
### CLASSIFY UNMATCHED SYNAPSES        ###
##########################################

if (nrow(to_classify) > 0) {

  #####################################
  ### PRE-COMPUTE ALPHA SHAPES      ###
  #####################################

  message("Pre-computing alpha shapes (this takes ~5-10 minutes)...")
  t_alpha_start <- Sys.time()

  # 1. Large regional volumes
  message("  Building regional volume alpha shapes...")
  volume_surfs <- list(
    neck = banc_neck_connective.surf,
    brain = banc_brain_neuropil.surf,
    optic_lobes = as.hxsurf(as.mesh3d(subset(banc_brain_neuropils.surf, "optic"))),
    sez = as.hxsurf(as.mesh3d(subset(banc_brain_neuropils.surf, "GNG|CAN|FLA|AMMC|SAD|PRW"))),
    central_brain = as.hxsurf(as.mesh3d(subset(banc_brain_neuropils.surf, "midbrain"))),
    vnc = banc_vnc_neuropil.surf
  )

  volume_ashapes <- lapply(names(volume_surfs), function(vol_name) {
    message(sprintf("    %s...", vol_name))
    alphashape3d::ashape3d(nat::xyzmatrix(volume_surfs[[vol_name]]),
                           alpha = alpha_value)
  })
  names(volume_ashapes) <- names(volume_surfs)

  # 2. Individual neuropils from brain and VNC
  message("  Building individual neuropil alpha shapes...")
  brain_nps <- sort(banc_brain_neuropils.surf$RegionList)
  vnc_nps <- sort(banc_vnc_neuropils.surf$RegionList)

  # Build neuropil -> region mapping
  np_region_map <- c(
    setNames(ifelse(grepl("^LO|^ME|^AME|^LOP", brain_nps), "optic_lobes",
             ifelse(grepl("^CAN|^GNG|^FLA|^AMMC|^SAD|^PRW", brain_nps), "suboesophageal_zone",
                    "central_brain")),
             brain_nps),
    setNames(rep("vnc", length(vnc_nps)), vnc_nps)
  )

  all_nps <- c(brain_nps, vnc_nps)
  all_np_surfs <- c(
    lapply(brain_nps, function(np) subset(banc_brain_neuropils.surf, np)),
    lapply(vnc_nps, function(np) subset(banc_vnc_neuropils.surf, np))
  )
  names(all_np_surfs) <- all_nps

  n_nps <- length(all_nps)
  np_ashapes <- vector("list", n_nps)
  names(np_ashapes) <- all_nps
  for (i in seq_len(n_nps)) {
    np_name <- all_nps[i]
    if (i == 1 || i %% 10 == 0 || i == n_nps) {
      message(sprintf("    [%d/%d] %s...", i, n_nps, np_name))
    }
    np_ashapes[[np_name]] <- alphashape3d::ashape3d(
      nat::xyzmatrix(all_np_surfs[[np_name]]),
      alpha = alpha_value)
  }

  t_alpha_end <- Sys.time()
  message(sprintf("Alpha shapes pre-computed in %s",
                  format(round(difftime(t_alpha_end, t_alpha_start, units = "mins"), 1))))
  message(sprintf("  %d regional volumes, %d individual neuropils",
                  length(volume_ashapes), length(np_ashapes)))

  ##########################################
  ### DEFINE CLASSIFICATION FUNCTION     ###
  ##########################################

  classify_chunk <- function(chunk_df) {
    # Initialize output columns
    chunk_df$neuropil <- ""
    chunk_df$region <- ""
    chunk_df$side <- ""

    points <- nat::xyzmatrix(chunk_df)

    # Determine left/right side
    lrdiffs <- bancr:::banc_lr_position(points, units = "nm")
    chunk_df$side <- ifelse(lrdiffs > 0, "right", "left")

    # Test against regional volumes
    for (vol_name in names(volume_ashapes)) {
      inside <- alphashape3d::inashape3d(points = points,
                                         as3d = volume_ashapes[[vol_name]],
                                         indexAlpha = "ALL")
      if (sum(inside)) {
        chunk_df$region[which(inside)] <- vol_name
      }
    }

    # Test against individual neuropils
    for (np_name in names(np_ashapes)) {
      inside <- alphashape3d::inashape3d(points = points,
                                         as3d = np_ashapes[[np_name]],
                                         indexAlpha = "ALL")
      if (sum(inside)) {
        # A synapse can be in multiple neuropils — append with comma
        idx <- which(inside)
        chunk_df$neuropil[idx] <- sapply(chunk_df$neuropil[idx], function(x) {
          paste(unique(unlist(strsplit(paste(x, np_name, sep = ","), split = ","))),
                collapse = ",")
        })
        # Only set region for points not already assigned by regional volume test
        unassigned <- idx[chunk_df$region[idx] == ""]
        if (length(unassigned)) chunk_df$region[unassigned] <- np_region_map[np_name]
      }
    }

    # Mark unassigned
    chunk_df$neuropil[chunk_df$neuropil == ""] <- "outside"
    chunk_df$region[chunk_df$region == ""] <- "outside"

    # Clean leading commas
    chunk_df$neuropil <- gsub("^,", "", chunk_df$neuropil)
    chunk_df$region <- gsub("^,", "", chunk_df$region)

    chunk_df
  }

  ##########################################
  ### PROCESS IN PARALLEL               ###
  ##########################################

  n_to_classify <- nrow(to_classify)
  n_chunks <- ceiling(n_to_classify / chunk_size)

  message(sprintf("Processing %s unmatched synapses in %d chunks of %s",
                  format(n_to_classify, big.mark = ","),
                  n_chunks,
                  format(chunk_size, big.mark = ",")))

  # Split into chunks
  chunk_indices <- split(seq_len(n_to_classify),
                         ceiling(seq_len(n_to_classify) / chunk_size))

  t_process_start <- Sys.time()
  chunks_done <- 0

  # Check for already-completed chunks (resumable)
  existing_chunks <- list.files(temp_dir, pattern = "\\.parquet$")
  n_existing <- length(existing_chunks)
  if (n_existing > 0) {
    message(sprintf("Found %d existing chunks — resuming from where we left off", n_existing))
  }

  # Process chunks in batches of n_cores for clean progress reporting
  # Workers inherit pre-computed alpha shapes via copy-on-write (fork-based)
  batch_starts <- seq(1, n_chunks, by = n_cores)

  for (bs in batch_starts) {
    batch_ci <- bs:min(bs + n_cores - 1, n_chunks)

    # Skip batches where all chunks already exist
    batch_files <- file.path(temp_dir, sprintf("chunk_%04d.parquet", batch_ci))
    if (all(file.exists(batch_files))) {
      chunks_done <- chunks_done + length(batch_ci)
      next
    }

    parallel::mclapply(batch_ci, function(ci) {
      chunk_file <- file.path(temp_dir, sprintf("chunk_%04d.parquet", ci))
      if (file.exists(chunk_file)) return(chunk_file)

      idx <- chunk_indices[[ci]]
      chunk_df <- to_classify[idx, , drop = FALSE]
      chunk_df <- classify_chunk(chunk_df)
      arrow::write_parquet(chunk_df, chunk_file)
      chunk_file
    }, mc.cores = n_cores)

    chunks_done <- chunks_done + length(batch_ci)
    elapsed <- as.numeric(difftime(Sys.time(), t_process_start, units = "mins"))
    rate <- chunks_done / max(elapsed, 0.01)
    eta <- (n_chunks - chunks_done) / rate

    message(sprintf("  Progress: %d/%d chunks (%.0f%%) | elapsed: %.1f min | ETA: %.1f min",
                    chunks_done, n_chunks,
                    100 * chunks_done / n_chunks,
                    elapsed, eta))
  }

  t_process_end <- Sys.time()
  message(sprintf("Chunk processing complete in %s",
                  format(round(difftime(t_process_end, t_process_start, units = "mins"), 1))))

  ##########################################
  ### COMBINE NEWLY CLASSIFIED           ###
  ##########################################

  message("Combining newly classified chunks...")
  t_combine_start <- Sys.time()
  chunk_files <- list.files(temp_dir, pattern = "\\.parquet$", full.names = TRUE)
  chunk_files <- sort(chunk_files)

  newly_classified <- dplyr::bind_rows(lapply(chunk_files, arrow::read_parquet))
  message(sprintf("Combined %s newly classified rows from %d chunks [%s]",
                  format(nrow(newly_classified), big.mark = ","),
                  length(chunk_files),
                  format(round(difftime(Sys.time(), t_combine_start, units = "secs"), 0))))

  # Run pointsnearby_banc on remaining "outside" synapses
  n_outside <- sum(grepl("outside", newly_classified$region) | grepl("outside", newly_classified$neuropil))
  if (n_outside > 0) {
    message(sprintf("Running nearest-mesh fallback on %s outside synapses (%.1f%%)...",
                    format(n_outside, big.mark = ","),
                    100 * n_outside / nrow(newly_classified)))
    t_nearby_start <- Sys.time()
    newly_classified <- pointsnearby_banc(newly_classified)
    message(sprintf("Nearest-mesh fallback complete [%s]",
                    format(round(difftime(Sys.time(), t_nearby_start, units = "mins"), 1))))
  }

  # Combine with lookup-matched rows
  if (!is.null(matched) && nrow(matched) > 0) {
    combined <- dplyr::bind_rows(matched, newly_classified)
    message(sprintf("Combined: %s lookup-matched + %s newly classified = %s total",
                    format(nrow(matched), big.mark = ","),
                    format(nrow(newly_classified), big.mark = ","),
                    format(nrow(combined), big.mark = ",")))
  } else {
    combined <- newly_classified
  }

  ##########################################
  ### UPDATE LOOKUP                      ###
  ##########################################

  message("Updating synapse neuropil lookup...")
  new_lookup_entries <- newly_classified[, c("id", "neuropil", "region", "side")]
  new_lookup_entries$id <- as.character(new_lookup_entries$id)

  if (file.exists(lookup_file)) {
    old_lookup <- arrow::read_parquet(lookup_file)
    old_lookup$id <- as.character(old_lookup$id)
    # New entries take precedence
    old_only <- old_lookup[!old_lookup$id %in% new_lookup_entries$id, ]
    updated_lookup <- dplyr::bind_rows(new_lookup_entries, old_only)
  } else {
    updated_lookup <- new_lookup_entries
  }

  arrow::write_parquet(updated_lookup, lookup_file)
  message(sprintf("  Lookup updated: %s total entries (+%s new)",
                  format(nrow(updated_lookup), big.mark = ","),
                  format(nrow(new_lookup_entries), big.mark = ",")))

  rm(newly_classified, matched, new_lookup_entries)
  gc()

} else {
  message("All synapses matched from lookup — skipping alpha shape computation")
  combined <- matched
}

##########################################
### CLEAN NEUROPIL NAMES              ###
##########################################

message("Cleaning neuropil names...")
combined <- combined %>%
  dplyr::mutate(
    neuropil = gsub("ITO_optic_|ITO_midbrain_|COURT_vnc_", "", neuropil),
    neuropil = gsub("^,", "", neuropil),
    region = gsub("^,", "", region)
  )

##########################################
### SAVE                              ###
##########################################

write_connectome_data(combined, output_parquet, format = "parquet")
message(sprintf("Saved: %s (%s rows)",
                output_parquet,
                format(nrow(combined), big.mark = ",")))

# Report summary
message("\nNeuropil assignment summary:")
region_summary <- table(combined$region, useNA = "ifany")
print(region_summary)

##########################################
### COMPLETION METRICS                ###
##########################################

message("\nCalculating completion metrics...")

# Load the processed neuropil parquet so this section can run independently
if (file.exists(output_parquet)) {
  message(sprintf("Loading processed synapses from %s...", output_parquet))
  combined <- arrow::read_parquet(output_parquet)
} else {
  stop("Output parquet not found — run neuropil classification first: ", output_parquet)
}

# Get neuron IDs from versioned meta (produced by banc-calculate-connectivity.R)
meta_file <- file.path(banc.connectivity.save.path,
                        sprintf("banc_%s_meta.feather", banc.version))
if (file.exists(meta_file)) {
  banc.meta <- arrow::read_feather(meta_file)
  neuron.ids <- unique(as.character(banc.meta$root_id))

  # Classify pre/post root_id as neuron vs fragment
  combined <- combined %>%
    dplyr::mutate(
      pre_status = ifelse(pre_root_id %in% neuron.ids, "neuron", "fragment"),
      post_status = ifelse(post_root_id %in% neuron.ids, "neuron", "fragment")
    )

  # 1. Gross capture rates: pre_status x post_status
  summary_gross <- combined %>%
    dplyr::count(pre_status, post_status) %>%
    dplyr::mutate(prop = round(n / sum(n), 4))
  readr::write_csv(summary_gross,
    file.path(banc.connectivity.save.path,
              sprintf("banc_%s_gross_capture_rates.csv", banc.version)))
  message("\n  Gross capture rates:")
  print(as.data.frame(summary_gross))

  # 2. Inside/outside mesh capture rates
  summary_inout <- combined %>%
    dplyr::mutate(in_mesh = ifelse(region == "outside", "outside", "inside")) %>%
    dplyr::group_by(in_mesh) %>%
    dplyr::count(pre_status, post_status, in_mesh) %>%
    dplyr::mutate(prop = round(n / sum(n), 4)) %>%
    dplyr::ungroup()
  readr::write_csv(summary_inout,
    file.path(banc.connectivity.save.path,
              sprintf("banc_%s_inout_capture_rates.csv", banc.version)))

  # 3. By region and side
  summary_region <- combined %>%
    dplyr::group_by(region, side) %>%
    dplyr::count(pre_status, post_status, side, region) %>%
    dplyr::mutate(prop = round(n / sum(n), 4)) %>%
    dplyr::ungroup()
  readr::write_csv(summary_region,
    file.path(banc.connectivity.save.path,
              sprintf("banc_%s_region_capture_rates.csv", banc.version)))

  # 4. By neuropil
  summary_neuropil <- combined %>%
    dplyr::group_by(region, side, neuropil) %>%
    dplyr::count(pre_status, post_status, side, region, neuropil) %>%
    dplyr::mutate(prop = round(n / sum(n), 4)) %>%
    dplyr::ungroup()
  readr::write_csv(summary_neuropil,
    file.path(banc.connectivity.save.path,
              sprintf("banc_%s_neuropil_capture_rates.csv", banc.version)))

  # Overall summary
  n_total <- nrow(combined)
  n_pre_neuron <- sum(combined$pre_status == "neuron")
  n_post_neuron <- sum(combined$post_status == "neuron")
  n_both_neuron <- sum(combined$pre_status == "neuron" & combined$post_status == "neuron")
  n_outside <- sum(combined$region == "outside")

  message(sprintf("\n  Total synapses: %s", format(n_total, big.mark = ",")))
  message(sprintf("  Pre neuron: %s (%.1f%%)",
                  format(n_pre_neuron, big.mark = ","), 100 * n_pre_neuron / n_total))
  message(sprintf("  Post neuron: %s (%.1f%%)",
                  format(n_post_neuron, big.mark = ","), 100 * n_post_neuron / n_total))
  message(sprintf("  Both neuron: %s (%.1f%%)",
                  format(n_both_neuron, big.mark = ","), 100 * n_both_neuron / n_total))
  message(sprintf("  Outside mesh: %s (%.1f%%)",
                  format(n_outside, big.mark = ","), 100 * n_outside / n_total))
  message(sprintf("  Saved 4 completion metrics CSVs to %s", banc.connectivity.save.path))
} else {
  message("  Meta feather not found — skipping completion metrics")
  message(sprintf("  (Expected: %s)", meta_file))
}

# Cleanup temp chunks
message("\nCleaning up temporary chunk files...")
unlink(temp_dir, recursive = TRUE)

message(sprintf("### banc: neuropil inclusion calculation complete [total: %s] ###",
                format(round(difftime(Sys.time(), t_total_start, units = "mins"), 1))))

})
