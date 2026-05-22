#' banc-calculate-synapses — Per-neuron synapse / mitochondria metrics.
#'
#' For each `root_888` neuron: if a split CSV exists (from
#' `banc-calculate-split.R`), use it (includes axon/dendrite compartment
#' labels); otherwise extract from the v888 synapse parquet. Tracks
#' `source` ("split" or "table") so re-runs automatically upgrade
#' table-derived neurons when new split CSVs appear.
#'
#' @section Reads:
#'   - `<banc.l2split.save.path>/synapses/<root_id>.csv` (preferred)
#'   - `<banc.connectivity.save.path>/banc_<ver>_synapses_<src>.parquet` (fallback)
#'   - SeaTable `banc_meta`
#'
#' @section Writes:
#'   - `<banc.save.path>/banc_detections.feather`
#'     (cols: input_connections, output_connections,
#'      input_side_index, output_side_index,
#'      mitochondria, mitochondria_volume, source)
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_metrics.sh`,
#'   `o2/production/o2_banc_synapses.sh`,
#'   `o2/production/o2_banc_update.sh`
#'
#' @section Used by:
#'   `banc/share/banc-data.R` (joins into banc_888_metrics.feather)

###########################################################
### Calculate synapse metrics for BANC neurons (v888)
###
### For each root_888 neuron from seatable:
###   1. If split CSV exists (from banc-calculate-split.R),
###      use it (includes axon/dendrite compartment labels)
###   2. Otherwise, extract from the v888 synapse table
###
### On re-run, neurons with source="table" are re-processed
### when a split CSV has since become available.
###
### Metrics:
###   input_connections, output_connections
###   input_side_index, output_side_index
###   mitochondria, mitochondria_volume
###
### Output: banc_detections.feather
###########################################################
source("banc/banc-startup.R")
library(data.table)

local({

message("### banc: calculating synapse metrics ###")
t_start <- Sys.time()
bancr::choose_banc()

#######################
### CONFIGURATION   ###
#######################

# banc.version set in banc-startup.R
save.path <- file.path("/n/data1/hms/neurobio/wilson/connectomes/banc",
                        paste0("banc_", banc.version))
synapse_size_threshold <- 5

# Split synapse CSVs (from banc-calculate-split.R)
split_syn_dir <- file.path(banc.l2split.save.path, "synapses")

# Output
feather_file <- file.path(banc.save.path, "banc_detections.feather")

###########################
### Get target neurons  ###
###########################

bc <- banctable_query("SELECT _id, status, root_id FROM banc_meta") %>%
  banc_filter_neurons() %>%
  dplyr::mutate(root_id = as.character(root_id))

target_ids <- unique(bc$root_id[!is.na(bc$root_id) & bc$root_id != ""])
if (exists("banc.test.ids", envir = .GlobalEnv))
  target_ids <- intersect(target_ids, as.character(banc.test.ids))
message(sprintf("Target neurons (root_%s): %d", banc.version, length(target_ids)))

###########################
### Load existing state ###
###########################

if (file.exists(feather_file)) {
  existing <- arrow::read_feather(feather_file) %>%
    dplyr::mutate(root_id = as.character(root_id)) %>%
    dplyr::mutate(dplyr::across(
      dplyr::any_of(c("mitochondria", "mitochondria_volume",
                       "input_connections", "output_connections",
                       "input_side_index", "output_side_index")),
      as.numeric)) %>%
    dplyr::distinct(root_id, .keep_all = TRUE)
} else {
  existing <- data.frame(root_id = character(0))
}

# Ensure source column for tracking split vs table origin
if (!"source" %in% names(existing)) existing$source <- "table"

###########################
### Identify work       ###
###########################

# Available split CSVs
if (dir.exists(split_syn_dir)) {
  split_files <- list.files(split_syn_dir, pattern = "\\.csv$", full.names = FALSE)
  available_split_ids <- intersect(tools::file_path_sans_ext(split_files), target_ids)
} else {
  available_split_ids <- character(0)
}

# New neurons not in feather
new_ids <- setdiff(target_ids, existing$root_id)

# Existing table-derived neurons upgradeable to split
upgradeable_ids <- existing$root_id[existing$source == "table" &
                                     existing$root_id %in% available_split_ids]

# Work sets
to_process_split <- intersect(union(new_ids, upgradeable_ids), available_split_ids)
to_process_table <- setdiff(new_ids, available_split_ids)

n_total <- length(to_process_split) + length(to_process_table)
message(sprintf("  New: %d | Upgradeable (table->split): %d | Total: %d",
                length(new_ids), length(upgradeable_ids), n_total))

if (n_total == 0) {
  message("All synapse metrics up to date.")
  return(invisible())
}

###################################
### Process split neurons       ###
###################################

split_results <- data.table()

if (length(to_process_split) > 0) {
  message(sprintf("Processing %d neurons from split CSVs...", length(to_process_split)))
  pb <- utils::txtProgressBar(max = length(to_process_split), style = 3)

  split_list <- lapply(seq_along(to_process_split), function(i) {
    utils::setTxtProgressBar(pb, i)
    rid <- to_process_split[i]
    f <- file.path(split_syn_dir, paste0(rid, ".csv"))
    if (!file.exists(f)) return(NULL)

    syns <- tryCatch(
      fread(f, select = c("prepost", "X", "Y", "Z")),
      error = function(e) {
        tryCatch({
          dt <- fread(f, select = c("prepost", "x", "y", "z"))
          setnames(dt, c("x", "y", "z"), c("X", "Y", "Z"))
          dt
        }, error = function(e2) NULL)
      }
    )
    if (is.null(syns) || nrow(syns) == 0) return(NULL)

    n_in <- sum(syns$prepost == 1L, na.rm = TRUE)
    n_out <- sum(syns$prepost == 0L, na.rm = TRUE)

    # Laterality
    pos_nm <- banc_raw2nm(as.matrix(syns[, .(X, Y, Z)]))
    lr <- bancr:::banc_lr_position(pos_nm, units = "nm")
    sides <- ifelse(lr > 0, "right", "left")

    in_mask <- syns$prepost == 1L
    out_mask <- syns$prepost == 0L

    data.table(
      root_id = rid,
      input_connections = n_in,
      output_connections = n_out,
      input_side_index = if (n_in == 0) NA_real_ else
        round((sum(sides[in_mask] == "right") - sum(sides[in_mask] == "left")) / n_in, 4),
      output_side_index = if (n_out == 0) NA_real_ else
        round((sum(sides[out_mask] == "right") - sum(sides[out_mask] == "left")) / n_out, 4),
      source = "split"
    )
  })

  close(pb)
  split_results <- rbindlist(Filter(Negate(is.null), split_list))
  message(sprintf("  Split results: %d neurons", nrow(split_results)))
}

###################################
### Process table neurons       ###
###################################

table_results <- data.table()

if (length(to_process_table) > 0) {
  message(sprintf("Loading v%s synapse table...", banc.version))

  # Prefer parquet (faster), fall back to raw CSV
  syn_candidates <- c(
    file.path(banc.connectivity.save.path,
              sprintf("banc_%s_synapses.parquet", banc.version)),
    file.path(save.path,
              sprintf("banc_%s_synapses_enriched.parquet", banc.version)),
    file.path(banc.save.path, paste0("v", banc.version),
              "synapses_v2_human_readable.csv"),
    bancsynapses
  )
  syn_path <- NULL
  for (sp in syn_candidates) {
    if (file.exists(sp)) { syn_path <- sp; break }
  }
  if (is.null(syn_path))
    stop("Could not find v", banc.version, " synapse table. Checked:\n  ",
         paste(syn_candidates, collapse = "\n  "))
  message(sprintf("  Using: %s", syn_path))

  if (grepl("\\.parquet$", syn_path)) {
    # Parquet: columns are id, pre_root_id, post_root_id, X, Y, Z, size
    all_syns <- as.data.table(arrow::read_parquet(syn_path,
      col_select = c("id", "pre_root_id", "post_root_id", "X", "Y", "Z", "size")))
  } else {
    # Raw CSV: 15 columns, skip header row, select by position
    # 1=id 8=ctr_x 9=ctr_y 10=ctr_z 11=size 13=pre_root_id 15=post_root_id
    all_syns <- fread(syn_path, skip = 1L, header = FALSE,
                      select = c(1L, 8L, 9L, 10L, 11L, 13L, 15L))
    setnames(all_syns, c("id", "X", "Y", "Z", "size", "pre_root_id", "post_root_id"))
  }

  all_syns[, `:=`(pre_root_id = as.character(pre_root_id),
                  post_root_id = as.character(post_root_id))]

  # Filter to relevant neurons, size threshold, no autapses
  all_syns <- all_syns[
    size > synapse_size_threshold &
    pre_root_id != post_root_id &
    (pre_root_id %in% to_process_table | post_root_id %in% to_process_table)
  ]
  message(sprintf("  Filtered: %s synapses", format(nrow(all_syns), big.mark = ",")))

  # Compute laterality for all synapses at once (vectorized)
  pos_nm <- banc_raw2nm(as.matrix(all_syns[, .(X, Y, Z)]))
  all_syns[, side := ifelse(bancr:::banc_lr_position(pos_nm, units = "nm") > 0,
                             "right", "left")]

  # Input metrics (neuron is post_root_id)
  input_metrics <- all_syns[post_root_id %in% to_process_table,
    .(input_connections = .N,
      input_side_index = round((sum(side == "right") - sum(side == "left")) / .N, 4)),
    by = .(root_id = post_root_id)]

  # Output metrics (neuron is pre_root_id)
  output_metrics <- all_syns[pre_root_id %in% to_process_table,
    .(output_connections = .N,
      output_side_index = round((sum(side == "right") - sum(side == "left")) / .N, 4)),
    by = .(root_id = pre_root_id)]

  # Merge
  table_results <- data.table(root_id = to_process_table)
  table_results <- merge(table_results, input_metrics, by = "root_id", all.x = TRUE)
  table_results <- merge(table_results, output_metrics, by = "root_id", all.x = TRUE)
  table_results[is.na(input_connections), input_connections := 0L]
  table_results[is.na(output_connections), output_connections := 0L]
  table_results[, source := "table"]

  rm(all_syns, pos_nm)
  gc()
  message(sprintf("  Table results: %d neurons", nrow(table_results)))
}

###################################
### Mitochondria (batch)        ###
###################################

all_to_process <- c(to_process_split, to_process_table)
message(sprintf("Querying mitochondria for %d neurons...", length(all_to_process)))

mito_summary <- tryCatch({
  mitos <- banc_mitochondria(all_to_process)
  as.data.table(mitos)[, .(
    mitochondria = .N,
    mitochondria_volume = sum(volume, na.rm = TRUE)
  ), by = .(root_id = as.character(pt_root_id))]
}, error = function(e) {
  message("  Mitochondria query failed: ", e$message)
  data.table(root_id = character(0),
             mitochondria = numeric(0),
             mitochondria_volume = numeric(0))
})
message(sprintf("  Mitochondria data for %d neurons", nrow(mito_summary)))

###########################
### Combine results     ###
###########################

new_data <- rbindlist(list(split_results, table_results), fill = TRUE)
new_data <- merge(new_data, mito_summary, by = "root_id", all.x = TRUE)
new_data[is.na(mitochondria), `:=`(mitochondria = 0, mitochondria_volume = 0)]
new_data[, root_id := as.character(root_id)]

# Replace NaN with NA
num_cols <- c("input_connections", "output_connections",
              "input_side_index", "output_side_index",
              "mitochondria", "mitochondria_volume")
for (col in intersect(num_cols, names(new_data))) {
  new_data[is.nan(get(col)), (col) := NA_real_]
}

###########################
### Update feather      ###
###########################

if (nrow(existing) > 0) {
  unchanged <- as.data.table(existing)[!root_id %in% new_data$root_id]
  # Align columns
  all_cols <- union(names(new_data), names(unchanged))
  for (col in setdiff(all_cols, names(unchanged))) unchanged[, (col) := NA]
  for (col in setdiff(all_cols, names(new_data))) new_data[, (col) := NA]
  final <- rbindlist(list(new_data, unchanged), use.names = TRUE, fill = TRUE)
} else {
  final <- new_data
}

final <- as.data.frame(final)
final <- final[!duplicated(final$root_id), ]

arrow::write_feather(final, feather_file)
message(sprintf("### banc: synapse metrics updated: %d new/upgraded, %d total [%s] ###",
                nrow(new_data), nrow(final),
                format(round(difftime(Sys.time(), t_start, units = "mins"), 1))))

###########################
### Diagnostic plot     ###
###########################

tryCatch({
  syn_df <- final %>%
    dplyr::filter(!is.na(input_connections) | !is.na(output_connections)) %>%
    dplyr::select(root_id, input_connections, output_connections, source) %>%
    tidyr::pivot_longer(cols = c(input_connections, output_connections),
                        names_to = "type", values_to = "count") %>%
    dplyr::mutate(type = gsub("_connections", "", type))

  has_synapses <- sum(!is.na(final$input_connections))
  n_split <- sum(final$source == "split", na.rm = TRUE)
  n_table <- sum(final$source == "table", na.rm = TRUE)

  g <- ggplot2::ggplot(syn_df, ggplot2::aes(x = count, fill = type)) +
    ggplot2::geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
    ggplot2::scale_x_log10() +
    ggplot2::labs(
      title = sprintf("Synapse count distribution (%s)", Sys.Date()),
      subtitle = sprintf("%d/%d neurons (%.1f%%) | split: %d, table: %d",
                         has_synapses, nrow(final),
                         100 * has_synapses / nrow(final),
                         n_split, n_table),
      x = "Synapse count (log scale)", y = "Neuron count", fill = "Type") +
    ggplot2::scale_fill_manual(values = c("input" = "#0072B2", "output" = "#D55E00")) +
    ggplot2::theme_minimal()
  ggplot2::ggsave(plot = g,
                  filename = "inst/images/banc_synapse_distribution.png",
                  width = 10, height = 6, dpi = 300, bg = "transparent")
}, error = function(e) message("  Plot failed: ", e$message))

})
