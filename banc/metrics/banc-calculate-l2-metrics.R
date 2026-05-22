#' banc-calculate-l2-metrics — Merge per-neuron L2 + split metrics into a single feather.
#'
#' Extracted columns: `l2_nodes`, `l2_cable_length_um`, `pd_width`, `segregation_index`.
#'
#' @section Reads:
#'   - per-neuron L2 + split metric CSVs (from `banc-l2.R` / `banc-calculate-split.R`)
#'
#' @section Writes:
#'   - `banc_l2_metrics.feather`
#'
#' @section Invoked by:
#'   production v888 rebuild chain.
#'
#' @section Used by:
#'   `banc/share/banc-data.R` (joins into `banc_888_metrics.feather`).

###########################################################
### Read L2 skeleton metrics and split metrics CSVs
### and merge into banc_l2_metrics.feather
###
### Reads per-neuron metrics CSVs from banc.metrics.save.path
### (written by banc-l2.R) and split metrics CSVs from
### banc.l2split.save.path/metrics/ (written by split pipeline)
###
### Extracts: l2_nodes, l2_cable_length_um, pd_width,
###           segregation_index
###
### Saves results to banc_l2_metrics.feather
###########################################################
source("banc/banc-startup.R")

local({

message("### banc: reading L2 and split metrics ###")
t_start <- Sys.time()
bancr::choose_banc()

###########################
### Read current state  ###
###########################

bc <- banctable_query("SELECT _id, status, root_id from banc_meta") %>%
  banc_filter_neurons()
bc.root.ids <- unique(bc$root_id)
# Honor test mode if active
if (exists("banc.test.ids", envir = .GlobalEnv))
  bc.root.ids <- intersect(bc.root.ids, banc.test.ids)

# Read existing L2 metrics feather
feather_file <- file.path(banc.save.path, "banc_l2_metrics.feather")
if (file.exists(feather_file)) {
  banc.l2metrics <- arrow::read_feather(feather_file) %>%
    dplyr::distinct(root_id, .keep_all = TRUE)
} else {
  banc.l2metrics <- data.frame(root_id = character(0))
}

# Determine which neurons have metrics files but not yet in feather
has_l2 <- banc.l2metrics %>%
  dplyr::filter(!is.na(l2_nodes)) %>%
  dplyr::pull(root_id)

###########################
### Read L2 metrics     ###
###########################

# Find metrics files for neurons not yet integrated
metrics_files <- file.path(banc.metrics.save.path, paste0(bc.root.ids, ".csv"))
metrics_files <- metrics_files[file.exists(metrics_files)]
metrics_ids <- gsub("\\.csv$", "", basename(metrics_files))
new_metrics_files <- metrics_files[!metrics_ids %in% has_l2]

l2_data <- NULL
if (length(new_metrics_files) > 0) {
  message(sprintf("Reading %d L2 metrics files...", length(new_metrics_files)))
  n_mfiles <- length(new_metrics_files)
  by.query.metrics <- foreach::foreach(mfile = new_metrics_files, i = seq_along(new_metrics_files)) %do% {
    if (i == 1 || i %% 500 == 0 || i == n_mfiles) {
      message(sprintf("  L2 metrics: %d/%d (%.0f%%)", i, n_mfiles, 100 * i / n_mfiles))
    }
    tryCatch({
      mdf <- readr::read_csv(mfile, col_types = banc.col.types, show_col_types = FALSE)
      numeric_columns <- sapply(mdf, is.numeric)
      mdf[numeric_columns] <- lapply(mdf[numeric_columns], function(x) { x[!is.finite(x)] <- NA; x })
      mdf[numeric_columns] <- round(mdf[numeric_columns], digits = 2)
      mdf[, c("l2_n_trees", "l2_endpoints", "l2_branchpoints", "l2_segments", "l2_root")] <- NULL
      mdf
    }, error = function(e) {
      message(sprintf("  Error reading %s: %s", basename(mfile), e$message))
      NULL
    })
  }
  by.query.metrics <- by.query.metrics[unlist(lapply(by.query.metrics, is.data.frame))]

  if (length(by.query.metrics) > 0) {
    l2_data <- do.call(plyr::rbind.fill, by.query.metrics) %>%
      dplyr::mutate(l2_cable_length_um = round(l2_cable_length / 1000, 2)) %>%
      dplyr::distinct(root_id, .keep_all = TRUE) %>%
      dplyr::arrange(dplyr::desc(l2_nodes))
    message(sprintf("Read L2 metrics for %d neurons", nrow(l2_data)))
  }
}

###########################
### Read split metrics  ###
###########################

split.metrics.folder <- file.path(banc.l2split.save.path, "metrics")
has_split <- banc.l2metrics %>%
  dplyr::filter(!is.na(pd_width) | !is.na(segregation_index)) %>%
  dplyr::pull(root_id)

split_data <- NULL
if (dir.exists(split.metrics.folder)) {
  split_files <- list.files(split.metrics.folder, pattern = "\\.csv$", full.names = TRUE)
  if (length(split_files) > 0) {
    n_sfiles <- length(split_files)
    message(sprintf("Reading %d split metrics files...", n_sfiles))
    by.query.split <- foreach::foreach(mfile = split_files, i = seq_along(split_files)) %do% {
      if (i == 1 || i %% 500 == 0 || i == n_sfiles) {
        message(sprintf("  Split metrics: %d/%d (%.0f%%)", i, n_sfiles, 100 * i / n_sfiles))
      }
      tryCatch({
        mdf <- readr::read_csv(mfile, col_types = banc.col.types, show_col_types = FALSE)
        numeric_columns <- sapply(mdf, is.numeric)
        mdf[numeric_columns] <- round(mdf[numeric_columns], digits = 2)
        # Remove summary columns that conflict with L2 metrics
        for (col in c("n_trees", "endpoints", "branchpoints", "segments", "root")) {
          if (col %in% colnames(mdf)) mdf[[col]] <- NULL
        }
        mdf
      }, error = function(e) {
        message(sprintf("  Error reading %s: %s", basename(mfile), e$message))
        NULL
      })
    }
    by.query.split <- by.query.split[unlist(lapply(by.query.split, is.data.frame))]

    if (length(by.query.split) > 0) {
      split_data <- do.call(plyr::rbind.fill, by.query.split) %>%
        dplyr::mutate(cable_length_um = round(as.numeric(cable_length) / 1000, 2)) %>%
        dplyr::distinct(root_id, .keep_all = TRUE) %>%
        dplyr::arrange(dplyr::desc(nodes))
      # Remove columns that would conflict
      for (col in c("nucleus_id", "nucleus_position_nm", "root_position_nm")) {
        if (col %in% colnames(split_data)) split_data[[col]] <- NULL
      }
      # Filter to current root_ids and those not already processed
      split_data <- split_data %>%
        dplyr::filter(root_id %in% bc.root.ids & !root_id %in% has_split)
      message(sprintf("Read split metrics for %d neurons", nrow(split_data)))
    }
  }
}

###########################
### Merge and update    ###
###########################

if (is.null(l2_data) && is.null(split_data)) {
  message("No new metrics to integrate. Nothing to do.")
  return(invisible())
}

# Combine L2 and split metrics
if (!is.null(l2_data) && !is.null(split_data)) {
  new_data <- dplyr::full_join(l2_data, split_data, by = "root_id")
} else if (!is.null(l2_data)) {
  new_data <- l2_data
} else {
  new_data <- split_data
}
new_data <- new_data %>% dplyr::distinct(root_id, .keep_all = TRUE)

# Update feather
if (nrow(banc.l2metrics) > 0) {
  unchanged <- banc.l2metrics %>% dplyr::filter(!root_id %in% new_data$root_id)
  final <- dplyr::bind_rows(new_data, unchanged) %>%
    dplyr::distinct(root_id, .keep_all = TRUE)
} else {
  final <- new_data
}

arrow::write_feather(final, feather_file)
message(sprintf("### banc: L2 + split metrics updated for %d neurons [%s] ###",
                nrow(new_data),
                format(round(difftime(Sys.time(), t_start, units = "mins"), 1))))

###########################
### Diagnostic plot     ###
###########################

tryCatch({
  has_l2_final <- sum(!is.na(final$l2_cable_length_um))
  cable_df <- final %>% dplyr::filter(!is.na(l2_cable_length_um), l2_cable_length_um > 0)
  g <- ggplot2::ggplot(cable_df, ggplot2::aes(x = l2_cable_length_um)) +
    ggplot2::geom_histogram(bins = 50, fill = "#0072B2", alpha = 0.8) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(title = sprintf("L2 cable length distribution (%s)", Sys.Date()),
                  subtitle = sprintf("%d/%d neurons with L2 metrics (%.1f%%)",
                                     has_l2_final, nrow(final),
                                     100 * has_l2_final / nrow(final)),
                  x = expression("Cable length ("*mu*"m, log scale)"),
                  y = "Neuron count") +
    ggplot2::theme_minimal()
  ggplot2::ggsave(plot = g,
                  filename = "inst/images/banc_l2_cable_length_distribution.png",
                  width = 10, height = 6, dpi = 300, bg = "transparent")
}, error = function(e) message("  Plot failed: ", e$message))

})
