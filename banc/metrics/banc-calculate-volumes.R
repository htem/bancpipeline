#' banc-calculate-volumes — Per-neuron volume via the CAVE L2 cache.
#'
#' Sums L2-chunk `size_nm3` per neuron; incremental.
#'
#' @section Reads:
#'   - SeaTable `banc_meta`, CAVE L2 cache
#'
#' @section Writes:
#'   - `banc_volumes.feather`
#'
#' @section Invoked by:
#'   production v888 rebuild chain.
#'
#' @section Used by:
#'   `banc/share/banc-data.R` (joins into `banc_888_metrics.feather`).

###########################################################
### Calculate neuron volumes via L2 cache for BANC neurons
###
### Uses banc_neuron_volume() to sum L2 chunk size_nm3
### values from CAVE. Caches results per neuron.
###
### Saves results to banc_volumes.feather
###########################################################
source("banc/banc-startup.R")

local({

message("### banc: calculating neuron volumes ###")
t_start <- Sys.time()
bancr::choose_banc()

###########################
### Read current state  ###
###########################

bc <- banctable_query_cached() %>%
  banc_filter_neurons() %>%
  dplyr::mutate(root_id = as.character(root_id))

# Read existing volumes feather
feather_file <- file.path(banc.save.path, "banc_volumes.feather")
if (file.exists(feather_file)) {
  banc.volumes <- arrow::read_feather(feather_file) %>%
    dplyr::mutate(root_id = as.character(root_id)) %>%
    dplyr::distinct(root_id, .keep_all = TRUE)
} else {
  banc.volumes <- data.frame(root_id = character(0))
}

# Determine which neurons need volume calculation
has_volume <- banc.volumes %>%
  dplyr::filter(!is.na(volume_nm3)) %>%
  dplyr::pull(root_id)

to_process <- bc %>%
  dplyr::filter(!root_id %in% has_volume)
# Honor test mode if active
if (exists("banc.test.ids", envir = .GlobalEnv))
  to_process <- to_process %>% dplyr::filter(root_id %in% as.character(banc.test.ids))

message(sprintf("%d neurons need volume calculation", nrow(to_process)))

if (nrow(to_process) == 0) {
  message("All volumes are up to date. Nothing to do.")
  return(invisible())
}

###########################
### Calculate volumes   ###
###########################

# Create volume cache directory
volumes.save.path <- file.path(banc.save.path, "volumes")
dir.create(volumes.save.path, showWarnings = FALSE, recursive = TRUE)

volume_list <- pbapply::pblapply(seq_len(nrow(to_process)), function(i) {
  rid <- to_process$root_id[i]
  vol_file <- file.path(volumes.save.path, paste0(rid, ".csv"))

  tryCatch({
    if (file.exists(vol_file)) {
      vol <- readr::read_csv(vol_file, show_col_types = FALSE, progress = FALSE,
                              col_types = readr::cols(root_id = "c"))
    } else {
      vol <- banc_neuron_volume(rid, OmitFailures = FALSE)
      vol$root_id <- as.character(vol$root_id)
      readr::write_csv(vol, vol_file)
    }
    vol
  }, error = function(e) {
    message(sprintf("  Error for %s: %s", rid, e$message))
    NULL
  })
})
names(volume_list) <- to_process$root_id
volume_list <- Filter(Negate(is.null), volume_list)

if (length(volume_list) == 0) {
  message("No volumes calculated. Nothing to save.")
  return(invisible())
}

new_data <- dplyr::bind_rows(volume_list) %>%
  dplyr::distinct(root_id, .keep_all = TRUE)
message(sprintf("Calculated volumes for %d neurons", nrow(new_data)))

###########################
### Update feather      ###
###########################

if (nrow(banc.volumes) > 0) {
  unchanged <- banc.volumes %>% dplyr::filter(!root_id %in% new_data$root_id)
  final <- dplyr::bind_rows(new_data, unchanged) %>%
    dplyr::distinct(root_id, .keep_all = TRUE)
} else {
  final <- new_data
}

arrow::write_feather(final, feather_file)
message(sprintf("### banc: volumes updated for %d neurons [%s] ###",
                nrow(new_data),
                format(round(difftime(Sys.time(), t_start, units = "mins"), 1))))

###########################
### Diagnostic plot     ###
###########################

tryCatch({
  has_vol <- sum(!is.na(final$volume_nm3))
  vol_df <- final %>% dplyr::filter(!is.na(volume_nm3), volume_nm3 > 0)
  g <- ggplot2::ggplot(vol_df, ggplot2::aes(x = volume_nm3 / 1e9)) +
    ggplot2::geom_histogram(bins = 50, fill = "#0072B2", alpha = 0.8) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(title = sprintf("Neuron volume distribution (%s)", Sys.Date()),
                  subtitle = sprintf("%d/%d neurons with volume data (%.1f%%)",
                                     has_vol, nrow(final),
                                     100 * has_vol / nrow(final)),
                  x = expression("Volume ("*mu*"m"^3*", log scale)"),
                  y = "Neuron count") +
    ggplot2::theme_minimal()
  ggplot2::ggsave(plot = g,
                  filename = "inst/images/banc_volume_distribution.png",
                  width = 10, height = 6, dpi = 300, bg = "transparent")
}, error = function(e) message("  Plot failed: ", e$message))

})
