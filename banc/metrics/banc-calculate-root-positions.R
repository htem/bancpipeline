#' banc-calculate-root-positions — Compute root / soma position per BANC neuron.
#'
#' Prefers nucleus position if available; otherwise picks the furthest L2
#' skeleton leaf outside neuropil + neck connective. Incremental.
#'
#' @section Reads:
#'   - SeaTable `banc_meta`, L2 SWCs
#'
#' @section Writes:
#'   - `banc_root_positions.feather`
#'
#' @section Invoked by:
#'   production v888 rebuild chain.

###########################################################
### Calculate root/soma positions for BANC neurons
###
### For each neuron missing root_position_nm:
###   - If it has a nucleus: use nucleus_position_nm
###   - Otherwise: estimate from L2 skeleton leaf endpoints
###     (pick furthest leaf outside neuropil + neck connective)
###
### Saves results to banc_root_positions.feather
###########################################################
source("banc/banc-startup.R")

local({

message("### banc: calculating root positions ###")
t_start <- Sys.time()
bancr::choose_banc()

###########################
### Read current state  ###
###########################

bc <- banctable_query("SELECT _id, status, proofread, root_id, supervoxel_id, position, root_position, root_position_nm, nucleus_id, nucleus_position, nucleus_position_nm from banc_meta") %>%
  banc_filter_neurons() %>%
  dplyr::mutate(root_id = as.character(root_id),
                nucleus_id = as.character(nucleus_id))

# Read existing root positions feather
feather_file <- file.path(banc.save.path, "banc_root_positions.feather")
if (file.exists(feather_file)) {
  banc.positions <- arrow::read_feather(feather_file) %>%
    dplyr::mutate(root_id = as.character(root_id),
                  nucleus_id = as.character(nucleus_id)) %>%
    dplyr::distinct(root_id, .keep_all = TRUE)
} else {
  banc.positions <- data.frame(root_id = character(0))
}

# Helper: detect suspiciously large coordinate strings
too_many_digits <- function(x) {
  sapply(strsplit(x, ","), function(parts) {
    any(nchar(trimws(parts)) > 7)
  })
}

# Adopt valid SeaTable root_position into the feather —
# users sometimes manually insert good positions (raw coords), so always trust SeaTable
seatable_valid <- bc %>%
  dplyr::filter(!is.na(root_position),
                root_position != "") %>%
  dplyr::select(root_id, root_position_seatable = root_position)

# Filter out rows where root_position looks like nm values (not raw voxels)
bad_pos <- too_many_digits(seatable_valid$root_position_seatable)
if (any(bad_pos)) {
  message(sprintf("  Excluding %d root_position values that appear to be nm, not voxels",
                  sum(bad_pos)))
  seatable_valid <- seatable_valid[!bad_pos, ]
}

if (nrow(seatable_valid)) {
  # Convert all valid SeaTable raw coords to nm
  seatable_valid <- tryCatch({
    raw_coords <- xyzmatrix(seatable_valid$root_position_seatable)
    nm_coords <- banc_raw2nm(raw_coords)
    seatable_valid$root_position_nm_seatable <- apply(
      nm_coords, 1, function(r) paste(r, collapse = ","))
    seatable_valid
  }, error = function(e) {
    message("  SeaTable coord conversion failed: ", e$message)
    seatable_valid$root_position_nm_seatable <- NA_character_
    seatable_valid
  })

  if (nrow(banc.positions)) {
    # Update existing feather entries with SeaTable values
    banc.positions <- banc.positions %>%
      dplyr::left_join(seatable_valid, by = "root_id") %>%
      dplyr::mutate(
        root_position = ifelse(!is.na(root_position_seatable),
                               root_position_seatable, root_position),
        root_position_nm = ifelse(!is.na(root_position_nm_seatable),
                                   root_position_nm_seatable, root_position_nm)
      ) %>%
      dplyr::select(-root_position_seatable, -root_position_nm_seatable)
  }

  # Add neurons with valid SeaTable positions not yet in the feather
  new_from_seatable <- seatable_valid %>%
    dplyr::filter(!root_id %in% banc.positions$root_id,
                  !is.na(root_position_nm_seatable)) %>%
    dplyr::transmute(root_id,
                     root_position = root_position_seatable,
                     root_position_nm = root_position_nm_seatable,
                     nucleus_id = NA_character_)
  if (nrow(new_from_seatable)) {
    banc.positions <- dplyr::bind_rows(banc.positions, new_from_seatable)
  }

  n_adopted <- sum(banc.positions$root_id %in% seatable_valid$root_id)
  message(sprintf("Adopted %d valid root positions from SeaTable (%d new to feather)",
                  n_adopted, nrow(new_from_seatable)))
}

# Neurons that still need a root position: not in feather OR no valid position
has_root <- banc.positions %>%
  dplyr::filter(!is.na(root_position_nm), root_position_nm != "") %>%
  dplyr::pull(root_id)

to_process <- bc %>%
  dplyr::filter(!root_id %in% has_root)
# Honor test mode if active
if (exists("banc.test.ids", envir = .GlobalEnv))
  to_process <- to_process %>% dplyr::filter(root_id %in% as.character(banc.test.ids))

message(sprintf("%d neurons need root position calculation", nrow(to_process)))

if (nrow(to_process) == 0) {
  # Still save the feather in case SeaTable positions were adopted above
  if (nrow(banc.positions)) arrow::write_feather(banc.positions, feather_file)
  message("All root positions are up to date. Nothing to do.")
  return(invisible())
}

###########################
### Calculate positions ###
###########################

# Read L2 SWC files for neurons to process
nfiles <- paste0(file.path(banc.l2swc.save.path, unique(to_process$root_id)), ".swc")
nfiles <- nfiles[file.exists(nfiles)]
nfiles.ids <- gsub(".swc", "", basename(nfiles))

# Get soma/root points from nucleus positions where available
soma.positions <- to_process %>%
  dplyr::mutate(root_position_nm = ifelse(is.na(root_position_nm), nucleus_position_nm, root_position_nm)) %>%
  dplyr::distinct(root_id, supervoxel_id, nucleus_id, root_position_nm) %>%
  dplyr::filter(!is.na(root_position_nm))

if (nrow(soma.positions)) {
  soma.positions[, c("X", "Y", "Z")] <- nat::xyzmatrix(soma.positions$root_position_nm)
  soma.positions <- dplyr::distinct(soma.positions, supervoxel_id, root_id, .keep_all = TRUE)
} else {
  soma.positions <- data.frame()
}

# For neurons without nucleus, estimate from L2 skeleton leaf endpoints
banc.missing.nuclei <- subset(to_process, is.na(root_position_nm) & is.na(nucleus_position_nm))
if (nrow(banc.missing.nuclei)) {
  not.root.files <- nfiles[nfiles.ids %in% banc.missing.nuclei$root_id]
  message(sprintf("Estimating root positions for %d neurons without nuclei", length(not.root.files)))

  with_progress({
    p <- progressor(steps = length(not.root.files))
    missing.roots <- foreach(nrfile = not.root.files,
                             .combine = rbind,
                             .packages = c("nat", "bancr"),
                             .errorhandling = 'pass') %do% {
      p(sprintf("Processing ID: %s", nrfile))
      tryCatch({
        neuron <- read.neuron(nrfile)
        id <- gsub(".swc", "", basename(nrfile))
        leaves <- nat::endpoints(neuron)
        npoints1 <- nat::xyzmatrix(neuron)[leaves, ]
        if (nrow(npoints1) == 0) return(NULL)
        npoints <- npoints1

        # Exclude points inside neuropil
        pin <- nat::pointsinside(x = npoints, surf = bancr::banc_neuropil.surf)
        npoints2 <- data.frame(npoints[!pin, ])
        if (nrow(npoints2) > 0 && sum(!pin) > 2) npoints <- npoints2

        # Exclude points in neck connective
        pin <- nat::pointsinside(x = npoints, surf = bancr::banc_neck_connective.surf)
        npoints3 <- data.frame(npoints[!pin, ])
        if (nrow(npoints3) > 0 && sum(!pin) > 2) npoints <- npoints3

        if (!nrow(npoints)) return(NULL)
        if (is.null(nrow(npoints))) npoints <- matrix(npoints, ncol = 3)

        npoints <- as.data.frame(npoints)
        npoints$nucleus_id <- "0"
        npoints$root_id <- id
        npoints$root_position_nm <- apply(npoints, 1, function(x) paste(x[c("X", "Y", "Z")], collapse = ","))
        data.frame(npoints[1, ])
      }, error = function(e) {
        message(sprintf("  Error for %s: %s", basename(nrfile), e$message))
        NULL
      })
    }
  })
  soma.positions <- plyr::rbind.fill(soma.positions, missing.roots)
}

if (nrow(soma.positions) == 0) {
  message("No root positions calculated. Nothing to save.")
  return(invisible())
}

# Convert nm to raw coords
soma.positions$root_position <- apply(banc_nm2raw(xyzmatrix(soma.positions$root_position_nm)), 1, bancr:::paste_coords)
soma.positions$root_position <- gsub("\\(|\\)", "", soma.positions$root_position)

###########################
### Update feather      ###
###########################

new_data <- soma.positions %>%
  dplyr::select(root_id, root_position, root_position_nm, nucleus_id) %>%
  dplyr::distinct(root_id, .keep_all = TRUE)

# Calculate side from root position
new_data$side <- NA_character_
has_pos <- !is.na(new_data$root_position_nm) & new_data$root_position_nm != ""
if (any(has_pos)) {
  tryCatch({
    roots <- nat::xyzmatrix(new_data$root_position_nm[has_pos])
    lrdiffs <- bancr:::banc_lr_position(roots, units = "nm")
    new_data$side[has_pos] <- ifelse(lrdiffs > 0, "right", "left")
  }, error = function(e) message("  Side calculation failed: ", e$message))
}

# Merge with existing data (only own columns)
if (nrow(banc.positions) > 0) {
  unchanged <- banc.positions %>% dplyr::filter(!root_id %in% new_data$root_id)
  final <- dplyr::bind_rows(new_data, unchanged) %>%
    dplyr::distinct(root_id, .keep_all = TRUE)
} else {
  final <- new_data
}

arrow::write_feather(final, feather_file)
message(sprintf("### banc: root positions updated for %d neurons, total %d in feather [%s] ###",
                nrow(new_data), nrow(final),
                format(round(difftime(Sys.time(), t_start, units = "mins"), 1))))

###########################
### Diagnostic plot     ###
###########################

tryCatch({
  plot_df <- data.frame(
    category = c("Has root position", "Missing root position"),
    count = c(sum(!is.na(final$root_position_nm)), sum(is.na(final$root_position_nm)))
  )
  g <- ggplot2::ggplot(plot_df, ggplot2::aes(x = category, y = count, fill = category)) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label = count), vjust = -0.3) +
    ggplot2::labs(title = sprintf("Root position coverage (%s)", Sys.Date()),
                  subtitle = sprintf("%d/%d neurons (%.1f%%)",
                                     plot_df$count[1], sum(plot_df$count),
                                     100 * plot_df$count[1] / sum(plot_df$count)),
                  x = NULL, y = "Neuron count") +
    ggplot2::scale_fill_manual(values = c("Has root position" = "#0072B2",
                                           "Missing root position" = "#D55E00")) +
    ggplot2::theme_minimal()
  ggplot2::ggsave(plot = g,
                  filename = "inst/images/banc_root_position_coverage.png",
                  width = 8, height = 6, dpi = 300, bg = "transparent")
}, error = function(e) message("  Plot failed: ", e$message))

})
