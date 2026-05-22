#' banc-calculate-regions — Assign region of innervation per BANC neuron.
#'
#' Tests L2-skeleton neuropil containment to assign each neuron one of
#' `ventral_nerve_cord`, `central_brain`, `optic_lobe`, `brain`, `rind`.
#' Incremental.
#'
#' @section Reads:
#'   - SeaTable `banc_meta`, L2 SWCs
#'
#' @section Writes:
#'   - `banc_regions.feather`
#'
#' @section Invoked by:
#'   production v888 rebuild chain.

###########################################################
### Calculate region of innervation for BANC neurons
###
### Reads L2 SWC files and checks neuropil containment
### to assign one of:
###   ventral_nerve_cord, central_brain, optic_lobe, brain, rind
###
### Saves results to banc_regions.feather
###########################################################
source("banc/banc-startup.R")

local({

message("### banc: calculating regions of innervation ###")
t_start <- Sys.time()
bancr::choose_banc()

###########################
### Read current state  ###
###########################

bc <- banctable_query("SELECT _id, status, root_id, region, root_region from banc_meta") %>%
  banc_filter_neurons() %>%
  dplyr::mutate(root_id = as.character(root_id))

# Read existing regions feather
feather_file <- file.path(banc.save.path, "banc_regions.feather")
if (file.exists(feather_file)) {
  banc.regions <- arrow::read_feather(feather_file) %>%
    dplyr::mutate(root_id = as.character(root_id)) %>%
    dplyr::distinct(root_id, .keep_all = TRUE)
} else {
  banc.regions <- data.frame(root_id = character(0), region = character(0))
}

# Determine which neurons need region calculation
# Missing region in both seatable and feather
has_region_feather <- banc.regions %>%
  dplyr::filter(!is.na(region) & !region %in% c("", "0")) %>%
  dplyr::pull(root_id)
has_region_st <- bc %>%
  dplyr::filter(!is.na(region) & !region %in% c("", "0")) %>%
  dplyr::pull(root_id)
has_region <- unique(c(has_region_feather, has_region_st))

to_process <- bc %>%
  dplyr::filter(!root_id %in% has_region)
# Honor test mode if active
if (exists("banc.test.ids", envir = .GlobalEnv))
  to_process <- to_process %>% dplyr::filter(root_id %in% as.character(banc.test.ids))

message(sprintf("%d neurons need region calculation", nrow(to_process)))

if (nrow(to_process) == 0) {
  # Even when nothing to compute, backfill the feather with seatable regions
  # so the combined feather has complete data
  st_regions <- bc %>%
    dplyr::filter(!is.na(region) & !region %in% c("", "0") & !root_id %in% has_region_feather) %>%
    dplyr::select(root_id, region) %>%
    dplyr::mutate(region = dplyr::case_when(
      region == "vnc" ~ "ventral_nerve_cord",
      region == "midbrain" ~ "central_brain",
      region == "optic" ~ "optic_lobe",
      TRUE ~ region
    ))
  if (nrow(st_regions) > 0) {
    message(sprintf("Backfilling %d seatable regions into feather", nrow(st_regions)))
    # Put st_regions first so seatable values override existing NA values
    final <- dplyr::bind_rows(st_regions, banc.regions) %>%
      dplyr::distinct(root_id, .keep_all = TRUE)
    arrow::write_feather(final, feather_file)
  } else {
    message("All regions are up to date. Nothing to do.")
  }
  return(invisible())
}

###########################
### Ensure L2 SWCs exist ##
###########################
# Any to_process neuron missing an L2 SWC blocks region calculation. Rather
# than silently skipping (the long-standing 35-day-stale-feather bug), invoke
# banc-l2.R for the missing IDs. banc-l2.R already supports a `banc.test.ids`
# global to scope which root_ids to fetch — set that before sourcing so the
# subroutine only does the work we need.

.swc_paths_for <- function(rids) {
  paste0(file.path(banc.l2swc.save.path, unique(rids)), ".swc")
}

.initial_nfiles <- .swc_paths_for(to_process$root_id)
.missing_swc_ids <- to_process$root_id[!file.exists(.initial_nfiles)]

if (length(.missing_swc_ids) > 0) {
  message(sprintf("  %d / %d to_process neurons missing L2 SWC — invoking banc-l2.R...",
                  length(.missing_swc_ids), nrow(to_process)))
  # banc-l2.R uses banctable_query_cached() — clear cache so newly-added
  # SeaTable rows are visible.
  if (exists("banctable_cache_clear")) {
    tryCatch(banctable_cache_clear(), error = function(e) NULL)
  }
  assign("banc.test.ids", .missing_swc_ids, envir = .GlobalEnv)
  on.exit({
    if (exists("banc.test.ids", envir = .GlobalEnv))
      rm("banc.test.ids", envir = .GlobalEnv)
  }, add = TRUE)
  # Source into an isolated env so banc-l2.R's locals (nfiles, banc.ids, ...)
  # don't clobber ours. Side effects we care about (SWC files written to
  # banc.l2swc.save.path) are unaffected by sourcing scope.
  tryCatch(
    source("banc/metrics/banc-l2.R",
           local = new.env(parent = .GlobalEnv)),
    error = function(e) message(sprintf("  banc-l2.R subroutine failed: %s", e$message))
  )
}

###########################
### Calculate regions   ###
###########################

# Re-scan SWC presence after banc-l2.R run
nfiles <- .swc_paths_for(to_process$root_id)
nfiles <- nfiles[file.exists(nfiles)]
.still_missing <- setdiff(to_process$root_id,
                          gsub("\\.swc$", "", basename(nfiles)))

if (length(.still_missing) > 0) {
  message(sprintf("  WARN: %d to_process neurons still lack L2 SWC after fetch (likely no CAVE L2 representation):",
                  length(.still_missing)))
  .show_n <- min(10L, length(.still_missing))
  message(sprintf("    %s%s",
                  paste(head(.still_missing, .show_n), collapse = ", "),
                  if (length(.still_missing) > .show_n) sprintf(", ... (%d more)",
                                                                 length(.still_missing) - .show_n) else ""))
}

if (length(nfiles) == 0) {
  message("No L2 SWC files available for region calculation; running SeaTable -> feather backfill anyway and exiting.")
  # Fall through to the backfill branch (mirrors lines 56-76 above) so any
  # manually-set SeaTable regions still land in the feather even when no
  # new SWCs were available this run.
  st_regions <- bc %>%
    dplyr::filter(!is.na(region) & !region %in% c("", "0") &
                    !root_id %in% has_region_feather) %>%
    dplyr::select(root_id, region) %>%
    dplyr::mutate(region = dplyr::case_when(
      region == "vnc" ~ "ventral_nerve_cord",
      region == "midbrain" ~ "central_brain",
      region == "optic" ~ "optic_lobe",
      TRUE ~ region
    ))
  if (nrow(st_regions) > 0) {
    message(sprintf("  Backfilling %d SeaTable regions into feather", nrow(st_regions)))
    final <- dplyr::bind_rows(st_regions, banc.regions) %>%
      dplyr::distinct(root_id, .keep_all = TRUE)
    arrow::write_feather(final, feather_file)
  }
  return(invisible())
}

# Load surfaces
optic_lobes <- subset(banc_brain_neuropils.surf, "optic")
central_brain <- subset(banc_brain_neuropils.surf, "midbrain|central_brain")

message(sprintf("Calculating regions for %d neurons", length(nfiles)))
with_progress({
  p <- progressor(steps = length(nfiles))
  results <- foreach(nfile = nfiles,
                     .combine = rbind,
                     .packages = c("nat", "bancr"),
                     .errorhandling = 'pass') %do% {
    tryCatch({
      p(sprintf("Processing ID: %s", nfile))
      neuron <- read.neuron(nfile)
      id <- gsub(".swc", "", basename(nfile))
      points <- nat::xyzmatrix(neuron)

      in.vnc <- any(pointsinside(x = points, surf = banc_vnc_neuropil.surf))
      in.brain <- any(pointsinside(x = points, surf = central_brain))
      in.optic_lobes <- any(pointsinside(x = points, surf = optic_lobes))
      in.neck <- any(pointsinside(x = points, surf = banc_neck_connective.surf))
      in.greater.brain <- any(pointsinside(x = points, surf = banc_brain_neuropil.surf))
      in.cortex <- any(pointsinside(x = points, surf = banc.surf))

      # "neck_connective" is no longer a valid region value — ascending/
      # descending neurons are identified via super_class and the cervical
      # connective lives as a separate "tract" entry in SeaTable.
      # Neurons that span brain+VNC fall through to whichever volume their
      # soma sits in (or NA if neither).
      if (in.vnc) {
        region <- "ventral_nerve_cord"
      } else if (in.optic_lobes) {
        region <- "optic_lobe"
      } else if (in.brain) {
        region <- "central_brain"
      } else if (in.greater.brain) {
        region <- "brain"
      } else if (in.cortex) {
        region <- "rind"
      } else {
        region <- NA
      }
      data.frame(root_id = id, region = region)
    },
    error = function(e) {
      message(sprintf("  Error for %s: %s", basename(nfile), e$message))
      data.frame(root_id = gsub(".swc", "", basename(nfile)), region = NA)
    })
  }
})

results <- results %>%
  dplyr::filter(!is.na(region)) %>%
  dplyr::distinct(root_id, .keep_all = TRUE)

if (nrow(results) == 0) {
  message("No regions calculated. Nothing to save.")
  return(invisible())
}

###########################
### Update feather      ###
###########################

# Also pull in any seatable regions for neurons not in our results
# (seatable regions are manually curated and should be preserved)
st_regions <- bc %>%
  dplyr::filter(!is.na(region) & !region %in% c("", "0") & !root_id %in% results$root_id) %>%
  dplyr::select(root_id, region) %>%
  # Normalize legacy names
  dplyr::mutate(region = dplyr::case_when(
    region == "vnc" ~ "ventral_nerve_cord",
    region == "midbrain" ~ "central_brain",
    region == "optic" ~ "optic_lobe",
    TRUE ~ region
  ))
new_data <- dplyr::bind_rows(results, st_regions)

# Merge with existing feather (only root_id + region columns)
if (nrow(banc.regions) > 0) {
  # Update region for existing root_ids
  banc.regions$region[match(new_data$root_id, banc.regions$root_id, nomatch = 0)] <-
    new_data$region[match(banc.regions$root_id, new_data$root_id, nomatch = 0) > 0]

  # Add new root_ids not yet in feather
  new_rids <- new_data %>% dplyr::filter(!root_id %in% banc.regions$root_id)
  if (nrow(new_rids) > 0) {
    banc.regions <- dplyr::bind_rows(banc.regions, new_rids)
  }
  final <- banc.regions %>% dplyr::distinct(root_id, .keep_all = TRUE)
} else {
  final <- new_data
}

arrow::write_feather(final, feather_file)
message(sprintf("### banc: regions updated for %d neurons [%s] ###",
                nrow(results),
                format(round(difftime(Sys.time(), t_start, units = "mins"), 1))))

# Summary of what was processed vs skipped this run
.computed <- nrow(results)
.skipped_no_swc <- if (exists(".still_missing")) length(.still_missing) else 0L
message(sprintf("  Summary: %d to_process | %d processed (region computed) | %d skipped (no L2 SWC after fetch)",
                nrow(to_process), .computed, .skipped_no_swc))

###########################
### Diagnostic plot     ###
###########################

tryCatch({
  region_counts <- final %>%
    dplyr::mutate(region = ifelse(is.na(region) | region == "", "unassigned", region)) %>%
    dplyr::count(region) %>%
    dplyr::arrange(dplyr::desc(n))
  g <- ggplot2::ggplot(region_counts, ggplot2::aes(x = stats::reorder(region, n), y = n, fill = region)) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label = n), hjust = -0.2, size = 3) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = sprintf("Region distribution (%s)", Sys.Date()),
                  subtitle = sprintf("%d neurons total", nrow(final)),
                  x = NULL, y = "Neuron count") +
    ggplot2::theme_minimal()
  ggplot2::ggsave(plot = g,
                  filename = "inst/images/banc_region_distribution.png",
                  width = 10, height = 6, dpi = 300, bg = "transparent")
}, error = function(e) message("  Plot failed: ", e$message))

})
