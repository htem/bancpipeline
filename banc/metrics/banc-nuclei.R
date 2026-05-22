#' banc-nuclei — CAVE nuclei → nm coords + L/R side + containing neuropil.
#'
#' Feeds `root_position` fallbacks in `banc-calculate-root-positions.R`.
#'
#' @section Reads:
#'   - CAVE `banc_nuclei()` (both tables)
#'
#' @section Writes:
#'   - `banc_nuclei.csv` / `.feather`

######################
### Process nuclei ###
######################
source("banc/banc-startup.R")

# Direct us to the BANC dataset
bancr::choose_banc()

# Get nuclei positions
banc.nuclei <- banc_nuclei(rawcoords = TRUE, table = "both")

# Prepare all data frames
banc.nuclei.mod <- banc.nuclei %>%
  dplyr::mutate(nucleus_position_nm = nat::xyzmatrix2str(bancr::banc_raw2nm(nucleus_position),
                                                          format = "%.0f, %.0f, %.0f")) %>%
  dplyr::mutate(nucleus_id = as.character(nucleus_id)) %>%
  dplyr::ungroup() %>%
  dplyr::rename(nucleus_supervoxel_id = pt_supervoxel_id) %>%
  dplyr::select(nucleus_id, nucleus_supervoxel_id, root_id, nucleus_position, nucleus_position_nm) %>%
  dplyr::distinct(nucleus_id, .keep_all = TRUE) %>%
  dplyr::select(root_id, nucleus_id, nucleus_supervoxel_id, nucleus_position, nucleus_position_nm)

# Determine whether neuron is on the left or right of BANC
nuclei.points <- nat::xyzmatrix(banc.nuclei.mod$nucleus_position_nm)
lrdiffs <- bancr:::banc_lr_position(nuclei.points,units = "nm")
sides <- ifelse(lrdiffs>0,"right","left")
banc.nuclei.mod$side <- sides

# Determine the neuropil region for the ID
banc.nuclei.mod[,c("x","y","z")] <- nat::xyzmatrix(banc.nuclei.mod$nucleus_position_nm)
banc.nuclei.mod$region <- NA
banc.nuclei.mod$neuropil <- NA
nuclei.regions <- pointsnearby_banc(banc.nuclei.mod, id = "nucleus_id")
nuclei.regions <- nuclei.regions %>%
  dplyr::mutate(neuropil = gsub("outside_|_L$|_R$","",neuropil)) %>%
  dplyr::select(-x, -y, -z, -region, -nucleus_position_nm) %>%
  dplyr::rename(pt_root_id=root_id, 
                pt_supervoxel_id=nucleus_supervoxel_id,
                pt_position=nucleus_position)

# Save
readr::write_csv(nuclei.regions, file.path(banc.meta.save.path,"banc_nuclei_with_near_neuropils.csv"))
system(sprintf("gsutil cp %s/%s gs://lee-lab_brain-and-nerve-cord-fly-connectome/meta", banc.meta.save.path, "banc_nuclei_with_near_neuropils.csv"))

