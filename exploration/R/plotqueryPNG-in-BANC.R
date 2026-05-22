### -----------------------------------------------------------------
### Legacy exploratory analysis by a member of the Wilson Lab.
### Historical record only — not part of the bancpipeline release.
### -----------------------------------------------------------------

# Script to load query neurons from BANC (using full folder name) 
#   and overlay matched FAFB/MANC neurons in Neuroglancer

#############################################################
####################### Load libraries ######################
#############################################################

library(bancr)
library(tidyverse)
library(arrow)

#############################################################
####################### Load datasets #######################
#############################################################

# define location of datasets to load
fafb.meta.path <- "//research.files.med.harvard.edu/Neurobio/wilsonlab/banc/connectivity/flywire_783_meta.feather"
manc.meta.path <- "//research.files.med.harvard.edu/Neurobio/wilsonlab/banc/connectivity/manc_1.2.1_meta.feather"
banc.meta.path <- "//research.files.med.harvard.edu/Neurobio/wilsonlab/banc/connectivity/banc_meta.feather"

# Get the meta data, cell types, etc.
fafb.meta <- arrow::read_feather(fafb.meta.path)
manc.meta <- arrow::read_feather(manc.meta.path)
banc.meta <- arrow::read_feather(banc.meta.path)

######################################################################
####################### Make plotting function #######################
######################################################################

process_neuron_data <- function(input_string, fafb.meta, manc.meta) {
  # Extract cell_type
  if (identical(grep("hit_cell_type_", input_string), integer(0))) {
    cell_type2use <- strsplit(input_string, "cell_type_")[[1]]
    cell_type2use <- strsplit(cell_type2use[2], "_nblast")[[1]]
    cell_type2use <- cell_type2use[1]
  } else {
    cell_type2use <- strsplit(input_string, "hit_cell_type_")[[1]]
    cell_type2use <- cell_type2use[2]
  }

  # Extrac root ids
  split_string <- input_string
  if (grepl("_nucleus_id_", split_string, fixed = TRUE)) {
    split_string <- strsplit(split_string, "_nucleus_id_")[[1]]
    split_string <- split_string[1]
  }
  split_string <- strsplit(split_string, "root_id_")[[1]]
  split_string <- strsplit(split_string[2], "_supervoxel_id_")[[1]]
  split_string <- split_string[1]
  
  if (identical(grep("_", split_string), integer(0))) {
    root_ids <- split_string
  } else{
    split_string <- strsplit(split_string, "_")[[1]]
    root_ids <- c(split_string)
  }

  # Find neuron in fafb
  fafb.id2use <- find_neuron(cell_type2use, fafb.meta, "fafb")
  
  # Find neuron in manc
  manc.id2use <- find_neuron(cell_type2use, manc.meta, "manc")
  
  # plot all neurons
  bancsee(banc_ids = root_ids,
          fafb_ids = c(fafb.id2use$id),
          hemibrain_ids = NULL,
          manc_ids = c(manc.id2use$id),
          open = TRUE)
}

find_neuron <- function(cell_type2use, meta_data, data_source) {
  id2use <- meta_data %>%
    dplyr::filter(str_detect(cell_type2use, paste0("^", cell_type, "$")))
  
  cell_type <- id2use$cell_type
  id <- if(data_source == "fafb") id2use$root_783 else id2use$bodyid
  
  if (identical(id, character(0)) && grepl("_", cell_type2use)) {
    # try using anything before an underscore "_"
    cell_type2use_split <- strsplit(cell_type2use, "_")[[1]]
    id2use <- meta_data %>%
      dplyr::filter(str_detect(cell_type2use_split[1], paste0("^", cell_type, "$")))
    cell_type <- id2use$cell_type
    id <- if(data_source == "fafb") id2use$root_783 else id2use$bodyid
  } else if (identical(id, character(0)) && !grepl("_", cell_type2use)) {
    # try adding underscore "_"
    cell_type2use_split <- paste(c(cell_type2use, "_"), collapse = "")
    id2use <- meta_data %>%
      dplyr::filter(grepl(cell_type2use_split[1], cell_type))
    cell_type <- id2use$cell_type
    id <- if(data_source == "fafb") id2use$root_783 else id2use$bodyid
  }
  
  if (identical(id, character(0))) {
    id <- NULL
    cell_type <- NULL
  } else {
    print(paste(toupper(data_source), 'cell type'))
    print(cell_type)
  }
  
  return(list(id = id, cell_type = cell_type))
}

#############################################################
####################### Get neurons to plot #################
#############################################################

input_string <- '1_cell_type_DNg36_nblast_score_0.66_root_id_720575941481446397_720575941524828516_nucleus_id_NA_NA_hit_bodyid_25642'
process_neuron_data(input_string, fafb.meta, manc.meta) 