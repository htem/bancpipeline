### -----------------------------------------------------------------
### Legacy exploratory analysis by a member of the Wilson Lab.
### Historical record only — not part of the bancpipeline release.
### -----------------------------------------------------------------

# Load libraries
library(bancr)
library(tidyverse)
library(arrow)

# Data paths
# O2 path
fafb.meta.path <- '/n/data1/hms/neurobio/wilson/banc/connectivity/flywire_783_meta.feather'
fafb.edgelist.path <- '/n/data1/hms/neurobio/wilson/banc/connectivity/flywire_783_edgelist.feather'
# Files1 path
#fafb.meta.path <- "/Volumes/neurobio/wilsonlab/banc/connectivity/flywire_783_meta.feather"
#fafb.edgelist.path <- "/Volumes/neurobio/wilsonlab/banc/connectivity/flywire_783_edgelist.feather"

# Get the meta data, cell types, etc.
fw.meta <- arrow::read_feather(fafb.meta.path) %>%
  dplyr::group_by(cell_type) %>%
  dplyr::arrange(dplyr::desc(top_p)) %>%
  dplyr::mutate(top_nt=top_nt[1]) %>%
  dplyr::ungroup()

# Let's make a plot of tangential cell input onto all hDelta neurons
fw.tangential <- fw.meta %>%
  dplyr::filter(grepl("CB.FB8|CB.FB7|CB.FB9",cell_type),
                cell_class %in% c("central_complex","CX"))
fw.hdelta <- fw.meta %>%
  dplyr::filter(grepl("hDeltaD",cell_type),
                cell_class  %in% c("central_complex","CX"))

# First, collect the IDs we're interested in
tangential_ids <- fw.tangential$root_783
hdelta_ids <- fw.hdelta$root_783

# Read the edgelist
fw.elist <- arrow::read_feather(fafb.edgelist.path)

# Filter the edgelist for tangential -> hDelta connections
fw.elist.th <- fw.elist %>%
  dplyr::filter(count >= 15,
                pre %in% tangential_ids,
                post %in% hdelta_ids)

#for future analyses cutoff edge list (>10 or so)