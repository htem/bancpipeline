######################
### Assign regions ###
######################
source("banc/banc-startup.R")

#################################
### Read synapses by neuropil ###
#################################

# Get BANC synapses
banc.syns.nps <- arrow::read_parquet(
  file.path(banc.connectivity.save.path, "banc_synapses.parquet")) %>%
  dplyr::filter(post_status!="glia",
                pre_status!="glia"
  ) %>%
  dplyr::select(id, neuropil, region, side)

##############################
### Write to synapse files ###
##############################



#############################
### Calculate SEZ neurons ###
#############################


###################################
### Calculate abdominal neurons ###
###################################





