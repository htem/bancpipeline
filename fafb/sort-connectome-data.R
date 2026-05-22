#########################################
### ORGANISE CONNECTOME DATA MATCHING ###
#########################################
source("banc/banc-startup.R")
library(DBI)
library(RSQLite)
library(fs)

####################
### Flywire data ###
####################

# Read from legacy source SQLite
original_file <- "/n/data1/hms/neurobio/wilson/fafbz/hemibrainr_data/flywire_neurons/flywire_data.sqlite"
con <- dbConnect(RSQLite::SQLite(), original_file)
fw.meta.old <- dplyr::tbl(con, "flywire_meta") %>% dplyr::collect()
fw.synapses <- dplyr::tbl(con, "flywire_synapses") %>% dplyr::collect()
fw.edgelist <- dplyr::tbl(con, "flywire_edgelist") %>% dplyr::collect()
dbDisconnect(con)

# Update meta data
fw.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"flywire_meta.csv"),
                                            col_types = hemibrainr:::sql_col_types))
shared.cols <- setdiff(colnames(fw.meta.old),colnames(fw.meta))
meta.new <- dplyr::left_join(fw.meta,fw.meta.old[,c(shared.cols,"root_id")],by=c("root_783"="root_id"))

# Write as feather/parquet
arrow::write_feather(meta.new, file.path(banc.connectivity.save.path, "flywire_783_meta.feather"))
write_connectome_data(fw.synapses, file.path(banc.connectivity.save.path, "flywire_783_synapses.parquet"), format = "parquet")
arrow::write_feather(fw.edgelist, file.path(banc.connectivity.save.path, "flywire_783_edgelist.feather"))

######################
### hemibrain data ###
######################

# Read from legacy source SQLite
original_file <- "/n/data1/hms/neurobio/wilson/fafbz/hemibrainr_data/hemibrain_neurons/hemibrainr_data.sqlite"
con <- dbConnect(RSQLite::SQLite(), original_file)
hb.meta <- dplyr::tbl(con, "hemibrain_all_neurons_metrics_polypre_centrifugal_synapses") %>% dplyr::collect()
hb.synapses <- dplyr::tbl(con, "hemibrain_all_neurons_synapses_polypre_centrifugal_synapses") %>% dplyr::collect()
hb.edgelist <- dplyr::tbl(con, "hemibrain_all_neurons_edgelist_polypre_centrifugal_synapses") %>% dplyr::collect()
dbDisconnect(con)

# Write as feather/parquet
arrow::write_feather(hb.meta, file.path(banc.connectivity.save.path, "hemibrain_v.1.2.1_meta.feather"))
write_connectome_data(hb.synapses, file.path(banc.connectivity.save.path, "hemibrain_v.1.2.1_synapses.parquet"), format = "parquet")
arrow::write_feather(hb.edgelist, file.path(banc.connectivity.save.path, "hemibrain_v.1.2.1_edgelist.feather"))






