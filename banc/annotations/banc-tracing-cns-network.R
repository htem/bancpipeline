#' banc-tracing-cns-network — CNS-cluster / CNS-network curation + snapshot restore.
#'
#' Patches CNS cluster + CNS-network assignments from tracing sheets and
#' restores prior assignments from snapshots when needed.
#'
#' @section Reads:
#'   - `BANC-project/data/cns_cluster/spectral_clustering_assignments_banc_591.csv`
#'   - SeaTable `banc_meta`
#'
#' @section Writes:
#'   - SeaTable `banc_meta`: cols `cns_cluster`, `cns_network`
#'
#' @section Notes:
#'   - Split out 2026-05-21 from `banc/utilities/banc-tracing.R`.

###############################################################################
### BANC tracing: CNS-cluster and CNS-network updates + snapshot restore
###
### Patches CNS cluster + CNS-network assignments from tracing sheets, and restores prior assignments from snapshots when needed.
###
### Originally part of banc/utilities/banc-tracing.R; split out 2026-05-21
### into per-task scripts under banc/annotations/ because the original was
### never meant to be run end-to-end — each section is a one-off curation
### task. Sections may share state; if so, run them in order.
###############################################################################
source("banc/banc-startup.R")

###################
### CNS cluster ###
###################

# Get UMAP coordinates
cns.clusters <- readr::read_csv("/Users/papers/BANC-project/data/cns_cluster/spectral_clustering_assignments_banc_591.csv", 
                            col_types = banc.col.types) 
cns.clusters <- bancr::banc_updateids(cns.clusters)
bc.current <- banctable_query("SELECT _id, root_id, supervoxel_id, super_class, cns_cluster from banc_meta")
bc.update <- bc.current %>%
  dplyr::select(-cns_cluster) %>%
  dplyr::left_join(
    cns.clusters %>%
      dplyr::mutate(spectral_cluster = as.integer(spectral_cluster)+1,
                    cns_cluster = paste0("CNS_",str_pad(spectral_cluster,width=2,pad="0"))) %>%
      dplyr::select(root_id, cns_cluster) %>%
      dplyr::arrange(cns_cluster) %>%
      dplyr::distinct(root_id, .keep_all = TRUE),
    by = "root_id"
)  %>%
  dplyr::mutate(cns_cluster = dplyr::case_when(
    super_class%in%c("glia","sensory", "sensory_ascending","motor","visceral_circulatory","not_a_neuron","visual_centrifugal","visual_projection") ~ NA,
    TRUE ~ cns_cluster)) %>%
  dplyr::select(`_id`, root_id, cns_cluster)
bc.update <- as.data.frame(bc.update)
bc.update[is.na(bc.update)] <- ''
banctable_update_rows(base = 'banc_meta', 
                      table = "banc_meta", 
                      df = bc.update, 
                      append_allowed = FALSE, 
                      chunksize = 1000)

####################
### CNS networks ###
####################

# Get UMAP coordinates
cns.clusters <- readr::read_csv("/Users/papers/BANC-project/data/cns_network/spectral_clustering_min_connection_strength_1_banc_version_626_cluster_count_13_cluster_seed_10_embedding_seed_3.csv", 
                                col_types = banc.col.types) 
cns.clusters <- bancr::banc_updateids(cns.clusters)
bc.current <- banctable_query("SELECT _id, root_id, supervoxel_id, super_class, cns_network from banc_meta")
bc.update <- bc.current %>%
  dplyr::select(-cns_network) %>%
  dplyr::left_join(
    cns.clusters %>%
      dplyr::mutate(spectral_cluster = as.integer(spectral_cluster),
                    cns_network = paste0("CNS_",str_pad(spectral_cluster,width=2,pad="0"))) %>%
      dplyr::select(root_id, cns_network) %>%
      dplyr::arrange(cns_network) %>%
      dplyr::distinct(root_id, .keep_all = TRUE),
    by = "root_id"
  )  %>%
  dplyr::mutate(cns_network = dplyr::case_when(
    super_class%in%c("glia","sensory", "sensory_ascending","motor","visceral_circulatory","not_a_neuron") ~ NA,
    TRUE ~ cns_network)) %>%
  dplyr::select(`_id`, root_id, cns_network) %>%
  dplyr::mutate(cns_network = dplyr::case_when(
    grepl("flange",cns_network) ~ "flange median bundle",
    TRUE ~ cns_network
  ))
bc.update <- as.data.frame(bc.update)
bc.update[is.na(bc.update)] <- ''
banctable_update_rows(base = 'banc_meta', 
                      table = "banc_meta", 
                      df = bc.update, 
                      append_allowed = FALSE, 
                      chunksize = 1000)

###############
### Restore ###
###############

bc.current <- banctable_query("SELECT _id, nucleus_position, nucleus_position_nm, root_id, root_position, root_position_nm, root_region, status from banc_meta")
bc.old <- readr::read_csv("/Volumes/neurobio/wilsonlab/banc/meta/snapshots/2026-02-01_12-25_banc_seatable.csv",
                          col_types = banc.col.types)
bc.restore <- bc.old %>%
  dplyr::select(`_id`, root_position, root_position_nm, root_region, status) %>%
  dplyr::mutate(root_position = dplyr::case_when(
    grepl("\\.",root_position) ~ NA,
    TRUE ~ root_position
  )) %>%
  dplyr::mutate(root_position_nm = dplyr::case_when(
    grepl("\\.",root_position_nm) ~ NA,
    TRUE ~ root_position_nm
  )) 
bc.restore$root_position <- apply(banc_nm2raw(xyzmatrix(bc.restore$root_position_nm)), 1, bancr:::paste_coords)
bc.restore$root_position <- gsub("\\(|\\)", "", bc.restore$root_position)  
bc.restore$root_position_nm <- apply(banc_raw2nm(xyzmatrix(bc.restore$root_position)), 1, bancr:::paste_coords)
bc.restore$root_position_nm <- gsub("\\(|\\)", "", bc.restore$root_position_nm)
bc.restore$root_position_nm[bc.restore$root_position_nm=="NA, NA, NA"] <- NA
bc.restore$root_position[bc.restore$root_position=="NA, NA, NA"] <- NA
bc.restore <- bc.restore %>%
  dplyr::mutate(root_region = dplyr::case_when(
    is.na(root_position) ~ NA,
    TRUE ~ root_region
  )) %>%
  dplyr::filter(
    `_id` %in% bc.current$`_id`
  )
bc.restore.2 <- bc.current %>%
  dplyr::filter(!`_id` %in% bc.restore$`_id`) %>%
  dplyr::mutate(root_position = dplyr::case_when(
    !is.na(nucleus_position) ~ NA,
    TRUE ~ nucleus_position
  )) %>%
  dplyr::mutate(root_position_nm = dplyr::case_when(
    !is.na(nucleus_position_nm) ~ NA,
    TRUE ~ nucleus_position_nm
  )) %>%
  dplyr::select(`_id`, root_position, root_position_nm, root_region, status)
bc.restore.2$root_region <- NA
bc.restored <- rbind(bc.restore,bc.restore.2)
bc.restored <-as.data.frame(bc.restored)
bc.restored[is.na(bc.restored)] <- ''
banctable_update_rows(base='banc_meta', 
                      table = "banc_meta", 
                      df = bc.restored, 
                      append_allowed = FALSE, 
                      chunksize = 1000)


