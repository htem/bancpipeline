#' banc-pn-by-connectivity — Name projection neurons by connectivity to known glomeruli.
#'
#' @section Reads:
#'   - SeaTable `banc_meta`, CAVE upstream connectivity
#'
#' @section Writes:
#'   - SeaTable `banc_meta`: col `cell_type` (PN-only, conditional)

###############################################################################
### BANC annotation: PN naming by connectivity to known glomeruli
###
### Aggregates per-PN connectivity, maps glomerular labels, derives 'from connectivity' names for unmatched PNs and reconciles them against existing cell types.
###
### Originally part of banc/utilities/banc-tracing.R; split out 2026-05-21
### into per-task scripts under banc/annotations/ because the original was
### never meant to be run end-to-end — each section is a one-off curation
### task. Sections may share state; if so, run them in order.
###############################################################################
source("banc/banc-startup.R")

###############################
### Aggregate by PN celltype ###
###############################

# Aggregate synapse counts per sensory neuron × PN cell type
orn_targets <- conn.top %>%
  dplyr::group_by(pre_pt_root_id, post_cell_type) %>%
  dplyr::summarise(total_count = sum(count), .groups = "drop")

########################################
### Define top targets with tolerance ###
########################################

tol <- 0.3  # 30% tolerance

orn_best <- orn_targets %>%
  dplyr::group_by(pre_pt_root_id) %>%
  dplyr::mutate(
    max_count = max(total_count, na.rm = TRUE),
    thresh    = (1 - tol) * max_count,
    is_top    = total_count >= thresh
  ) %>%
  dplyr::filter(is_top) %>%           # keep only "top band" targets
  dplyr::mutate(n_top = dplyr::n()) %>%  # how many within that band?
  dplyr::ungroup()

#####################################
### Map PNs to glomerular labels  ###
#####################################

orn_best <- orn_best %>%
  dplyr::mutate(
    glomerulus = dplyr::case_when(
      post_cell_type == "VP1l+_lvPN"              ~ "VP1l",  # HRN
      post_cell_type == "VP1d_il2PN"              ~ "VP1d",  # HRN
      post_cell_type == "VP1m+_lvPN"              ~ "VP1m",  # TRN
      post_cell_type == "VP3+_l2PN"               ~ "VP3",   # TRN
      post_cell_type == "VP5+_l2PN,VP5+VP2_l2PN"  ~ "VP5",   # HRN
      post_cell_type == "VP4_vPN"          ~ "VP4",  # HRN
      post_cell_type == "VP2+_adPN"               ~ "VP2",   # TRN
      TRUE ~ sub("_(.*)$", "", post_cell_type)
    )
  )

###########################################
### Set sensory prefix from glomerulus  ###
###########################################

orn_best <- orn_best %>%
  dplyr::mutate(
    sensory_prefix = dplyr::case_when(
      glomerulus %in% c("VP1l", "VP1d", "VP5", "VP4") ~ "HRN",
      glomerulus %in% c("VP1m", "VP3", "VP2")          ~ "TRN",
      TRUE                                             ~ "ORN"
    )
  )

#########################################
### Build "from connectivity" names   ###
#########################################

orn_names_from_conn <- orn_best %>%
  dplyr::mutate(
    conn_cell_type = dplyr::if_else(
      n_top == 1,
      paste0(sensory_prefix, "_", glomerulus),  # clear single best glomerulus
      sensory_prefix                             # ambiguous: only ORN/TRN/HRN
    )
  ) %>%
  dplyr::distinct(pre_pt_root_id, conn_cell_type, n_top)

#########################################
### Compare with existing cell types  ###
#########################################

orn_compare <- orn_names_from_conn %>%
  dplyr::left_join(
    bc %>%
      dplyr::distinct(root_id, .keep_all = TRUE) %>%
      dplyr::select(pre_pt_root_id = root_id, pre_cell_type = cell_type),
    by = "pre_pt_root_id"
  ) %>%
  # VM6 special handling (keep existing label if both are VM6 variants)
  dplyr::mutate(
    conn_cell_type = dplyr::case_when(
      grepl("VM6", pre_cell_type) & grepl("VM6", conn_cell_type) ~ pre_cell_type,
      TRUE ~ conn_cell_type
    )
  ) %>%
  dplyr::mutate(agree = pre_cell_type == conn_cell_type) %>%
  dplyr::select(pre_pt_root_id, pre_cell_type, conn_cell_type, n_top, agree)

# Example: disagreements where there is a clear top band (n_top == 1)
orn_disagreements <- orn_compare %>%
  dplyr::filter(!agree, n_top == 1)
head(orn_disagreements)

# Update seatable
bc.restored<-bc %>%
  dplyr::distinct(`_id`,root_id) %>%
  dplyr::left_join(orn_disagreements %>%
                     dplyr::select(root_id=pre_pt_root_id,
                                   cell_type=conn_cell_type),
                   by="root_id") %>%
  dplyr::filter(!is.na(cell_type))
bc.restored <-as.data.frame(bc.restored)
bc.restored[is.na(bc.restored)] <- ''
banctable_update_rows(base='banc_meta', 
                      table = "banc_meta", 
                      df = bc.restored, 
                      append_allowed = FALSE, 
                      chunksize = 1000)

