#' banc-tracing-status — Status updates + roughly-proofread sensory marking (curation).
#'
#' Per-task curation script split out from the old
#' `banc/utilities/banc-tracing.R`. Updates status from external tracing
#' sheets, applies super_class corrections, cleans up orphans, marks
#' sensory neurons roughly-proofread, and integrates the
#' `brain_sensory_rescue` gsheet.
#'
#' @section Reads:
#'   - SeaTable `banc_meta`
#'   - Google Sheets referenced by `banc.keys` (e.g.
#'     `gsheet_banc_tracing_neck_fixing`, `gsheet_brain_sensory_rescue`)
#'
#' @section Writes:
#'   - SeaTable `banc_meta`: cols `status`, `super_class`,
#'     `roughly_proofread` (block-dependent)
#'
#' @section Notes:
#'   - Manual; sections may share state, run in order. Many push calls are
#'     commented out by default.

###############################################################################
### BANC tracing: status updates + roughly-proofread sensories
###
### Updates from external tracing sheets, super_class corrections, orphan cleanup, roughly-proofread marking, and the brain_sensory_rescue gsheet integration.
###
### Originally part of banc/utilities/banc-tracing.R; split out 2026-05-21
### into per-task scripts under banc/annotations/ because the original was
### never meant to be run end-to-end — each section is a one-off curation
### task. Sections may share state; if so, run them in order.
###############################################################################
source("banc/banc-startup.R")

########################################
### Update neurons we have looked at ###
########################################

# # Read BANC meta seatable
# bancr:::banctable_updateids()
# bc <- banctable_query(sql = "select _id, root_id, super_class, status from banc_meta")
# 
# # Read data from 6-11 neck fixing task
# data <- googlesheets4::read_sheet(banc.keys$gsheet_banc_tracing_neck_fixing, skip = 3)
# tracing.done <- data %>%
#   dplyr::select(reviewed_root_id_latest, `Matching Review`, `Proofreading Status`) %>%
#   dplyr::filter(`Proofreading Status` %in% c("Done","No need for proofreading"))
# tracing.done.ids <- banc_latestid(tracing.done$reviewed_root_id_latest)
# tracing.done.ids <- setdiff(tracing.done.ids,c("0",""," "))
# tracing.done.ids <- unique(tracing.done.ids)
# glia <- data %>%
#   dplyr::select(reviewed_root_id_latest, `Matching Review`, `Proofreading Status`) %>%
#   dplyr::filter(`Proofreading Status` %in% c("Glial Cell"))
# glia.ids <- banc_latestid(glia$reviewed_root_id_latest)
# glia.ids <- setdiff(glia.ids,c("0",""," "))
# glia.ids <- unique(glia.ids)
# 
# # Take note of resolved issues
# bc.new <- bc %>%
#   dplyr::filter(root_id %in% c(tracing.done.ids,glia.ids)) %>%
#   dplyr::rowwise() %>%
#   dplyr::mutate(status = dplyr::case_when(
#     root_id %in% tracing.done.ids ~ subtract_status(status,"TRACING_ISSUE"),
#     root_id %in% glia.ids ~ subtract_status(status,"TRACING_ISSUE"),
#     TRUE ~ status
#   )) %>%
#   dplyr::mutate(status = dplyr::case_when(
#     root_id %in% tracing.done.ids ~ append_status(status,"TRACING_ISSUE_RESOLVED"),
#     root_id %in% glia.ids ~ append_status(status,"GLIA"),
#     TRUE ~ status
#   )) %>%
#   dplyr::mutate(super_class = dplyr::case_when(
#     #root_id %in% glia.ids ~ "glia",
#     TRUE ~ super_class
#   ))
# 
# # Update
# banctable_update_rows(base='banc_meta', 
#                       table = "banc_meta", 
#                       df = as.data.frame(bc.new), 
#                       append_allowed = FALSE, 
#                       chunksize = 1000)


#######################################
### Correct super_class ###
#######################################

bc <- banctable_query(sql = "select _id, root_id, super_class from banc_meta")
bc.update <- bc %>%
  dplyr::filter(super_class=="0") %>%
  dplyr::mutate(super_class="")
banctable_update_rows(base='banc_meta', 
                      table = "banc_meta", 
                      df = bc.update, 
                      append_allowed = FALSE, 
                      chunksize = 1000)

#######################################
### Request: data for NBLAST matrix ###
#######################################

# Get seatbale inventory
bc <- banctable_query(sql = "select cell_type, root_id, root_626, fafb_match, fafb_nblast_match, fafb_nblast, manc_match, manc_nblast_match, manc_nblast from banc_meta")
bc.chosen <- bc %>%
  dplyr::filter((!is.na(fafb_match)&(fafb_match==fafb_nblast_match))|
                (!is.na(manc_match)&(manc_match==manc_nblast_match))) %>%
  dplyr::select(cell_type, root_id, root_626, fafb_match, fafb_nblast, manc_match, manc_nblast)
banc.tracing.save.path <- "tracing/tracing_issues/"
current_datetime <- Sys.time()
datetime_string <- format(current_datetime, "%Y-%m-%d")
readr::write_csv(bc.chosen, file.path(banc.tracing.save.path,paste0(datetime_string,"_banc_to_fafb_and_manc_best_matches.csv")))

#####################################
### Remove old orphan annotations ###
#####################################

# Get seatbale inventory
bc <- banctable_query(sql = "select _id, root_id, position, nucleus_id, proofread, roughly_proofread, status, region, super_class, cell_class, cell_type, l2_nodes from banc_meta")
bc.valid <- bc %>%
  dplyr::filter(grepl("TRACING_ISSUE_TADPOLE", status),
                l2_nodes>50) %>%
  dplyr::distinct(root_id, .keep_all = TRUE) %>%
  dplyr::arrange(l2_nodes)

banctable_update_rows(base='banc_meta', 
                      table = "banc_meta", 
                      df = bc.valid %>%
                        dplyr::mutate(status = subtract_status(status,"TRACING_ISSUE_TADPOLE")) %>%
                        dplyr::distinct(`_id`,root_id,status), 
                      append_allowed = FALSE, 
                      chunksize = 1000)

bc.invalid <- bc %>%
  dplyr::filter(grepl("TRACING_ISSUE_TADPOLE", status),
                l2_nodes<50,
                !is.na(cell_type),
                !is.na(nucleus_id)) %>%
  dplyr::distinct(root_id, .keep_all = TRUE) %>%
  dplyr::arrange(l2_nodes)

banctable_update_rows(base='banc_meta', 
                      table = "banc_meta", 
                      df = bc.invalid %>%
                        dplyr::mutate(status = append_status(status,"FAFB_PNG_MATCH_WRONG, MANC_PNG_MATCH_WRONG, HEMIBRAIN_PNG_MATCH_WRONG")) %>%
                        dplyr::distinct(`_id`,root_id,status), 
                      append_allowed = FALSE, 
                      chunksize = 1000)

#################################
### Marking roughly proofread ###
#################################

# Get seatbale inventory
bc <- banctable_query(sql = "select _id, root_id, position, proofread, roughly_proofread, status, region, super_class, cell_class, cell_type, l2_nodes from banc_meta")
bc.rough <- bc %>%
  dplyr::filter(proofread != "TRUE",
                !grepl("NOT_A_NEURON|DEBRIS|GLIA|TRACHEA|MERGE_MONSTER|MERGE", status),
                !grepl("glia|trachea|not_a_neuron|unknown",super_class),
                !grepl("glia|trachea|not_a_neuron|unknown|astrocyte|nerve",cell_type),
                !is.na(cell_type),
                position!="0",
                roughly_proofread != "TRUE",
                !is.na(position),
                l2_nodes>50) %>%
  dplyr::distinct(root_id, .keep_all = TRUE) %>%
  dplyr::arrange(l2_nodes)

# Annotate 
update <- bancr:::banc_annotate_proofreading_notes(positions = nat::xyzmatrix(bc.rough$position), 
                                                   user_id = 355,
                                                   units = "raw",
                                                   label = "roughly proofread", 
                                                   datastack_name = "brain_and_nerve_cord",
                                                   use_admin_creds = TRUE)
banctable_update_rows(base='banc_meta', 
                      table = "banc_meta", 
                      df = bc.rough %>%
                        dplyr::mutate(roughly_proofread="TRUE", status = subtract_status(status,"TRACING_ISSUE_TADPOLE")) %>%
                        dplyr::distinct(`_id`,root_id,roughly_proofread,status), 
                      append_allowed = FALSE, 
                      chunksize = 1000)

###############################################################
### Update roughly proofread and large fragment annotations ###
###############################################################

# Get proofreading notes
banc.proofreading.notes <- banc_proofreading_notes(live=2)

banc.proofreading.rough <- banc.proofreading.notes %>%
  dplyr::filter(tag %in% c("roughly proofread","large fragment", "arbor damaged")) %>%
  dplyr::mutate(pt_root_id = as.character(pt_root_id)) %>%
  dplyr::pull(pt_root_id)
rough <- readr::read_csv(file = "data/others/roughly_proofread_from_princeton_2025.csv", 
                         col_types = banc.col.types)
rough$root_id_updates <- bancr::banc_latestid(rough$root_id)
rough.df <- rough %>%
  dplyr::filter(!root_id_updates %in% bc$root_id,
                !root_id_updates %in%banc.proofreading.rough) %>%
  dplyr::distinct(root_id_updates, .keep_all = TRUE)
rough.df$position <- NA
n <- nrow(rough.df)
pb <- txtProgressBar(min = 0, max = n, style = 3)
for (i in seq_len(n)) {
  id <- rough.df$root_id_updates[i]
  res <- try({
    l2   <- bancr::banc_read_l2skel(id, rawcoords = TRUE)
    strh <- nat::strahler_order(l2[[1]])
    pos  <- nat::xyzmatrix(l2[[1]]$d)[which.max(strh$points), ]
    pos.raw <- nat:::xyzmatrix2str(bancr:::banc_nm2raw(pos))
    pos.raw
  }, silent = TRUE)
  if (!inherits(res, "try-error")) {
    rough.df$position[i] <- res
  } else {
    rough.df$position[i] <- NA
  }
  setTxtProgressBar(pb, i)
}
close(pb)

# roughly proofread
rough.positions.rough <- rough.df %>%
  dplyr::filter(tag=="roughly_proofread", !is.na(position)) %>%
  dplyr::distinct(position) %>%
  dplyr::pull(position)
rough.positions.damaged <- rough.df %>%
  dplyr::filter(tag=="large_fragment", !is.na(position)) %>%
  dplyr::distinct(position) %>%
  dplyr::pull(position)

# Annotate 
update <- bancr:::banc_annotate_proofreading_notes(positions = nat::xyzmatrix(rough.positions.rough), 
                                                   user_id = 355,
                                                   units = "raw",
                                                   label = "roughly proofread", 
                                                   datastack_name = "brain_and_nerve_cord",
                                                   use_admin_creds = TRUE)
update <- bancr:::banc_annotate_proofreading_notes(positions = nat::xyzmatrix(rough.positions.damaged), 
                                                   user_id = 355,
                                                   units = "raw",
                                                   label = "arbor is damaged", 
                                                   datastack_name = "brain_and_nerve_cord",
                                                   use_admin_creds = TRUE)
# bancr::banctable_append_rows(base='banc_meta',
#                              table = "banc_meta",
#                              df = rough.df %>%
#                                dplyr::distinct(root_id=root_id_updates,position) %>%
#                                dplyr::mutate(roughly_proofread="TRUE",
#                                              proofread="FALSE"),
#                              chunksize = 1000)  

##########################################
### Update roughly proofread sensories ###
##########################################

# Get proofreading notes
banc.proofreading.notes <- banc_proofreading_notes(live=2)
bc <- banctable_query(sql = "SELECT root_id, position, proofread, roughly_proofread, status from banc_meta")
banc.proofreading.rough <- banc.proofreading.notes %>%
  dplyr::filter(tag %in% c("roughly proofread","large fragment", "arbor damaged")) %>%
  dplyr::mutate(pt_root_id = as.character(pt_root_id)) %>%
  dplyr::pull(pt_root_id)
gsheet_id <- banc.keys$gsheet_banc_tracing_brain_sensory_rescue
if (is.null(gsheet_id) || !nzchar(gsheet_id))
  stop("banc.keys$gsheet_banc_tracing_brain_sensory_rescue not set. ",
       "Add the row to data/private/keys.csv — see CLAUDE.md for the canonical name.")
url <- sprintf("https://docs.google.com/spreadsheets/d/%s/edit", gsheet_id)
rough <- googlesheets4::read_sheet(url, sheet = "brain_sensory_rescue")
rough <- rough %>%
  dplyr::filter(!is.na(correctedPosition)) %>%
  dplyr::select(root_id,position=correctedPosition,is_sensory,tag,notes)
rough$root_id <- bancr::banc_latestid(rough$root_id)
rough.df <- rough %>%
  dplyr::filter(is_sensory == "TRUE",!is.na(root_id),root_id!="0") %>%
  dplyr::filter(!root_id %in% bc$root_id,
                !root_id %in%banc.proofreading.rough) %>%
  dplyr::distinct(root_id, .keep_all = TRUE) %>%
  dplyr::mutate(cell_type = ifelse(grepl("OSN",tag),"ORN",NA),
                super_class="sensory")

# roughly proofread
rough.positions.rough <- rough.df %>%
  dplyr::filter(!grepl("contralateral_missing|damaged",tag), !is.na(position)) %>%
  dplyr::distinct(position) %>%
  dplyr::pull(position)
rough.positions.damaged <- rough.df %>%
  dplyr::filter(grepl("contralateral_missing|damaged",tag), !is.na(position)) %>%
  dplyr::distinct(position) %>%
  dplyr::pull(position)

# Annotate 
update <- bancr:::banc_annotate_proofreading_notes(positions = nat::xyzmatrix(rough.positions.rough), 
                                                   user_id = 355,
                                                   units = "raw",
                                                   label = "roughly proofread", 
                                                   datastack_name = "brain_and_nerve_cord",
                                                   use_admin_creds = TRUE)
update <- bancr:::banc_annotate_proofreading_notes(positions = nat::xyzmatrix(rough.positions.damaged), 
                                                   user_id = 355,
                                                   units = "raw",
                                                   label = "arbor is damaged", 
                                                   datastack_name = "brain_and_nerve_cord",
                                                   use_admin_creds = TRUE)
bancr::banctable_append_rows(base='banc_meta',
                             table = "banc_meta",
                             df = rough.df %>%
                               dplyr::distinct(root_id,position,root_position=position,cell_type,super_class) %>%
                               dplyr::mutate(roughly_proofread="TRUE",
                                             proofread="FALSE",
                                             status="SENSORY_SEARCH_2"),
                             chunksize = 1000)

