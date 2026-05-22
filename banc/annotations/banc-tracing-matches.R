#' banc-tracing-matches — Cross-dataset match updates + Hampel-lab integration (curation).
#'
#' Pulls match decisions from the Princeton gsheet (key
#' `gsheet_princeton_matches`) and integrates the Hampel-lab tracing sheet.
#' Originally part of `banc/utilities/banc-tracing.R`; split out 2026-05-21
#' into per-task scripts. Sections may share state — run in order.
#'
#' @section Reads:
#'   - Google Sheet referenced by `banc.keys$gsheet_princeton_matches`
#'   - Hampel-lab tracing sheet (via `banc.keys`)
#'   - SeaTable `banc_meta`
#'
#' @section Writes:
#'   - SeaTable `banc_meta`: match columns (cell-block dependent)
#'
#' @section Notes:
#'   - Manual / curation script. SeaTable push calls live inside the
#'     individual blocks; many are commented out by default.
#'   - Requires `data/private/keys.csv` populated.

###############################################################################
### BANC tracing: cross-dataset match updates + Hampel-lab integration
###
### Pulls match decisions from the Princeton sheet (gsheet_princeton_matches in keys.csv) and integrates the Hampel-lab tracing sheet.
###
### Originally part of banc/utilities/banc-tracing.R; split out 2026-05-21
### into per-task scripts under banc/annotations/ because the original was
### never meant to be run end-to-end — each section is a one-off curation
### task. Sections may share state; if so, run them in order.
###############################################################################
source("banc/banc-startup.R")


########################
#### Update matches ####
########################

# # Read princeton google sheet
# data <- googlesheets4::read_sheet(banc.keys$gsheet_princeton_matches)
# data <- subset(data, result == "good_match")
# data$banc_nblast_match <- bancr::banc_updateids(data$banc_nblast_match )
# data$manc_id <- as.character(data$manc_id)
# 
# # Read BANC meta seatable
# bancr:::banctable_updateids()
# bc <- banctable_query(sql = "select _id, root_id, cell_type_source, notes, manc_match from banc_meta")
# 
# # Update
# bc.new <- data %>%
#   dplyr::rename(notes_princeton = notes) %>%
#   dplyr::left_join(bc, by = c("banc_nblast_match"="root_id")) %>%
#   dplyr::mutate(cell_type_source = dplyr::case_when(
#     is.na(manc_match) ~ "Princeton",
#     TRUE ~ cell_type_source
#   )) %>%
#   dplyr::mutate(manc_match = dplyr::case_when(
#     is.na(manc_match) ~ manc_id,
#     TRUE ~ manc_match
#   )) %>%
#   dplyr::rowwise() %>%
#   dplyr::mutate(notes = dplyr::case_when(
#     !is.na(notes_princeton) ~ paste(na.omit(c(notes,notes_princeton)), collapse = "; "),
#     TRUE ~ notes
#   )) %>%
#   dplyr::distinct(`_id`, root_id=banc_nblast_match, cell_type_source, notes, manc_match)
# 
# # Update
# banctable_update_rows(base='banc_meta', 
#                       table = "banc_meta", 
#                       df = as.data.frame(bc.new), 
#                       append_allowed = FALSE, 
#                       chunksize = 1000)
#   

##################################
### Integrate hampel lab sheet ###
##################################

# Get data
data <- readxl::read_xlsx("data/others/from_hampel_for_upload_october.xlsx") %>%
  dplyr::select(root_id, cell_function, other_names, side, hemilineage, position, cell_type_source, nerve)

# Update
data <- banc_updateids(data)

# Get BANC data
bc <- banctable_query(sql = "select _id, root_id, cell_function, other_names, side, hemilineage, position, cell_type_source, nerve from banc_meta")
bc[bc=="NA"] <- ""

# Fixes
fix_cell_type_source <- function(cell_type_source){
  a = sort(unique(unlist(strsplit(tolower(cell_type_source),split=","))))
  a = a[!a%in%c("",",")]
  a = gsub("^ ","",a)
  a[a=="seeds/hampel lab,seeds/hampel lab"] <- "seeds/hampel lab"
  a[a=="seeds/hampel lab,seeds/hampel lab,wilson lab"] <- "wilson lab"
  paste(a,collapse=",")
}

# Create update
bc.new <- data %>%
  dplyr::distinct(root_id, .keep_all = TRUE) %>%
  dplyr::left_join(bc, by = "root_id") %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    position = dplyr::case_when(
      is.na(position.y) ~ position.x,
      TRUE ~ position.y
    ),
    side = dplyr::case_when(
      is.na(side.y) ~ side.x,
      TRUE ~ side.x
    ),
    cell_function = dplyr::case_when(
        is.na(cell_function.y)|grepl("auto",cell_function.y) ~ cell_function.x,
        is.na(cell_function.x) ~ cell_function.y,
        !is.na(cell_function.y) ~ paste(unique(c(cell_function.y,cell_function.x)),collapse=", "),
        TRUE ~ cell_function.x
      ),
    hemilineage = dplyr::case_when(
      is.na(hemilineage.y)|grepl("auto",hemilineage.y) ~ hemilineage.x,
      is.na(hemilineage.x) ~ hemilineage.y,
      !is.na(hemilineage.y) ~ paste(unique(c(hemilineage.y,hemilineage.x)),collapse=", "),
      TRUE ~ hemilineage.x
    ),
    nerve = dplyr::case_when(
      is.na(nerve.y)|grepl("auto",nerve.y) ~ nerve.x,
      is.na(nerve.x) ~ nerve.y,
      !is.na(nerve.y) ~ paste(unique(c(nerve.y,nerve.x)),collapse=", "),
      TRUE ~ nerve.x
    ),
    cell_type_source = dplyr::case_when(
      is.na(cell_type_source.y)|grepl("auto",cell_type_source.y) ~ cell_type_source.x,
      is.na(cell_type_source.x) ~ cell_type_source.y,
      !is.na(cell_type_source.y) ~ paste(unique(c(cell_type_source.y,cell_type_source.x)),collapse=", "),
      TRUE ~ cell_type_source.x
    ),
    other_names = dplyr::case_when(
      is.na(other_names.y)|grepl("auto",other_names.y) ~ other_names.x,
      is.na(other_names.x) ~ other_names.y,
      !is.na(other_names.y) ~ paste(unique(c(other_names.y,other_names.x)),collapse=", "),
      TRUE ~ other_names.x
    ),
    cell_type_source = tolower(cell_type_source)) %>%
  dplyr::ungroup() %>%
  dplyr::select(-contains(".x"), -contains(".y")) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(cell_type_source = fix_cell_type_source(cell_type_source))

# Update
bc.update <- subset(bc.new, !is.na(bc.new$`_id`))
bc.nrows <- subset(bc.new, is.na(bc.new$`_id`))
bc.nrows$`_id` <- NULL
bc.nrows[is.na(bc.nrows)] <- ''
bc.update[is.na(bc.update)] <- ''
banctable_update_rows(base='banc_meta', 
                      table = "banc_meta", 
                      df = as.data.frame(bc.update)[,c("root_id",
                                                       "_id", 
                                                       "position",
                                                       "hemilineage", 
                                                       "cell_type_source", 
                                                       "other_names"
                                                       )], 
                      append_allowed = FALSE, 
                      chunksize = 1000)
if(nrow(bc.nrows)){
  banctable_append_rows(base='banc_meta',
                        table = "banc_meta",
                        df = bc.nrows,
                        chunksize = 1000)  
}

# Add NGL links to Hampel interest neurons
bc.hampel.orig <- banctable_query(sql = "select _id, nucleus_id, status, root_id, banc_match, fafb_match, ngl_link, cell_type_source, other_names, cell_type, hemilineage from banc_meta")
bc.hampel <- bc.hampel.orig %>%
  dplyr::filter(!grepl("DELETE|GLIA|NOT_A_NEURON",status),
                grepl('hampel|seeds',cell_type_source)|grepl("BM_",cell_type)|grepl("LB23",hemilineage)) %>%
  dplyr::mutate(update = TRUE) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(ngl_link = dplyr::case_when(
    !is.na(banc_match) ~ bancsee(banc_ids = c(root_id,banc_match),
                                 fafb_ids = c(fafb_match),
                                 nuclei_ids = c(nucleus_id)),
    is.na(banc_match) ~ bancsee(banc_ids = c(root_id),
                                fafb_ids = c(fafb_match),
                                nuclei_ids = c(nucleus_id)),
    TRUE ~ ngl_link
  )) %>%
  dplyr::distinct(`_id`, root_id, ngl_link) %>%
  dplyr::filter(!is.na(ngl_link))
bc.update <- bc.hampel
bc.update[is.na(bc.update)] <- ''
banctable_update_rows(base='banc_meta', 
                      table = "banc_meta", 
                      df = as.data.frame(bc.update), 
                      append_allowed = FALSE, 
                      chunksize = 1000)

# Get BANC data
bc.hampel.orig <- banctable_query(sql = "select _id, root_id, supervoxel_id, position, other_names, cell_type_source from banc_meta")
bc.hampel <- bc.hampel.orig %>%
  dplyr::filter(grepl('hampel|seeds',cell_type_source),
                !is.na(other_names),
                !grepl("\\?",other_names)) %>%
  dplyr::select(pt_root_id = root_id, pt_supervoxel_id = supervoxel_id, pt_position = position, tag = other_names) %>%
  dplyr::mutate(tag2 = "neuron identity", user_id = '125') %>%
  separate_rows(tag, sep = ",\\s*") %>%
  dplyr::filter(!grepl("\\?|lab",tag)) %>%
  dplyr::mutate(tag = gsub("^ ","",tag))

# Save
readr::write_csv(bc.hampel, file.path(banc.meta.save.path,"seeds_hampel_lab_annotations.csv"))
system(sprintf("gsutil cp %s/%s gs://lee-lab_brain-and-nerve-cord-fly-connectome/meta", banc.meta.save.path, "seeds_hampel_lab_annotations.csv"))


