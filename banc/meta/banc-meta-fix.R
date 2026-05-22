#' banc-meta-fix — One-off spot-fixes for BANC metadata issues.
#'
#' Ad-hoc fixer script for the franken / flywire / manc metas, run by hand
#' when a class of metadata bugs surfaces. Not part of any automated run.
#'
#' @section Reads:
#'   - `<banc.meta.save.path>/flywire_meta.csv`
#'   - `<banc.meta.save.path>/manc_meta.csv`
#'   - `<banc.path>/frankenbrain_v.1.4_meta.csv`
#'   - SeaTable `banc_meta`
#'
#' @section Writes:
#'   - back to the meta CSVs above (per-fix-block; consult the script)
#'
#' @section Notes:
#'   - Manual; do not invoke from production.

######################################
### BANC SPOT FIX META DATA ISSUES ###
######################################
source("banc/banc-startup.R")
redo <- FALSE
version <- banc.nblast.version

#################
### Read data ###
#################

# get meta
hms.path <- "/Volumes/neurobio/wilsonlab/banc/meta"
banc.path <- "/Users/papers/BANC-project/data/meta"
fw.meta <- suppressWarnings(readr::read_csv(file = file.path(hms.path,"flywire_meta.csv"), 
                                            col_types = banc.col.types))
mc.meta <- suppressWarnings(readr::read_csv(file = file.path(hms.path,"manc_meta.csv"), 
                                            col_types = banc.col.types))
franken.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.path,"frankenbrain_v.1.4_meta.csv"), 
                                            col_types = banc.col.types))

# get updates from BANCTABLE
bc <- banctable_query()
bc.neck <- bc %>%
  dplyr::filter(!grepl("DELETE|GLIA|NOT_A_NEURON", status),
                grepl("ascending|descending", super_class))

########################
### Edits to Franken ###
########################

# Remove erroneous MANC/FAFB neuron entries
fafb.remove <- na.omit(franken.meta$fafb_id[franken.meta$dataset=="cross-matched"])
manc.remove <- na.omit(franken.meta$manc_id[franken.meta$dataset=="cross-matched"])
franken.meta$banc_id <- bancr::banc_updateids(franken.meta$banc_id)
franken.meta.updated <- franken.meta %>%
  dplyr::filter(!id %in% fafb.remove, !id %in% manc.remove)

# Update missing FAFB meta data
fw.meta <- fw.meta[!duplicated(fw.meta$root_783),]
in.fafb <- match(fw.meta$root_783, franken.meta.updated$fafb_id)
some.nas <- which(!is.na(in.fafb))
in.franken <- in.fafb[some.nas]
fafb.chosen.cols <- intersect(colnames(fw.meta),colnames(franken.meta.updated))
fafb.chosen.cols <- setdiff(fafb.chosen.cols,c("seed","cell_function"))
update <- fw.meta[some.nas, fafb.chosen.cols]
franken.meta.updated[in.franken,fafb.chosen.cols] <- update
franken.meta.updated$nerve[franken.meta.updated$nerve=="CvC"] <- "CV"  

# Further processing
franken.meta.updated <- franken.meta.updated %>%
  dplyr::mutate(dataset = ifelse(!is.na(banc_id),"cross-matched",dataset)) %>%
  dplyr::mutate(cell_function = dplyr::case_when(
    super_class %in% c("sensory_ascending") ~ "sensory_ascending",
    super_class %in% c("ascending") ~ "ascending",
    super_class %in% c("sensory_descending") ~ "sensory_descending",
    super_class %in% c("descending") ~ "descending",
    TRUE ~ gsub("_neuron","",cell_function)
  )) %>%
  # region="neck_connective" no longer set here — ascending/descending neurons
  # are identified via super_class; cervical_connective is a tract entry.
  # dplyr::mutate(region = dplyr::case_when(
  #   super_class %in% c("sensory_ascending","ascending","sensory_descending","descending") ~ "neck_connective",
  #   TRUE ~ region
  # )) %>%

  dplyr::rowwise() %>%
  dplyr::mutate(seed = dplyr::case_when(
    super_class %in% c("sensory_ascending","ascending","sensory_descending","descending") ~ paste0(super_class,"_",side),
    (flow%in%"vnc_afferent") & !is.na(cell_function) & !is.na(nerve) ~ paste(cell_function,snakecase::to_snake_case(nerve),collapse="_",sep="_"),
    (flow%in%"vnc_afferent") & !is.na(cell_function) & !is.na(target) ~ paste(cell_function,snakecase::to_snake_case(target),side,collapse="_",sep="_"),
    (flow%in%"vnc_afferent") & !is.na(cell_function) & is.na(target) ~ paste(cell_function,side,collapse="_",sep="_"),
    (super_class%in%c("ascending","descending","descending_sensory","ascending_sensory")) & !is.na(side) & !is.na(super_class) ~ paste(super_class,side,collapse="_",sep="_"),
    TRUE ~ gsub("__","_",seed)
  ))  %>%
  dplyr::mutate(seed = gsub("_R$|_r$","_right",seed),
                seed = gsub("_L$|_l$","_left",seed),
                seed = gsub("_R_|_r_","_right_",seed),
                seed = gsub("_L_|_l_","_left_",seed)) %>%
  dplyr::ungroup()

# Consolidate duplicate IDs with special handling for match columns
franken.meta.consolidated <- franken.meta.updated %>%
  dplyr::group_by(id) %>%
  dplyr::summarise(
    # Special handling for id columns - concatenate unique non-NA values
    dplyr::across(c(fafb_id, manc_id), 
           ~paste(unique(na.omit(.)), collapse = ",")),
    # Standard handling for all other columns - first non-NA value
    dplyr::across(!c(fafb_id, manc_id), 
           ~if(all(is.na(.))) {
             NA
           } else {
             first(na.omit(.))
           })) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(bc.neck[,c("root_id","an_type","dn_type","side")], by = c("banc_id"="root_id")) %>%
  dplyr::mutate(side = dplyr::case_when(
    is.na(side.y)  ~ side.x,
    is.na(side.x)  ~ side.y,
    TRUE ~ side.x
  )) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(cell_sub_class = dplyr::case_when(
    super_class == "ascending" ~ paste0("AN_",an_type),
    super_class == "descending" ~ paste0("DN_",dn_type),
    TRUE ~ cell_sub_class
  )) %>%
  dplyr::mutate(seed = dplyr::case_when(
    super_class == "ascending" ~ paste0("AN_",an_type,"_",side),
    super_class == "descending" ~ paste0("DN_",dn_type,"_",side),
    TRUE ~ gsub("__","_",seed)
  )) %>%
  dplyr::select(-side.x, -side.y, -an_type, -dn_type) %>%
  dplyr::mutate(dataset = dplyr::case_when(
    !is.na(banc_id)&banc_id!=""&!is.na(fafb_id)&fafb_id!=""&!is.na(manc_id)&manc_id!="" ~ "BANC-FAFB-MANC",
    !is.na(banc_id)&banc_id!=""&!is.na(fafb_id)&fafb_id!="" ~ "BANC-FAFB",
    !is.na(banc_id)&banc_id!=""&!is.na(manc_id)&manc_id!="" ~ "BANC-MANC",
    !is.na(fafb_id)&fafb_id!=""&!is.na(manc_id)&manc_id!="" ~ "FAFB-MANC",
    !is.na(fafb_id)&fafb_id!="" ~ "FAFB",
    !is.na(manc_id)&manc_id!="" ~ "MANC",
    !is.na(banc_id)&banc_id!="" ~ "BANC",
    TRUE ~ dataset
  ))

########################
### Write to sources ###
########################
  
readr::write_csv(x=franken.meta.consolidated,
                 file=file.path(banc.path,"frankenbrain_v.1.4.1_meta.csv"))
arrow::write_feather(franken.meta.consolidated,
                     "/Volumes/neurobio/wilsonlab/banc/connectivity/frankenbrain_v.1.4_meta.feather")
arrow::write_feather(franken.meta.consolidated,
                     "/Users/abates/HMS Dropbox/Alexander Bates/neuroanat/connectomes/frankenbrain_v.1.4_meta.feather")
  
# nb <- dplyr::tbl(con, "neck_bridge") %>%
#   dplyr::collect()
  
