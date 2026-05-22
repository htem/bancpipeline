#' banc-tracing-png-tags — Add tracing-issue tags driven by PNG review outputs.
#'
#' Walks the per-neuron PNG review folders, extracts root_ids / supervoxel_ids
#' from filenames, and applies tracing-issue tags + emits Neuroglancer links
#' for problem neurons.
#'
#' @section Reads:
#'   - PNG match folders under `<banc.nblast.{dataset}.save.path>/`
#'   - SeaTable `banc_meta`
#'
#' @section Writes:
#'   - SeaTable `banc_meta`: col `status` (TRACING_ISSUE tags)
#'   - per-batch NGL link CSVs
#'
#' @section Notes:
#'   - Split out 2026-05-21 from `banc/utilities/banc-tracing.R`.

###############################################################################
### BANC tracing: PNG-derived issue tags + Neuroglancer links
###
### Adds tracing-issue tags driven by PNG review outputs, and produces Neuroglancer links for problem neurons.
###
### Originally part of banc/utilities/banc-tracing.R; split out 2026-05-21
### into per-task scripts under banc/annotations/ because the original was
### never meant to be run end-to-end — each section is a one-off curation
### task. Sections may share state; if so, run them in order.
###############################################################################
source("banc/banc-startup.R")

################################################
### Add tracing issue tag based on PNG files ###
################################################

# Extract updated ID # CB334
extract_id <- function(path, update = FALSE) { 
  # Try to extract supervoxel_id first
  sv_match <- regmatches(path, regexpr("supervoxel_id_[0-9]+", path))
  if(!update){
    root_match <- regmatches(path, regexpr("root_id_[0-9]+", path))
    result <- as.character(sub("root_id_", "", root_match))
  }else{
    if (length(sv_match) > 0 && !is.na(sv_match)) {
      result <-bancr::banc_rootid(as.character(sub("supervoxel_id_", "", sv_match)))
    } else {
      # Fallback to first root_id
      root_match <- regmatches(path, regexpr("root_id_[0-9]+", path))
      if (length(root_match) > 0 && !is.na(root_match)) {
        result <-bancr::banc_latestid(as.character(sub("root_id_", "", root_match)))
      } else {
        result <- NA  # Return NA if neither found
      }
    }    
  }
  result
}

# SLP235/7 | AN09B018 | AN12B004-CB0148 | AN01A055-CB0264
# unzip -v /Volumes/neurobio/wilsonlab/banc/matching/fafb/correct/good.zip -d 3_good/
# rsync -av --ignore-existing --exclude='.DS_Store' /Users/abates/Desktop/matching/fafb/correct/good/ /Volumes/neurobio/wilsonlab/banc/matching/fafb/correct/3_good/

# Remove matched neurons
good.paths <- c("/Users/abates/Desktop/matching/fafb/correct/orphans/",
               "/Users/abates/Desktop/matching/manc/correct/orphans",
               "/Users/abates/Desktop/matching/hemibrain/correct/orphans",
               "/Users/abates/Desktop/matching/fafb/correct/bad",
               "/Users/abates/Desktop/matching/fafb/correct/glia",
               "/Users/abates/Desktop/matching/manc/correct/glia",
               "/Users/abates/Desktop/matching/hemibrain/correct/glia",
               "/Users/abates/Desktop/matching/fafb/correct/")
ids.good <- c()
for(path in good.paths){
  pngs.good <- unique(basename(list.files(path, recursive = TRUE)))
  ids.good <- c(ids.good,unlist(unname(pbapply::pbsapply(pngs.good,extract_id,update=FALSE))) ) 
}
path <- "/Users/abates/Desktop/matching/fafb/todo/"
pngs.current <- list.files(path, recursive = TRUE, full.names = TRUE)
ids.current <- pbapply::pbsapply(pngs.current,extract_id,update=FALSE)
pngs.del <- pngs.current[ids.current%in%ids.good]
sum(file.remove(pngs.del))
bancr:::express_lane("/Users/abates/Desktop/matching/fafb/todo/","^1")
pngs.current <- list.files(path, recursive = TRUE, full.names = TRUE)
ids.current <- pbapply::pbsapply(pngs.current,extract_id,update=FALSE)
pngs.del <- pngs.current[ids.current%in%ids.good]
sum(file.remove(pngs.del))

# Remove matching folders where all hits are optic lobe intrinsic or photoreceptors
optic_types <- unique(c(
  fafb.meta$cell_type[fafb.meta$super_class == "optic_lobe_intrinsic"],
  fafb.meta$cell_type[grepl("photoreceptor", fafb.meta$super_class, ignore.case = TRUE)]
))
optic_types <- optic_types[!is.na(optic_types)]
todo_path <- "/Volumes/neurobio/wilsonlab/banc/matching/fafb/images/todo/"
todo_dirs <- list.dirs(todo_path, recursive = TRUE, full.names = TRUE)
# Keep only leaf directories (matching folders, not super_class parents)
todo_dirs <- todo_dirs[!todo_dirs %in% dirname(todo_dirs[todo_dirs != todo_path])]
todo_dirs <- setdiff(todo_dirs, todo_path)
for (d in todo_dirs) {
  pngs <- list.files(d, pattern = "\\.png$", full.names = FALSE)
  if (length(pngs) == 0) next
  hit_types <- gsub(".*hit_cell_type_(.+)\\.png$", "\\1", pngs)
  if (all(hit_types %in% optic_types)) {
    message("Removing optic-only folder: ", d)
    unlink(d, recursive = TRUE)
  }
}

# Update bad IDs
path <- "/Users/abates/Desktop/matching/fafb/correct/bad"
pngs <- unique(basename(list.files(path, recursive = TRUE)))
pngs.new <- pbapply::pbsapply(pngs,extract_id)
pngs.new <- setdiff(unique(pngs.new),"0")
pngs.new <- na.omit(pngs.new)

# Update easy.fix IDs
path <- "/Users/abates/Desktop/matching/fafb/correct/broken"
pngs <- unique(basename(list.files(path, recursive = TRUE)))
easy.fix <- pbapply::pbsapply(pngs,extract_id)
easy.fix <- setdiff(unique(easy.fix),"0")
easy.fix <- na.omit(easy.fix)

# Update tadpoles
path <- "/Users/abates/Desktop/matching/fafb/correct/orphans/"
pngs <- unique(basename(list.files(path, recursive = TRUE)))
orphan.fix <- pbapply::pbsapply(pngs,extract_id)
orphan.fix <- setdiff(unique(orphan.fix),"0")
orphan.fix <- na.omit(orphan.fix)

# Read BANC meta seatable 
#bancr:::banctable_updateids()
bc.orig <- banctable_query(sql = "select _id, root_id, status, super_class, l2_nodes from banc_meta") 
bc.status <- bc.orig %>%
  dplyr::filter(root_id %in% c(easy.fix,orphan.fix)) %>%
  dplyr::rowwise() %>%
  #dplyr::mutate(status = subtract_status(status,"TRACING_ISSUE_RESOLVED")) %>%
  dplyr::mutate(status = dplyr::case_when(
    #grepl("sensory",super_class) ~ append_status(status,"TRACING_ISSUE_SENSORY"),
    root_id %in% orphan.fix ~ append_status(status,"TRACING_ISSUE_TADPOLE"),
    #root_id %in% easy.fix ~ append_status(status,"TRACING_ISSUE_MISSING_AXON"),
    TRUE ~ append_status(status,"TRACING_ISSUE_3")
  )) %>%
  dplyr::ungroup() %>%
  # dplyr::rowwise() %>%
  # dplyr::mutate(update = dplyr::case_when(
  #   is.na(ngl_link) ~ ngl_link,
  #   root_id %in% tryCatch(c(ngl_segments(ngl_link)), error = function(e)) ~ FALSE,
  #   TRUE ~ TRUE
  # )) %>%
  # dplyr::mutate(ngl_link = dplyr::case_when(
  #   !update ~ ngl_link,
  #   !is.na(root_id) ~ bancsee(banc_ids = c(root_id),
  #                                fafb_ids = c(fafb_nblast_match),
  #                                manc_ids = c(manc_nblast_match),
  #                                nuclei_ids = c(nucleus_id)),
  #   TRUE ~ ngl_link
  # )) %>%
  as.data.frame()

# Update
if(nrow(bc.status)){
  bc.update <- bc.status %>%
    dplyr::select(`_id`, root_id, status)
  banctable_update_rows(base='banc_meta', 
                        table = "banc_meta", 
                        df = as.data.frame(bc.update), 
                        append_allowed = FALSE, 
                        chunksize = 1000)    
}

# Remove tadpole from large neurons
bc.status <- bc.orig %>%
  dplyr::filter(grepl("TRACING_ISSUE_TADPOLE",status),
                l2_nodes >= 100) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(status =subtract_status(status,"TRACING_ISSUE_TADPOLE")) %>%
  dplyr::mutate(status =subtract_status(status,"TOO_SMALL")) %>%
  dplyr::ungroup() %>%
  as.data.frame()

# Update
if(nrow(bc.status)){
  bc.update <- bc.status %>%
    dplyr::select(`_id`, root_id, status)
  banctable_update_rows(base='banc_meta', 
                        table = "banc_meta", 
                        df = bc.update, 
                        append_allowed = FALSE, 
                        chunksize = 1000)    
}

#####################################################
### Provide neuroglancer link for problem neurons ###
#####################################################

# Read BANC meta seatable
#bancr:::banctable_updateids()
bc <- banctable_query(sql = "select _id, root_id, region, l2_nodes, nucleus_id, banc_match, banc_nblast_match, fafb_match, fafb_nblast_match, manc_match, manc_nblast_match, hemibrain_match, hemibrain_nblast_match, status, ngl_link from banc_meta")

# Add ngl link for neurons in need of tracing
bc[bc=="no_match"] <- NA
bc.new <- bc %>%
  dplyr::filter(grepl("TRACING_ISSUE|TOO_SMALL",status), 
                !grepl("TRACING_ISSUE_RESOLVED|NOT_A_NEURON|GLIA|TRACHEA|MERGE_MONSTER|SENT_TO_PRINCETON|SENT_TO_AELYSIA",status),
                !is.na(root_id)) %>%
  dplyr::mutate(banc_match = 
                  dplyr::case_when(
                    is.na(banc_match) ~ banc_nblast_match,
                    TRUE ~ banc_match
                  ),
                fafb_match = 
                  dplyr::case_when(
                    is.na(fafb_match) ~ fafb_nblast_match,
                    TRUE ~ fafb_match
                  ),
                manc_match = 
                  dplyr::case_when(
                    is.na(manc_match) ~ manc_nblast_match,
                    TRUE ~ manc_match
                  ),
                hemibrain_match = 
                  dplyr::case_when(
                    is.na(hemibrain_match) ~ hemibrain_nblast_match,
                    TRUE ~ hemibrain_match
                  )) %>%
  dplyr::mutate(hemibrain_match = gsub("m"," ",hemibrain_match)) %>%
  dplyr::mutate(update = TRUE) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(update = dplyr::case_when(
    is.na(ngl_link) ~ ngl_link,
    root_id %in% tryCatch(c(ngl_segments(ngl_link)), error = function(e) character(0)) ~ FALSE,
    TRUE ~ TRUE
  )) %>%
  dplyr::mutate(ngl_link = dplyr::case_when(
    !update ~ ngl_link,
    !is.na(banc_match) ~ bancsee(banc_ids = c(root_id,banc_match),
                                                      fafb_ids = c(fafb_match),
                                                      manc_ids = c(manc_match),
                                                      nuclei_ids = c(nucleus_id)),
    is.na(banc_match) ~ bancsee(banc_ids = c(root_id),
                                 fafb_ids = c(fafb_match),
                                 manc_ids = c(manc_match),
                                 nuclei_ids = c(nucleus_id)),
  TRUE ~ ngl_link
  )) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(status = append_status(status,"SENT_TO_PRINCETON_2"))

# Update
bc.update <- bc.new %>%
  dplyr::select(`_id`, root_id, ngl_link)
banctable_update_rows(base='banc_meta', 
                      table = "banc_meta", 
                      df = as.data.frame(bc.update), 
                      append_allowed = FALSE, 
                      chunksize = 1000)

# Write as .csv file
bc.tracing2 <- bc.new %>%
  dplyr::filter(grepl("TRACING_ISSUE_MISSING_AXON",status)) %>%
  dplyr::mutate(l2_nodes = as.numeric(l2_nodes)) %>%
  dplyr::arrange(region, l2_nodes) %>%
  dplyr::select(root_id, 
                nucleus_id,
                banc_match, 
                fafb_match, 
                hemibrain_match, 
                manc_match,
                region,
                l2_nodes,
                status,
                ngl_link) %>%
  dplyr::mutate(notes="",
                annotator="")
banc.tracing.save.path <- "tracing/tracing_issues/"
current_datetime <- Sys.time()
datetime_string <- format(current_datetime, "%Y-%m-%d")
readr::write_csv(bc.tracing, file.path(banc.tracing2.save.path,paste0(datetime_string,"_banc_tracing_issues_missing_axon.csv")))

# Update
bc.update2 <- bc.tracing2 %>%
  dplyr::select(`_id`, root_id, status)
banctable_update_rows(base='banc_meta', 
                      table = "banc_meta", 
                      df = as.data.frame(bc.update2), 
                      append_allowed = FALSE, 
                      chunksize = 1000)

# Write as .csv file
bc.tracing3 <- bc.new %>%
  dplyr::filter(!is.na(nucleus_id),
                !proofread,
                grepl("midbrain|central_brain",region),
                !grepl("TRACING_ISSUE_MISSING_AXON",status)) %>%
  dplyr::mutate(l2_nodes = as.numeric(l2_nodes)) %>%
  dplyr::arrange(region, l2_nodes) %>%
  dplyr::select(root_id, 
                nucleus_id,
                banc_match, 
                fafb_match, 
                hemibrain_match, 
                manc_match,
                region,
                l2_nodes,
                status,
                ngl_link) %>%
  dplyr::mutate(notes="",
                annotator="")
banc.tracing.save.path <- "tracing/tracing_issues/"
current_datetime <- Sys.time()
datetime_string <- format(current_datetime, "%Y-%m-%d")
readr::write_csv(bc.tracing3, file.path(banc.tracing.save.path,paste0(datetime_string,"_banc_tracing_issues_midbrain_too_small.csv")))

# Update
bc.update3 <- bc.tracing %>%
  dplyr::select(`_id`, root_id, status)
banctable_update_rows(base='banc_meta', 
                      table = "banc_meta", 
                      df = as.data.frame(bc.update3), 
                      append_allowed = FALSE, 
                      chunksize = 1000)
