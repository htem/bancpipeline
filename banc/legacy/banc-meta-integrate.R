##############################################
### Integrate BANC meta data with seatable ###
##############################################
source("banc/banc-startup.R")
message("##### Building seatable meta data #####")

# Update BANC meta seatable
bancr:::banctable_updateids()

# Save snapshot
bc.orig <- banctable_query()
dir.create(file.path(banc.meta.save.path,"snapshots"), showWarnings = FALSE)
current_datetime <- Sys.time()
datetime_string <- format(current_datetime, "%Y-%m-%d_%H-%M")
readr::write_csv(bc.orig, file.path(file.path(banc.meta.save.path,"snapshots"),paste0(datetime_string,"_banc_seatable.csv")))

#####################
### get meta data ###
#####################

# Get save folder
tracing.folder <- "tracing"
dir.create(tracing.folder, showWarnings = FALSE)

# Read meta
banc.ids <- readr::read_csv(file=file.path(banc.meta.save.path,"banc_ids.csv"),
                            col_types = banc.col.types, 
                            show_col_types = FALSE) %>%
  dplyr::filter(!is.na(root_id)) %>%
  dplyr::distinct(root_id, .keep_all = TRUE)
# fw.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"flywire_meta.csv"), 
#                                             col_types = hemibrainr:::sql_col_types))
# mc.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"manc_meta.csv"), 
#                                             col_types = hemibrainr:::sql_col_types))
franken.meta <- franken_meta()
fw.meta <- franken.meta %>%
  dplyr::filter(!is.na(fafb_id)) %>%
  dplyr::rename(root_783=fafb_id)
mc.meta <- franken.meta %>%
  dplyr::filter(!is.na(manc_id)) %>%
  dplyr::rename(bodyid=manc_id)
hb.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"hemibrain_meta.csv"), 
                                            col_types = hemibrainr:::sql_col_types))
fc.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"fanc_meta.csv"),
                                            col_types = hemibrainr:::sql_col_types))

malecns.meta <- readr::read_csv(file.path(banc.meta.save.path,"malecns_09_meta.csv"), col_types = banc.col.types)
mcns.meta <- malecns.meta %>%
  dplyr::filter(!is.na(malecns_09_id)) %>%
  dplyr::rename(bodyid=malecns_09_id)

# Metric columns
metric.cols <- c(#"supervoxel_id", "position", 
  #"super_class", "cell_class", #"cell_type", # "nerve",
  "proofread", 
  "nucleus_id", 
  "nucleus_supervoxel_id", 
  "nucleus_position_nm", 
  "l2_nodes", 
  "l2_cable_length_um",
  "input_connections", "output_connections",
  "input_side_index", "output_side_index",
  "banc_nblast_match", 
  "banc_nblast_match_supervoxel_id", 
  "banc_nblast",
  "fafb_nblast_match", 
  "fafb_nblast",
  "manc_nblast_match", 
  "manc_nblast",
  "fanc_nblast_match", 
  "fanc_nblast",
  "hemibrain_nblast_match",
  "hemibrain_nblast",
  "malecns_nblast_match",
  "malecns_nblast")
drop.cols <- c("_locked", "_locked_by", 
               "_archived", "_creator", 
               "_ctime", "_last_modifier", "_mtime")

######################################
### combine with latest CAVE pulls ###
######################################

# Perform the join and update
bc1 <- banctable_query()
bc1[,drop.cols] <- NULL
bc <- bc1 %>%
  dplyr:::rename(root_id.x=root_id, position.x=position) %>%
  dplyr::full_join(banc.ids %>% 
                     dplyr:::select(root_id, supervoxel_id, position) %>% 
                     dplyr:::rename(root_id.y=root_id, position.y=position), 
                   by = "supervoxel_id") %>%
  dplyr::mutate(
    root_id = dplyr::coalesce(root_id.y, root_id.x, root_id.y),
    position = dplyr::coalesce(position.y, position.x, position.y)
  ) %>%
  dplyr::select(-any_of(c("root_id.x", "root_id.y", "position.x", "position.y"))) %>%
  dplyr::distinct(root_id, .keep_all = TRUE) %>%
  dplyr::distinct(supervoxel_id, .keep_all = TRUE) %>%
  dplyr::distinct(position, .keep_all = TRUE)  %>%
  dplyr::filter(!is.na(root_id),!is.na(supervoxel_id))

# Read BANC meta local
banc.meta <- readr::read_csv(file=file.path(banc.meta.save.path,"banc_meta.csv"),
                             col_types = banc.col.types, 
                             show_col_types = FALSE) %>%
  dplyr::arrange(cell_type,desc(fafb_nblast),desc(manc_nblast),desc(malecns_nblast),desc(banc_nblast))
banc.meta <- banc.meta %>%
  dplyr::filter(!is.na(root_id), !is.na(supervoxel_id), root_id!="0")
banc.meta <- banc_updateids(banc.meta) 
readr::write_csv(banc.meta,file=file.path(banc.meta.save.path,"banc_meta.csv"))
banc.meta <- banc.meta[!duplicated(banc.meta$root_id),]
metric.cols <- intersect(metric.cols,colnames(banc.meta))

# Wipe columns that we want to totally auto-update
bc[,colnames(bc)%in%metric.cols] <- NULL

# Join
bc$fafb_png_match <- NULL
bc$manc_png_match <- NULL
bc$fanc_png_match <- NULL
bc$hemibrain_png_match <- NULL
bc$malecns_png_match <- NULL
bc$fafb_nblast_match <- NULL
bc$manc_nblast_match <- NULL
bc$fanc_nblast_match <- NULL
bc$hemibrain_png_match <- NULL
bc$malecns_nblast_match <- NULL
bc.new <- dplyr::full_join(bc, 
                           banc.meta[,c("root_id",metric.cols)], 
                           by="root_id")
bc.new[bc.new=="NA"] <- NA
bc.new[bc.new==""] <- NA
bc.new[bc.new==" "] <- NA
bc.new[] <- lapply(bc.new, function(x) {
  if(is.character(x)) {
    x[grepl("^auto", x, ignore.case = TRUE)] <- NA
  }
  return(x)
})
bc.new$status[bc.new$status=="NA"]<-""

#############################################
### Conditionally update matching columns ###
#############################################

# Human reviewed matches
mirror.matches.df <- readr::read_csv(file=file.path(banc.meta.save.path,"banc_mirror_reviewed_matches.csv"),
                                     col_types = banc.col.types, 
                                     show_col_types = FALSE) %>%
  dplyr::filter(!grepl("\\_",match_id)) %>%
  dplyr::rename(root_id=pt_root_id,
                banc_png_match=match_root_id,
                banc_png_match_supervoxel_id=match_supervoxel_id)  %>%
  dplyr::filter(valid=='t')
manc.matches.df <- readr::read_csv(file=file.path(banc.meta.save.path,"banc_manc_reviewed_matches.csv"),
                                   col_types = banc.col.types, 
                                   show_col_types = FALSE) %>%
  dplyr::filter(!grepl("\\_",match_id)) %>%
  dplyr::rename(root_id=pt_root_id,
                manc_png_match=match_id) %>%
  dplyr::filter(valid=='t')
fafb.matches.df <- readr::read_csv(file=file.path(banc.meta.save.path,"banc_fafb_reviewed_matches.csv"),
                                   col_types = banc.col.types, 
                                   show_col_types = FALSE) %>%
  dplyr::filter(!grepl("\\_",match_id)) %>%
  dplyr::rename(root_id=pt_root_id,
                fafb_png_match=match_id) %>%
  dplyr::filter(valid=='t')
hemibrain.matches.df <- readr::read_csv(file=file.path(banc.meta.save.path,"banc_hemibrain_reviewed_matches.csv"),
                                        col_types = banc.col.types, 
                                        show_col_types = FALSE) %>%
  dplyr::filter(!grepl("\\_",match_id)) %>%
  dplyr::rename(root_id=pt_root_id,
                hemibrain_png_match=match_id) %>%
  dplyr::filter(valid=='t')
fanc.matches.df <- readr::read_csv(file=file.path(banc.meta.save.path,"banc_fanc_reviewed_matches.csv"),
                                   col_types = banc.col.types,
                                   show_col_types = FALSE) %>%
  dplyr::filter(!grepl("\\_",match_id)) %>%
  dplyr::rename(root_id=pt_root_id,
                fanc_png_match=match_id) %>%
  dplyr::filter(valid=='t')
malecns.matches.df <- readr::read_csv(file=file.path(banc.meta.save.path,"banc_malecns_reviewed_matches.csv"),
                                      col_types = banc.col.types,
                                      show_col_types = FALSE) %>%
  dplyr::filter(!grepl("\\_",match_id)) %>%
  dplyr::rename(root_id=pt_root_id,
                malecns_png_match=match_id) %>%
  dplyr::filter(valid=='t')

# Matches adjusted for side
banc.meta.fafb.nb <- arrow::read_feather(file.path(banc.meta.save.path,"banc_fafb_783_nblast.feather")) %>%
  dplyr::filter(valid=='t',!grepl("\\_",match_id)) %>%
  dplyr::rename(root_id=pt_root_id,
                fafb_match=match_id)
banc.meta.manc.nb <- arrow::read_feather(file.path(banc.meta.save.path,"banc_manc_v1.2.1_nblast.feather")) %>%
  dplyr::filter(valid=='t',!grepl("\\_",match_id)) %>%
  dplyr::rename(root_id=pt_root_id,
                manc_match=match_id)
banc.meta.hemibrain.nb <- arrow::read_feather(file.path(banc.meta.save.path,"banc_hemibrain_v1.2.1_nblast.feather")) %>%
  dplyr::filter(valid=='t',!grepl("\\_",match_id))%>%
  dplyr::rename(root_id=pt_root_id,
                hemibrain_match=match_id)
banc.meta.fanc.nb <- arrow::read_feather(file.path(banc.meta.save.path,"banc_fanc_1116_nblast.feather")) %>%
  dplyr::filter(valid=='t') %>%
  dplyr::rename(root_id=pt_root_id,
                fanc_match=match_id)
banc.meta.malecns.nb <- arrow::read_feather(file.path(banc.meta.save.path,"banc_malecns_v0.9_nblast.feather")) %>%
  dplyr::filter(valid=='t',!grepl("\\_",match_id)) %>%
  dplyr::rename(root_id=pt_root_id,
                malecns_match=match_id)
banc.meta.mirror.nb <- arrow::read_feather(file.path(banc.meta.save.path,"banc_mirror_nblast.feather")) %>%
  dplyr::filter(valid=='t') %>%
  dplyr::rename(root_id=pt_root_id,
                banc_match=match_root_id,
                banc_match_supervoxel_id=match_supervoxel_id)

# Combine together
matches.df <- plyr::rbind.fill(mirror.matches.df,
                               fafb.matches.df,
                               manc.matches.df,
                               fanc.matches.df,
                               hemibrain.matches.df,
                               malecns.matches.df,
                               banc.meta.fafb.nb,
                               banc.meta.manc.nb,
                               banc.meta.fanc.nb,
                               banc.meta.hemibrain.nb,
                               banc.meta.malecns.nb)

# Update match columns
update.cols <- c("banc_match","banc_match_supervoxel_id","fafb_match","manc_match", "hemibrain_match", "fanc_match", "malecns_match",
                 "banc_png_match","banc_png_match_supervoxel_id","fafb_png_match","manc_png_match","hemibrain_png_match","fanc_png_match","malecns_png_match")
bc.new[,c("banc_png_match","banc_png_match_supervoxel_id","fafb_png_match","manc_png_match","hemibrain_png_match","fanc_png_match","malecns_png_match")]<-""
conflicts <- data.frame()
for(update.col in update.cols){
  dataset.df <- matches.df[!is.na(matches.df[[update.col]]),]
  for(i in 1:nrow(bc.new)){
    entry <- bc.new[i,update.col]
    new.entry <- dataset.df[match(bc.new$root_id[i], dataset.df$root_id),update.col]
    if(!length(new.entry)){
      next
    }
    if(is.na(new.entry)|new.entry==''){
      next
    }
    if(!is.na(entry)&!is.na(new.entry)&entry==new.entry){
      if(is.na(bc.new[i,"cell_type_source"])){
        bc.new[i,"cell_type_source"] <- paste(sort(unique(c(bc.new[i,"cell_type_source"],"wilson lab"))),collapse=",")
      }
    }
    if(is.na(entry)|entry%in%c("NA",""," ","unknown","issue","ISSUE","ISSUE!")){
      if(!is.na(new.entry)){
        #message("adding ", update.col,": ", new.entry," for ", bc.new[i,"root_id"])
        bc.new[i,update.col] <- new.entry
        bc.new[i,"cell_type_source"] <- paste(sort(unique(c(bc.new[i,"cell_type_source"],"wilson lab"))),collapse=",")
      }
    }else{
      if(entry!=new.entry){
        message("conflict ", update.col,": current ", entry, " versus suggested ", new.entry," for ", bc.new[i,"root_id"])
        cnflts <- data.frame(`_id`= bc.new[i,"_id"], status = bc.new[i,"status"], root_id = bc.new[i,"root_id"], side = bc.new[i,"side"], 
                             type = update.col, current = entry, suggested = new.entry)
        conflicts <- rbind(conflicts,cnflts)
      }
    }
  }
}

# # Assess conflicts
# if(nrow(conflicts)){
#   try({
#     
#     # Join conflicts with fw.meta for 'fafb_match' type
#     fafb_conflicts <- conflicts %>%
#       filter(type == 'fafb_match') %>%
#       dplyr::rename(banc_side = side) %>%
#       left_join(fw.meta[,c("root_783","side","cell_type")], by = c("current"="root_783")) %>%
#       dplyr::rename(current_cell_type = cell_type, current_side = side) %>%
#       left_join(fw.meta[,c("root_783","side","cell_type")], by = c("suggested"="root_783")) %>%
#       dplyr::rename(suggested_cell_type = cell_type, suggested_side = side)
#     
#     # Side conflicts we can resolve
#     fafb_conflicts.side <- fafb_conflicts %>%
#       dplyr::filter(current_cell_type==suggested_cell_type,
#                     banc_side==suggested_side)
#     bc.new[match(fafb_conflicts.side$root_id,bc.new$root_id),"fafb_match"] <- fafb_conflicts.side$suggested
#     
#     # Matching conflicts
#     fafb_conflicts.issues <- fafb_conflicts %>%
#       dplyr::filter(current_cell_type!=suggested_cell_type) %>%
#       dplyr::rowwise() %>%
#       dplyr::mutate(ngl_link = bancsee(
#         banc_ids = root_id,
#         fafb_ids = c(current,suggested)))
#     
#     # Join conflicts with hb.meta for 'hemibrain_match' type
#     # hb_conflicts <- conflicts %>%
#     #   filter(type == 'hemibrain_match') %>%
#     #   dplyr::rename(banc_side = side) %>%
#     #   left_join(hb.meta[,c("bodyid", "side", "cell_type")], by = c("current"="bodyid")) %>%
#     #   dplyr::rename(current_cell_type = cell_type, current_side = side) %>%
#     #   left_join(hb.meta[,c("bodyid", "side","cell_type")], by = c("suggested"="bodyid")) %>%
#     #   dplyr::rename(suggested_cell_type = cell_type, suggested_side = side)
#     
#     # Join conflicts with mc.meta for 'manc_match' type
#     manc_conflicts <- conflicts %>%
#       filter(type == 'manc_match') %>%
#       dplyr::rename(banc_side = side) %>%
#       left_join(mc.meta[,c("bodyid", "side","cell_type")], by = c("current"="bodyid")) %>%
#       dplyr::rename(current_cell_type = cell_type, current_side = side) %>%
#       left_join(mc.meta[,c("bodyid", "side","cell_type")], by = c("suggested"="bodyid")) %>%
#       dplyr::rename(suggested_cell_type = cell_type, suggested_side = side)
#     
#     # Side conflicts we can resolve
#     manc_conflicts.side <- manc_conflicts %>%
#       dplyr::filter(current_cell_type==suggested_cell_type,
#                     banc_side==suggested_side)
#     bc.new[match(manc_conflicts.side$root_id,bc.new$root_id),"manc_match"] <- manc_conflicts.side$suggested
#     
#     # Matching conflicts
#     manc_conflicts.issues <- manc_conflicts %>%
#       dplyr::filter(current_cell_type!=suggested_cell_type) %>%
#       dplyr::rowwise() %>%
#       dplyr::mutate(ngl_link = bancsee(
#         banc_ids = root_id,
#         manc_ids = c(current,suggested)))
#     
#     # examine final product
#     final_conflicts <- bind_rows(
#       fafb_conflicts.issues,
#       # hb_conflicts,
#       manc_conflicts.issues) %>%
#       dplyr::mutate(fafb_match = "",
#                     manc_match = "") %>%
#       dplyr::arrange(status, current_cell_type, root_id, suggested_cell_type)
#     
#     # Save
#     current_datetime <- Sys.time()
#     datetime_string <- format(current_datetime, "%Y-%m-%d")
#     write.csv(x=final_conflicts,
#               file=file.path(tracing.folder,paste0(datetime_string,"_matching_conflicts.csv")))
#   })
# }  

# Update match columns
update.cols <- c("region", "cell_class", "nerve", "nucleus_id") # "super_class"
for(update.col in update.cols){
  bc.meta.df <- banc.meta[!is.na(banc.meta[[update.col]]),]
  for(i in 1:nrow(bc.new)){
    entry <- bc.new[i,update.col]
    new.entry <- bc.meta.df[match(bc.new$root_id[i], bc.meta.df$root_id),][[update.col]]
    if(is.na(new.entry)|new.entry==''){
      next
    }
    if(is.na(entry)|entry%in%c("NA",""," ","unknown","issue","ISSUE","ISSUE!")){
      message("adding ", update.col,": ", new.entry," for ", bc.new[i,"root_id"])
      bc.new[i,update.col] <- new.entry
    }
  }
}

##############################################################
### Update matches if cell type labels given but not match ###
##############################################################

# # Make some matches by proxy
# bc.cts <- unique(bc.new$cell_type)
# bc.cts <-bc.cts[!is.na(bc.cts)]
# fw.proxy.matches <- fw.meta %>%
#   dplyr::filter(cell_type %in% bc.cts) %>%
#   dplyr::rename(fafb_match=root_783) %>%
#   dplyr::distinct(fafb_match, .keep_all = TRUE) %>%
#   dplyr::select(fafb_match, side, cell_type) %>%
#   dplyr::left_join(bc.new[,c("root_id","side","cell_type")], 
#                    by = c("side","cell_type"),
#                    relationship = "many-to-many") %>%
#   dplyr::filter(!is.na(root_id)) %>%
#   dplyr::mutate(source = "proxy")
# mc.proxy.matches <- mc.meta %>%
#   dplyr::rename(manc_match=bodyid) %>%
#   dplyr::filter(cell_type %in% bc.cts) %>%
#   dplyr::mutate(side = ifelse(side=="unknown",root_side,side)) %>%
#   dplyr::filter(side %in% c("RHS","LHS","right","left")) %>%
#   dplyr::mutate(side = ifelse(side=="RHS","right",side)) %>%
#   dplyr::mutate(side = ifelse(side=="LHS","left",side)) %>%
#   dplyr::distinct(manc_match, side, cell_type) %>%
#   dplyr::left_join(bc.new[,c("root_id","side","cell_type")], 
#                    by = c("side","cell_type"),
#                    relationship = "many-to-many") %>%
#   dplyr::filter(!is.na(root_id)) %>%
#   dplyr::mutate(source = "proxy")
# hb.proxy.matches <- hb.meta %>%
#   dplyr::rename(hemibrain_match=bodyid) %>%
#   dplyr::arrange(dplyr::desc(post), dplyr::desc(pre)) %>%
#   dplyr::filter(cell_type %in% bc.cts) %>%
#   dplyr::distinct(hemibrain_match,cell_type) %>%
#   dplyr::left_join(bc.new[,c("root_id","cell_type")], 
#                    by = c("cell_type"),
#                    relationship = "many-to-many") %>%
#   dplyr::filter(!is.na(root_id)) %>%
#   dplyr::mutate(source = "proxy")
# 
# # Combine
# ct.community.matches$source <- "community"
# proxy.matches.df <- plyr::rbind.fill(ct.community.matches,
#                                      fw.proxy.matches,
#                                      mc.proxy.matches,
#                                      hb.proxy.matches)
# 
# # Update match columns
# update.cols <- c("cell_type_source", "cell_class", "super_class", "cell_type","banc_match","fafb_match","manc_match", "hemibrain_match")
# for(update.col in update.cols){
#   dataset.df <- proxy.matches.df[!is.na(proxy.matches.df[[update.col]]),]
#   for(i in 1:nrow(bc.new)){
#     entry <- bc.new[i,update.col]
#     new.entry <- dataset.df[match(bc.new$root_id[i], dataset.df$root_id),update.col]
#     if(!length(new.entry)){
#       next
#     }
#     if(is.na(new.entry)|new.entry==''){
#       next
#     }
#     if(!is.na(entry)&(entry==new.entry)){
#       if(update.col=="cell_type"&(is.na(bc.new[i,"cell_type_source"])|bc.new[i,"cell_type_source"]=='')){
#         bc.new[i,"cell_type_source"] <- dataset.df[match(bc.new$root_id[i], dataset.df$root_id),"source"]
#       }
#     }
#     if(is.na(entry)|entry%in%c("NA",""," ","unknown","issue","ISSUE","ISSUE!")){
#       if(!is.na(new.entry)){
#         message("adding proxy ", update.col,": ", new.entry," for ", bc.new[i,"root_id"])
#         bc.new[i,update.col] <- new.entry
#         if(update.col=="cell_type"&(is.na(bc.new[i,"cell_type_source"])|bc.new[i,"cell_type_source"]=='')){
#           bc.new[i,"cell_type_source"] <- dataset.df[match(bc.new$root_id[i], dataset.df$root_id),"source"]
#         }
#       }
#     }else{
#       # if(entry!=new.entry){
#       #   message("proxy conflict ", update.col,": old ", entry, " versus new ", new.entry," for ", bc.new[i,"root_id"])
#       # }
#     }
#   }
# }

###########################################################
### Update cell info if matches given but not cell info ###
###########################################################
bc.new <- bc.new %>%
  dplyr::mutate(fafb_match = dplyr::case_when(
    is.na(fafb_match)&!fafb_png_match%in%c("no_match","","0") ~ fafb_png_match,
    TRUE ~ fafb_match
  )) %>%
  dplyr::mutate(manc_match = dplyr::case_when(
    is.na(manc_match)&!manc_png_match%in%c("no_match","","0") ~ manc_png_match,
    TRUE ~ manc_match
  )) %>%
  dplyr::mutate(hemibrain_match = dplyr::case_when(
    is.na(hemibrain_match)&!hemibrain_png_match%in%c("no_match","","0") ~ hemibrain_png_match,
    TRUE ~ hemibrain_match
  )) %>%
  dplyr::mutate(fanc_match = dplyr::case_when(
    is.na(fanc_match)&fanc_png_match%in%c("no_match","","0") ~ fanc_png_match,
    TRUE ~ fanc_match
  )) %>%
  dplyr::mutate(banc_match = dplyr::case_when(
    is.na(banc_match)&!banc_png_match%in%c("no_match","","0") ~ banc_png_match,
    TRUE ~ banc_match
  )) %>%
  dplyr::mutate(banc_match_supervoxel_id = dplyr::case_when(
    is.na(banc_match_supervoxel_id)&!banc_png_match_supervoxel_id%in%c("no_match","","0") ~ banc_png_match_supervoxel_id,
    TRUE ~ banc_match_supervoxel_id
  ))

# Integrate FAFB NBLAST results
banc.meta.fafb <- banc.meta
banc.meta.fafb$fafb_nblast_match <- NULL 
banc.meta.fafb$fafb_nblast <- NULL 
fafb.cols <- colnames(banc.meta.fafb)[grepl("^fafb",colnames(banc.meta.fafb))]
banc.meta.fafb <- banc.meta[,c("root_id",fafb.cols)]
banc.meta.fafb <- banc.meta.fafb[!duplicated(banc.meta.fafb$root_id),]
bc.new.fafb <- bc.new %>%
  dplyr::left_join(banc.meta.fafb, by = "root_id") %>%
  dplyr::mutate(cell_type = dplyr::case_when(
    !is.na(cell_type) ~ cell_type,
    !is.na(fafb_match)&(is.na(cell_type)|grepl("auto",cell_type)|cell_type=="") ~ fw.meta$cell_type[match(fafb_match,fw.meta$root_783)],      
    TRUE ~ cell_type
  ),
  fafb_cell_type = dplyr::case_when(
    !is.na(fafb_match)&fafb_match!="" ~ fw.meta$cell_type[match(fafb_match,fw.meta$root_783)],      
    is.na(fafb_cell_type) ~ fafb_auto_cell_type,
    TRUE ~ fafb_cell_type
  ),
  auto_type = is.na(cell_type)|grepl("auto",cell_type)
  ) %>%
  dplyr::mutate(
    # hemilineage = dplyr::case_when(
    #   !grepl("auto",hemilineage)&!is.na(hemilineage) ~ hemilineage,
    #   !auto_type ~ fw.meta$hemilineage[match(cell_type,fw.meta$cell_type)],
    #   is.na(hemilineage) ~ fafb_ito_lee_hemilineage,
    #   TRUE ~ hemilineage
    # ),
    # nerve = dplyr::case_when(
    #   !grepl("auto",nerve)&!is.na(nerve) ~ nerve,
    #   !auto_type ~ fw.meta$nerve[match(cell_type,fw.meta$cell_type)],
    #   is.na(nerve) ~ fafb_nerve,
    #   TRUE ~ nerve
    # ),
    # fafb_super_class = dplyr::case_when(
    #   !auto_type ~ fw.meta$super_class[match(fafb_cell_type,fw.meta$cell_type)],
    #   TRUE ~ fafb_super_class.y
    # ),
    # fafb_cell_class = dplyr::case_when(
    #   !auto_type ~ fw.meta$cell_class[match(fafb_cell_type,fw.meta$cell_type)],
    #   TRUE ~ fafb_cell_class.y
    # ),
    top_nt = dplyr::case_when(
      !auto_type ~ fw.meta$top_nt[match(cell_type,fw.meta$cell_type)],
      is.na(fafb_top_nt) ~ fafb_top_nt,
      TRUE ~ top_nt
    )
    # known_nt = dplyr::case_when(
    #   !auto_type ~ fw.meta$known_nt[match(cell_type,fw.meta$cell_type)],
    #   TRUE ~ known_nt
    # ),
    # known_nt_source = dplyr::case_when(
    #   !auto_type ~ fw.meta$known_nt_source[match(cell_type,fw.meta$cell_type)],
    #   TRUE ~ known_nt_source
    # ),
    # cell_function = dplyr::case_when(
    #   !grepl("auto",cell_function)&!is.na(cell_function) ~ cell_function,
    #   !auto_type ~ fw.meta$cell_function[match(cell_type,fw.meta$cell_type)],
    #   is.na(cell_function) ~ fafb_cell_function,
    #   TRUE ~ cell_function
    # )
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(status = dplyr::case_when(
    is.na(cell_type)|is.na(fafb_match)|grepl("neck",region) ~ status,
    cell_type != fw.meta$cell_type[match(fafb_match,fw.meta$root_783)] ~ append_status(status,"FAFB_ALT_MATCH"),
    TRUE ~ status
  )) %>%
  dplyr::ungroup() %>%
  dplyr::select(-any_of(fafb.cols), -auto_type)
bc.new <- bc.new.fafb

# Integrate MANC NBLAST results
banc.meta.manc <- banc.meta %>%
  dplyr::distinct(root_id, .keep_all = TRUE) %>%
  dplyr::select(root_id, dplyr::starts_with("manc")) %>%
  dplyr::select(-manc_nblast_match, -manc_nblast) %>%
  as.data.frame()
manc.cols <- colnames(banc.meta.manc)[grepl("^manc",colnames(banc.meta.manc))]
bc.new.manc <- bc.new %>%
  dplyr::left_join(banc.meta.manc, by = "root_id") %>%
  dplyr::mutate(cell_type = dplyr::case_when(
    !is.na(cell_type) ~ cell_type,
    !is.na(manc_match)&(is.na(cell_type)|grepl("auto",cell_type)|cell_type=="") ~ mc.meta$cell_type[match(manc_match,mc.meta$bodyid)],      
    TRUE ~ cell_type
  ),
  manc_cell_type = dplyr::case_when(
    !is.na(manc_match)&manc_match!="" ~ mc.meta$cell_type[match(manc_match,mc.meta$bodyid)],      
    is.na(manc_cell_type) ~ manc_auto_cell_type,
    TRUE ~ manc_cell_type
  ),
  auto_type = is.na(cell_type)|grepl("auto",cell_type)
  ) %>%
  dplyr::mutate(
    #hemilineage = dplyr::case_when(
    #   !grepl("auto",hemilineage)&!is.na(hemilineage) ~ hemilineage,
    #   !auto_type ~ mc.meta$hemilineage[match(cell_type,mc.meta$cell_type)],
    #   is.na(hemilineage) ~ manc_hemilineage,
    #   TRUE ~ hemilineage
    # ),
    # nerve = dplyr::case_when(
    #   !grepl("auto",nerve)&!is.na(nerve) ~ nerve,
    #   !auto_type ~ mc.meta$nerve[match(cell_type,mc.meta$type)],
    #   is.na(nerve) ~ manc_nerve,
    #   TRUE ~ nerve
    # ),
    # manc_super_class = dplyr::case_when(
    #   !auto_type ~ mc.meta$super_class[match(manc_cell_type,mc.meta$cell_type)],
    #   TRUE ~ manc_super_class
    # ),
    # manc_cell_class = dplyr::case_when(
    #   !auto_type ~ mc.meta$cell_class[match(manc_cell_type,mc.meta$cell_type)],
    #   TRUE ~ manc_cell_class.y
    # ),
    top_nt = dplyr::case_when(
      !grepl("auto",top_nt)&!is.na(top_nt) ~ top_nt,
      !auto_type ~ mc.meta$top_nt[match(cell_type,mc.meta$cell_type)],
      is.na(manc_top_nt) ~ manc_top_nt,
      TRUE ~ top_nt
    )
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(status = dplyr::case_when(
    is.na(cell_type)|is.na(manc_match)|grepl("neck",region) ~ status,
    cell_type != mc.meta$cell_type[match(manc_match,mc.meta$bodyid)] ~ append_status(status,"MANC_ALT_MATCH"),
    TRUE ~ status
  )) %>%
  dplyr::ungroup() %>%
  dplyr::select(-any_of(manc.cols), -auto_type)
bc.new <- bc.new.manc

# Integrate maleCNS NBLAST results
banc.meta.malecns <- banc.meta %>%
  dplyr::distinct(root_id, .keep_all = TRUE) %>%
  dplyr::select(root_id, dplyr::starts_with("malecns")) %>%
  dplyr::select(-malecns_nblast_match, -malecns_nblast) %>%
  as.data.frame()
malecns.cols <- colnames(banc.meta.malecns)[grepl("^malecns",colnames(banc.meta.malecns))]
bc.new.malecns <- bc.new %>%
  dplyr::left_join(banc.meta.malecns, by = "root_id") %>%
  dplyr::mutate(cell_type = dplyr::case_when(
    !is.na(cell_type) ~ cell_type,
    !is.na(malecns_match)&(is.na(cell_type)|grepl("auto",cell_type)|cell_type=="") ~ mcns.meta$cell_type[match(malecns_match,mcns.meta$bodyid)],
    TRUE ~ cell_type
  ),
  malecns_cell_type = dplyr::case_when(
    !is.na(malecns_match)&malecns_match!="" ~ mcns.meta$cell_type[match(malecns_match,mcns.meta$bodyid)],
    is.na(malecns_cell_type) ~ malecns_auto_cell_type,
    TRUE ~ malecns_cell_type
  ),
  auto_type = is.na(cell_type)|grepl("auto",cell_type)
  ) %>%
  dplyr::mutate(
    top_nt = dplyr::case_when(
      !grepl("auto",top_nt)&!is.na(top_nt) ~ top_nt,
      !auto_type ~ mcns.meta$top_nt[match(cell_type,mcns.meta$cell_type)],
      is.na(malecns_top_nt) ~ malecns_top_nt,
      TRUE ~ top_nt
    )
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(status = dplyr::case_when(
    is.na(cell_type)|is.na(malecns_match)|grepl("neck",region) ~ status,
    cell_type != mcns.meta$cell_type[match(malecns_match,mcns.meta$bodyid)] ~ append_status(status,"MALECNS_ALT_MATCH"),
    TRUE ~ status
  )) %>%
  dplyr::ungroup() %>%
  dplyr::select(-any_of(malecns.cols), -auto_type)
bc.new <- bc.new.malecns

# Integrate hemibrain NBLAST results
banc.meta.hemibrain <- banc.meta %>%
  dplyr::distinct(root_id, .keep_all = TRUE) %>%
  dplyr::select(root_id, dplyr::starts_with("hemibrain")) %>%
  dplyr::select(-hemibrain_nblast_match, -hemibrain_nblast) %>%
  dplyr::filter(!is.na(hemibrain_auto_cell_type)) %>%
  as.data.frame()
hemibrain.cols <- colnames(banc.meta.hemibrain)[grepl("^hemibrain",colnames(banc.meta.hemibrain))]
bc.new.hemibrain <- bc.new %>%
  dplyr::left_join(banc.meta.hemibrain, by = "root_id") %>%
  dplyr::mutate(
    cell_type = dplyr::case_when(
      !is.na(cell_type) ~ cell_type,
      !is.na(hemibrain_match)&(is.na(cell_type)|grepl("auto",cell_type)|cell_type=="") ~ hb.meta$cell_type[match(hemibrain_match,hb.meta$bodyid)],
      TRUE ~ cell_type
    ),
    hemibrain_cell_type = dplyr::case_when(
      !is.na(hemibrain_match)&hemibrain_match!="" ~ hb.meta$cell_type[match(hemibrain_match,hb.meta$bodyid)],      
      is.na(hemibrain_cell_type) ~ hemibrain_auto_cell_type,
      TRUE ~ hemibrain_cell_type
    ),
    auto_type = is.na(cell_type)|grepl("auto",cell_type)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(-all_of(hemibrain.cols), -auto_type)
bc.new <- bc.new.hemibrain

# Integrate FANC NBLAST results
banc.meta.fanc <- banc.meta %>%
  dplyr::distinct(root_id, .keep_all = TRUE) %>%
  dplyr::select(root_id, dplyr::starts_with("fanc")) %>%
  dplyr::select(-fanc_nblast_match, -fanc_nblast) %>%
  as.data.frame()
fanc.cols <- colnames(banc.meta.fanc)[grepl("^fanc",colnames(banc.meta.fanc))]
bc.new.fanc <- bc.new %>%
  dplyr::left_join(banc.meta.fanc, by = "root_id") %>%
  dplyr::mutate(cell_type = dplyr::case_when(
    !is.na(cell_type) ~ cell_type,
    !is.na(fanc_match)&(is.na(cell_type)|grepl("auto",cell_type)) ~ fc.meta$cell_type[match(fanc_match,fc.meta$cell_id)],      
    TRUE ~ cell_type
  ),
  fanc_cell_type = dplyr::case_when(
    !is.na(fanc_match)&fanc_match!="" ~ fc.meta$cell_type[match(fanc_match,fc.meta$cell_id)],  
    is.na(fanc_cell_type) ~ fanc_auto_cell_type,
    TRUE ~ fanc_cell_type
  ),
  auto_type = is.na(cell_type)|grepl("auto",cell_type)
  ) %>%
  # dplyr::mutate(
  #   # hemilineage = dplyr::case_when(
  #   #   !grepl("auto",hemilineage)&!is.na(hemilineage) ~ hemilineage,
  #   #   !auto_type ~ fc.meta$hemilineage[match(fanc_match,fc.meta$cell_id)],
  #   #   is.na(hemilineage) ~ fanc_ito_lee_hemilineage,
  #   #   TRUE ~ hemilineage
  #   # ),
  #   # nerve = dplyr::case_when(
  #   #   !grepl("auto",nerve)&!is.na(nerve) ~ nerve,
  #   #   !auto_type ~ fc.meta$nerve[match(fanc_match,fc.meta$cell_id)],
  #   #   is.na(nerve) ~ fanc_nerve,
  #   #   TRUE ~ nerve
  #   # )
  #   # ),
  #   # fanc_cell_class = dplyr::case_when(
  #   #   !grepl("auto",cell_class)&!is.na(cell_class) ~ cell_class,
  #   #   !auto_type ~ fc.meta$cell_class[match(fanc_match,fc.meta$cell_id)],
  #   #   is.na(cell_class) ~ fanc_cell_class,
  #   #   TRUE ~ cell_class
  #   # )
  # ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(status = dplyr::case_when(
    is.na(cell_type)|is.na(fanc_match)|grepl("neck",region) ~ status,
    cell_type != fc.meta$cell_type[match(fanc_match,fc.meta$cell_id)] ~ append_status(status,"FANC_ALT_MATCH"),
    TRUE ~ status
  )) %>%
  dplyr::ungroup() %>%
  dplyr::select(-any_of(fanc.cols), -auto_type)
bc.new <- bc.new.fanc

######################################
### Collect manual matching issues ###
######################################

# # Detect correct matches
# matches1 <- list.files(file.path(banc.mirror.correct.match.path,"1_perfect"), pattern = "png")
# matches2 <- list.files(file.path(banc.mirror.correct.match.path,"2_confident"), pattern = "png")
# matches3 <- list.files(file.path(banc.mirror.correct.match.path,"3_good"), pattern = "png")
# matches4 <- list.files(file.path(banc.mirror.correct.match.path,"4_likely"), pattern = "png")
# matches5 <- list.files(file.path(banc.mirror.correct.match.path,"5_possible"), pattern = "png")
# matches6 <- list.files(file.path(banc.mirror.correct.match.path,"6_no_match"), pattern = "png")
# blue_wrong <- list.files(file.path(banc.mirror.correct.match.path,"blue_wrong"), pattern = "png")
# red_wrong <- list.files(file.path(banc.mirror.correct.match.path,"red_wrong"), pattern = "png")
# investigate <- list.files(file.path(banc.mirror.correct.match.path,"investigate"), pattern = "png")
# mirror.matches <- unique(c(matches1,
#                            matches2,
#                            matches3,
#                            matches4,
#                            matches5))
# for(file in unique(c(investigate))){
#   root_id <- gsub(".*_root_id_|_nucleus_id_.*","",file)
#   root_id <- banc_updateids(root_id)
#   match.pos <- match(root_id,bc.new$root_id)
#   if(is.na(match.pos)){
#     next
#   }
#   if(!grepl("INVESTIGATE_RESOLVED",bc.new[match.pos,'status'])){
#     bc.new[match.pos,] <- bancr:::banc_update_status(bc.new[match.pos,],
#                                                      update="INVESTIGATE")
#   }
# }
# for(file in red_wrong){
#   try({
#     root_id <- gsub(".*_hit_id_|_hit_nucleus_id_.*","",file)
#     #root_id <- banc_updateids(root_id)
#     match.pos <- match(root_id,bc.new$root_id)
#     if(is.na(match.pos)){
#       next
#     }
#     if(!grepl("TRACING_ISSUE_RESOLVED",bc.new[match.pos,'status'])){
#       bc.new[match.pos,] <- bancr:::banc_update_status(bc.new[match.pos,],
#                                                        update="TRACING_ISSUE")
#     }
#   })
# }
# for(file in blue_wrong){
#   try({
#     root_id <- gsub(".*_root_id_|_nucleus_id_.*","",file)
#     #root_id <- banc_updateids(root_id)
#     match.pos <- match(root_id,bc.new$root_id)
#     if(is.na(match.pos)){
#       next
#     }
#     if(!grepl("TRACING_ISSUE_RESOLVED",bc.new[match.pos,'status'])){
#       bc.new[match.pos,] <- bancr:::banc_update_status(bc.new[match.pos,],
#                                                        update="TRACING_ISSUE")
#     }
#   })
# }
# # for(file in matches6){
# #   try({
# #     root_id <- gsub(".*_root_id_|_nucleus_id_.*","",file)
# #     #root_id <- banc_updateids(root_id)
# #     match.pos <- match(root_id,bc.new$root_id)
# #     if(is.na(match.pos)){
# #       next
# #     }
# #     if(is.na(bc.new[match.pos,'banc_png_match'])|bc.new[match.pos,'banc_png_match']==""){
# #       bc.new[match.pos,"banc_png_match"] <- "no_match"
# #     }
# #   })
# # }
# 
# # Detect correct FAFB matches
# matches1 <- list.files(file.path(banc.fafb.correct.match.path,"1_perfect"), pattern = "png")
# matches2 <- list.files(file.path(banc.fafb.correct.match.path,"2_confident"), pattern = "png")
# matches3 <- list.files(file.path(banc.fafb.correct.match.path,"3_good"), pattern = "png")
# matches4 <- list.files(file.path(banc.fafb.correct.match.path,"4_likely"), pattern = "png")
# matches5 <- list.files(file.path(banc.fafb.correct.match.path,"5_possible"), pattern = "png")
# matches6 <- list.files(file.path(banc.fafb.correct.match.path,"6_no_match"), pattern = "png")
# blue_wrong <- list.files(file.path(banc.fafb.correct.match.path,"blue_wrong"), pattern = "png")
# green_wrong <- list.files(file.path(banc.fafb.correct.match.path,"green_wrong"), pattern = "png")
# red_wrong <- list.files(file.path(banc.fafb.correct.match.path,"red_wrong"), pattern = "png")
# investigate <- list.files(file.path(banc.fafb.correct.match.path,"investigate"), pattern = "png")
# fafb.matches <- unique(c(matches1,
#                          matches2,
#                          matches3,
#                          matches4))
# for(file in unique(c(investigate))){
#   try({
#     root_id <- unlist(strsplit(gsub(".*_root_id_|_nucleus_id_.*","",file),split = "_"))
#     root_id <- banc_updateids(root_id)
#     match.pos <- match(root_id,bc.new$root_id)
#     match.pos <- match.pos[!is.na(match.pos)]
#     if(!length(match.pos)){
#       next
#     }
#     if(!grepl("INVESTIGATE_RESOLVED",bc.new[match.pos,'status'])){
#       bc.new[match.pos,] <- bancr:::banc_update_status(bc.new[match.pos,],
#                                                        update="INVESTIGATE")
#     }
#   })
# }
# for(file in c(blue_wrong,green_wrong)){
#   try({
#     root_id <- unlist(strsplit(gsub(".*_root_id_|_nucleus_id_.*","",file),split = "_"))
#     #root_id <- banc_updateids(root_id)
#     match.pos <- match(root_id,bc.new$root_id)
#     match.pos <- match.pos[!is.na(match.pos)]
#     if(!length(match.pos)){
#       next
#     }
#     if(!grepl("TRACING_ISSUE_RESOLVED",bc.new[match.pos,'status'])){
#       bc.new[match.pos,] <- bancr:::banc_update_status(bc.new[match.pos,],
#                                                        update="TRACING_ISSUE")
#     }
#   })
# }
# for(file in matches6){
#   root_id <- unlist(strsplit(gsub(".*_root_id_|_nucleus_id_.*","",file),split = "_"))
#   #root_id <- banc_updateids(root_id)
#   match.pos <- match(root_id,bc.new$root_id)
#   match.pos <- match.pos[!is.na(match.pos)]
#   if(!length(match.pos)){
#     next
#   }
#   for(match.p in match.pos){
#     if(is.na(bc.new[match.p,'fafb_png_match'])|bc.new[match.p,'fafb_png_match']==""){
#       bc.new[match.p,"fafb_png_match"] <- "no_match"
#     } 
#   }
# }
# 
# # Detect correct MANC matches
# matches1 <- list.files(file.path(banc.manc.correct.match.path,"1_perfect"), pattern = "png")
# matches2 <- list.files(file.path(banc.manc.correct.match.path,"2_confident"), pattern = "png")
# matches3 <- list.files(file.path(banc.manc.correct.match.path,"3_good"), pattern = "png")
# matches4 <- list.files(file.path(banc.manc.correct.match.path,"4_likely"), pattern = "png")
# matches5 <- list.files(file.path(banc.manc.correct.match.path,"5_possible"), pattern = "png")
# matches6 <- list.files(file.path(banc.manc.correct.match.path,"6_no_match"), pattern = "png")
# blue_wrong <- list.files(file.path(banc.manc.correct.match.path,"blue_wrong"), pattern = "png")
# green_wrong <- list.files(file.path(banc.manc.correct.match.path,"green_wrong"), pattern = "png")
# red_wrong <- list.files(file.path(banc.manc.correct.match.path,"red_wrong"), pattern = "png")
# investigate <- list.files(file.path(banc.manc.correct.match.path,"investigate"), pattern = "png")
# manc.matches <- unique(c(matches1,
#                          matches2,
#                          matches3,
#                          matches4))
# for(file in unique(c(investigate))){
#   try({
#     
#   })
#   root_id <- unlist(strsplit(gsub(".*_root_id_|_nucleus_id_.*","",file),split = "_"))
#   root_id <- banc_updateids(root_id)
#   match.pos <- match(root_id,bc.new$root_id)
#   match.pos <- match.pos[!is.na(match.pos)]
#   if(!length(match.pos)){
#     next
#   }
#   if(!grepl("INVESTIGATE_RESOLVED",bc.new[match.pos,'status'])){
#     bc.new[match.pos,] <- bancr:::banc_update_status(bc.new[match.pos,],
#                                                      update="INVESTIGATE")
#   }
# }
# for(file in c(blue_wrong,green_wrong)){
#   root_id <- unlist(strsplit(gsub(".*_root_id_|_nucleus_id_.*","",file),split = "_"))
#   #root_id <- banc_updateids(root_id)
#   match.pos <- match(root_id,bc.new$root_id)
#   match.pos <- match.pos[!is.na(match.pos)]
#   if(!length(match.pos)){
#     next
#   }
#   if(!grepl("TRACING_ISSUE_RESOLVED",bc.new[match.pos,'status'])){
#     bc.new[match.pos,] <- bancr:::banc_update_status(bc.new[match.pos,],
#                                                      update="TRACING_ISSUE")
#   }
# }
# # for(file in matches6){
# #   root_id <- unlist(strsplit(gsub(".*_root_id_|_nucleus_id_.*","",file),split = "_"))
# #   #root_id <- banc_updateids(root_id)
# #   match.pos <- match(root_id,bc.new$root_id)
# #   match.pos <- match.pos[!is.na(match.pos)]
# #   if(!length(match.pos)){
# #     next
# #   }
# #   for(match.p in match.pos){
# #     if(is.na(bc.new[match.p,'manc_png_match'])|bc.new[match.p,'manc_png_match']==""){
# #       bc.new[match.p,"manc_png_match"] <- "no_match"
# #     } 
# #   }
# # }
# 
# # Detect correct hemibrain matches
# matches1 <- list.files(file.path(banc.hemibrain.correct.match.path,"1_perfect"), pattern = "png")
# matches2 <- list.files(file.path(banc.hemibrain.correct.match.path,"2_confident"), pattern = "png")
# matches3 <- list.files(file.path(banc.hemibrain.correct.match.path,"3_good"), pattern = "png")
# matches4 <- list.files(file.path(banc.hemibrain.correct.match.path,"4_likely"), pattern = "png")
# matches5 <- list.files(file.path(banc.hemibrain.correct.match.path,"5_possible"), pattern = "png")
# matches6 <- list.files(file.path(banc.hemibrain.correct.match.path,"6_no_match"), pattern = "png")
# blue_wrong <- list.files(file.path(banc.hemibrain.correct.match.path,"blue_wrong"), pattern = "png")
# green_wrong <- list.files(file.path(banc.hemibrain.correct.match.path,"green_wrong"), pattern = "png")
# red_wrong <- list.files(file.path(banc.hemibrain.correct.match.path,"red_wrong"), pattern = "png")
# investigate <- list.files(file.path(banc.hemibrain.correct.match.path,"investigate"), pattern = "png")
# hemibrain.matches <- unique(c(matches1,
#                               matches2,
#                               matches3,
#                               matches4))
# for(file in unique(c(investigate))){
#   try({
#     
#   })
#   root_id <- unlist(strsplit(gsub(".*_root_id_|_nucleus_id_.*","",file),split = "_"))
#   root_id <- banc_updateids(root_id)
#   match.pos <- match(root_id,bc.new$root_id)
#   match.pos <- match.pos[!is.na(match.pos)]
#   if(!length(match.pos)){
#     next
#   }
#   if(!grepl("INVESTIGATE_RESOLVED",bc.new[match.pos,'status'])){
#     bc.new[match.pos,] <- bancr:::banc_update_status(bc.new[match.pos,],
#                                                      update="INVESTIGATE")
#   }
# }
# for(file in c(blue_wrong,green_wrong)){
#   root_id <- unlist(strsplit(gsub(".*_root_id_|_nucleus_id_.*","",file),split = "_"))
#   #root_id <- banc_updateids(root_id)
#   match.pos <- match(root_id,bc.new$root_id)
#   match.pos <- match.pos[!is.na(match.pos)]
#   if(!length(match.pos)){
#     next
#   }
#   if(!grepl("TRACING_ISSUE_RESOLVED",bc.new[match.pos,'status'])){
#     bc.new[match.pos,] <- bancr:::banc_update_status(bc.new[match.pos,],
#                                                      update="TRACING_ISSUE")
#   }
# }
# # for(file in matches6){
# #   root_id <- unlist(strsplit(gsub(".*_root_id_|_nucleus_id_.*","",file),split = "_"))
# #   #root_id <- banc_updateids(root_id)
# #   match.pos <- match(root_id,bc.new$root_id)
# #   match.pos <- match.pos[!is.na(match.pos)]
# #   if(!length(match.pos)){s
# #     next
# #   }
# #   for(match.p in match.pos){
# #     if(is.na(bc.new[match.p,'hemibrain_png_match'])|bc.new[match.p,'hemibrain_png_match']==""){
# #       bc.new[match.p,"hemibrain_png_match"] <- "no_match"
# #     } 
# #   }
# # }

##############
### Issues ###
##############

bc.new[bc.new=="NA,NA,NA"] <- ""
duplicates <- unique(bc.new$root_id[duplicated(bc.new$root_id)])
bc.new <- bc.new %>%
  ### cell name uniformity ###
  dplyr::mutate(cell_type = dplyr::case_when(
    !is.na(cell_type) ~ cell_type,
    !is.na(manc_cell_type)&grepl("neck",region)&grepl("ascending|sensory",super_class)&!grepl("auto",manc_cell_type) ~ manc_cell_type,
    !is.na(fafb_cell_type)&grepl("neck",region)&grepl("descending",super_class)&!grepl("auto",fafb_cell_type) ~ fafb_cell_type,
    !is.na(fafb_cell_type)&!grepl("auto",fafb_cell_type) ~ fafb_cell_type,
    !is.na(manc_cell_type)&!grepl("auto",manc_cell_type) ~ manc_cell_type,
    !is.na(hemibrain_cell_type)&!grepl("auto",hemibrain_cell_type) ~ hemibrain_cell_type,
    !is.na(fanc_cell_type)&!grepl("auto",fanc_cell_type) ~ fanc_cell_type,
    !is.na(fafb_cell_type) ~ fafb_cell_type,
    !is.na(manc_cell_type) ~ manc_cell_type,
    !is.na(hemibrain_cell_type) ~ hemibrain_cell_type,
    !is.na(fanc_cell_type) ~ fanc_cell_type,
    TRUE ~ cell_type
  )) %>%
  ### Status wrangling ###
  dplyr::mutate(root_position_nm = dplyr::case_when(
    is.na(root_position_nm)|root_position_nm=="" ~ nucleus_position_nm,
    TRUE ~ root_position_nm
  )) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(status = dplyr::case_when(
    # Check if the id exists in both root_id and banc_match
    root_id %in% banc_match | banc_match %in% root_id ~ 
      dplyr::case_when(
        # If it exists in both, check if the cell types match
        cell_type != cell_type[match(root_id, banc_match)] |
          cell_type != cell_type[match(banc_match, root_id)] ~ append_status(status,"LR_TYPE_CONFLICT"),
        # If none of the above conditions are met
        TRUE ~ status
      ),
    # If the id doesn't exist in both columns
    TRUE ~ status
  )) %>%
  dplyr::mutate(status = dplyr::case_when(
    is.na(banc_match)|banc_match%in%c("NA",""," ") ~ status,
    side == .data$side[match(banc_match,.data$root_id)] ~ append_status(status,"SIDE_CONFLICT"),
    side != .data$side[match(banc_match,.data$root_id)] ~ subtract_status(status,"SIDE_CONFLICT"),
    TRUE ~ status
  )) %>%
  dplyr::mutate(status = dplyr::case_when(
    is.na(root_position_nm)|root_position_nm=="" ~ append_status(status,"UNROOTED"),
    !is.na(root_position_nm)|root_position_nm!=""&grepl("UNROOTED",status) ~ subtract_status(status,"UNROOTED"),
    TRUE ~ status
  )) %>%
  dplyr::mutate(status = dplyr::case_when(
    grepl("TRACING_ISSUE_RESOLVED",status) ~ subtract_status(status,"TRACING_ISSUE"),
    TRUE ~ status
  )) %>%
  dplyr::mutate(status = dplyr::case_when(
    !is.na(banc_match)|banc_match!='' ~ subtract_status(status,"NO_MIRROR_MATCH"),
    TRUE ~ status
  )) %>%
  dplyr::mutate(status = dplyr::case_when(
    !is.na(fafb_match)|fafb_match!='' ~ subtract_status(status,"NO_FAFB_MATCH"),
    TRUE ~ status
  )) %>%
  dplyr::mutate(status = dplyr::case_when(
    !is.na(fanc_match)|fanc_match!='' ~ subtract_status(status,"NO_FANC_MATCH"),
    TRUE ~ status
  )) %>%
  dplyr::mutate(status = dplyr::case_when(
    !is.na(manc_match)|manc_match!='' ~ subtract_status(status,"NO_MANC_MATCH"),
    TRUE ~ status
  )) %>%
  dplyr::mutate(status = dplyr::case_when(
    !is.na(hemibrain_match)|hemibrain_match!='' ~ subtract_status(status,"NO_HEMIBRAIN_MATCH"),
    TRUE ~ status
  )) %>%
  dplyr::mutate(status = dplyr::case_when(
    grepl("CONFLICT_RESOLVED",status) ~ subtract_status(status,"FAFB_TYPE_CONFLICT"),
    TRUE ~ status
  )) %>%
  dplyr::mutate(status = dplyr::case_when(
    grepl("CONFLICT_RESOLVED",status) ~ subtract_status(status,"MANC_TYPE_CONFLICT"),
    TRUE ~ status
  )) %>%
  dplyr::mutate(status = dplyr::case_when(
    grepl("CONFLICT_RESOLVED",status) ~ subtract_status(status,"HEMIBRAIN_TYPE_CONFLICT"),
    TRUE ~ status
  )) %>%
  dplyr::mutate(status = dplyr::case_when(
    grepl("CONFLICT_RESOLVED",status) ~ subtract_status(status,"FANC_TYPE_CONFLICT"),
    TRUE ~ status
  )) %>%
  # dplyr::mutate(status = dplyr::case_when(
  #   !root_id %in% duplicates ~ subtract_status(status,"DUPLICATED"),
  #   root_id %in% duplicates ~ append_status(status,"DUPLICATED"),
  #   TRUE ~ status
  # )) %>%
  dplyr::mutate(status = dplyr::case_when(
    grepl("NOT_TOO_SMALL|TOO_SMALL_RESOLVED",status) ~ subtract_status(status,"TOO_SMALL"),
    is.na(l2_cable_length_um) ~ status,
    l2_cable_length_um > 50 ~ subtract_status(status,"TOO_SMALL"),
    l2_cable_length_um <= 50 ~ append_status(status,"TOO_SMALL"),
    TRUE ~ status
  )) %>%
  dplyr::mutate(status = dplyr::case_when(
    is.na(fafb_match)|is.na(side) ~ status,
    grepl("ascending|sensory_ascending|ascending, sensory neuron",cell_class) ~ subtract_status(status,"FAFB_MATCH_WRONG_SIDE"),
    side == fw.meta[match(fafb_match,fw.meta$root_783),"side"]  ~ subtract_status(status,"FAFB_MATCH_WRONG_SIDE"),
    side != fw.meta[match(fafb_match,fw.meta$root_783),"side"]  ~ append_status(status,"FAFB_MATCH_WRONG_SIDE"),
    TRUE ~ status
  )) %>%
  dplyr::mutate(status = dplyr::case_when(
    is.na(manc_match)|is.na(side) ~ status,
    grepl("descending",cell_class) ~ subtract_status(status,"MANC_MATCH_WRONG_SIDE"),
    side == mc.meta[match(manc_match,mc.meta$bodyid),"side"]  ~ subtract_status(status,"MANC_MATCH_WRONG_SIDE"),
    side != mc.meta[match(manc_match,mc.meta$bodyid),"side"]  ~ append_status(status,"MANC_MATCH_WRONG_SIDE"),
    TRUE ~ status
  )) %>%
  dplyr::ungroup()
# dplyr::mutate(ngl_link = dplyr::case_when(
#   status!=""&region%in%c("neck_connective")&(is.na(ngl_link)|ngl_link=="") ~ bancsee(banc_ids = c(root_id,banc_match),
#                                                       fafb_ids = c(fafb_match),
#                                                       manc_ids = c(manc_match),
#                                                       nuclei_ids = c(nucleus_id)),
#   TRUE ~ ngl_link
# )) %>%

# Update sides
root_no_side <- (!is.na(bc.new$root_position_nm)&bc.new$root_position_nm!="")&(bc.new$side==""|is.na(bc.new$side))
root_no_side[is.na(root_no_side)] <- FALSE
if(sum(root_no_side)){
  roots <- nat::xyzmatrix(bc.new$root_position_nm[root_no_side])
  lrdiffs <- bancr:::banc_lr_position(roots,units = "nm")
  sides <- ifelse(lrdiffs>0,"right","left")
  bc.new$side[root_no_side] <- sides  
}

# Correct nerves
# banc.peripheral.nerves <- banc_peripheral_nerves() %>%
#   dplyr::mutate(pt_root_id=as.character(pt_root_id)) %>%
#   dplyr::rename(nerve=tag) %>%
#   dplyr::distinct(pt_root_id, .keep_all = TRUE)
# bc.new <- bc.new %>%
#   dplyr::left_join(banc.peripheral.nerves[,c("pt_root_id","nerve")], 
#                    by = c("root_id"="pt_root_id")) %>%
#   dplyr::mutate(nerve = dplyr::case_when(
#     is.na(nerve.x) & !is.na(nerve.y) ~ nerve.y,
#     !is.na(nerve.x) & is.na(nerve.y) ~ nerve.x,
#     !is.na(nerve.x) & !is.na(nerve.y) ~ nerve.y,
#     TRUE ~ nerve.x
#   )) %>%
#   dplyr::select(-nerve.x,nerve.y)

####################
### maleCNS type ###
####################

# malecns meta
malecns.meta <- readr::read_csv(file.path(banc.meta.save.path,"malecns_09_meta.csv"), col_types = banc.col.types) %>%
  dplyr::select(
    malecns_09_id,
    cell_type,
    mcns_fafb_cell_type = fafb_cell_type,
    mcns_hemibrain_cell_type = hemibrain_cell_type,
    mcns_manc_cell_type = manc_cell_type,
    mcns_manc_match = manc_match
  ) %>%
  dplyr::mutate(cell_type=dplyr::coalesce(cell_type,mcns_fafb_cell_type,mcns_hemibrain_cell_type,mcns_manc_cell_type)) %>%
  dplyr::filter(!is.na(cell_type))
bc.new$malecns_cell_type <- NA
bc.new <- bc.new %>%
  dplyr::left_join(malecns.meta %>%
                     dplyr::filter(!is.na(mcns_fafb_cell_type)) %>%
                     dplyr::distinct(mcns_fafb_cell_type, .keep_all = TRUE) %>%
                     dplyr::select(malecns_cell_type_from_fafb = cell_type,mcns_fafb_cell_type),
                   by = c("fafb_cell_type"="mcns_fafb_cell_type")) %>%
  dplyr::left_join(malecns.meta %>%
                     dplyr::filter(!is.na(mcns_hemibrain_cell_type)) %>%
                     dplyr::distinct(mcns_hemibrain_cell_type, .keep_all = TRUE) %>%
                     dplyr::select(malecns_cell_type_from_hemibrain = cell_type,mcns_hemibrain_cell_type),
                   by = c("hemibrain_cell_type"="mcns_hemibrain_cell_type")) %>%
  dplyr::left_join(malecns.meta %>%
                     dplyr::filter(!is.na(mcns_manc_cell_type)) %>%
                     dplyr::distinct(mcns_manc_cell_type, .keep_all = TRUE) %>%
                     dplyr::select(malecns_cell_type_from_manc = cell_type,mcns_manc_cell_type),
                   by = c("manc_cell_type"="mcns_manc_cell_type")) %>%
  dplyr::left_join(malecns.meta %>%
                     dplyr::filter(!is.na(mcns_manc_match)) %>%
                     dplyr::distinct(mcns_manc_match, .keep_all = TRUE) %>%
                     dplyr::select(malecns_cell_type_from_manc_id = cell_type,mcns_manc_match),
                   by = c("manc_match"="mcns_manc_match")) %>%
  dplyr::mutate(malecns_cell_type=dplyr::coalesce(malecns_cell_type_from_fafb,malecns_cell_type_from_manc,malecns_cell_type_from_manc_id,malecns_cell_type_from_hemibrain)) %>%
  dplyr::mutate(malecns_cell_type = dplyr::case_when(
    cell_type %in% !!malecns.meta$cell_type ~ cell_type,
    TRUE ~ malecns_cell_type
  ))

#######################
### Update seatable ###
#######################

# Some clean up
bc.new$status<-gsub("NA\\,|\\,NA","",bc.new$status)
bc.new$l2_cable_length_um <- as.numeric(bc.new$l2_cable_length_um)
bc.new$l2_cable_length_um[is.na(bc.new$l2_cable_length_um )] <- 0

# Update
bc.new$banc_nblast<-round(as.numeric(bc.new$banc_nblast),3)
bc.update <- bc.new %>%
  dplyr::filter(root_id != "0") %>%
  dplyr::mutate(across(where(is.character), ~replace_na(., "")),
                across(where(is.factor), ~as.factor(replace_na(as.character(.), ""))),
                across(where(is.numeric), ~replace_na(., 0)),
                across(where(is.logical), ~.)) %>%
  dplyr::select(# fully automated          
    `_id`, root_id, position, supervoxel_id, # identity 
    nucleus_id, nucleus_position, nucleus_position_nm, nucleus_supervoxel_id, # nucleus
    # root_region,
    proofread, # metrics
    #input_connections, output_connections, input_side_index, output_side_index,
    #l2_nodes, l2_cable_length_um, # metric2
    # region, 
    side, #manc_super_class, fafb_super_class,# fafb_cell_class, manc_cell_class# anatomy                            
    # top_nt, top_nt_conf,  # useful meta 
    cell_type,
    malecns_cell_type, malecns_png_match, malecns_nblast_match, malecns_nblast, # maleCNS matching data
    manc_cell_type, manc_png_match, manc_nblast_match, manc_nblast, # MANC matching data
    fafb_cell_type, fafb_png_match, fafb_nblast_match, fafb_nblast, # FAFB matching data
    hemibrain_cell_type, hemibrain_png_match, hemibrain_nblast_match, hemibrain_nblast, # hemibrain matching data
    fanc_cell_type, fanc_png_match, fanc_nblast_match, fanc_nblast, # FANC matching data
    banc_match_supervoxel_id, banc_png_match, banc_png_match_supervoxel_id, banc_nblast_match_supervoxel_id, banc_nblast_match, banc_nblast, # BANC matching data
    # semi-automated
    # other_names, cell_type_source, # scrounge CAVE info
    # root_position, root_position_nm,
    # status, hemilineage, nerve, 
    # cell_type, # useful active meta, super_class, cell_sub_class       
    banc_match, fafb_match, manc_match, fanc_match, hemibrain_match, malecns_match # 1:1 neuron matches
  ) %>%
  dplyr::mutate(cell_type = ifelse(grepl("auto",cell_type),NA,cell_type))# NOT: super_class, cell_class, cell_sub_class

# Has anything changes in the meantime? If so, do not touch that row.
bc2 <- banctable_query()
bc2[,drop.cols] <- NULL
updated.cols <- intersect(colnames(bc2),colnames(bc.update))
bc2 <- bc2[,updated.cols]
bc2 <- bc2[,intersect(colnames(bc2),colnames(bc1))]
bckeep <- bc2 %>%
  dplyr::semi_join(bc1, by = names(bc2))
bckeep <- bckeep$`_id`
bckeep2 <- bc.orig[,intersect(colnames(bc.orig),colnames(bc.new))] %>%
  dplyr::anti_join(bc.new, by = c("status","fafb_match","fanc_match","manc_match","hemibrain_match","malecns_match"))
bckeep2 <- bckeep2$`_id`
bckeep <- intersect(bckeep,bckeep2)

# # Use only the columns we need from bc2: _id and its current cell_type
# bc2_ct <- bc2 %>%
#   dplyr::select(`_id`, cell_type_bc2 = cell_type)
# 
# # Join, compare cell_type, and set cell_type_source = "bates" where needed
# bc.update <- bc.update %>%
#   dplyr::left_join(bc2_ct, by = "_id") %>%
#   dplyr::mutate(
#     cell_type_source = dplyr::case_when(
#       is.na(cell_type_bc2) ~ append_status(cell_type_source,"bates"),
#       is.na(cell_type) ~ cell_type_source,
#       cell_type_bc2 != cell_type ~ append_status(cell_type_source,"bates"),
#       TRUE ~ cell_type_source
#     )
#   ) %>%
#   dplyr::select(-cell_type_bc2)

# Free space by moving optic stuff to big data storage
# banctable_move_to_bigdata(where = "`region` = 'optic'")

# Update already recorded
bc.update.present <- bc.update %>%
  dplyr::filter(`_id`%in%bckeep, 
                `_id`!="")%>%
  as.data.frame()
bc.update.present[is.na(bc.update.present)] <- ""
bc.update.present <- bc.update.present %>%
  dplyr::anti_join(bc2, by = updated.cols)
if(nrow(bc.update.present)){
  message("Modifying ", nrow(bc.update.present), " extant neurons in banctable")
  banctable_update_rows(base='banc_meta', 
                        table = "banc_meta", 
                        df = bc.update.present, 
                        append_allowed = FALSE, 
                        chunksize = 1000)  
}

# Append new rows
new.rows <- bc.update %>%
  dplyr::filter(bc.update$`_id`=="",
                !root_id %in% bc.update.present$root_id,
                !supervoxel_id %in% bc.update.present$supervoxel_id,
                !root_id %in% bc2$root_id,
                !supervoxel_id %in% bc2$supervoxel_id) %>%
  dplyr::distinct(root_id, .keep_all = TRUE)
new.rows <- as.data.frame(new.rows)
new.rows$`_id` <- NULL
new.rows[is.na(new.rows)] <- ""
new.rows$status <- NULL
if(nrow(new.rows)){
  message("Adding ", nrow(new.rows), " new neurons to banctable")
  banctable_append_rows(base='banc_meta',
                        table = "banc_meta",
                        df = new.rows,
                        chunksize = 1000)  
}

# Correct nerve
# # banc.peripheral.nerves <- banc_peripheral_nerves()
# n.ids <- unique(banc.peripheral.nerves$pt_root_id)
# banc.peripheral.nerves$nerve = banc.peripheral.nerves$tag
# banc.peripheral.nerves$pt_root_id = as.character(banc.peripheral.nerves$pt_root_id)
# new.rows <- new.rows %>%
#   dplyr::filter(root_id %in% n.ids) %>%
#   dplyr::select(-nerve) %>%
#   dplyr::left_join(banc.peripheral.nerves[,c("pt_root_id","nerve")],
#                    by=c("root_id"="pt_root_id")) %>%
#   dplyr::mutate(super_class = "sensory") %>%
#   dplyr::mutate(cell_class = dplyr::case_when(
#     !is.na(cell_class) ~ "sensory neuron",
#     TRUE ~ cell_class
#   )) %>%
#   dplyr::rowwise() %>%
#   dplyr::mutate(root_position_nm = dplyr::case_when(
#     !is.na(position) ~ paste(bancr::banc_raw2nm(position),collapse=","),
#     TRUE ~ NA
#   )) %>%
#   dplyr::ungroup() %>%
#   as.data.frame()
# roots <- nat::xyzmatrix(new.rows$root_position_nm)
# lrdiffs <- bancr:::banc_lr_position(roots,units = "nm")
# sides <- ifelse(lrdiffs>0,"right","left")
# new.rows$side <- sides

# Write
readr::write_csv(bc.update, file.path(banc.meta.save.path,"banc_seatable.csv"))

# Notes:
# Maybe incorrect: MANC 27125 for 720575941542610437
