#!/usr/bin/env Rscript
source("banc/banc-startup.R")
version <- "746"

# File paths
save.path <- '/n/data1/hms/neurobio/wilson/connectomes/banc'
dir.create(save.path)

####################
### GET SYNAPSES ###
####################

# csv info
desired_columns <- c('id', 'size', 'pre_root_id', 'post_root_id', 'ctr_x', 'ctr_y', 'ctr_z')
col_types <- cols(
  id = col_character(),
  size = col_double(),
  pre_root_id = col_character(),
  post_root_id = col_character(),
  ctr_x = col_double(),
  ctr_y = col_double(),
  ctr_z = col_double(),
  .default = col_double()
)
column_names <- c('id', 'pre_x', 'pre_y', 'pre_z', 'post_x', 'post_y', 'post_z',
                  'ctr_x', 'ctr_y', 'ctr_z', 'size', 'pre_supervoxel_id',
                  'pre_root_id', 'post_supervoxel_id', 'post_root_id')

# # Get the synapse file remotely, you have it locally
version.path <- file.path(banc.save.path,"v746")
bancsynapses <- file.path(version.path,"synapses_v2_human_readable.csv")
if(!file.exists(bancsynapses)){
  dir.create(version.path)
  system(sprintf("gsutil cp gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/v746/synapses_v2_human_readable.csv.gz %s ",
                 file.path(version.path,"synapses_v2_human_readable.csv.gz")))
  system(sprintf("gunzip %s",file.path(version.path,"synapses_v2_human_readable.csv.gz")))
}

# Get BANC meta
banc.meta <- banctable_query() 
franken.meta <- franken_meta()
banc.meta$banc_746_id <- banc_rootid(banc.meta$supervoxel_id, version = "746")

############
### META ###
############

## meta data: full meta data for useful labels, both anatomy (e.g. cell type, cell class) and function (i.e. body parts innervated, sparse known cell functions)
banc.meta <- banc.meta %>%
  dplyr::arrange(cell_type, cell_sub_class, cell_class, super_class, hemilineage, nerve) %>%
  dplyr::filter(!grepl("glia|trachea|not_a_neuron",super_class)&!grepl("glia|trachea|not_a_neuron",cell_class)&!grepl("GLIA|TRACHEA|DEBRIS|NOT_A_NEURON",status)) %>%
  dplyr::distinct(banc_746_id, 
                  .keep_all = TRUE) %>%
  dplyr::select(banc_746_id,
                supervoxel_id,
                region,
                side,
                hemilineage,
                nerve,
                flow,
                super_class,
                cell_class,
                cell_sub_class,
                cell_type,
                neurotransmitter_predicted,
                neurotransmitter_score,
                body_part_sensory,
                body_part_effector,
                cell_function,
                cell_function_detailed,
                status) 
# %>%
#   dplyr::left_join(franken.meta %>%
#                      dplyr::filter(!is.na(cell_type)) %>%
#                      dplyr::distinct(cell_type,.keep_all = TRUE) %>%
#                      dplyr::select(cell_type, 
#                                      cell_function,
#                                      cell_function_detailed),
#                    by = "cell_type")
root.ids <- banc.meta$banc_746_id
arrow::write_feather(banc.meta, file.path(save.path,"banc_746_meta.feather"))

################
### SYNAPSES ###
################

# csv info
desired_columns <- c('id', 'size', 'pre_root_id', 'post_root_id', 'ctr_x', 'ctr_y', 'ctr_z')
col_types <- cols(
  id = col_character(),
  size = col_double(),
  pre_root_id = col_character(),
  post_root_id = col_character(),
  ctr_x = col_double(),
  ctr_y = col_double(),
  ctr_z = col_double(),
  .default = col_double()
)
column_names <- c('id', 'pre_x', 'pre_y', 'pre_z', 'post_x', 'post_y', 'post_z',
                  'ctr_x', 'ctr_y', 'ctr_z', 'size', 'pre_supervoxel_id',
                  'pre_root_id', 'post_supervoxel_id', 'post_root_id')

# Get BANC synapses
banc.syns <- vroom::vroom(bancsynapses,
                          col_names = column_names,
                          col_select = dplyr::all_of(desired_columns),
                          col_types = col_types,
                          skip = 0) %>%
  dplyr::rename(x=ctr_x,
                y=ctr_y,
                z=ctr_z,
                pre = pre_root_id,
                post = post_root_id) %>%
  dplyr::filter(pre != post,
                pre %in% !!root.ids | post %in% !!root.ids) %>%
  dplyr::distinct(id, .keep_all  = TRUE) %>%
  tibble::as_tibble()

# banc NPs
syn.ids <- unique(banc.syns$id)
banc.nps <- vroom::vroom(file.path(banc.connectivity.save.path,"banc_synapses_to_neuropils_v2.csv"),
                         col_names = c("id", "X", "Y", "Z", "side", "region", "neuropil"),
                         col_select = dplyr::all_of(c("id","side", "region", "neuropil")),
                         col_types = cols(
                           id = col_character(),
                           X = col_double(),
                           Y = col_double(),
                           Z = col_double(),
                           .default = col_character()
                         ),
                         skip = 1) %>%
  dplyr::distinct(id, .keep_all  = TRUE) %>%
  dplyr::filter(id %in% !!syn.ids) %>%
  tibble::as_tibble()

# # Join
# banc.syns <- banc.syns %>%
#   dplyr::left_join(banc.nps, by = "id")
chunk_size <- 2000000L
n <- nrow(banc.syns)
n_chunks <- ceiling(n / chunk_size)
indexes <- split(seq_len(n), ceiling(seq_len(n) / chunk_size))
# function to process one chunk index vector
process_chunk <- function(idxs) {
  chunk <- banc.syns[idxs, , drop = FALSE]               # slice rows
  sids <- unique(chunk$id)
  bnps <- banc.nps %>%
    dplyr::filter(id %in% !!sids) 
  # perform left join for this chunk only
  res <- dplyr::left_join(chunk, bnps, by = "id")
  # free chunk and run gc
  rm(chunk); gc()
  res
}
# run over chunks and row-bind results
banc.syns <- purrr::map_dfr(indexes, process_chunk, .progress = TRUE)
gc()

# Nt predictions
bancsynapses.nts <- "/n/data3_vast/hms/neurobio/htem2/users/kd193/banc_neurotransmitter_pred/banc_nt_inference/results/banc_nt_pred_v2/banc_nt_prediction_w_sizethresh_5_11102025.parquet"
banc.nt <- arrow::read_parquet(bancsynapses.nts)
banc.syns <- dplyr::left_join(banc.syns,
                              banc.nt  %>%
                                dplyr::mutate(id = as.character(id)) %>%
                                dplyr::rename(syn_top_nt = predicted_nt, syn_top_p = probability),
                              by="id")

# Save
pre.ids.missing <- unique(setdiff(root.ids,banc.syns$pre))
post.ids.missing <- unique(setdiff(root.ids,banc.syns$post))
arrow::write_parquet(banc.syns, file.path(save.path,"banc_746_synapses.parquet"),
                     version = "2.6", 
                     compression = "snappy",
                     compression_level = NULL,
                     chunk_size = 100000,
                     use_dictionary = TRUE,
                     allow_truncated_timestamps = FALSE)

################
### EDGELIST ###
################

## simple edgelists: proofread neurons and connections between them
banc.elist.simp <- banc.syns %>%
  dplyr::group_by(post) %>%
  dplyr::mutate(total_input = dplyr::n()) %>%
  dplyr::group_by(pre, post) %>%
  dplyr::mutate(count = dplyr::n(),
                norm = count/total_input) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(pre,post,count,norm,total_input)
pre.ids.missing <- unique(setdiff(root.ids,banc.elist.simp$pre))
post.ids.missing <- unique(setdiff(root.ids,banc.elist.simp$post))
arrow::write_feather(banc.elist.simp, file.path(save.path,"banc_746_simple_edgelist.feather"))

#################
### SKELETONS ###
#################

## skeletons: including compartmentalised skeletons (optional: axons vs. dendrites)

### In BANC space
banc.in.banc.space <- banc.l2swc.save.path
wd <- getwd()
# out_zip <- file.path(save.path, "banc_banc_space_l2_swc.zip")
# oldwd <- setwd(banc.in.banc.space); on.exit(setwd(oldwd), add = TRUE)
# files <- list.files(".", recursive = TRUE, all.files = TRUE, no.. = TRUE)
# zip::zipr(zipfile = out_zip, files = files)
# cat("Wrote:", out_zip, "\n")
# setwd(wd)
new.folder <- file.path(save.path, "banc_banc_space_l2_swc")
dir.create(new.folder)
files <- list.files(banc.in.banc.space, full.names = TRUE)
file.copy(from = files, to = new.folder, overwrite = FALSE)
files <- list.files(new.folder,full.names = TRUE)
toremove <- files[!gsub("\\.swc","",basename(files))%in%root.ids]
file.remove(toremove)

################
### CUT OUTS ###
################

# Go by cut-out
cut.outs <- c("mushroom_body","antennal_lobe","central_complex","optic","suboesophageal_zone","front_leg","abdominal_neuromere")
for(cut.out in cut.outs){
  
  # Save path
  cut.out.good <- snakecase::to_snake_case(cut.out)
  save.path.cut.out <- file.path(save.path,cut.out.good)
  dir.create(save.path.cut.out)
  if(cut.out=="optic"){
    cut.out = "optic_lobe_intrinsic|optic_lobe_sensory|visual_project_visual_centrifugal|optic_lobe_bilateral"
    cut.out = "^LO|^LOP|^AME|^ME"
  }
  if(cut.out=="antennal_lobe"){
    cut.out = "antennal_lobe|olfactory_receptor|thermosensory_receptor|hygrosensory_receptor|CSD"
  }
  if(cut.out=="suboesophageal_zone"){
    cut.out = "^FLA|^SEZ|^GNG|^SAD|^AMMC|^PRW"
  }
  if(cut.out=="front_leg"){
    cut.out = "^LegNp\\(T1\\)|T1|^ProNM-T1|^LNp_T1"
  }
  if(cut.out=="abdominal_neuromere"){
    cut.out = "^ANm|^ABDNM"
  }
  
  # Custom cuts
  if(cut.out=="mushroom_body"){
    cut.out = "mushroom_body|kenyon_cell|APL|DPM|LHMB1|OA-VPM3"
    kc.ids <- banc.meta %>%
      dplyr::filter(cell_class=="kenyon_cell",side=="right") %>%
      dplyr::pull(banc_746_id)
    kc.elist.simp <- banc.elist.simp %>%
      dplyr::filter(pre%in%kc.ids|post%in%kc.ids) %>%
      dplyr::mutate(pre=ifelse(pre%in%kc.ids,"KC",pre),
                    post=ifelse(post%in%kc.ids,"KC",post)) %>%
      dplyr::group_by(pre,post) %>%
      dplyr::mutate(count=sum(count)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(count>=100)
    chosen.ids <- setdiff(unique(c(kc.elist.simp$pre,kc.elist.simp$post)),kc.ids)
    chosen.ids <- setdiff(chosen.ids,"KC")
    banc.meta.cutout <- banc.meta %>%
      dplyr::filter(side=="right" & (grepl(cut.out,super_class)|
                                       grepl(cut.out,cell_class)|
                                       grepl(cut.out,cell_sub_class)|
                                       grepl(cut.out,cell_type))|
                      banc_746_id %in% chosen.ids)
  }else if(cut.out.good%in%c("optic","front_leg")){
    banc.chosen.pre <- banc.syns %>%
      dplyr::filter(grepl(cut.out,neuropil),side=="right") %>%
      dplyr::group_by(pre) %>%
      dplyr::mutate(count = dplyr::n()) %>%
      dplyr::ungroup() %>%
      dplyr::filter(count>=100) %>%
      dplyr::pull(pre)
    banc.chosen.post <- banc.syns %>%
      dplyr::filter(grepl(cut.out,neuropil),side=="right") %>%
      dplyr::group_by(post) %>%
      dplyr::mutate(count = dplyr::n()) %>%
      dplyr::ungroup() %>%
      dplyr::filter(count>=100) %>%
      dplyr::pull(post)
    chosen.ids <- unique(c(banc.chosen.pre,banc.chosen.post))
    banc.meta.cutout <- banc.meta %>%
      dplyr::filter(banc_746_id %in% chosen.ids)
  }else if(cut.out.good%in%c("suboesophageal_zone","abdominal_neuromere")){
    banc.chosen.pre <- banc.syns %>%
      dplyr::filter(grepl(cut.out,neuropil)) %>%
      dplyr::group_by(pre) %>%
      dplyr::mutate(count = dplyr::n()) %>%
      dplyr::ungroup() %>%
      dplyr::filter(count>=100) %>%
      dplyr::pull(pre)
    banc.chosen.post <- banc.syns %>%
      dplyr::filter(grepl(cut.out,neuropil)) %>%
      dplyr::group_by(post) %>%
      dplyr::mutate(count = dplyr::n()) %>%
      dplyr::ungroup() %>%
      dplyr::filter(count>=100) %>%
      dplyr::pull(post)
    chosen.ids <- unique(c(banc.chosen.pre,banc.chosen.post))
    banc.meta.cutout <- banc.meta %>%
      dplyr::filter(banc_746_id %in% chosen.ids)
  }else{
    banc.meta.cutout <- banc.meta %>%
      dplyr::filter(grepl(cut.out,super_class)|
                      grepl(cut.out,cell_class)|
                      grepl(cut.out,cell_sub_class)|
                      grepl(cut.out,cell_type))
  }
  root.ids.cut.out <- na.omit(unique(banc.meta.cutout$banc_746_id))
  
  # Edgelist
  banc.elist.simp.cut.out <- banc.elist.simp %>%
    dplyr::filter(pre %in% root.ids.cut.out & post %in% root.ids.cut.out)
  
  # Synapses
  banc.syns.cut.out <- banc.syns %>%
    dplyr::filter(pre %in% root.ids.cut.out | post %in% root.ids.cut.out)
  
  # Save
  arrow::write_feather(banc.meta.cutout, file.path(save.path.cut.out,sprintf("banc_746_%s_meta.feather",cut.out.good)))
  arrow::write_feather(banc.elist.simp.cut.out, file.path(save.path.cut.out,sprintf("banc_746_%s_simple_edgelist.feather",cut.out.good)))
  arrow::write_feather(banc.syns.cut.out, file.path(save.path.cut.out,sprintf("banc_746_%s_synapses.feather",cut.out.good)))
}

###########
### OBJ ###
###########

save.path.obj <- file.path(save.path,"obj")
dir.create(save.path.obj)
save.path.obj.np <- file.path(save.path.obj,"neuropils")
dir.create(save.path.obj.np)
Rvcg::vcgObjWrite(as.mesh3d(bancr::banc.surf), filename = file.path(save.path.obj,"banc_volume_nm.obj"))
Rvcg::vcgObjWrite(as.mesh3d(bancr::banc_vnc_neuropil.surf), filename = file.path(save.path.obj,"banc_vnc_neuropil_nm.obj"))
Rvcg::vcgObjWrite(as.mesh3d(bancr::banc_brain_neuropil.surf), filename = file.path(save.path.obj,"banc_brain_neuropil_nm.obj"))
regs <- bancr::banc_vnc_neuropils.surf$RegionList
for(reg in regs){
  Rvcg::vcgObjWrite(as.mesh3d(subset(banc_vnc_neuropils.surf,reg)), filename = file.path(save.path.obj.np,sprintf("banc_neuropil_%s_nm.obj",reg)))
}
regs <- bancr::banc_brain_neuropils.surf$RegionList
for(reg in regs){
  Rvcg::vcgObjWrite(as.mesh3d(subset(banc_brain_neuropils.surf,reg)), filename = file.path(save.path.obj.np,sprintf("banc_neuropil_%s_nm.obj",reg)))
}

##############
### BUCKET ###
##############

# Send to google bucket
system("gsutil -m rsync -r /n/data1/hms/neurobio/wilson/connectomes/banc gs://sjcabs_2025_data/banc")












