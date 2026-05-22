#!/usr/bin/env Rscript
source("banc/banc-startup.R")

# File paths
malecns.path <- '/n/data1/hms/neurobio/wilson/malecns'
save.path <- '/n/data1/hms/neurobio/wilson/connectomes/malecns'
dir.create(malecns.path)
dir.create(save.path)

# Direct Google Cloud Storage URLs (decoded from the Proofpoint links)
urls <- c(
  "https://storage.googleapis.com/flyem-male-cns/v0.9/connectome-data/flat-connectome/connectome-weights-male-cns-v0.9-minconf-0.5.feather",
  "https://storage.googleapis.com/flyem-male-cns/v0.9/connectome-data/flat-connectome/body-annotations-male-cns-v0.9-minconf-0.5.feather",
  "https://storage.googleapis.com/flyem-male-cns/v0.9/connectome-data/flat-connectome/body-neurotransmitters-male-cns-v0.9.feather",
  "https://storage.googleapis.com/flyem-male-cns/v0.9/connectome-data/flat-connectome/body-stats-male-cns-v0.9-minconf-0.5.feather",
  "https://storage.googleapis.com/flyem-male-cns/v0.9/connectome-data/flat-connectome/syn-points-male-cns-v0.9-minconf-0.5.feather",
  "https://storage.googleapis.com/flyem-male-cns/v0.9/connectome-data/flat-connectome/syn-partners-male-cns-v0.9-minconf-0.5.feather",
  "https://storage.googleapis.com/flyem-male-cns/v0.9/connectome-data/flat-connectome/tbar-neurotransmitters-male-cns-v0.9.feather")

# Local filenames
local_files <- base::file.path(malecns.path, base::basename(urls))

# Download each file
# base::options(timeout = base::max(3600, base::getOption("timeout")))
# for (i in base::seq_along(urls)) {
#   if (base::file.exists(local_files[i])) {
#     base::message("Removing partial file: ", local_files[i])
#     base::file.remove(local_files[i])
#   }
#   base::message("Downloading ", urls[i])
#   utils::download.file(
#     url      = urls[i],
#     destfile = local_files[i],
#     mode     = "wb",
#     quiet    = FALSE
#   )
# }

# malecns meta
# malecns.meta.orig <- readr::read_csv(file.path(banc.meta.save.path,"malecns_09_meta.csv"), col_types = banc.col.types) %>%
#   dplyr::select(
#     malecns_09_id,
#     region,
#     side,
#     hemilineage,
#     nerve,
#     flow,
#     super_class,
#     cell_class,
#     cell_sub_class,
#     cell_type,
#     fafb_cell_type,
#     hemibrain_cell_type,
#     manc_cell_type,
#     neurotransmitter_predicted,
#     neurotransmitter_score,
#     cell_function,
#     cell_function_detailed,
#     body_part_sensory,
#     body_part_effector,
#     dimorphism,
#     optic_lobe_hex_1,
#     optic_lobe_hex_2,
#     status
#   ) %>%
#   dplyr::arrange(
#     malecns_09_id,
#     region,
#     side,
#     hemilineage,
#     nerve,
#     flow,
#     super_class,
#     cell_class,
#     cell_sub_class,
#     cell_type,
#     fafb_cell_type,
#     hemibrain_cell_type,
#     manc_cell_type,
#     neurotransmitter_predicted,
#     neurotransmitter_score,
#     cell_function,
#     cell_function_detailed,
#     body_part_sensory,
#     body_part_effector,
#     dimorphism,
#     optic_lobe_hex_1,
#     optic_lobe_hex_2,
#     status
#   ) %>%
#   dplyr::mutate(neurotransmitter_score = round(neurotransmitter_score,digits = 6),
#                 cell_type=dplyr::coalesce(cell_type,fafb_cell_type,hemibrain_cell_type,manc_cell_type)) %>%
#   dplyr::arrange(cell_type)

# Initialise seatable
# malecns.meta.orig$neurotransmitter_score[is.na(malecns.meta.orig$neurotransmitter_score)] <- 0
#bancr::banctable_append_rows(df = malecns.meta.orig, table = "malecns", base = "cns_meta", bigdata = TRUE)

# Read seatable
malecns.meta <- franken_meta(sql = "SELECT * FROM malecns", base = "cns_meta")
bodyids <- unique(malecns.meta$malecns_09_id)

# Save
arrow::write_feather(malecns.meta, file.path(save.path,"malecns_09_meta.feather"))

# # maleCNS edgelist
# malecns.edgelist.simp <- arrow::read_feather(file.path(malecns.path,"connectome-weights-male-cns-v0.9-minconf-0.5.feather"))
# malecns.edgelist.simp <- malecns.edgelist.simp %>%
#   dplyr::rename(pre=body_pre,post=body_post,count=weight) %>%
#   dplyr::group_by(post) %>%
#   dplyr::mutate(total_inputs = sum(count)) %>%
#   dplyr::ungroup() %>%
#   dplyr::mutate(norm=count/total_inputs)
# arrow::write_feather(malecns.meta, file.path(save.path,"malecns_09_simple_edgelist.feather"))

# maleCNS synapses
malecns.syns <- arrow::read_feather(file.path(malecns.path,"syn-points-male-cns-v0.9-minconf-0.5.feather"))
malecns.presyns <- malecns.syns %>%
  dplyr::filter(body %in% bodyids, kind=="PreSyn") %>%
  dplyr::select(x,
                y,
                z,
                pre_label=compartment) %>%
  dplyr::mutate(pre_label = dplyr::case_when(
    pre_label=="linker" ~ "primary_dendrite",
    pre_label=="cell-body-fiber" ~ "primary_neurite",
    TRUE ~ pre_label
  )) 
malecns.postsyns <- malecns.syns %>%
  dplyr::filter(body %in% bodyids, kind=="PostSyn") %>%
  dplyr::select(x,
                y,
                z,
                post_label=compartment) %>%
  dplyr::mutate(post_label = dplyr::case_when(
    post_label=="linker" ~ "primary_dendrite",
    post_label=="cell-body-fiber" ~ "primary_neurite",
    TRUE ~ post_label
  ))
rm('malecns.syns')
gc()

malecns.partners <- arrow::read_feather(file.path(malecns.path,"syn-partners-male-cns-v0.9-minconf-0.5.feather"))
malecns.partners <- malecns.partners %>%
  dplyr::filter(body_post %in% bodyids | body_pre %in% bodyids) %>%
  dplyr::select(x=x_pre,
                y=y_pre,
                z=z_pre,
                x_post,
                y_post,
                z_post,
                confidence = conf_pre,
                pre = body_pre,
                post = body_post,
                neuropil = primary_post) %>%
  dplyr::mutate(confidence = round(confidence,6),
                neuropil = as.character(neuropil),
                prepost = 0,
                neuropil = as.character(neuropil),
                side = dplyr::case_when(
                  grepl("\\(R\\)",neuropil) ~ "right",
                  grepl("\\(L\\)",neuropil) ~ "left",
                  TRUE ~ NA
                ),
                neuropil = gsub("\\(R\\)|\\(L\\)","",neuropil)
  )

malecns.nts <- arrow::read_feather(file.path(malecns.path,"tbar-neurotransmitters-male-cns-v0.9.feather")) %>%
  dplyr::select(
    id = point_id,
    #confidence = conf,
    x, y, z,
    region = major,
    acetylcholine = nt_acetylcholine_prob,
    gaba = nt_gaba_prob,
    glutamate = nt_glutamate_prob,
    dopamine = nt_dopamine_prob,
    serotonin = nt_serotonin_prob,
    octopamine = nt_octopamine_prob,
    histamine = nt_histamine_prob
  ) %>%
  # round numeric columns in one call
  dplyr::mutate(across(c(acetylcholine:histamine), ~ round(.x, 6)))
gc()

# build numeric matrix of neurotransmitter probs
nt_mat <- as.matrix(dplyr::select(malecns.nts, acetylcholine:histamine))

# compute per-row top probability, preserving NA when all cols are NA
syn_top_p <- apply(nt_mat, 1, function(r) if (all(is.na(r))) NA_real_ else max(r, na.rm = TRUE))

# for the name of the top neurotransmitter:
# replace NA with -Inf so max.col works (but we'll mask rows that were all-NA afterwards)
nt_mat_for_max <- nt_mat
nt_mat_for_max[is.na(nt_mat_for_max)] <- -Inf
idx <- max.col(nt_mat_for_max, ties.method = "first")    # index of first max per row
syn_top_nt <- colnames(nt_mat)[idx]
# mask names where there was no non-NA value
syn_top_nt[is.na(syn_top_p)] <- NA_character_

# attach back to the tibble
malecns.nts <- malecns.nts %>%
  dplyr::mutate(syn_top_p  = syn_top_p,
                syn_top_nt = syn_top_nt)
gc()

# malecns.synapses <- malecns.partners %>%
#   dplyr::left_join(malecns.nts,by=c("x","y","z")) %>%
#   dplyr::left_join(malecns.presyns,by=c("x","y","z")) %>%
#   dplyr::left_join(malecns.postsyns,by=c("x","y","z"))
chunk_size <- 2000000L
base_df <- malecns.partners
n <- nrow(base_df)
n_chunks <- ceiling(n / chunk_size)
indexes <- split(seq_len(n), ceiling(seq_len(n) / chunk_size))
process_chunk <- function(idxs) {
  # slice rows for this chunk
  chunk <- base_df[idxs, , drop = FALSE]
  
  # NT predictions and pre compartment join on presynaptic (x, y, z)
  # Post compartment joins on postsynaptic (x_post, y_post, z_post)
  res <- chunk %>%
    left_join(malecns.nts,      by = c("x", "y", "z")) %>%
    left_join(malecns.presyns,  by = c("x", "y", "z")) %>%
    left_join(malecns.postsyns, by = c("x_post" = "x", "y_post" = "y", "z_post" = "z"))
  
  rm(chunk)
  gc()
  res
}
malecns.synapses <- map_dfr(indexes, process_chunk, .progress = TRUE)
gc()

# Save
malecns.synapses <- malecns.synapses %>%
  dplyr::select(-x_post, -y_post, -z_post) %>%
  dplyr::mutate(pre=as.character(pre),
                post=as.character(post))
arrow::write_parquet(x = malecns.synapses, 
                     sink = file.path(save.path,"malecns_09_synapses.parquet"),
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
malecns.elist.simp <- malecns.synapses %>%
  dplyr::group_by(post) %>%
  dplyr::mutate(total_input = dplyr::n()) %>%
  dplyr::group_by(pre, post) %>%
  dplyr::mutate(count = dplyr::n(),
                norm = count/total_input) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(pre,post,count,norm,total_input)
pre.ids.missing <- unique(setdiff(bodyids,malecns.elist.simp$pre))
post.ids.missing <- unique(setdiff(bodyids,malecns.elist.simp$post))
arrow::write_feather(malecns.elist.simp, file.path(save.path,"malecns_09_simple_edgelist.feather"))

## skeletons: including compartmentalised skeletons (optional: axons vs. dendrites)
malecns.elist <- malecns.synapses %>%
  dplyr::group_by(post, post_label) %>%
  dplyr::mutate(compartment_input = dplyr::n()) %>%
  dplyr::group_by(pre, pre_label, post, post_label) %>%
  dplyr::mutate(count = dplyr::n(),
                norm = count/compartment_input,
                post_label = ifelse(is.na(post_label),"unknown",post_label),
                pre_label = ifelse(is.na(pre_label),"unknown",pre_label)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(pre, pre_label, post, post_label, count, norm, compartment_input)
pre.ids.missing <- unique(setdiff(bodyids,malecns.elist$pre))
post.ids.missing <- unique(setdiff(bodyids,malecns.elist$post))
arrow::write_feather(malecns.elist, file.path(save.path,"malecns_09_split_edgelist.feather"))

#################
### SKELETONS ###
#################

## skeletons: including compartmentalised skeletons (optional: axons vs. dendrites)

### In malecns space
malecns.in.malecns.space <- "/n/data1/hms/neurobio/wilson/banc/matching/malecns/JRCFIB2022M/skeletons-swc"
# wd <- getwd()
# skels <- malevnc::malecns_read_neurons(bodyids,unit="raw")
# write.neurons(skels,malecns.in.malecns.space)
# out_zip <- file.path(save.path, "malecns_malecns_space_swc.zip")
# oldwd <- setwd(malecns.in.malecns.space); on.exit(setwd(oldwd), add = TRUE)
# files <- list.files(".", recursive = TRUE, all.files = TRUE, no.. = TRUE)
# zip::zipr(zipfile = out_zip, files = files)
# cat("Wrote:", out_zip, "\n")
# setwd(wd)
new.folder <- file.path(save.path, "malecns_malecns_space_swc")
dir.create(new.folder)
files <- list.files(malecns.in.malecns.space, full.names = TRUE)
file.copy(from = files, to = new.folder, overwrite = FALSE)

### In BANC space
#malecns.in.banc.space <- "/n/data1/hms/neurobio/wilson/banc/matching/malecns/banc_space_split/swc"
version <- banc.nblast.malecns.version
malecns.in.banc.space <- file.path(banc.nblast.malecns.swc.save.path, version)
# out_zip <- file.path(save.path, "malecns_banc_space_swc.zip")
# oldwd <- setwd(malecns.in.banc.space); on.exit(setwd(oldwd), add = TRUE)
# files <- list.files(".", recursive = TRUE, all.files = TRUE, no.. = TRUE)
# if (!requireNamespace("zip", quietly = TRUE)) install.packages("zip")
# zip::zipr(zipfile = out_zip, files = files)
# cat("Wrote:", out_zip, "\n")
# setwd(wd)
new.folder <- file.path(save.path, "malecns_banc_space_swc")
dir.create(new.folder)
files <- list.files(malecns.in.banc.space, full.names = TRUE)
file.copy(from = files, to = new.folder, overwrite = TRUE)

################
### CUT OUTS ###
################

# Region cut-outs disabled 2026-05-02: we no longer host per-region
# subsets on GCS. The whole-CNS meta + edgelist + synapse files
# above are sufficient; downstream consumers can filter as needed.
# Original logic preserved in git history if it's ever needed back.

# cut.outs <- c("mushroom_body","antennal_lobe","central_complex","optic","suboesophageal_zone","front_leg","abdominal_neuromere","optic_lobe_hex_08","optic_lobe_hex_04")
# for(cut.out in cut.outs){
#
#   # Save path
#   cut.out.good <- snakecase::to_snake_case(cut.out)
#   save.path.cut.out <- file.path(save.path,cut.out.good)
#   dir.create(save.path.cut.out)
#   if(cut.out=="optic"){
#     cut.out = "optic_lobe_intrinsic|optic_lobe_sensory|visual_project_visual_centrifugal|optic_lobe_bilateral"
#     cut.out = "^LO|^LOP|^AME|^ME"
#   }
#   if(cut.out=="antennal_lobe"){
#     cut.out = "antennal_lobe|olfactory_receptor|thermosensory_receptor|hygrosensory_receptor|CSD"
#   }
#   if(cut.out=="suboesophageal_zone"){
#     cut.out = "^FLA|^SEZ|^GNG|^SAD|^AMMC|^PRW"
#   }
#   if(cut.out=="front_leg"){
#     cut.out = "^LegNp\\(T1\\)|T1|^ProNM-T1|^LNp_T1"
#   }
#   if(cut.out=="abdominal_neuromere"){
#     cut.out = "^ANm|^ABDNM|^ADNM"
#   }
#
#   # Custom cuts
#   if(cut.out=="mushroom_body"){
#     cut.out = "mushroom_body|kenyon_cell|APL|DPM|LHMB1|OA-VPM3"
#     kc.ids <- malecns.meta %>%
#       dplyr::filter(cell_class=="kenyon_cell",side=="right") %>%
#       dplyr::pull(malecns_09_id)
#     kc.elist.simp <- malecns.elist.simp %>%
#       dplyr::filter(pre%in%kc.ids|post%in%kc.ids) %>%
#       dplyr::mutate(pre=ifelse(pre%in%kc.ids,"KC",pre),
#                     post=ifelse(post%in%kc.ids,"KC",post)) %>%
#       dplyr::group_by(pre,post) %>%
#       dplyr::mutate(count=sum(count)) %>%
#       dplyr::ungroup() %>%
#       dplyr::filter(count>=100)
#     chosen.ids <- setdiff(unique(c(kc.elist.simp$pre,kc.elist.simp$post)),kc.ids)
#     chosen.ids <- setdiff(chosen.ids,"KC")
#     malecns.meta.cutout <- malecns.meta %>%
#       dplyr::filter(side=="right" & (grepl(cut.out,super_class)|
#                                        grepl(cut.out,cell_class)|
#                                        grepl(cut.out,cell_sub_class)|
#                                        grepl(cut.out,cell_type))|
#                       malecns_09_id %in% chosen.ids)
#   }else if(cut.out.good%in%c("optic","front_leg")){
#     malecns.chosen.pre <- malecns.synapses %>%
#       dplyr::filter(grepl(cut.out,neuropil),side=="right") %>%
#       dplyr::group_by(pre) %>%
#       dplyr::mutate(count = dplyr::n()) %>%
#       dplyr::ungroup() %>%
#       dplyr::filter(count>=100) %>%
#       dplyr::pull(pre)
#     malecns.chosen.post <- malecns.synapses %>%
#       dplyr::filter(grepl(cut.out,neuropil),side=="right") %>%
#       dplyr::group_by(post) %>%
#       dplyr::mutate(count = dplyr::n()) %>%
#       dplyr::ungroup() %>%
#       dplyr::filter(count>=100) %>%
#       dplyr::pull(post)
#     chosen.ids <- unique(c(malecns.chosen.pre,malecns.chosen.post))
#     malecns.meta.cutout <- malecns.meta %>%
#       dplyr::filter(malecns_09_id %in% chosen.ids)
#   }else if(cut.out.good%in%c("suboesophageal_zone","abdominal_neuromere")){
#     malecns.chosen.pre <- malecns.synapses %>%
#       dplyr::filter(grepl(cut.out,neuropil)) %>%
#       dplyr::group_by(pre) %>%
#       dplyr::mutate(count = dplyr::n()) %>%
#       dplyr::ungroup() %>%
#       dplyr::filter(count>=100) %>%
#       dplyr::pull(pre)
#     malecns.chosen.post <- malecns.synapses %>%
#       dplyr::filter(grepl(cut.out,neuropil)) %>%
#       dplyr::group_by(post) %>%
#       dplyr::mutate(count = dplyr::n()) %>%
#       dplyr::ungroup() %>%
#       dplyr::filter(count>=100) %>%
#       dplyr::pull(post)
#     chosen.ids <- unique(c(malecns.chosen.pre,malecns.chosen.post))
#     malecns.meta.cutout <- malecns.meta %>%
#       dplyr::filter(malecns_09_id %in% chosen.ids)
#   }else if(cut.out.good=="optic_lobe_hex_08"){
#     malecns.meta.cutout1 <- malecns.meta %>%
#       dplyr::filter(optic_lobe_hex_1=="8"&optic_lobe_hex_2=="8")
#     downs1 <- malecns.elist.simp %>%
#       dplyr::filter(post %in% !!malecns.meta.cutout1$malecns_09_id,
#                     count>5)
#     malecns.meta.cutout <- malecns.meta %>%
#       dplyr::filter(malecns_09_id%in%!!downs1$pre|malecns_09_id%in%!!downs1$post)
#   }else if(cut.out.good=="optic_lobe_hex_04"){
#     malecns.meta.cutout1 <- malecns.meta %>%
#       dplyr::filter(optic_lobe_hex_1=="4"&optic_lobe_hex_2=="4")
#     downs1 <- malecns.elist.simp %>%
#       dplyr::filter(post %in% !!malecns.meta.cutout1$malecns_09_id,
#                     count>5)
#     malecns.meta.cutout <- malecns.meta %>%
#       dplyr::filter(malecns_09_id%in%!!downs1$pre|malecns_09_id%in%!!downs1$post)
#   }else{
#     malecns.meta.cutout <- malecns.meta %>%
#       dplyr::filter(grepl(cut.out,super_class)|
#                       grepl(cut.out,cell_class)|
#                       grepl(cut.out,cell_sub_class)|
#                       grepl(cut.out,cell_type))
#   }
#   bodyids.cut.out <- na.omit(unique(malecns.meta.cutout$malecns_09_id))
#
#   # Edgelist
#   malecns.elist.simp.cut.out <- malecns.elist.simp %>%
#     dplyr::filter(pre %in% bodyids.cut.out & post %in% bodyids.cut.out)
#
#   # Synapses
#   malecns.syns.cut.out <- malecns.synapses %>%
#     dplyr::filter(pre %in% bodyids.cut.out | post %in% bodyids.cut.out)
#
#   # Save
#   arrow::write_feather(malecns.meta.cutout, file.path(save.path.cut.out,sprintf("malecns_09_%s_meta.feather",cut.out.good)))
#   arrow::write_feather(malecns.elist.simp.cut.out, file.path(save.path.cut.out,sprintf("malecns_09_%s_simple_edgelist.feather",cut.out.good)))
#   arrow::write_feather(malecns.syns.cut.out, file.path(save.path.cut.out,sprintf("malecns_09_%s_synapses.feather",cut.out.good)))
# }

###########
### OBJ ###
###########
library(malecns)
save.path.obj <- file.path(save.path,"obj")
dir.create(save.path.obj)
save.path.obj.np <- file.path(save.path.obj,"neuropils")
dir.create(save.path.obj.np)
Rvcg::vcgObjWrite(as.mesh3d(malecns_shell.surf/c(8,8,8)), filename = file.path(save.path.obj,"malecns_brain_volume_raw.obj"))
Rvcg::vcgObjWrite(as.mesh3d(malecns.surf/c(8,8,8)), filename = file.path(save.path.obj,"malecns_brain_neuropil_raw.obj"))
Rvcg::vcgObjWrite(as.mesh3d(malecnsvnc_shell.surf/c(8,8,8)), filename = file.path(save.path.obj,"malecns_vnc_volume_raw.obj"))
Rvcg::vcgObjWrite(as.mesh3d(malecnsvnc.surf/c(8,8,8)), filename = file.path(save.path.obj,"malecns_vnc_neuropil_raw.obj"))
rois <- neuprintr::neuprint_ROIs(dataset = "male-cns:v0.9")
rois <- rois[!grepl("_col_|_layer_",rois)]
for(roi in rois){
  try({
    mesh=neuprint_ROI_mesh(roi,dataset="male-cns:v0.9")
    Rvcg::vcgObjWrite(as.mesh3d(mesh), filename = file.path(save.path.obj.np,sprintf("malecns_neuropil_%s_raw.obj",roi)))
  })
}

##############
### BUCKET ###
##############

# Send to google bucket. Target nests under malecns_09/ to match the
# layout used by banc_888/, fafb_783/, manc_121/, hemibrain_121/.
# Excludes historical region-cutout subdirs which still exist on the
# local disk from earlier runs but are no longer hosted on GCS.
cutout.skip <- "^(mushroom_body|antennal_lobe|central_complex|optic|suboesophageal_zone|front_leg|abdominal_neuromere|optic_lobe_hex_04|optic_lobe_hex_08)/"
system(sprintf(
  "gsutil -m rsync -r -x '%s' /n/data1/hms/neurobio/wilson/connectomes/malecns gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/malecns_09",
  cutout.skip))
























