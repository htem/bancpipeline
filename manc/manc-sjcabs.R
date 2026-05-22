#!/usr/bin/env Rscript
source("banc/banc-startup.R")

# File paths
split.master.folder <- banc.nblast.manc.split.save.path
split.folder <- file.path(split.master.folder,"swc")
synapses.folder <- file.path(split.master.folder,"synapses")
metrics.folder <- file.path(split.master.folder,"metrics")
save.path <- '/n/data1/hms/neurobio/wilson/connectomes/manc'
dir.create(save.path)
manc.data.prefix <- "manc_1.2.1"

# Essential data:

## meta data: full meta data for useful labels, both anatomy (e.g. cell type, cell class) and function (i.e. body parts innervated, sparse known cell functions)
manc.meta <- franken_meta() %>%
  dplyr::filter(grepl("MANC",dataset), !is.na(manc_id)) %>%
  dplyr::distinct(manc_121_id = manc_id, 
                  .keep_all = TRUE) %>%
  dplyr::select(manc_121_id,
                region,
                side,
                hemilineage,
                nerve,
                flow,
                super_class,
                cell_class,
                cell_sub_class,
                cell_type,
                neurotransmitter_predicted = top_nt,
                cell_function,
                cell_function_detailed,
                body_part_sensory,
                body_part_effector)
root.ids <- manc.meta$manc_121_id
arrow::write_feather(manc.meta, file.path(save.path,"manc_121_meta.feather"))

## synapse tables: from CAVE, with a neuropil column + axon-dendrite split
manc.synapses <- arrow::read_parquet(
  file.path(banc.connectivity.save.path, paste0(manc.data.prefix, "_synapses.parquet"))) %>%
  dplyr::select(connector_id,
                x,y,z,
                confidence,
                prepost,
                pre = pre_id,
                post = post_id,
                neuropil = inside,
                #scores, cleft_scores,
                syn_top_p = conf_nt_p,
                syn_top_nt = conf_nt,
                # gaba, acetylcholine, glutamate, octopamine, serotonin, dopamine,
                # strahler_order,
                # geodesic_distance,
                # geodesic_distance_norm,
                pre_label=label,
                status) %>%
  dplyr::filter(pre %in% root.ids | post %in% root.ids) %>%
  dplyr::mutate(pre_label = snakecase::to_snake_case(pre_label))
manc.postsynapses <- arrow::read_parquet(
  file.path(banc.connectivity.save.path, paste0(manc.data.prefix, "_synapses.parquet"))) %>%
  dplyr::filter(prepost==1) %>%
  dplyr::filter(post_id %in% root.ids) %>%
  dplyr::select(connector_id,
                post=post_id,
                post_label = label) %>%
  dplyr::mutate(post_label = snakecase::to_snake_case(post_label))
manc.synapses <- manc.synapses %>%
  dplyr::left_join(manc.postsynapses,by=c("connector_id","post")) %>%
  dplyr::mutate(pre_label=ifelse(is.na(pre_label),"unknown",pre_label),
                post_label=ifelse(is.na(post_label),"unknown",post_label))

# Assign neuropils

# Put MANC synapses back into MANC space
xyz <- nat::xyzmatrix(manc.synapses)
xyz <- bancr::banc_to_JRC2018F(xyz, 
                               region="vnc", 
                               method="tpsreg", 
                               banc.units = "nm", 
                               inverse = FALSE)
xyz <- nat.templatebrains::xform_brain(xyz, 
                                       reference = "MANC", 
                                       sample = "JRCVNC2018F")
nat::xyzmatrix(manc.synapses) <- xyz

# Save
pre.ids.missing <- unique(setdiff(root.ids,manc.synapses$pre))
post.ids.missing <- unique(setdiff(root.ids,manc.synapses$post))
arrow::write_parquet(manc.synapses, file.path(save.path,"manc_121_synapses.parquet"),
                     version = "2.6", 
                     compression = "snappy",
                     compression_level = NULL,
                     chunk_size = 100000,
                     use_dictionary = TRUE,
                     allow_truncated_timestamps = FALSE)

## simple edgelists: proofread neurons and connections between them
manc.elist.simp <- arrow::read_feather(
  file.path(banc.connectivity.save.path, paste0(manc.data.prefix, "_edgelist.feather"))) %>%
  dplyr::distinct(pre,post,count,norm,total_input=post_count)
pre.ids.missing <- unique(setdiff(root.ids,manc.elist.simp$pre))
post.ids.missing <- unique(setdiff(root.ids,manc.elist.simp$post))
arrow::write_feather(manc.elist.simp, file.path(save.path,"manc_121_simple_edgelist.feather"))

## split edgelists: proofread neurons and connections between them
manc.elist <- arrow::read_feather(
  file.path(banc.connectivity.save.path, paste0(manc.data.prefix, "_edgelist.feather")))
pre.ids.missing <- unique(setdiff(root.ids,manc.elist$pre))
post.ids.missing <- unique(setdiff(root.ids,manc.elist$post))
arrow::write_feather(manc.elist, file.path(save.path,"manc_121_split_edgelist.feather"))

## skeletons: including compartmentalised skeletons (optional: axons vs. dendrites)

### In manc space
manc.in.manc.space <- "/n/data1/hms/neurobio/wilson/banc/matching/manc/manc_manc_microns_space/"
# wd <- getwd()
# skels <- malevnc::manc_read_neurons(root.ids,unit="microns")
# write.neurons(skels,manc.in.manc.space)
# out_zip <- file.path(save.path, "manc_manc_space_swc.zip")
# oldwd <- setwd(manc.in.manc.space); on.exit(setwd(oldwd), add = TRUE)
# files <- list.files(".", recursive = TRUE, all.files = TRUE, no.. = TRUE)
# zip::zipr(zipfile = out_zip, files = files)
# cat("Wrote:", out_zip, "\n")
# setwd(wd)
new.folder <- file.path(save.path, "manc_manc_space_swc")
dir.create(new.folder)
files <- list.files(manc.in.manc.space, full.names = TRUE)
file.copy(from = files, to = new.folder, overwrite = FALSE)

### In BANC space
manc.in.banc.space <- "/n/data1/hms/neurobio/wilson/banc/matching/manc/banc_space_split/swc"
# out_zip <- file.path(save.path, "manc_banc_space_split_swc.zip")
# oldwd <- setwd(manc.in.banc.space); on.exit(setwd(oldwd), add = TRUE)
# files <- list.files(".", recursive = TRUE, all.files = TRUE, no.. = TRUE)
# if (!requireNamespace("zip", quietly = TRUE)) install.packages("zip")
# zip::zipr(zipfile = out_zip, files = files)
# cat("Wrote:", out_zip, "\n")
# new.folder <- file.path(save.path, "manc_banc_space_split_swc")
# dir.create(new.folder)
files <- list.files(manc.in.banc.space, full.names = TRUE)
file.copy(from = files, to = new.folder, overwrite = FALSE)

###########
### OBJ ###
###########

# MANC.surf is in microns
save.path.obj <- file.path(save.path,"obj")
dir.create(save.path.obj)
save.path.obj.np <- file.path(save.path.obj,"neuropils")
dir.create(save.path.obj.np)
Rvcg::vcgObjWrite(as.mesh3d(MANC.tissue.surf), filename = file.path(save.path.obj,"manc_volume_microns.obj"))
Rvcg::vcgObjWrite(as.mesh3d(MANC.surf), filename = file.path(save.path.obj,"manc_neuropil_microns.obj"))
regs <- JRCFIBVNC2020MNP.surf$RegionList
for(reg in regs){
  Rvcg::vcgObjWrite(as.mesh3d(subset(JRCFIBVNC2020MNP.surf,reg)), filename = file.path(save.path.obj.np,sprintf("manc_neuropil_%s_microns.obj",reg)))
}


################
### CUT OUTS ###
################

# Go by cut-out
cut.outs <- c("front_leg","abdominal_neuromere")
for(cut.out in cut.outs){
  
  # Save path
  cut.out.good <- snakecase::to_snake_case(cut.out)
  save.path.cut.out <- file.path(save.path,cut.out.good)
  dir.create(save.path.cut.out)
  if(cut.out=="front_leg"){
    cut.out = "^LegNp\\(T1\\)|T1|^ProNM-T1|^LNp_T1"
  }
  if(cut.out=="abdominal_neuromere"){
    cut.out = "^ANm|^ABDNM"
  }
  
  # Custom cuts
  if(cut.out.good%in%c("front_leg","abdominal_neuromere")){
    manc.chosen.pre <- manc.synapses %>%
      dplyr::filter(grepl(cut.out,neuropil),side=="right") %>%
      dplyr::group_by(pre) %>%
      dplyr::mutate(count = dplyr::n()) %>%
      dplyr::ungroup() %>%
      dplyr::filter(count>=100) %>%
      dplyr::pull(pre)
    manc.chosen.post <- manc.synapses %>%
      dplyr::filter(grepl(cut.out,neuropil),side=="right") %>%
      dplyr::group_by(post) %>%
      dplyr::mutate(count = dplyr::n()) %>%
      dplyr::ungroup() %>%
      dplyr::filter(count>=100) %>%
      dplyr::pull(post)
    chosen.ids <- unique(c(manc.chosen.pre,manc.chosen.post))
    manc.meta.cutout <- manc.meta %>%
      dplyr::filter(manc_121_id %in% chosen.ids)
  }else{
    manc.meta.cutout <- manc.meta %>%
      dplyr::filter(grepl(cut.out,super_class)|
                      grepl(cut.out,cell_class)|
                      grepl(cut.out,cell_sub_class)|
                      grepl(cut.out,cell_type))
  }
  bodyids.cut.out <- na.omit(unique(manc.meta.cutout$manc_121_id))
  
  # Edgelist
  manc.elist.simp.cut.out <- manc.elist.simp %>%
    dplyr::filter(pre %in% bodyids.cut.out & post %in% bodyids.cut.out)
  
  # Synapses
  manc.syns.cut.out <- manc.synapses %>%
    dplyr::filter(pre %in% bodyids.cut.out | post %in% bodyids.cut.out)
  
  # Save
  arrow::write_feather(manc.meta.cutout, file.path(save.path.cut.out,sprintf("manc_121_%s_meta.feather",cut.out.good)))
  arrow::write_feather(manc.elist.simp.cut.out, file.path(save.path.cut.out,sprintf("manc_121_%s_simple_edgelist.feather",cut.out.good)))
  arrow::write_feather(manc.syns.cut.out, file.path(save.path.cut.out,sprintf("manc_121_%s_synapses.feather",cut.out.good)))
}

##############
### BUCKET ###
##############

# Send to google bucket
system("gsutil -m rsync -r /n/data1/hms/neurobio/wilson/connectomes/manc gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data")




