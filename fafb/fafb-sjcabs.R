#!/usr/bin/env Rscript

#' fafb-sjcabs — Assemble + publish the FAFB v783 SJCABS-style data bundle.
#'
#' Pulls FAFB v783 metadata, edgelist, NBLAST scores and synapse / skeleton
#' artefacts; aggregates into a per-version SQLite + companion feathers
#' under save.path; then rsyncs to the lab GCS bucket
#' (`compiled_data/fafb_783/`). Mirrors `manc/manc-sjcabs.R` and
#' `malecns/malecns-sjcabs.R`.
#'
#' @section Reads:
#'   - `franken_meta()` (FAFB rows), `hemibrainr` snapshots, FAFB SQLite
#'   - SeaTable `banc_meta` for FAFB cross-match columns
#'
#' @section Writes:
#'   - `/n/data1/hms/neurobio/wilson/connectomes/fafb/fafb_783_*.{feather,sqlite}`
#'   - GCS: `gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/fafb_783/`
#'
#' @section Notes:
#'   - Wrangling logic was originally developed in `fafbpipeline`
#'     (`flyconnectome/fafbpipeline:flywire/flywire-sjcabs.R`). This is the
#'     GCS-publish-side mirror, copied 2026-05-22 so bancpipeline owns the
#'     publish step. `banc/banc-startup.R` covers the library load that the
#'     original `flywire-startup.R` provided.
#'   - Run from bancpipeline repo root: `Rscript fafb/fafb-sjcabs.R`.
source("banc/banc-startup.R")
library(bancr)

# Essential data:

# master save path:
flywire.save.path <- "/n/data1/hms/neurobio/wilson/fafbz//hemibrainr_data/flywire_neurons/783"
sql.file <- "flywire_783_data.sqlite"
save.path <- '/n/data1/hms/neurobio/wilson/connectomes/fafb'
dir.create(save.path)

## meta data: full meta data for useful labels, both anatomy (e.g. cell type, cell class) and function (i.e. body parts innervated, sparse known cell functions)
fafb.meta <- franken_meta() %>%
  dplyr::filter(grepl("FAFB",dataset), !is.na(fafb_id)) %>%
  dplyr::distinct(fafb_783_id = fafb_id, 
                  .keep_all = TRUE) %>%
  dplyr::select(fafb_783_id,
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
                neurotransmitter_score = FAFB_top_nt_conf,
                cell_function,
                cell_function_detailed,
                body_part_sensory,
                body_part_effector,
                status = FAFB_status,
                sexually_dimorphic = dimorphism,
                soma_dcv_density,
                cell_dcv_density_um,
                cell_dcv_density_um3)
root.ids <- fafb.meta$fafb_783_id
arrow::write_feather(fafb.meta, file.path(save.path,"fafb_783_meta.feather"))

## synapse tables: from CAVE, with a neuropil column + axon-dendrite split
con <- DBI::dbConnect(RSQLite::SQLite(), file.path(flywire.save.path,sql.file))
fw.presynapses <- dplyr::tbl(con, "presynapses") %>%
  dplyr::select(id=offset,
                connector_id,
                x,y,z,
                confidence = cleft_scores, 
                syn_top_nt,
                syn_top_p, 
                gaba, acetylcholine, glutamate, octopamine, serotonin, dopamine,
                prepost,
                pre = pre_id,
                post = post_id,
                neuropil = inside,
                # strahler_order,
                # geodesic_distance,
                # geodesic_distance_norm,
                # status,
                pre_label = label) %>%
  dplyr::filter(pre %in% root.ids | post %in% root.ids,
                pre != post,
                prepost == 0) %>%
  dplyr::collect()
fw.postsynapses <- dplyr::tbl(con, "postsynapses") %>%
  dplyr::filter(post_id %in% root.ids,
                prepost==1) %>%
  dplyr::select(id=offset,
                post_label = label) %>%
  dplyr::collect()
DBI::dbDisconnect(con)
fafb.synapses <- dplyr::left_join(fw.presynapses,fw.postsynapses,by="id") %>%
  dplyr::mutate(pre_label = snakecase::to_snake_case(pre_label),
                post_label = snakecase::to_snake_case(post_label),
                side = dplyr::case_when(
                  grepl("_R",neuropil) ~ "right",
                  grepl("_L",neuropil) ~ "left",
                  TRUE ~ NA
                ),
                neuropil = gsub("_R|_L","",neuropil),
                pre_label = ifelse(is.na(pre_label),"unknown",pre_label),
                post_label = ifelse(is.na(post_label),"unknown",post_label)) %>%
  dplyr::distinct(id, .keep_all = TRUE)
pre.ids.missing <- unique(setdiff(root.ids,fafb.synapses$pre))
post.ids.missing <- unique(setdiff(root.ids,fafb.synapses$post))
arrow::write_parquet(fafb.synapses, file.path(save.path,"fafb_783_synapses.parquet"),
                     version = "2.6", 
                     compression = "snappy",
                     compression_level = NULL,
                     chunk_size = 100000,
                     use_dictionary = TRUE,
                     allow_truncated_timestamps = FALSE)

## simple edgelists: proofread neurons and connections between them
con <- DBI::dbConnect(RSQLite::SQLite(), file.path(flywire.save.path,sql.file))
fafb.elist.simp <- dplyr::tbl(con, "edgelist_simple") %>%
  dplyr::collect() %>%
  dplyr::distinct(pre,post,count,norm,total_input=post_count)
DBI::dbDisconnect(con)
pre.ids.missing <- unique(setdiff(root.ids,fafb.elist.simp$pre))
post.ids.missing <- unique(setdiff(root.ids,fafb.elist.simp$post))
arrow::write_feather(fafb.elist.simp, file.path(save.path,"fafb_783_simple_edgelist.feather"))

## skeletons: including compartmentalised skeletons (optional: axons vs. dendrites)
fafb.elist <- fafb.synapses %>%
  dplyr::group_by(post, post_label) %>%
  dplyr::mutate(compartment_input = dplyr::n()) %>%
  dplyr::group_by(pre, pre_label, post, post_label) %>%
  dplyr::mutate(count = dplyr::n(),
                norm = count/compartment_input) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(pre, pre_label, post, post_label, count, norm, compartment_input)
pre.ids.missing <- unique(setdiff(root.ids,fafb.elist$pre))
post.ids.missing <- unique(setdiff(root.ids,fafb.elist$post))
arrow::write_feather(fafb.elist, file.path(save.path,"fafb_783_split_edgelist.feather"))

# ## split edgelists: proofread neurons and connections between them
# con <- DBI::dbConnect(RSQLite::SQLite(), file.path(flywire.save.path,sql.file))
# fafb.elist <- dplyr::tbl(con, "edgelist") %>%
#   dplyr::collect()
# DBI::dbDisconnect(con)
# pre.ids.missing <- unique(setdiff(root.ids,fafb.elist$pre))
# post.ids.missing <- unique(setdiff(root.ids,fafb.elist$post))
# arrow::write_feather(fafb.elist, file.path(save.path,"fafb_783_split_edgelist.feather"))

############
### DCVS ###
############

## DCV cell tables
con <- DBI::dbConnect(RSQLite::SQLite(), file.path(flywire.save.path,sql.file))
fw.cell.dcvs <- dplyr::tbl(con, "cell_dcvs") %>%
  dplyr::filter(root_783 %in% root.ids) %>%
  dplyr::select(id,
                x,y,z,
                size,
                supervoxel_id,
                root_783,
                neuropil) %>%
  dplyr::collect()
DBI::dbDisconnect(con)
arrow::write_feather(fw.cell.dcvs, file.path(save.path,"fafb_783_cell_dcv_detection.feather"))
rm('fw.cell.dcvs')

## DCV soma tables
con <- DBI::dbConnect(RSQLite::SQLite(), file.path(flywire.save.path,sql.file))
fw.soma.dcvs <- dplyr::tbl(con, "soma_dcvs") %>%
  dplyr::filter(root_783 %in% root.ids) %>%
  dplyr::collect()
DBI::dbDisconnect(con)
arrow::write_feather(fw.soma.dcvs, file.path(save.path,"fafb_783_soma_dcv_detection.feather"))
rm('fw.soma.dcvs')

#################
### SKELETONS ###
#################

### In FAFB space
fafb.in.fafb.space <- "/n/data1/hms/neurobio/wilson/fafbz//hemibrainr_data/flywire_neurons/783/split"
# wd <- getwd()
# out_zip <- file.path(save.path, "fafb_fafb_space_swc.zip")
# oldwd <- setwd(fafb.in.fafb.space); on.exit(setwd(oldwd), add = TRUE)
# files <- list.files(".", recursive = TRUE, all.files = TRUE, no.. = TRUE)
# zip::zipr(zipfile = out_zip, files = files)
# cat("Wrote:", out_zip, "\n")
# setwd(wd)
new.folder <- file.path(save.path, "fafb_fafb_space_swc")
dir.create(new.folder)
files <- list.files(fafb.in.fafb.space, full.names = TRUE)
file.copy(from = files, to = new.folder, overwrite = FALSE)

### In BANC space
fafb.in.banc.space <- "/n/data1/hms/neurobio/wilson/banc/matching/fafb/fafb_banc_space_swc/"
# out_zip <- file.path(save.path, "fafb_banc_space_swc.zip")
# oldwd <- setwd(fafb.in.banc.space); on.exit(setwd(oldwd), add = TRUE)
# files <- list.files(".", recursive = TRUE, all.files = TRUE, no.. = TRUE)
# zip::zipr(zipfile = out_zip, files = files)
# cat("Wrote:", out_zip, "\n")
# setwd(wd)
new.folder <- file.path(save.path, "fafb_banc_space_swc")
dir.create(new.folder)
files <- list.files(fafb.in.banc.space, full.names = TRUE)
file.copy(from = files, to = new.folder, overwrite = FALSE)

# src=/n/data1/hms/neurobio/wilson/fafbz//hemibrainr_data/flywire_neurons/783/split/
# dst=/n/data1/hms/neurobio/wilson/connectomes/fafb/fafb_fafb_space_swc/
# rsync -av "$src" "$dst"

# src=/n/data1/hms/neurobio/wilson/banc/matching/fafb/fafb_banc_space_swc/
# dst=/n/data1/hms/neurobio/wilson/connectomes/fafb/fafb_banc_space_swc/
# rsync -av "$src" "$dst"


################
### CUT OUTS ###
################

# Go by cut-out
cut.outs <- c("mushroom_body","antennal_lobe","central_complex","optic","suboesophageal_zone")
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
  
  # Custom cuts
  if(cut.out=="mushroom_body"){
    cut.out = "mushroom_body|kenyon_cell|APL|DPM|LHMB1|OA-VPM3"
    kc.ids <- fafb.meta %>%
      dplyr::filter(cell_class=="kenyon_cell",side=="right") %>%
      dplyr::pull(fafb_783_id)
    kc.elist.simp <- fafb.elist.simp %>%
      dplyr::filter(pre%in%kc.ids|post%in%kc.ids) %>%
      dplyr::mutate(pre=ifelse(pre%in%kc.ids,"KC",pre),
                    post=ifelse(post%in%kc.ids,"KC",post)) %>%
      dplyr::group_by(pre,post) %>%
      dplyr::mutate(count=sum(count)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(count>=100)
    chosen.ids <- setdiff(unique(c(kc.elist.simp$pre,kc.elist.simp$post)),kc.ids)
    chosen.ids <- setdiff(chosen.ids,"KC")
    fafb.meta.cutout <- fafb.meta %>%
      dplyr::filter(side=="right" & (grepl(cut.out,super_class)|
                                       grepl(cut.out,cell_class)|
                                       grepl(cut.out,cell_sub_class)|
                                       grepl(cut.out,cell_type))|
                      fafb_783_id %in% chosen.ids)
  }else if(cut.out.good%in%c("suboesophageal_zone","optic")){
    fafb.chosen.pre <- fafb.synapses %>%
      dplyr::filter(grepl(cut.out,neuropil),side=="right") %>%
      dplyr::group_by(pre) %>%
      dplyr::mutate(count = dplyr::n()) %>%
      dplyr::ungroup() %>%
      dplyr::filter(count>=100) %>%
      dplyr::pull(pre)
    fafb.chosen.post <- fafb.synapses %>%
      dplyr::filter(grepl(cut.out,neuropil),side=="right") %>%
      dplyr::group_by(post) %>%
      dplyr::mutate(count = dplyr::n()) %>%
      dplyr::ungroup() %>%
      dplyr::filter(count>=100) %>%
      dplyr::pull(post)
    chosen.ids <- unique(c(fafb.chosen.pre,fafb.chosen.post))
    fafb.meta.cutout <- fafb.meta %>%
      dplyr::filter(fafb_783_id %in% chosen.ids)
  }else{
    fafb.meta.cutout <- fafb.meta %>%
      dplyr::filter(grepl(cut.out,super_class)|
                      grepl(cut.out,cell_class)|
                      grepl(cut.out,cell_sub_class)|
                      grepl(cut.out,cell_type))
  }
  root.ids.cut.out <- na.omit(unique(fafb.meta.cutout$fafb_783_id))
  
  # Edgelist
  fafb.elist.simp.cut.out <- fafb.elist.simp %>%
    dplyr::filter(pre %in% root.ids.cut.out & post %in% root.ids.cut.out)
  
  # Synapses
  fafb.syns.cut.out <- fafb.synapses %>%
    dplyr::filter(pre %in% root.ids.cut.out | post %in% root.ids.cut.out)
  
  # Save
  arrow::write_feather(fafb.meta.cutout, file.path(save.path.cut.out,sprintf("fafb_783_%s_meta.feather",cut.out.good)))
  arrow::write_feather(fafb.elist.simp.cut.out, file.path(save.path.cut.out,sprintf("fafb_783_%s_simple_edgelist.feather",cut.out.good)))
  arrow::write_feather(fafb.syns.cut.out, file.path(save.path.cut.out,sprintf("fafb_783_%s_synapses.feather",cut.out.good)))
}

###########
### OBJ ###
###########

save.path.obj <- file.path(save.path,"obj")
dir.create(save.path.obj)
save.path.obj.np <- file.path(save.path.obj,"neuropils")
dir.create(save.path.obj.np)
Rvcg::vcgObjWrite(as.mesh3d(FAFB14NP.surf), filename = file.path(save.path.obj,"fafb14_volume_raw.obj"))
regs <- FAFB14NP.surf$RegionList
for(reg in regs){
  Rvcg::vcgObjWrite(as.mesh3d(subset(FAFB14NP.surf,reg)), filename = file.path(save.path.obj.np,sprintf("fafb14_neuropil_%s_raw.obj",reg)))
}

##############
### BUCKET ###
##############

# Send to google bucket
system("gsutil -m rsync -r /n/data1/hms/neurobio/wilson/connectomes/fafb gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/fafb_783")














