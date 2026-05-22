#!/usr/bin/env Rscript

#' hemibrain-sjcabs — Assemble + publish the hemibrain v1.2.1 SJCABS data bundle.
#'
#' Pulls hemibrain v1.2.1 metadata, edgelist, NBLAST scores and synapse /
#' skeleton artefacts; aggregates into a per-version SQLite + companion
#' feathers under save.path; then rsyncs to the lab GCS bucket
#' (`compiled_data/hemibrain_121/`). Mirrors `manc/manc-sjcabs.R` and
#' `malecns/malecns-sjcabs.R`.
#'
#' @section Reads:
#'   - `franken_meta()` (hemibrain rows), `hemibrainr` snapshots
#'   - SeaTable `banc_meta` for hemibrain cross-match columns
#'   - `<banc.connectivity.save.path>/banc_<ver>_*` for BANC↔hemibrain joins
#'
#' @section Writes:
#'   - `/n/data1/hms/neurobio/wilson/connectomes/hemibrain/hemibrain_1.2.1_*.{feather,sqlite}`
#'   - GCS: `gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/hemibrain_121/`
#'
#' @section Notes:
#'   - Wrangling logic was originally developed in `fafbpipeline`
#'     (`flyconnectome/fafbpipeline:hemibrain/hemibrain-sjcabs.R`). This is
#'     the GCS-publish-side mirror, copied 2026-05-22 so bancpipeline owns
#'     the publish step. `banc/banc-startup.R` covers the library load that
#'     the original `hemibrain-startup.R` provided.
#'   - Run from bancpipeline repo root:
#'     `Rscript hemibrain/hemibrain-sjcabs.R`.
source("banc/banc-startup.R")
library(bancr)

# File paths
banc.connectivity.save.path <- "/n/data1/hms/neurobio/wilson/banc/connectivity"
save.path <- '/n/data1/hms/neurobio/wilson/connectomes/hemibrain'
dir.create(save.path)
sql.file <- "hemibrain_1.2.1_data.sqlite"

# Essential data:

## meta data: full meta data for useful labels, both anatomy (e.g. cell type, cell class) and function (i.e. body parts innervated, sparse known cell functions)
franken.meta <- franken_meta() %>%
  dplyr::distinct(cell_type, .keep_all = TRUE) %>%
  dplyr::select(region,
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
banc.meta <- banctable_query() %>%
  dplyr::distinct(hemibrain_match, .keep_all = TRUE) %>%
  dplyr::select(root_id,
                  hemibrain_match,
                  hemibrain_cell_type,
                  region,
                  hemilineage,
                  nerve,
                  flow,
                  super_class,
                  cell_class,
                  cell_sub_class,
                  cell_type,
                  neurotransmitter_predicted,
                  cell_function,
                  cell_function_detailed,
                  body_part_sensory,
                  body_part_effector)
hb.meta <- fafbseg::flytable_query("SELECT * FROM hb_info",  base="hemibrain") %>%
  dplyr::mutate(type = ifelse(is.na(type_corrected),type,type_corrected)) %>%
  dplyr::select(hemibrain_121_id=bodyId,
                instance,
                cell_type=type,
                super_class,
                cell_class,
                neurotransmitter_predicted = top_nt,
                neurotransmitter_score = top_nt_p,
                status,
                cropped,
                root=somaLocation,
                hemilineage=ito_lee_hemilineage) 
hemibrain.meta <- hb.meta %>%
  dplyr::left_join(
    franken.meta,
    by = "cell_type",
    suffix = c("_hb", "_fr")
  ) %>%
  dplyr::left_join(
    banc.meta %>%
      dplyr::filter(!is.na(hemibrain_match)),
    by = c("hemibrain_121_id" = "hemibrain_match"),
    suffix = c("", "_bymatch")
  ) %>%
  dplyr::left_join(
    banc.meta %>%
      dplyr::filter(!is.na(hemibrain_cell_type)) %>%
      dplyr::distinct(hemibrain_cell_type, .keep_all = TRUE),
    by = c("cell_type" = "hemibrain_cell_type"),
    suffix = c("", "_bytype")
  ) %>%
  dplyr::mutate(
    # Location / anatomy fields
    region_u         = dplyr::coalesce(.data$region, .data$region_bymatch, .data$region_bytype),
    hemilineage_u    = dplyr::coalesce(.data$hemilineage_hb, .data$hemilineage, .data$hemilineage_fr, .data$hemilineage_bytype),
    nerve_u          = dplyr::coalesce(.data$nerve, .data$nerve_bymatch, .data$nerve_bytype),
    flow_u           = dplyr::coalesce(.data$flow, .data$flow_bymatch, .data$flow_bytype),
    
    # Taxonomy
    super_class_u    = dplyr::coalesce(.data$super_class_fr, .data$super_class, .data$super_class_bytype, .data$super_class_hb),
    cell_class_u     = dplyr::coalesce(.data$cell_class_fr, .data$cell_class, .data$cell_class_bytype, .data$cell_class_hb),
    cell_sub_class_u = dplyr::coalesce(.data$cell_sub_class, .data$cell_sub_class_bymatch, .data$cell_sub_class_bytype),
    
    # Cell type
    cell_type_u      = dplyr::coalesce(.data$cell_type_bymatch, .data$cell_type_bytype,.data$cell_type), 
    
    # NT + function
    nt_pred_u        = .data$neurotransmitter_predicted_hb,
    nt_pred_score    = .data$neurotransmitter_score,
    cf_u             = dplyr::coalesce(.data$cell_function, .data$cell_function_bymatch, .data$cell_function_bytype),
    cfd_u            = dplyr::coalesce(.data$cell_function_detailed, .data$cell_function_detailed_bymatch, .data$cell_function_detailed_bytype),
    bps_u            = dplyr::coalesce(.data$body_part_sensory, .data$body_part_sensory_bytype, .data$body_part_sensory_bymatch),
    bpe_u            = dplyr::coalesce(.data$body_part_effector, .data$body_part_effector_bymatch, .data$body_part_effector_bytype)
    ) %>%
  dplyr::transmute(
    hemibrain_121_id,
    instance,
    cell_type            = cell_type_u,
    region               = region_u,
    hemilineage          = hemilineage_u,
    nerve                = nerve_u,
    flow                 = flow_u,
    super_class          = super_class_u,
    cell_class           = cell_class_u,
    cell_sub_class       = cell_sub_class_u,
    neurotransmitter_predicted = nt_pred_u,
    neurotransmitter_score = nt_pred_score,
    cell_function        = cf_u,
    cell_function_detailed = cfd_u,
    body_part_sensory    = bps_u,
    body_part_effector   = bpe_u,
    status,
    cropped,
    root
  ) %>%
  dplyr::mutate(super_class = dplyr::case_when(
    super_class == "central_tbc" ~ "central_brain_intrinsic",
    TRUE ~ super_class
  )) %>%
  dplyr::mutate(cell_class = dplyr::case_when(
    cell_class == "ALIN" ~ "antennal_lobe_input_neuron",
    cell_class == "ALLN" ~ "antennal_lobe_local_neuron",
    cell_class == "ALON" ~ "antennal_lobe_output_neuron",
    cell_class == "ALPN" ~ "antennal_lobe_projection_neuron",
    cell_class == "DAN" ~ "mushroom_body_dopaminergic_neuron",
    cell_class == "CX" ~ "central_complex",
    TRUE ~ cell_class
  )) %>%
  dplyr::distinct()
root.ids <- hemibrain.meta$hemibrain_121_id
# hemibrain.meta[is.na(hemibrain.meta)] <- ""
# banctable_append_rows(df = hemibrain.meta, table = "hemibrain", base = "cns_meta", bigdata = TRUE)
arrow::write_feather(hemibrain.meta, file.path(save.path,"hemibrain_121_meta.feather"))

## synapse tables: from CAVE, with a neuropil column + axon-dendrite split
hemibrain.presynapses <- read_csv(file.path(hemibrain.save.path,"hemibrain_all_neurons_synapses_polypre_centrifugal_synapses.csv"), 
                               col_types = hemibrainr:::sql_col_types) %>%
  dplyr::filter(prepost==0) %>%
  dplyr::select(connector_id,
                x,y,z,
                confidence,
                #scores, cleft_scores,
                syn_top_p, 
                syn_top_nt,
                gaba, acetylcholine, glutamate, octopamine, serotonin, dopamine,
                prepost,
                pre = bodyid,
                post = partner,
                neuropil = inside,
                # strahler_order,
                # geodesic_distance,
                # geodesic_distance_norm,
                pre_label=label) %>%
  dplyr::mutate(side="right") %>%
  dplyr::filter(pre %in% root.ids | post %in% root.ids) %>%
  dplyr::collect() %>%
  dplyr::mutate(pre_label = snakecase::to_snake_case(pre_label))
hemibrain.postsynapses <- read_csv(file.path(hemibrain.save.path,"hemibrain_all_neurons_synapses_polypre_centrifugal_synapses.csv"), 
                                 col_types = hemibrainr:::sql_col_types) %>%
  dplyr::filter(prepost==1) %>%
  dplyr::filter(bodyid %in% root.ids) %>%
  dplyr::select(connector_id,
                post=bodyid,
                post_label = label) %>%
  dplyr::collect() %>%
  dplyr::mutate(post_label = snakecase::to_snake_case(post_label))
hemibrain.synapses <- hemibrain.presynapses %>%
  dplyr::left_join(hemibrain.postsynapses,by=c("connector_id","post")) %>%
  dplyr::mutate(pre_label=ifelse(is.na(pre_label),"unknown",pre_label),
                post_label=ifelse(is.na(post_label),"unknown",post_label),
                side = dplyr::case_when(
                  grepl("\\(R\\)",neuropil) ~ "right",
                  grepl("\\(L\\)",neuropil) ~ "left",
                  TRUE ~ NA
                ),
                neuropil = gsub("\\(R\\)|\\(L\\)","",neuropil))
pre.ids.missing <- unique(setdiff(root.ids,hemibrain.synapses$pre))
post.ids.missing <- unique(setdiff(root.ids,hemibrain.synapses$post))
arrow::write_parquet(hemibrain.synapses, file.path(save.path,"hemibrain_121_synapses.parquet"))

## split edgelists: proofread neurons and connections between them
hemibrain.elist <- read_csv(file.path(hemibrain.save.path,"hemibrain_all_neurons_edgelist_polypre_centrifugal_synapses.csv"), col_types = hemibrainr:::sql_col_types) %>%
  dplyr::mutate(post_count=as.numeric(post_count)) %>%
  dplyr::distinct(pre,
                  post,
                  post_label,
                  pre_label,
                  count,
                  norm,
                  connection,
                  total_input=post_count)
pre.ids.missing <- unique(setdiff(root.ids,hemibrain.elist$pre))
post.ids.missing <- unique(setdiff(root.ids,hemibrain.elist$post))
arrow::write_feather(hemibrain.elist, file.path(save.path,"hemibrain_121_split_edgelist.feather"))

## simple edgelists: proofread neurons and connections between them
hemibrain.elist.simp <- hemibrain.elist %>%
  dplyr::group_by(pre, post) %>%
  dplyr::mutate(count = sum(count),
                norm = count/sum(total_input)) %>%
  dplyr::ungroup() %>%
  dplyr::select(pre, post, count, norm)
pre.ids.missing <- unique(setdiff(root.ids,hemibrain.elist.simp$pre))
post.ids.missing <- unique(setdiff(root.ids,hemibrain.elist.simp$post))
arrow::write_feather(hemibrain.elist.simp, file.path(save.path,"hemibrain_121_simple_edgelist.feather"))

## skeletons: including compartmentalised skeletons (optional: axons vs. dendrites)

### In hemibrain space
hemibrain.in.banc.space <- "/n/data1/hms/neurobio/wilson/fafbz/hemibrainr_data//hemibrain_neurons/split"
# out_zip <- file.path(save.path, "hemibrain_hemibrain_raw_space_swc.zip")
# oldwd <- setwd(hemibrain.in.banc.space); on.exit(setwd(oldwd), add = TRUE)
# files <- list.files(".", recursive = TRUE, all.files = TRUE, no.. = TRUE)
# if (!requireNamespace("zip", quietly = TRUE)) install.packages("zip")
# zip::zipr(zipfile = out_zip, files = files)
# cat("Wrote:", out_zip, "\n")
new.folder <- file.path(save.path, "hemibrain_hemibrain_raw_space_swc")
dir.create(new.folder)
files <- list.files(fafb.in.banc.space, full.names = TRUE)
file.copy(from = files, to = new.folder, overwrite = FALSE)

### In BANC space
hemibrain.in.banc.space <- "/n/data1/hms/neurobio/wilson/banc/matching/hemibrain/hemibrain_banc_space_swc/elastix_tpsreg_240721"
# out_zip <- file.path(save.path, "hemibrain_banc_space_swc.zip")
# oldwd <- setwd(hemibrain.in.banc.space); on.exit(setwd(oldwd), add = TRUE)
# files <- list.files(".", recursive = TRUE, all.files = TRUE, no.. = TRUE)
# if (!requireNamespace("zip", quietly = TRUE)) install.packages("zip")
# zip::zipr(zipfile = out_zip, files = files)
# cat("Wrote:", out_zip, "\n")
new.folder <- file.path(save.path, "hemibrain_banc_space_swc")
dir.create(new.folder)
files <- list.files(fafb.in.banc.space, full.names = TRUE)
file.copy(from = files, to = new.folder, overwrite = FALSE)

# src=/n/data1/hms/neurobio/wilson/fafbz/hemibrainr_data//hemibrain_neurons/split
# dst=/n/data1/hms/neurobio/wilson/connectomes/hemibrain/hemibrain_hemibrain_raw_space_swc/
# rsync -av "$src" "$dst"

# src=/n/data1/hms/neurobio/wilson/banc/matching/hemibrain/hemibrain_banc_space_swc/elastix_tpsreg_240721
# dst=/n/data1/hms/neurobio/wilson/connectomes/hemibrain/hemibrain_banc_space_swc/
# rsync -av "$src" "$dst"

###########
### OBJ ###
###########

save.path.obj <- file.path(save.path,"obj")
dir.create(save.path.obj)
save.path.obj.np <- file.path(save.path,"neuropils")
dir.create(save.path.obj.np)
Rvcg::vcgObjWrite(as.mesh3d(hemibrainr::hemibrain.surf), filename = file.path(save.path.obj,"hemibrain_volume_raw.obj"))
Rvcg::vcgObjWrite(as.mesh3d(hemibrainr::hemibrain_microns.surf), filename = file.path(save.path.obj,"hemibrain_volume_microns.obj"))
regs <- hemibrain_al.surf$RegionList
for(reg in regs){
  Rvcg::vcgObjWrite(as.mesh3d(subset(hemibrain_al.surf,reg)), filename = file.path(save.path.obj.np,sprintf("hemibrain_antennal_lobe_glomerulus_%s_raw.obj",reg)))
}
for(reg in regs){
  Rvcg::vcgObjWrite(as.mesh3d(subset(hemibrain_al_microns.surf,reg)), filename = file.path(save.path.obj.np,sprintf("hemibrain_antennal_lobe_glomerulus_%s_microns.obj",reg)))
}

################
### CUT OUTS ###
################

# Go by cut-out
cut.outs <- c("mushroom_body","antennal_lobe","central_complex")
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
    kc.ids <- hemibrain.meta %>%
      dplyr::filter(cell_class=="kenyon_cell") %>%
      dplyr::pull(hemibrain_121_id)
    kc.elist.simp <- hemibrain.elist.simp %>%
      dplyr::filter(pre%in%kc.ids|post%in%kc.ids) %>%
      dplyr::mutate(pre=ifelse(pre%in%kc.ids,"KC",pre),
                    post=ifelse(post%in%kc.ids,"KC",post)) %>%
      dplyr::group_by(pre,post) %>%
      dplyr::mutate(count=sum(count)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(count>=100)
    chosen.ids <- setdiff(unique(c(kc.elist.simp$pre,kc.elist.simp$post)),kc.ids)
    chosen.ids <- setdiff(chosen.ids,"KC")
    hemibrain.meta.cutout <- hemibrain.meta %>%
      dplyr::filter((grepl(cut.out,super_class)|
                                       grepl(cut.out,cell_class)|
                                       grepl(cut.out,cell_sub_class)|
                                       grepl(cut.out,cell_type))|
                      hemibrain_121_id %in% chosen.ids)
  }else if(cut.out.good%in%c("suboesophageal_zone","optic")){
    hemibrain.chosen.pre <- hemibrain.synapses %>%
      dplyr::filter(grepl(cut.out,neuropil)) %>%
      dplyr::group_by(pre) %>%
      dplyr::mutate(count = dplyr::n()) %>%
      dplyr::ungroup() %>%
      dplyr::filter(count>=100) %>%
      dplyr::pull(pre)
    hemibrain.chosen.post <- hemibrain.synapses %>%
      dplyr::filter(grepl(cut.out,neuropil)) %>%
      dplyr::group_by(post) %>%
      dplyr::mutate(count = dplyr::n()) %>%
      dplyr::ungroup() %>%
      dplyr::filter(count>=100) %>%
      dplyr::pull(post)
    chosen.ids <- unique(c(hemibrain.chosen.pre,hemibrain.chosen.post))
    hemibrain.meta.cutout <- hemibrain.meta %>%
      dplyr::filter(hemibrain_121_id %in% chosen.ids)
  }else{
    hemibrain.meta.cutout <- hemibrain.meta %>%
      dplyr::filter(grepl(cut.out,super_class)|
                      grepl(cut.out,cell_class)|
                      grepl(cut.out,cell_sub_class)|
                      grepl(cut.out,cell_type))
  }
  root.ids.cut.out <- na.omit(unique(hemibrain.meta.cutout$hemibrain_121_id))
  
  # Edgelist
  hemibrain.elist.simp.cut.out <- hemibrain.elist.simp %>%
    dplyr::filter(pre %in% root.ids.cut.out & post %in% root.ids.cut.out)
  
  # Synapses
  hemibrain.syns.cut.out <- hemibrain.synapses %>%
    dplyr::filter(pre %in% root.ids.cut.out | post %in% root.ids.cut.out)
  
  # Save
  arrow::write_feather(hemibrain.meta.cutout, file.path(save.path.cut.out,sprintf("hemibrain_121_%s_meta.feather",cut.out.good)))
  arrow::write_feather(hemibrain.elist.simp.cut.out, file.path(save.path.cut.out,sprintf("hemibrain_121_%s_simple_edgelist.feather",cut.out.good)))
  arrow::write_feather(hemibrain.syns.cut.out, file.path(save.path.cut.out,sprintf("hemibrain_121_%s_synapses.feather",cut.out.good)))
}

##############
### BUCKET ###
##############

# Send to google bucket
system("gsutil -m rsync -r /n/data1/hms/neurobio/wilson/connectomes/hemibrain gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/hemibrain_121")





