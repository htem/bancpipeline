#' banc-sjcabs — Export per-neuron SWCs, region cut-outs, neuropil meshes for collaborator share.
#'
#' @section Reads:
#'   - `banc_<ver>_meta.feather`, `banc_<ver>_synapses_enriched.parquet`,
#'     `banc_<ver>_simple_edgelist.feather`, `banc_<ver>_edgelist_split.feather`
#'
#' @section Writes (under versioned save path):
#'   - `banc_banc_space_l2_swc/`, `<cut_out>/`, `obj/` (+ GCS mirror via `banc-sjcabs-upload.R`)

###########################################################
### Export BANC data for sharing (SJCABS)
###
### Reads compiled data files from banc-data.R, exports
### skeletons, cut-outs, and OBJ meshes, then syncs
### everything to the GCS bucket.
###
### Dependencies (must run first):
###   banc/share/banc-data.R → meta, synapses, edgelists
###
### Outputs (in save.path = /n/data1/hms/neurobio/wilson/connectomes/banc):
###   banc_{version}_meta.feather
###   banc_{version}_synapses_enriched.parquet
###   banc_{version}_simple_edgelist.feather
###   banc_{version}_edgelist_split.feather
###   banc_banc_space_l2_swc/   (L2 skeletons)
###   {cut_out}/                 (7 brain regions)
###   obj/                       (neuropil meshes)
###
### Usage: Rscript banc/share/banc-sjcabs.R
###########################################################
source("banc/banc-startup.R")

local({

message("### banc: exporting data for sharing ###")
t_start <- Sys.time()
bancr::choose_banc()

#######################
### CONFIGURATION   ###
#######################

# BANC export version — must match banc-data.R
banc.version <- "821"

# Export path
save.path <- '/n/data1/hms/neurobio/wilson/connectomes/banc'
dir.create(save.path, showWarnings = FALSE)

version_id_col <- paste0("banc_", banc.version, "_id")

##############################
### READ DATA OUTPUTS      ###
##############################

message("\n=== Reading compiled data ===")

# Meta
meta_file <- file.path(banc.connectivity.save.path,
                        sprintf("banc_%s_meta.feather", banc.version))
if (!file.exists(meta_file)) {
  message("Meta feather not found — run banc/share/banc-data.R first")
  return(invisible())
}
banc.meta <- arrow::read_feather(meta_file)
root.ids <- banc.meta[[version_id_col]]
message(sprintf("  Meta: %d neurons", nrow(banc.meta)))

# Copy data files to export path
data_files <- list(
  meta = sprintf("banc_%s_meta.feather", banc.version),
  synapses = sprintf("banc_%s_synapses_enriched.parquet", banc.version),
  edgelist_simple = sprintf("banc_%s_edgelist_simple.feather", banc.version),
  edgelist_split = sprintf("banc_%s_edgelist_split.feather", banc.version)
)
for (name in names(data_files)) {
  src <- file.path(banc.connectivity.save.path, data_files[[name]])
  dst <- file.path(save.path, data_files[[name]])
  if (file.exists(src)) {
    file.copy(src, dst, overwrite = TRUE)
    message(sprintf("  Copied %s to export path", data_files[[name]]))
  } else {
    message(sprintf("  %s not found — skipping", data_files[[name]]))
  }
}

# Read edgelist and synapses for cut-outs
edgelist_src <- file.path(banc.connectivity.save.path,
                           sprintf("banc_%s_edgelist_simple.feather", banc.version))
synapses_src <- file.path(banc.connectivity.save.path,
                           sprintf("banc_%s_synapses_enriched.parquet", banc.version))

if (file.exists(edgelist_src)) {
  banc.elist.simp <- arrow::read_feather(edgelist_src)
} else {
  banc.elist.simp <- data.frame()
}

if (file.exists(synapses_src)) {
  banc.syns <- arrow::read_parquet(synapses_src)
} else {
  banc.syns <- data.frame()
}

#################
### SKELETONS ###
#################

message("\n=== Skeletons ===")

new.folder <- file.path(save.path, "banc_banc_space_l2_swc")
dir.create(new.folder, showWarnings = FALSE)
files <- list.files(banc.l2swc.save.path, full.names = TRUE)
file.copy(from = files, to = new.folder, overwrite = FALSE)
# Remove skeletons not in this version's neuron set
files <- list.files(new.folder, full.names = TRUE)
toremove <- files[!gsub("\\.swc", "", basename(files)) %in% root.ids]
if (length(toremove)) file.remove(toremove)
message(sprintf("  Copied %d L2 skeletons",
                length(list.files(new.folder, pattern = "\\.swc$"))))

################
### CUT OUTS ###
################

message("\n=== Cut outs ===")

cut.outs <- c("mushroom_body", "antennal_lobe", "central_complex",
              "optic", "suboesophageal_zone", "front_leg", "abdominal_neuromere")

for (cut.out in cut.outs) {
  cut.out.good <- snakecase::to_snake_case(cut.out)
  save.path.cut.out <- file.path(save.path, cut.out.good)
  dir.create(save.path.cut.out, showWarnings = FALSE)

  if (cut.out == "optic") {
    cut.out <- "^LO|^LOP|^AME|^ME"
  }
  if (cut.out == "antennal_lobe") {
    cut.out <- "antennal_lobe|olfactory_receptor|thermosensory_receptor|hygrosensory_receptor|CSD"
  }
  if (cut.out == "suboesophageal_zone") {
    cut.out <- "^FLA|^SEZ|^GNG|^SAD|^AMMC|^PRW"
  }
  if (cut.out == "front_leg") {
    cut.out <- "^LegNp\\(T1\\)|T1|^ProNM-T1|^LNp_T1"
  }
  if (cut.out == "abdominal_neuromere") {
    cut.out <- "^ANm|^ABDNM"
  }

  # Select neurons for this cut-out
  if (cut.out == "mushroom_body") {
    cut.out <- "mushroom_body|kenyon_cell|APL|DPM|LHMB1|OA-VPM3"
    kc.ids <- banc.meta %>%
      dplyr::filter(cell_class == "kenyon_cell", side == "right") %>%
      dplyr::pull(.data[[version_id_col]])
    kc.elist.simp <- banc.elist.simp %>%
      dplyr::filter(pre %in% kc.ids | post %in% kc.ids) %>%
      dplyr::mutate(pre = ifelse(pre %in% kc.ids, "KC", pre),
                    post = ifelse(post %in% kc.ids, "KC", post)) %>%
      dplyr::group_by(pre, post) %>%
      dplyr::mutate(count = sum(count)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(count >= 100)
    chosen.ids <- setdiff(unique(c(kc.elist.simp$pre, kc.elist.simp$post)), kc.ids)
    chosen.ids <- setdiff(chosen.ids, "KC")
    banc.meta.cutout <- banc.meta %>%
      dplyr::filter(side == "right" & (grepl(cut.out, super_class) |
                                         grepl(cut.out, cell_class) |
                                         grepl(cut.out, cell_sub_class) |
                                         grepl(cut.out, cell_type)) |
                      .data[[version_id_col]] %in% chosen.ids)
  } else if (cut.out.good %in% c("optic", "front_leg")) {
    if (nrow(banc.syns)) {
      banc.chosen.pre <- banc.syns %>%
        dplyr::filter(grepl(cut.out, neuropil), side == "right") %>%
        dplyr::group_by(pre_root_id) %>%
        dplyr::mutate(count = dplyr::n()) %>%
        dplyr::ungroup() %>%
        dplyr::filter(count >= 100) %>%
        dplyr::pull(pre_root_id)
      banc.chosen.post <- banc.syns %>%
        dplyr::filter(grepl(cut.out, neuropil), side == "right") %>%
        dplyr::group_by(post_root_id) %>%
        dplyr::mutate(count = dplyr::n()) %>%
        dplyr::ungroup() %>%
        dplyr::filter(count >= 100) %>%
        dplyr::pull(post_root_id)
      chosen.ids <- unique(c(banc.chosen.pre, banc.chosen.post))
    } else {
      chosen.ids <- character(0)
    }
    banc.meta.cutout <- banc.meta %>%
      dplyr::filter(.data[[version_id_col]] %in% chosen.ids)
  } else if (cut.out.good %in% c("suboesophageal_zone", "abdominal_neuromere")) {
    if (nrow(banc.syns)) {
      banc.chosen.pre <- banc.syns %>%
        dplyr::filter(grepl(cut.out, neuropil)) %>%
        dplyr::group_by(pre_root_id) %>%
        dplyr::mutate(count = dplyr::n()) %>%
        dplyr::ungroup() %>%
        dplyr::filter(count >= 100) %>%
        dplyr::pull(pre_root_id)
      banc.chosen.post <- banc.syns %>%
        dplyr::filter(grepl(cut.out, neuropil)) %>%
        dplyr::group_by(post_root_id) %>%
        dplyr::mutate(count = dplyr::n()) %>%
        dplyr::ungroup() %>%
        dplyr::filter(count >= 100) %>%
        dplyr::pull(post_root_id)
      chosen.ids <- unique(c(banc.chosen.pre, banc.chosen.post))
    } else {
      chosen.ids <- character(0)
    }
    banc.meta.cutout <- banc.meta %>%
      dplyr::filter(.data[[version_id_col]] %in% chosen.ids)
  } else {
    banc.meta.cutout <- banc.meta %>%
      dplyr::filter(grepl(cut.out, super_class) |
                      grepl(cut.out, cell_class) |
                      grepl(cut.out, cell_sub_class) |
                      grepl(cut.out, cell_type))
  }

  root.ids.cut.out <- na.omit(unique(banc.meta.cutout[[version_id_col]]))

  # Subset edgelist and synapses
  banc.elist.simp.cut.out <- banc.elist.simp %>%
    dplyr::filter(pre %in% root.ids.cut.out & post %in% root.ids.cut.out)
  banc.syns.cut.out <- banc.syns %>%
    dplyr::filter(pre_root_id %in% root.ids.cut.out | post_root_id %in% root.ids.cut.out)

  # Save
  arrow::write_feather(banc.meta.cutout,
    file.path(save.path.cut.out, sprintf("banc_%s_%s_meta.feather", banc.version, cut.out.good)))
  arrow::write_feather(banc.elist.simp.cut.out,
    file.path(save.path.cut.out, sprintf("banc_%s_%s_simple_edgelist.feather", banc.version, cut.out.good)))
  arrow::write_feather(banc.syns.cut.out,
    file.path(save.path.cut.out, sprintf("banc_%s_%s_synapses.feather", banc.version, cut.out.good)))
  message(sprintf("  %s: %d neurons, %d connections, %d synapses",
                  cut.out.good, nrow(banc.meta.cutout),
                  nrow(banc.elist.simp.cut.out), nrow(banc.syns.cut.out)))
}

###########
### OBJ ###
###########

message("\n=== OBJ meshes ===")

save.path.obj <- file.path(save.path, "obj")
dir.create(save.path.obj, showWarnings = FALSE)
save.path.obj.np <- file.path(save.path.obj, "neuropils")
dir.create(save.path.obj.np, showWarnings = FALSE)

Rvcg::vcgObjWrite(as.mesh3d(bancr::banc.surf),
                   filename = file.path(save.path.obj, "banc_volume_nm.obj"))
Rvcg::vcgObjWrite(as.mesh3d(bancr::banc_vnc_neuropil.surf),
                   filename = file.path(save.path.obj, "banc_vnc_neuropil_nm.obj"))
Rvcg::vcgObjWrite(as.mesh3d(bancr::banc_brain_neuropil.surf),
                   filename = file.path(save.path.obj, "banc_brain_neuropil_nm.obj"))
regs <- bancr::banc_vnc_neuropils.surf$RegionList
for (reg in regs) {
  Rvcg::vcgObjWrite(as.mesh3d(subset(banc_vnc_neuropils.surf, reg)),
                     filename = file.path(save.path.obj.np, sprintf("banc_neuropil_%s_nm.obj", reg)))
}
regs <- bancr::banc_brain_neuropils.surf$RegionList
for (reg in regs) {
  Rvcg::vcgObjWrite(as.mesh3d(subset(banc_brain_neuropils.surf, reg)),
                     filename = file.path(save.path.obj.np, sprintf("banc_neuropil_%s_nm.obj", reg)))
}
message("  OBJ meshes exported")

##############
### BUCKET ###
##############

message("\n=== GCS bucket sync ===")
system("gsutil -m rsync -r /n/data1/hms/neurobio/wilson/connectomes/banc gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data")
message("  Sync complete")

message(sprintf("\n### banc: export complete [%s] ###",
                format(round(difftime(Sys.time(), t_start, units = "mins"), 1))))

})
