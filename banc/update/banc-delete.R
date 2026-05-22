#' banc-delete — Delete redundant / stale per-neuron data files.
#'
#' Sweeps the per-neuron file trees (OBJ, SWC, L2 SWC, split SWC, synapses,
#' images) for files whose root_id no longer appears in SeaTable, removes
#' them, and reports coverage statistics per tree.
#'
#' @section Reads:
#'   - SeaTable `banc_meta`
#'   - per-neuron file trees under `banc.{obj,swc,l2swc,split,l2split}.save.path`
#'
#' @section Writes:
#'   - Removes stale per-neuron files (.obj / .swc / .csv / .png)
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_v888_rebuild.sh`,
#'   `o2/production/o2_banc_update.sh`,
#'   `o2/production/o2_banc_nblast.sh`,
#'   `o2/oneshots/o2_banc_v850_rebuild.sh`,
#'   `o2/oneshots/o2_banc_v890_rebuild.sh`

##########################################
### Delete redundant and old BANC data ###
##########################################
source("banc/banc-startup.R")

# Report covergae function
report_coverage <- function(dest, 
                            ids = NULL, 
                            done = NULL,
                            pattern = NULL,
                            id.update = FALSE,
                            banner = "BANC PIPELINE COUNT",
                            dataset = "BANC neurons"){
  
  # Files already computed
  if(is.null(done)){
    n.files <- list_files_age_sorted(dest, pattern = pattern, recursive = TRUE)
    if(grepl("root_id",n.files[1])){
      done <- regmatches(basename(n.files), regexpr("(?<=root_id_)\\d+", basename(n.files), perl = TRUE))
    }else{
      done <- gsub(".*_root_id_|\\.csv|\\.png|\\.swc|\\.obj","",basename(n.files))
    }
    if(id.update){
      done <- banc_updateids(done) 
    }
  }
  if(is.null(ids)){
    banc.ids <- banctable_query()
    ids <- unique(banc.ids$root_id)
    if(id.update){
      ids <- banc_updateids(ids) 
    }
  }
 
  # IDs remaining to compute
  undone <- setdiff(ids, done)
  
  # Announce
  percent <- 100*length(undone)/length(ids)
  message(sprintf("**** %s **** %s/%s (%s percent) %s uncomputed in: %s **** %s ****", 
          banner, length(undone), length(ids), round(percent,2), dataset, dest, banner))
  percent/100
}

# Read IDs
banc.meta <- banctable_query() %>%
  banc_filter_neurons()
banc.ids.local <- readr::read_csv(file = file.path(banc.meta.save.path, "banc_ids.csv"),
                                   col_types = banc.col.types, show_col_types = FALSE)
banc.root.ids <- unique(c(unique(banc.ids.local$root_id), unique(banc.ids.local$root_626),
                          unique(banc.ids.local$root_850), unique(banc.ids.local$root_888)))

# Stable version IDs — files for these versions are never deleted
stable_version_ids <- unique(c(banc.ids.local$root_626, banc.ids.local$root_850,
                               banc.ids.local$root_888))
stable_version_ids <- stable_version_ids[!is.na(stable_version_ids) & stable_version_ids != ""]
message(sprintf("Protected stable version IDs: %d (v626 + v850 + v888)", length(stable_version_ids)))

# Remove .swc files no longer in banc.ids
message("##### Deleting redundant BANC .swc/.obj/NBLAST files ... #####")
for(save.dir in c(
  # file.path(banc.nblast.mirror.save.path,"results"),
                  file.path(banc.nblast.mirror.save.path,"images","todo"),
                  file.path(banc.nblast.manc.save.path,"results"),
                  file.path(banc.nblast.manc.save.path,"images","todo"),
                  file.path(banc.nblast.fafb.save.path,"results"),
                  file.path(banc.nblast.fafb.save.path,"images","todo"),
                  file.path(banc.nblast.hemibrain.save.path,"results_with_mirrored"),
                  file.path(banc.nblast.hemibrain.save.path,"images","todo"),
                  banc.obj.save.path,
                  #banc.swc.save.path,
                  #banc.metrics.save.path,
                  banc.synapses.save.path,
                  banc.l2swc.save.path)){
  
  # Report coverage
  res <- report_coverage(ids = banc.root.ids, dest = save.dir)
  
  # Remove old files no longer in banc.root.ids
  message("Deleting redundant files from: ", save.dir)
  extant.files <- list.files(save.dir, full.names = TRUE, pattern = "\\.obj$|\\.swc$|\\.csv$|\\.png$", recursive = TRUE)
  if(length(extant.files) && (grepl('png',extant.files[1])|grepl('root_id',extant.files[1]))){
    extant.ids <- regmatches(basename(extant.files), regexpr("(?<=root_id_)\\d+", basename(extant.files), perl = TRUE))
  }else{
    extant.ids <- gsub("\\.obj$|\\.swc$|\\.csv$|\\.png$","", basename(extant.files))
  }
  out.of.date <- banc_islatest(extant.ids)
  is_stable <- extant.ids %in% stable_version_ids
  files.to.delete <- extant.files[(!out.of.date | !extant.ids %in% banc.root.ids) & !is_stable]
  n_protected <- sum(is_stable & (!out.of.date | !extant.ids %in% banc.root.ids))
  message("removing ", length(files.to.delete), " outdated files: ", round(length(files.to.delete)/length(extant.files),4)*100,
          "% (", n_protected, " stable version files protected)")
  file.remove(files.to.delete)
}

# # Remove old match folders
# message("##### Deleting redundant BANC screening files ... #####")
# for(save.dir in c(file.path(banc.nblast.mirror.save.path,"images"))){
#   # Remove old files no longer in banc.root.ids
#   pngs <- list.files(save.dir, pattern = "\\.png$", recursive = TRUE, full.names = TRUE)
#   pngs.to.delete <- c()
#   for(png in pngs){
#     id <- gsub(".*_root_id_|_nucleus_.*|\\.png","",basename(png))
#     if(!id%in%banc.root.ids){
#       pngs.to.delete <- c(pngs.to.delete,png)
#     }
#   }
#   message("Removing ", length(pngs.to.delete), " outdated files")
#   file.remove(files.to.delete)
# }

# Clean up old job files (>30 days)
job_dir <- "jobs"
if (dir.exists(job_dir)) {
  job_files <- list.files(job_dir, full.names = TRUE, pattern = "\\.(out|err)$")
  if (length(job_files)) {
    ages <- difftime(Sys.time(), file.mtime(job_files), units = "days")
    old_files <- job_files[ages > 30]
    if (length(old_files)) {
      message(sprintf("Removing %d job files older than 30 days", length(old_files)))
      file.remove(old_files)
    }
  }
}

######################################
## Delete v746 versioned resources ####
######################################

# # Uncomment to remove v746-specific files once v850 rebuild is confirmed.
# message("##### Deleting v746 versioned resources #####")
#
# # v746 synapse/connectivity data directory
# v746_dir <- file.path(banc.save.path, "v746")
# if (dir.exists(v746_dir)) {
#   message(sprintf("  Removing v746 data directory: %s", v746_dir))
#   unlink(v746_dir, recursive = TRUE)
# }
#
# # v746 influence results
# v746_influence <- file.path(banc.influence.save.path, "banc_746")
# v746_influence_rev <- file.path(banc.influence.save.path, "banc_746_reverse")
# for (d in c(v746_influence, v746_influence_rev)) {
#   if (dir.exists(d)) {
#     message(sprintf("  Removing v746 influence directory: %s", d))
#     unlink(d, recursive = TRUE)
#   }
# }
#
# # v746 meta/edgelist feathers on shared storage
# banc_share <- "/n/data1/hms/neurobio/wilson/connectomes/banc"
# v746_files <- c(
#   file.path(banc_share, "banc_746_meta.feather"),
#   file.path(banc_share, "banc_746_simple_edgelist.feather"),
#   file.path(banc_share, "banc_746_synapses.parquet"),
#   file.path(banc_share, "banc_746_synapses_enriched.parquet"),
#   file.path(banc_share, "banc_746_edgelist_split.feather")
# )
# for (f in v746_files) {
#   if (file.exists(f)) {
#     message(sprintf("  Removing: %s", f))
#     file.remove(f)
#   }
# }
#
# # v746 GCS files (run manually or uncomment)
# # system2("gsutil", c("rm", "gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/banc_746_meta.feather"))
# # system2("gsutil", c("rm", "-r", "gs://lee-lab_brain-and-nerve-cord-fly-connectome/connectivity/banc_746/"))
#
# message("##### v746 cleanup complete #####")

# Announce
message("##### BANCpipeline: deleted outdated BANC data files #####")


