#' banc-sort-folders — Bidirectional rsync between O2 working and fileserve mirrors.
#'
#' `/n/data1/.../banc/` ↔ `/n/files/Neurobio/.../banc/` (csv / sqlite / png / obj / swc).
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_images.sh`.

##############################
### ORGANISE BANC MATCHING ###
##############################
source("banc/banc-startup.R")

############
### Main ###
############

# Remote paths
cat("Synchronising matching files between O2 and the fileserve")
A <- '/n/data1/hms/neurobio/wilson/banc/'
B <- '/n/files/Neurobio/wilsonlab/banc/'

### meta ###

# # Remote synchronization of meta data
#sync_files(path(B, "connectivity/"), path(A, "connectivity/"), extensions = c("sqlite"), move.old = FALSE)
remove_empty_dirs(path(A, "connectivity/"))
sync_files(path(A, "connectivity/"), path(B, "connectivity/"), extensions = c("csv", "sqlite"), move.old = FALSE)
remove_empty_dirs(path(A, "connectivity/"))
sync_files(path(A, "meta/"), path(B, "meta/"), extensions = c("csv", "png"), move.old = FALSE)
remove_empty_dirs(path(A, "meta/"))

### logs ###
sync_files("logs/", path(A, "logs/"), extensions = c("txt"), move.old = FALSE)
sync_files(path(A, "logs/"), path(B, "logs/"), extensions = c("txt"), move.old = FALSE)
remove_empty_dirs(path(A, "logs/"))
# sync_files(path(A, "l2split/"), path(B, "l2split/"), extensions = c("csv", "png"), move.old = FALSE)
# remove_empty_dirs(path(A, "l2split/"))
# sync_files(path(A, "matching/neck_connective"), path(B, "matching/neck_connective"), extensions = c("csv", "png"), move.old = FALSE)
# remove_empty_dirs(path(A, "matching/"))
# sync_files(path(B, "influence/"), path(A, "influence/"), extensions = c("csv", "sqlite"), move.old = FALSE)
# remove_empty_dirs(path(A, "influence/"))

# # From clipboard!
# files.to.move <- unlist(strsplit(a,split="\\n"))
# files.to.move <- gsub("",":",files.to.move)
# pngs <- list.files(file.path(A, 'matching/fafb/images', banc.nblast.version, 'todo/'), full.names = TRUE, recursive = TRUE)
# pngs.to.move <- pngs[basename(pngs)%in%unlist(files.to.move)]
# pngs.to.move <- pngs.to.move[!file.exists(file.path('/n/data1/hms/neurobio/wilson/banc/matching/fafb/correct/3_good',basename(pngs.to.move)))]
# if(length(pngs.to.move)){
#   file.rename(pngs.to.move, file.path('/n/data1/hms/neurobio/wilson/banc/matching/fafb/correct/3_good/',basename(pngs.to.move)))
# }

### mirror ###

# # #  Move specific groups
# bc <- banctable_query()
# dir.create('/n/data1/hms/neurobio/wilson/banc/matching/mirror/images/todo/neck_connective')
# nc.neck <- subset(bc, region == "neck_connective")
# nc.neck.ids <- unique(nc.neck$root_id)
# move_files_preserve_structure(source_dir = path(A, "matching/mirror/images/todo"),
#                               dest_dir = path(A, "matching/mirror/images/todo/neck_connective"),
#                               ids=nc.neck.ids)
# remove_empty_dirs(path(A, "matching/mirror/images/todo"))
# move_files_preserve_structure(source_dir = path(A, "matching/mirror/images/done"),
#                               dest_dir = path(A, "matching/mirror/images/done/neck_connective"),
#                               ids=nc.neck.ids)
# remove_empty_dirs(path(A, "matching/mirror/images/done"))

### fafb ###
sync_files(path(B, "matching/fafb/correct/"), path(A, "matching/fafb/correct/"), extensions = c("png"), move.old = FALSE, remote.source = TRUE)
sort_todos(base_images=path(A, "matching/fafb/images", banc.nblast.version),base_correct=path(A, "matching/fafb/correct/"), use.latest = FALSE, remove.old.roots = FALSE)
#try({sort_todos(base_images=path(A, "matching/fafb/images", banc.nblast.version),base_correct=path(A, "matching/fafb/correct/"), use.latest = TRUE, remove.old.roots = FALSE)})
remove_empty_dirs(path(A, "matching/fafb/images", banc.nblast.version))
sync_files(path(A, "matching/fafb/images", banc.nblast.version), path(B, "matching/fafb/images/"), extensions = c("png"), move.old = FALSE)

### manc ###
sync_files(path(B, "matching/manc/correct/"), path(A, "matching/manc/correct/"), extensions = c("png"), move.old = FALSE, remote.source = TRUE)
sort_todos(base_images=path(A, "matching/manc/images", banc.nblast.version),base_correct=path(A, "matching/manc/correct/"), use.latest = FALSE, remove.old.roots = FALSE)
#try({sort_todos(base_images=path(A, "matching/manc/images", banc.nblast.version),base_correct=path(A, "matching/manc/correct/"), use.latest = TRUE, remove.old.roots = FALSE)})
remove_empty_dirs(path(A, "matching/manc/images", banc.nblast.version))
sync_files(path(A, "matching/manc/images", banc.nblast.version), path(B, "matching/manc/images/"), extensions = c("png"), move.old = FALSE)

### fanc ###
sync_files(path(B, "matching/fanc/correct/"), path(A, "matching/fanc/correct/"), extensions = c("png"), move.old = FALSE, remote.source = TRUE)
sort_todos(base_images=path(A, "matching/fanc/images", banc.nblast.version),base_correct=path(A, "matching/fanc/correct/"), use.latest = FALSE, remove.old.roots = FALSE)
#try({sort_todos(base_images=path(A, "matching/fanc/images", banc.nblast.version),base_correct=path(A, "matching/fanc/correct/"), use.latest = TRUE, remove.old.roots = FALSE)})
remove_empty_dirs(path(A, "matching/fanc/images/"))
sync_files(path(A, "matching/fanc/images", banc.nblast.version), path(B, "matching/fanc/images/"), extensions = c("png"), move.old = FALSE)

# Update correct mirror-matches, remote sync
sync_files(path(B, "matching/mirror/correct/"), path(A, "matching/mirror/correct/"),  extensions = c("png"), move.old = FALSE)
sort_todos(base_images=path(A, "matching/mirror/images/"),base_correct=path(A, "matching/mirror/correct/"), use.latest = FALSE, remove.old.roots = FALSE)
remove_empty_dirs(path(A, "matching/mirror/images/"))
sync_files(path(A, "matching/mirror/images/"), path(B, "matching/mirror/images/"), extensions = c("png"), move.old = FALSE)

### hemibrain ###
sync_files(path(B, "matching/hemibrain/correct/"), path(A, "matching/hemibrain/correct/"), extensions = c("png"), move.old = FALSE, remote.source = TRUE)
sort_todos(base_images=path(A, "matching/hemibrain/images/"),base_correct=path(A, "matching/hemibrain/correct/"), use.latest = FALSE, remove.old.roots = FALSE)
try({sort_todos(base_images=path(A, "matching/hemibrain/images/"),base_correct=path(A, "matching/hemibrain/correct/"), use.latest = TRUE, remove.old.roots = FALSE)})
remove_empty_dirs(path(A, "matching/hemibrain/images/"))
sync_files(path(A, "matching/hemibrain/images/"), path(B, "matching/hemibrain/images/"), extensions = c("png"), move.old = FALSE)

### malecns ###
sync_files(path(B, "matching/malecns/correct/"), path(A, "matching/malecns/correct/"), extensions = c("png"), move.old = FALSE, remote.source = TRUE)
sort_todos(base_images=path(A, "matching/malecns/images", banc.nblast.malecns.version),base_correct=path(A, "matching/malecns/correct/"), use.latest = FALSE, remove.old.roots = FALSE)
remove_empty_dirs(path(A, "matching/malecns/images", banc.nblast.malecns.version))
sync_files(path(A, "matching/malecns/images", banc.nblast.malecns.version), path(B, "matching/malecns/images/"), extensions = c("png"), move.old = FALSE)

# Announce
cat("File transfer and movement completed.\n")

# # Fix folder names
# bc <- banctable_query() %>%
#   dplyr::filter(!is.na(root_id))
# folder_path <- path(A, "matching/manc/images", banc.nblast.version)
# 
# # Get all folder names
# folder_names <- list.dirs(folder_path, full.names = TRUE, recursive = TRUE)
# 
# # Filter folders containing 'nucleus_id_'
# nucleus_folders <- folder_names[grep("nucleus_id_", folder_names)]
# 
# # Process each nucleus folder
# for (folder in nucleus_folders) {
#   
#   # Extract the current nucleus ID
#   current_id <- str_extract(basename(folder), "(?<=nucleus_id_)\\d{16}")
#   if(is.na(current_id)){
#     current_id <- str_extract(basename(folder), "(?<=root_id_)\\d{16}")
#     new_id <- bc$supervoxel_id[bc$nucleus_id==current_id]
#   }else{
#     new_id <- bc$supervoxel_id[bc$nucleus_id==current_id]
#   }
# 
#   # Generate a new random ID
#   if(!length(new_id)){
#     next
#   }
#   if(is.na(new_id)){
#     next
#   }
#   if(new_id=="0"){
#     next
#   }
# 
#   # Create the new folder name
#   new_folder <- str_replace(folder, paste0("nucleus_id_", current_id), paste0("supervoxel_id_", new_id))
#   
#   # Rename the folder
#   file.rename(file.path(folder_path, folder), file.path(folder_path, new_folder))
#   
#   # Rename files inside folder
#   png_files <- list.files(folder, pattern = "\\.png$", full.names = TRUE)
#   png_files <- png_files[grep("nucleus_id_", png_files)]
#   for(png_file in png_files){
#     new_png_file <- str_replace(png_file, paste0("nucleus_id_", current_id), paste0("supervoxel_id_", new_id))
#     file.rename(file.path(folder_path, folder, new_png_file), file.path(folder_path, folder, new_png_file))
#   }
#   
#   # Print the change
#   cat(paste("Renamed:", folder, "to", new_folder, "\n"))
# }




