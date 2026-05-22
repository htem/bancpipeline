#' banc-lr-nblast-images — Render per-query left-right BANC mirror NBLAST review PNGs.
#'
#' Restricted to BANC neurons with a strong FAFB or MANC NBLAST match.
#'
#' @section Reads:
#'   - per-query mirror NBLAST CSVs, BANC meshes
#'
#' @section Writes:
#'   - `<banc.nblast.mirror.save.path>/images/todo/*.png`
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_images.sh`.

##############################
### Left-right BANC IMAGES ###
##############################
source("banc/banc-startup.R")
redo <- FALSE

# Read IDs
banc.meta <- banctable_query() %>%
  dplyr::filter(fafb_nblast >= 0.65|manc_nblast >= 0.65) %>%
  dplyr::mutate(cell_class = gsub("auto:","",cell_class))
banc.root.ids <- unique(banc.meta$root_id)
# Honor test mode if active
if (exists("banc.test.ids", envir = .GlobalEnv))
  banc.root.ids <- intersect(banc.root.ids, banc.test.ids)

# Direct us to the BANC dataset
dir.results <- file.path(banc.nblast.mirror.save.path,"results")
dir.create(dir.results, showWarnings = FALSE, recursive = TRUE)
dir.images <- file.path(banc.nblast.mirror.save.path,"images","todo")
dir.create(dir.images, showWarnings = FALSE, recursive = TRUE)

# Find NBLAST CSV files and match to banc.meta by root_id or supervoxel_id
nblast.files <- list.files(dir.results, pattern = "\\.csv", full.names = TRUE)
nblast.file.roots <- gsub(".*_root_id_|\\.csv", "", basename(nblast.files))
nblast.file.svids <- stringr::str_match(basename(nblast.files), "supervoxel_id_([0-9]+)_root_id")[, 2]

# Match by root_id first, fall back to supervoxel_id (stable across edits)
match.by.root <- match(nblast.file.roots, banc.meta$root_id)
match.by.sv <- match(nblast.file.svids, banc.meta$supervoxel_id)
nblast.banc.idx <- ifelse(!is.na(match.by.root), match.by.root, match.by.sv)
nblast.in.meta <- !is.na(nblast.banc.idx)

message(sprintf("Found %d NBLAST files, %d matching banc.meta (%d by root_id, %d by supervoxel_id)",
                length(nblast.files), sum(nblast.in.meta),
                sum(!is.na(match.by.root) & nblast.in.meta),
                sum(is.na(match.by.root) & !is.na(match.by.sv))))

# Map files to super_class
nblast.files.matched <- nblast.files[nblast.in.meta]
nblast.super.class <- banc.meta$super_class[nblast.banc.idx[nblast.in.meta]]

# No pre-computed image list — banc_nblast_images() checks per-super_class dir only

# Iterate over BANC super_class
super_classes <- sample(unique(na.omit(nblast.super.class)))
for (sp in super_classes) {

  # Save folder
  sp.good <- gsub(" |, |,|>|\\.|\\:", "_", sp)
  dir.images.ind <- file.path(dir.images, sp.good)
  dir.create(dir.images.ind, recursive = TRUE, showWarnings = FALSE)

  # Select NBLAST files for this super_class
  nblast.files.chosen <- nblast.files.matched[nblast.super.class == sp]

  # Remove already-completed images
  if (!redo) {
    nblast.files.done <- paste0(gsub("matrix_", "", list.files(dir.images.ind)), ".csv")
    nblast.files.chosen <- nblast.files.chosen[!basename(nblast.files.chosen) %in% nblast.files.done]
  }
  completed.images <- NULL

  if (!length(nblast.files.chosen)) {
    message("No hits to process for: ", sp)
    next
  }

  # Create 2D screening files
  message("Working on ~", length(nblast.files.chosen), " for super_class: ", sp)
  banc_nblast_images(nblast.files = sample(nblast.files.chosen),
                     dir.images = dir.images.ind,
                     completed.images = completed.images,
                     query.meta = banc.meta,
                     other.meta = banc.meta,
                     query.obj.save.path = banc.obj.save.path,
                     other.obj.save.path = banc.obj.save.path,
                     query.id = "root_id",
                     other.id = "root_id",
                     numCores = numCores,
                     query.mirror = FALSE,
                     redo = redo,
                     max.hits = 3,
                     volume = NULL
  )
}

# Announce
message("##### BANCpipeline: banc LR images updated #####")
img.results <- file.path(banc.nblast.mirror.save.path,"images")
imgs <- list.files(img.results, recursive = TRUE, full.names = TRUE)
message(sprintf("##### we have generated a total of %s images or ~%s neurons", length(imgs),floor(length(imgs)/5)))

# Organise neck connective task
move_files_preserve_structure <- function(source_dir, dest_dir, ids, pattern = "*.png") {

  # List all files in the source directory, including sub-directories
  files_to_move <- dir_ls(source_dir, recurse = TRUE, type = "file", glob = pattern)
  files_to_move.moved <- dir_ls(dest_dir, recurse = TRUE, type = "file", glob = pattern)
  files_to_move <- files_to_move[!basename(files_to_move)%in%basename(files_to_move.moved)]
  queries <- regmatches(files_to_move, regexpr("(?<=root_id_)\\d+", files_to_move, perl = TRUE))
  correct.queries <- banc_latestid(queries)
  correct.queries <- correct.queries[correct.queries!="0"]
  files_to_move <- files_to_move[correct.queries%in%ids]
  files_to_move <- files_to_move[!is.na(correct.queries)]

  # For each file
  for (file in files_to_move) {
    # Get the relative path
    rel_path <- path_rel(file, start = source_dir)

    # Construct the new path
    new_path <- path(dest_dir, rel_path)

    # Create the directory structure if it doesn't exist
    dir_create(path_dir(new_path), recurse = TRUE)

    # Move the file
    file_move(file, new_path)

    # Print progress
    cat("Moved:", rel_path, "\n")
  }
}
