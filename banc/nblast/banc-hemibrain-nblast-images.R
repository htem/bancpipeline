#' banc-hemibrain-nblast-images — Render per-query BANC↔hemibrain NBLAST review PNGs.
#'
#' @section Reads:
#'   - per-query NBLAST CSVs under `<banc.nblast.hemibrain.save.path>/results_with_mirrored/<ver>`
#'   - BANC + hemibrain BANC-space meshes
#'
#' @section Writes:
#'   - `<banc.nblast.hemibrain.save.path>/images/<ver>/todo/*.png`
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_images.sh`, `o2/production/o2_banc_nblast_dataset.sh`.

#############################
### BANC-Hemibrain IMAGES ###
#############################
source("banc/banc-startup.R")
redo <- FALSE
version <- banc.nblast.version

# Get save folder
tracing.folder <- file.path("tracing","hemibrain")
dir.create(tracing.folder, showWarnings = FALSE)

# Direct us to the BANC dataset
dir.results <- file.path(banc.nblast.hemibrain.save.path,"results_with_mirrored",version)
banc.nblast.hemibrain.obj.save.path.version <- file.path(banc.nblast.hemibrain.obj.save.path,version)
dir.create(dir.results, recursive = TRUE, showWarnings = FALSE)
dir.images <- file.path(banc.nblast.hemibrain.save.path,"images",version,"todo")
dir.create(dir.images, showWarnings = FALSE, recursive = TRUE)
dir.done <- file.path(banc.nblast.hemibrain.save.path,"images",version,"done")

# Read IDs
banc.meta.all <- banctable_query()
banc.meta <- banc.meta.all %>%
  dplyr::filter(!region %in% c("ventral_nerve_cord"),
                is.na(hemibrain_match),
                is.na(cell_type),
                !grepl("glia|trachea",super_class),
                !grepl("glia|trachea",cell_class))
banc.root.ids <- unique(banc.meta$root_id)
# Honor test mode if active
if (exists("banc.test.ids", envir = .GlobalEnv))
  banc.root.ids <- intersect(banc.root.ids, banc.test.ids)

# Get hemibrain meta data
hb.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"hemibrain_meta.csv"),
                                            col_types = hemibrainr:::sql_col_types)) %>%
  dplyr::mutate(supervoxel_id = NA) %>%
  dplyr::mutate(super_class = ifelse(is.na(super_class)|super_class=="","other",super_class)) %>%
  dplyr::mutate(cell_class = ifelse(is.na(cell_class)|cell_class=="",super_class,cell_class)) %>%
  dplyr::mutate(cell_class = ifelse(is.na(cell_class),"other",cell_class))

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

# Pre-compute done folder list ONCE (non-recursive — just directory names, not 100k+ PNGs)
if (!redo) {
  message("Scanning done folders ...")
  all.done.files <- list.files(dir.done)
  message("  Found ", length(all.done.files), " done entries")
} else {
  all.done.files <- character(0)
}

# Iterate over BANC super_class
super_classes <- sample(unique(na.omit(nblast.super.class)))
for (sp in super_classes) {

  # Save folder
  sp.good <- gsub(" |, |,|>|\\.|\\:", "_", sp)
  dir.images.ind <- file.path(dir.images, sp.good)
  dir.images.ct <- file.path(dir.images, paste0(sp.good, "_by_cell_type"))
  dir.create(dir.images.ind, recursive = TRUE, showWarnings = FALSE)
  dir.create(dir.images.ct, recursive = TRUE, showWarnings = FALSE)

  # Choose neuroglancer view
  if (grepl("central_complex", sp)) {
    url <- "https://spelunker.cave-explorer.org/#!middleauth+https://global.daf-apis.com/nglstate/api/v1/6580710816088064"
  } else if (grepl("MB|KC|DAN|MBIN|MBON", sp)) {
    url <- "https://spelunker.cave-explorer.org/#!middleauth+https://global.daf-apis.com/nglstate/api/v1/5198096897146880"
  } else {
    url <- "https://spelunker.cave-explorer.org/#!middleauth+https://global.daf-apis.com/nglstate/api/v1/5213772009439232"
  }

  # Select NBLAST files for this super_class
  nblast.files.chosen <- nblast.files.matched[nblast.super.class == sp]

  # Remove already-completed images
  if (!redo) {
    nblast.files.done <- c(paste0(gsub("matrix_", "", list.files(dir.images.ind)), ".csv"),
                           paste0(gsub("matrix_", "", all.done.files), ".csv"))
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
                     other.meta = hb.meta,
                     query.obj.save.path = banc.obj.save.path,
                     other.obj.save.path = banc.nblast.hemibrain.obj.save.path.version,
                     query.id = "root_id",
                     other.id = "bodyid",
                     numCores = 1,
                     query.mirror = FALSE,
                     region = "brain",
                     redo = redo,
                     max.hits = 3,
                     volume = banc_hemibrain.surf
  )
}

# Announce
summary <- process_directories(dir.images, recurse = 0)
rownames(summary) <- NULL
print(knitr::kable(summary))

# Purge images for neurons that now have a hemibrain_match
done.images <- list.files(dir.images, recursive = TRUE, full.names = TRUE, pattern = "png")
if (length(done.images)) {
  matched.ids <- na.omit(unique(banc.meta.all$root_id[!is.na(banc.meta.all$hemibrain_match)]))
  done.ids <- sapply(done.images, extract_root_id)
  done.images.delete <- done.images[done.ids %in% matched.ids]
  if (length(done.images.delete)) {
    del <- file.remove(done.images.delete)
    message(sprintf("Purged %d images for already-matched neurons", sum(del)))
  }
}
