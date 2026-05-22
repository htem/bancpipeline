#' banc-calculate-skeletons — Skeletonise OBJ meshes into radius-bearing SWCs.
#'
#' Uses `fafbseg::skeletor()` (Python skeletor). Proofread neurons only;
#' tracks failures in `banc_failed.csv`, retries on Saturdays.
#'
#' @section Reads:
#'   - OBJ meshes, SeaTable `banc_meta`
#'
#' @section Writes:
#'   - `<banc.swc.save.path>/<root_id>.swc`, `banc_failed.csv`
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_split.sh`, `o2/production/o2_banc_update.sh`.

###########################################################
### Calculate accurate BANC skeletons from OBJ meshes
###
### Uses fafbseg::skeletor() (Python skeletor package) to
### convert OBJ meshes into SWC skeleton files with radius
### information at each point.
###
### Pegged to a BANC version (default 850): only processes
### root_{version} neurons from seatable.
###
### Only processes proofread neurons with OBJ files that
### have not yet been skeletonised. Failed skeletonisations
### are tracked in banc_failed.csv and retried on Saturdays.
###
### Input:  OBJ meshes in banc.obj.save.path
### Output: SWC files in banc.swc.save.path
###########################################################
source("banc/banc-startup.R")

local({

message("### banc: calculating accurate skeletons ###")
t_start <- Sys.time()

# Must set this due to bug in skeletor R code
choose_segmentation("flywire31")

#######################
### CONFIGURATION   ###
#######################

# banc.version set in banc-startup.R

###########################
### Read current state  ###
###########################

banc.meta <- banctable_query() %>%
  dplyr::filter(!grepl("GLIA|TRACHEA|NOT_A_NEURON|DELETE|DEBRIS", status)) %>%
  dplyr::filter(proofread == "TRUE" | roughly_proofread == "TRUE")
banc.root.ids <- unique(banc.meta$root_id)
message(sprintf("Target neurons (root_%s): %d", banc.version, length(banc.root.ids)))
# Honor test mode if active
if (exists("banc.test.ids", envir = .GlobalEnv))
  banc.root.ids <- intersect(banc.root.ids, banc.test.ids)

# Get OBJ files for current neurons
obj.full.files <- list.files(banc.obj.save.path, pattern = ".obj$", full.names = TRUE)
obj.ids <- gsub(pattern = "\\.obj$", "", basename(obj.full.files))
obj.full.files <- obj.full.files[obj.ids %in% banc.root.ids]
obj.ids <- intersect(obj.ids, banc.root.ids)

# Find OBJs not yet converted to SWC
swc.files <- gsub(pattern = "\\.swc$", "", list.files(banc.swc.save.path, pattern = ".swc$", full.names = FALSE))
# Clean up path artifacts in filenames
swc.files <- gsub("datahmneurobiwilsobanob", "", swc.files)
todo.files <- sample(obj.full.files[!obj.ids %in% swc.files])
message(sprintf("There are %d of %d .obj files not converted to accurate skeletons",
                length(todo.files), length(obj.full.files)))

# Check for previously failed skeletonisations
failed <- tryCatch(readr::read_csv(file.path(banc.save.path, "banc_failed.csv"),
                                    col_types = readr::cols("c", "c")),
                   error = function(e) NULL)
if (!is.null(failed)) {
  fails <- as.character(subset(failed, file.missing == "swc")$root_id)
  todo.ids <- gsub("\\.obj", "", basename(todo.files))
  todo.good.files <- todo.files[!todo.ids %in% fails]
} else {
  todo.good.files <- todo.files
  fails <- character(0)
}

obj.not.failed.before <- length(todo.good.files)
message(sprintf("There are %d .obj that have not failed skeletonisation before", obj.not.failed.before))

# On Saturdays, retry past failures
if (weekdays(Sys.time()) == "Saturday") {
  message(sprintf("Saturday: retrying %d past failures", length(fails)))
  todo.good.files <- todo.files
}

if (length(todo.good.files) == 0) {
  message("No skeletonisation needed. Nothing to do.")
  return(invisible())
}

###########################
### Skeletonise         ###
###########################

message(sprintf("Attempting to skeletonise %d .obj files...", length(todo.good.files)))

multiplier <- 1000
upper <- ifelse((numCores * multiplier) < length(todo.good.files),
                numCores * multiplier, length(todo.good.files))
batches <- split(sample(todo.good.files),
                 round(seq(from = 1, to = upper, length.out = length(todo.good.files))))

# Register parallel backend
cl <- setup_parallel()

# Progress bar
iterations <- length(batches)
pb <- utils::txtProgressBar(max = iterations, style = 3)
progress <- function(n) utils::setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Run skeletonisation
by.query <- foreach::foreach(batch = seq_along(batches),
                             .combine = 'c',
                             .errorhandling = 'pass',
                             .options.snow = opts) %dopar% {
  neuron.ids <- batches[[batch]]
  for (neuron.id in neuron.ids) {
    tryCatch({
      skels <- fafbseg::skeletor(obj = neuron.id,
                                  clean = TRUE,
                                  method = "wavefront",
                                  save.obj = NULL,
                                  mesh3d = FALSE,
                                  waves = 1,
                                  k.soma.search = 100,
                                  radius.soma.search = 10000,
                                  heal = TRUE,
                                  reroot = TRUE,
                                  heal.threshold = 1000000,
                                  heal.k = 10L,
                                  reroot_method = "density",
                                  brain = bancr::banc_neuropil.surf,
                                  elapsed = 10000,
                                  resample = 1000)
      skels[, "id"] <- names(skels) <- basename(gsub("\\.obj|\\/datahmneurobiwilsobanob", "", names(skels)))
      nat::write.neurons(skels, dir = banc.swc.save.path, format = 'swc', Force = TRUE)
    }, error = function(e) {
      message(sprintf("  Skeletonisation error for %s: %s", basename(neuron.id), e$message))
    })
  }
  NULL
}

stop_parallel(cl)

# Report errors
for (i in seq_along(by.query)) {
  if (!is.null(by.query[[i]])) message(by.query[[i]])
}

###########################
### Track results       ###
###########################

# Identify newly created SWC files
swc.files.new <- gsub(pattern = "\\.swc$", "", list.files(banc.swc.save.path, pattern = ".swc$", full.names = FALSE))
swc.files.new <- gsub("datahmneurobiwilsobanob", "", swc.files.new)
swc.files.new <- setdiff(swc.files.new, swc.files)
new.swc <- length(swc.files.new)

# Record which IDs failed
successful <- gsub("\\.swc", "", list.files(banc.swc.save.path, pattern = ".swc$", full.names = FALSE))
successful <- gsub("datahmneurobiwilsobanob", "", successful)
fw.ids.tried <- gsub(".obj", "", basename(todo.good.files))

new_failed <- data.frame(stringsAsFactors = FALSE)
fail.obj <- setdiff(fw.ids.tried, obj.ids)
if (length(fail.obj)) {
  new_failed <- rbind(new_failed, data.frame(root_id = fail.obj, file.missing = "obj"))
}
fail.swc <- setdiff(fw.ids.tried, successful)
if (length(fail.swc)) {
  new_failed <- rbind(new_failed, data.frame(root_id = fail.swc, file.missing = "swc"))
}
if (nrow(new_failed)) {
  new_failed$root_id <- as.character(new_failed$root_id)
  write.csv(new_failed, file.path(banc.save.path, "banc_failed.csv"), row.names = FALSE)
}

# Clean up path artifacts in filenames
files_to_rename <- list.files(banc.swc.save.path, pattern = ".swc$", full.names = TRUE)
new_names <- gsub("datahmneurobiwilsobanob", "", files_to_rename)
file.rename(from = files_to_rename, to = new_names)

message(sprintf("### banc: skeletonised %d .obj files out of %d attempted [%s] ###",
                new.swc, length(todo.good.files),
                format(round(difftime(Sys.time(), t_start, units = "mins"), 1))))

})
