#' banc-calculate-split — Flow-centrality axon/dendrite splits per proofread neuron.
#'
#' Attaches synapses from the versioned parquet to the skeleton (detailed
#' preferred, L2 fallback) and runs `hemibrainr::flow_centrality`. `Label`
#' ∈ {0=unknown, 1=soma, 2=axon, 3=dendrite, 4=primary.dendrite, 7=primary.neurite}.
#'
#' @section Reads:
#'   - `banc_<ver>_synapses_<src>.parquet`, detailed + L2 SWCs
#'
#' @section Writes (under `<banc.l2split.save.path>`):
#'   - `swc/<root_id>.swc`, `synapses/<root_id>.csv`, `metrics/banc_metrics_N.csv`
#'
#' @section Invoked by:
#'   production v888 rebuild chain.
#'
#' @section Used by:
#'   `banc/share/banc-data.R` (split edgelist + enriched synapses).
#'
#' @section Paper:
#'   Methods §"Flow centrality / axon-dendrite split".

###########################################################
### Calculate BANC flow centrality splits (v888)
###
### For each proofread root_888 neuron with a skeleton
### (L2 or detailed), splits into axon/dendrite
### compartments using flow_centrality.
###
### Synapses are extracted from the v888 synapse table
### (not per-neuron CAVE queries), then attached to
### skeletons for compartment labelling.
###
### Prefers detailed skeletons (from banc.swc.save.path)
### when available; falls back to L2 skeletons
### (from banc.l2swc.save.path).
###
### Outputs (labelled by compartment):
###   - SWC files:    l2split/swc/ (with Label column)
###   - Synapse CSVs: l2split/synapses/ (with Label column)
###   - Metrics CSVs: l2split/metrics/
###   - Images:       l2split/images/
###########################################################
source("banc/banc-startup.R")
library(data.table)

local({

message("### banc: calculating flow centrality splits ###")
t_start <- Sys.time()

# CLI flags
.cli_args <- commandArgs(trailingOnly = TRUE)
skip_images <- "--skip-images" %in% .cli_args
if (skip_images) message("  --skip-images: image generation disabled")

#######################
### CONFIGURATION   ###
#######################

# banc.version set in banc-startup.R
versioned.save.path <- file.path("/n/data1/hms/neurobio/wilson/connectomes/banc",
                                  paste0("banc_", banc.version))
synapse_size_threshold <- 5

# Output folders
split.master.folder <- banc.l2split.save.path
split.folder <- file.path(split.master.folder, "swc")
images.folder <- file.path(split.master.folder, "images")
synapses.folder <- file.path(split.master.folder, "synapses")
metrics.folder <- file.path(split.master.folder, "metrics")
dir.create(split.folder, recursive = TRUE, showWarnings = FALSE)
dir.create(images.folder, recursive = TRUE, showWarnings = FALSE)
dir.create(synapses.folder, recursive = TRUE, showWarnings = FALSE)
dir.create(metrics.folder, recursive = TRUE, showWarnings = FALSE)

###########################
### Read current state  ###
###########################

bc <- banctable_query()
banc.meta <- bc %>%
  banc_filter_neurons()
banc.root.ids <- unique(banc.meta$root_id)
# Honor test mode if active
if (exists("banc.test.ids", envir = .GlobalEnv))
  banc.root.ids <- intersect(banc.root.ids, banc.test.ids)
sensories <- subset(banc.meta, grepl("sensory|afferent", super_class) |
                      grepl("sensory|afferent", cell_class))$root_id
motors <- subset(banc.meta, grepl("motor|efferent|endocrine|visceral", super_class) |
                   grepl("motor|efferent|endocrine|visceral", cell_class))$root_id

# Find available skeletons: prefer detailed (SWC) over L2
# For each neuron, use the detailed skeleton if it exists, otherwise L2
l2_swcs <- list.files(banc.l2swc.save.path, pattern = "\\.swc$", full.names = TRUE)
l2_ids <- gsub("\\.swc$", "", basename(l2_swcs))
detailed_swcs <- list.files(banc.swc.save.path, pattern = "\\.swc$", full.names = TRUE)
detailed_ids <- gsub("\\.swc$", "", basename(detailed_swcs))
# Clean path artifacts
detailed_ids <- gsub("datahmneurobiwilsobanob", "", detailed_ids)

# Build skeleton file lookup: prefer detailed over L2
skeleton_lookup <- data.frame(
  root_id = unique(c(detailed_ids, l2_ids)),
  stringsAsFactors = FALSE
)
skeleton_lookup$swc_path <- sapply(skeleton_lookup$root_id, function(rid) {
  if (rid %in% detailed_ids) {
    idx <- which(detailed_ids == rid)[1]
    detailed_swcs[idx]
  } else {
    idx <- which(l2_ids == rid)[1]
    l2_swcs[idx]
  }
})
skeleton_lookup$source <- ifelse(skeleton_lookup$root_id %in% detailed_ids, "detailed", "l2")

# Filter to current root_888 neurons
skeleton_lookup <- skeleton_lookup %>%
  dplyr::filter(root_id %in% banc.root.ids)

# Skip already-split neurons (both SWC and CSV must exist)
done.swc <- gsub("\\.swc$", "", list.files(split.folder))
done.csv <- gsub("\\.csv$", "", list.files(synapses.folder))
done.ids <- intersect(done.swc, done.csv)

skeleton_lookup <- skeleton_lookup %>%
  dplyr::filter(!root_id %in% done.ids)

swcs <- skeleton_lookup$swc_path
names(swcs) <- skeleton_lookup$root_id
swcs <- sample(swcs)

n_detailed <- sum(skeleton_lookup$source == "detailed")
n_l2 <- sum(skeleton_lookup$source == "l2")
message(sprintf("Attempting to split %d BANC neurons (%d detailed, %d L2), already done: %d",
                length(swcs), n_detailed, n_l2, length(done.ids)))

if (length(swcs) == 0) {
  message("No neurons to split. Nothing to do.")
  return(invisible())
}

###################################
### Load v888 synapse table     ###
###################################

message(sprintf("Loading v%s synapse table...", banc.version))

syn_candidates <- c(
  file.path(banc.connectivity.save.path,
            sprintf("banc_%s_synapses.parquet", banc.version)),
  file.path(versioned.save.path,
            sprintf("banc_%s_synapses_enriched.parquet", banc.version)),
  file.path(banc.save.path, paste0("v", banc.version),
            "synapses_v2_human_readable.csv"),
  bancsynapses
)
syn_path <- NULL
for (sp in syn_candidates) {
  if (file.exists(sp)) { syn_path <- sp; break }
}
if (is.null(syn_path))
  stop("Could not find v", banc.version, " synapse table. Checked:\n  ",
       paste(syn_candidates, collapse = "\n  "))
message(sprintf("  Using: %s", syn_path))

if (grepl("\\.parquet$", syn_path)) {
  all_syns <- as.data.table(arrow::read_parquet(syn_path,
    col_select = c("id", "pre_root_id", "post_root_id", "X", "Y", "Z", "size")))
} else {
  # Raw CSV: 15 columns, skip header, select by position
  all_syns <- fread(syn_path, skip = 1L, header = FALSE,
                    select = c(1L, 8L, 9L, 10L, 11L, 13L, 15L))
  setnames(all_syns, c("id", "X", "Y", "Z", "size", "pre_root_id", "post_root_id"))
}

all_syns[, `:=`(pre_root_id = as.character(pre_root_id),
                post_root_id = as.character(post_root_id),
                id = as.character(id))]

# Filter: size threshold, no autapses, at least one end in our neuron set
all_neuron_ids <- names(swcs)
all_syns <- all_syns[
  size > synapse_size_threshold &
  pre_root_id != post_root_id &
  (pre_root_id %in% all_neuron_ids | post_root_id %in% all_neuron_ids)
]
message(sprintf("  Filtered: %s synapses for %d neurons",
                format(nrow(all_syns), big.mark = ","), length(all_neuron_ids)))

# Read root positions for rerooting
banc.roots <- tryCatch(
  arrow::read_feather(file.path(banc.save.path, "banc_root_positions.feather")) %>%
    as.data.frame() %>%
    dplyr::distinct(root_id, root_position, root_position_nm),
  error = function(e) data.frame()
)
banc.confirmed.roots <- bc %>%
  dplyr::distinct(root_id, root_position, root_position_nm) %>%
  rbind(banc.roots) %>%
  dplyr::distinct(root_id, root_position, root_position_nm) %>%
  dplyr::filter(!is.na(root_position), !is.na(root_position_nm),
                !grepl("NA", root_position_nm), !grepl("NA", root_position))

###################################
### Pre-extract per-batch syns  ###
###################################

# Batch processing
multiplier <- 10
upper <- ifelse((numCores * multiplier) < length(swcs), numCores * multiplier, length(swcs))
batches <- split(swcs, round(seq(from = 1, to = upper, length.out = length(swcs))))

# Pre-extract synapses per batch to avoid exporting full table to workers
message("Pre-extracting synapses per batch...")
batch_synapse_list <- lapply(batches, function(batch_swcs) {
  batch_ids <- names(batch_swcs)
  all_syns[pre_root_id %in% batch_ids | post_root_id %in% batch_ids]
})

# Free the full table
rm(all_syns)
gc()

###########################
### Split neurons       ###
###########################

# Register parallel backend
cl <- setup_parallel()

# Progress
iterations <- length(batches)
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

by.query <- foreach::foreach(batch = seq_along(batches),
                             batch_syns = batch_synapse_list,
                             .verbose = numCores > 1,
                             .combine = 'c',
                             .errorhandling = 'pass',
                             .export = c("skip_images"),
                             .options.snow = opts) %dopar% {

  swc <- batches[[batch]]

  # Skip already processed in this run
  target1 <- file.path(split.folder, paste0(names(swc), ".swc"))
  target2 <- file.path(synapses.folder, paste0(names(swc), ".csv"))
  exists <- sapply(target1, file.exists) & sapply(target2, file.exists)
  swc <- swc[!exists]
  if (!length(swc)) return(NULL)

  ids <- names(swc)
  message("running on: ", paste(ids, collapse = ", "))

  # Read neurons
  neurons <- banc_read_swc(unname(swc), ids = ids, meta = as.data.frame(bc))

  # Resample to 100nm
  neurons <- nat:::resample.neuronlist(neurons, 100, OmitFailures = TRUE)

  # Reroot at soma (use root_888 IDs, no latestid needed for versioned data).
  # nat::read.neurons does not set $id on individual neurons — only the
  # neuronlist's `df` attribute — so banc_reroot.neuron would throw
  # "a root_id in roots must be given" on every neuron. The neuronlist method
  # handles this via add_field_seq before dispatching; we replicate that here
  # and keep OmitFailures=TRUE so one bad neuron doesn't abort the batch.
  neurons <- bancr:::add_field_seq(neurons, entries = names(neurons), field = "id")
  bc.neurons.rerooted <- nlapply(neurons, banc_reroot,
                                  roots = banc.confirmed.roots,
                                  estimate = FALSE, .progress = "none",
                                  OmitFailures = TRUE)

  # Add synapses from pre-loaded table (not CAVE)
  # Build per-neuron connectors from the batch synapse subset
  bc.neurons.syns <- bc.neurons.rerooted
  for (nid in names(bc.neurons.syns)) {
    n_syns <- batch_syns[pre_root_id == nid | post_root_id == nid]
    if (nrow(n_syns) == 0) {
      bc.neurons.syns[[nid]]$connectors <- data.frame()
      next
    }
    connectors <- data.frame(
      connector_id = n_syns$id,
      pre_id = n_syns$pre_root_id,
      post_id = n_syns$post_root_id,
      prepost = ifelse(n_syns$post_root_id == nid, 1L, 0L),
      size = n_syns$size,
      X = n_syns$X, Y = n_syns$Y, Z = n_syns$Z,
      stringsAsFactors = FALSE
    )
    bc.neurons.syns[[nid]] <- tryCatch(
      banc_add_synapses(bc.neurons.syns[[nid]],
                         connectors = connectors,
                         id = nid, update.id = FALSE),
      error = function(e) {
        message("  banc_add_synapses failed for ", nid, ": ", e$message)
        bc.neurons.syns[[nid]]
      }
    )
  }

  # Drop neurons with no synapses
  good.ids <- c()
  for (i in seq_along(bc.neurons.syns)) {
    if (!nrow(bc.neurons.syns[[i]]$connectors) || !length(bc.neurons.syns[[i]]$connectors)) {
      message("no synapses for: ", names(bc.neurons.syns)[i])
    } else {
      good.ids <- c(good.ids, i)
    }
  }
  bc.neurons.syns <- bc.neurons.syns[good.ids]
  if (!length(bc.neurons.syns)) return(NULL)

  # Remove bad synapses (pre-split)
  bc.neurons.syns.rem <- bc.neurons.syns
  somas <- sapply(bc.neurons.syns, function(x) !is.null(x$tags$soma))
  for (i in seq_along(bc.neurons.syns)) {
    rem <- remove_bad_synapses(bc.neurons.syns[i], meshes = NULL,
                                soma = somas[i],
                                min.nodes.from.soma = 150,
                                min.nodes.from.pnt = 10,
                                primary.branchpoint = 0.25,
                                OmitFailures = TRUE,
                                .parallel = FALSE, wipe = TRUE)
    bc.neurons.syns.rem[[i]] <- rem[[1]]
  }

  # Flow centrality split
  bc.neurons.flow <- hemibrainr::flow_centrality(bc.neurons.syns.rem,
                                                  mode = mode, polypre = polypre,
                                                  split = split, .parallel = FALSE,
                                                  OmitFailures = FALSE)

  # Remove bad synapses (post-split)
  bc.neurons.flow.rem <- bc.neurons.flow
  for (i in seq_along(bc.neurons.flow)) {
    bc.neuron.flow.rem <- remove_bad_synapses(bc.neurons.flow[i], meshes = NULL,
                                               soma = TRUE,
                                               min.nodes.from.soma = 50,
                                               min.nodes.from.pnt = 5,
                                               primary.branchpoint = 0.25,
                                               OmitFailures = TRUE,
                                               .parallel = FALSE, wipe = TRUE)
    tryCatch({
      bc.neurons.flow.rem[[i]] <- bc.neuron.flow.rem[[1]]
    }, error = function(e) {
      bc.neurons.flow.rem[[i]] <<- bc.neurons.flow[[i]]
    })
  }
  bc.neurons.flow <- bc.neurons.flow.rem
  bc.neurons.flow <- bc.neurons.flow[sapply(bc.neurons.flow, function(x) is.neuron(x))]

  # Assert labels for sensory/motor neurons
  ids.sensories <- intersect(names(bc.neurons.flow), sensories)
  for (is in ids.sensories) {
    bc.neurons.flow[[is]] <- hemibrainr:::add_Label(bc.neurons.flow[[is]], Label = 2)
  }
  ids.motors <- intersect(names(bc.neurons.flow), motors)
  for (is in ids.motors) {
    bc.neurons.flow[[is]] <- hemibrainr:::add_Label(bc.neurons.flow[[is]], Label = 3)
  }

  # Fix unlabeled nodes by nearest labeled node
  for (id in names(bc.neurons.flow)) {
    neuron <- bc.neurons.flow[[id]]
    syns.df <- neuron$connectors
    d <- neuron$d
    g <- nat::as.ngraph(neuron)
    labeled_nodes <- which(d$Label != 0)
    unlabeled_nodes <- which(d$Label == 0)
    if (length(unlabeled_nodes) && length(labeled_nodes) > 0) {
      dm <- igraph::distances(g, v = unlabeled_nodes, to = labeled_nodes)
      closest_labeled_idx <- apply(dm, 1, which.min)
      d$Label[unlabeled_nodes] <- d$Label[labeled_nodes[closest_labeled_idx]]
      neuron$d <- d
      if (nrow(syns.df)) {
        node_coords <- as.matrix(neuron$d[, c("X", "Y", "Z")])
        syn_coords <- as.matrix(syns.df[, c("x", "y", "z")])
        nearest_node_idx <- RANN::nn2(node_coords, syn_coords, k = 1)$nn.idx[, 1]
        syns.df$Label <- neuron$d$Label[nearest_node_idx]
        neuron$connectors <- syns.df
      }
    }
    bc.neurons.flow[[id]] <- neuron
  }

  # Assign Strahler order
  bc.neurons.split <- assign_strahler(bc.neurons.flow, .parallel = FALSE, OmitFailures = TRUE)

  # Carry over labels to connectors
  bc.neurons.split <- nat::nlapply(bc.neurons.split, hemibrainr:::carryover_labels,
                                    .parallel = FALSE, OmitFailures = TRUE, .progress = "none")
  bc.neurons.split[, "root_id"] <- names(bc.neurons.split)

  # Save labelled SWC files
  message("writing data: ")
  write.it <- capture.output(
    write.neurons(nl = bc.neurons.split, dir = split.folder,
                  files = names(bc.neurons.split), Force = TRUE,
                  format = "swc", include.data.frame = FALSE)
  )

  # Save labelled synapse CSVs (data.table::fwrite for speed)
  for (n in names(bc.neurons.split)) {
    csv <- bc.neurons.split[n][[1]]$connectors
    data.table::fwrite(x = csv, file = file.path(synapses.folder, paste0(n, ".csv")))
  }

  # Calculate segregation index
  bc.neurons.batch <- hemibrainr:::segregation_index.neuronlist(bc.neurons.split)

  # Scale to microns for metrics
  bc.neurons.split.microns <- hemibrainr::scale_neurons(bc.neurons.batch,
                                                         scaling = 1 / 1000,
                                                         .parallel = FALSE,
                                                         OmitFailures = TRUE)
  bc.neurons.split.microns[, ] <- bc.neurons.split.microns[, c("root_id", "nucleus_id", "root_position_nm")]
  mets <- hemibrain_compartment_metrics(bc.neurons.split.microns,
                                         OmitFailures = TRUE, .parallel = FALSE,
                                         delta = 5, resample = NULL, locality = FALSE)
  colnames(mets) <- snakecase::to_snake_case(colnames(mets))
  bc.metrics <- mets
  bc.metrics$root_id <- bc.metrics$id
  bc.metrics$id <- NULL
  bc.metrics$projection_score <- NULL
  bc.metrics <- bc.metrics[!duplicated(bc.metrics$root_id), ]

  # Projection scores
  proj.scores <- hemibrainr:::projection_score.neuronlist(bc.neurons.split.microns,
                                                           .parallel = FALSE,
                                                           OmitFailures = TRUE)
  colnames(proj.scores) <- c("root_id", "projection_score")
  bc.metrics <- dplyr::left_join(bc.metrics, proj.scores, by = "root_id")

  # Summary of skeleton structure
  summaries <- summary(bc.neurons.split.microns)[, c("nodes", "cable.length")]
  colnames(summaries) <- snakecase::to_snake_case(colnames(summaries))
  summaries$root_id <- rownames(summaries)
  summaries$nsoma <- NULL
  keep <- c("root_id", setdiff(colnames(bc.metrics), colnames(summaries)))
  bc.metrics <- dplyr::full_join(bc.metrics[, keep], summaries, by = "root_id")
  colnames(bc.metrics) <- snakecase::to_snake_case(colnames(bc.metrics))

  # Round and save metrics
  bc.metrics.sq <- round_dataframe(bc.metrics, digits = 4)
  message(sprintf("metrics calculated for %d neurons", nrow(bc.metrics.sq)))
  nmets <- length(list.files(metrics.folder)) + 1
  data.table::fwrite(bc.metrics.sq,
                     file = file.path(metrics.folder, sprintf("banc_metrics_%d.csv", nmets)))

  # Save comparison images. Skipped entirely with --skip-images, and per-neuron
  # if an image already exists. This is the per-neuron bottleneck:
  # banc_neuron_comparison_plot can take hours on large neurons and previously
  # stalled the entire batch.
  if (!skip_images) {
    for (id in names(bc.neurons.split)) {
      img_path <- file.path(images.folder, paste0(id, ".png"))
      if (file.exists(img_path)) next
      try({
        img <- banc_neuron_comparison_plot(
          filename = img_path,
          neuron1 = bc.neurons.split[id],
          neuron1.info = id)
      })
    }
  }

  if (length(bc.neurons.split) == 1) {
    message("completed: ", names(bc.neurons.split))
  } else {
    message("completed batch: ", batch, ": neuron count: ", length(bc.neurons.split))
  }
}

stop_parallel(cl)

# Report errors
for (i in seq_along(by.query)) {
  if (!is.null(by.query[[i]])) message(by.query[[i]])
}

message(sprintf("### banc: flow centrality splitting complete [%s] ###",
                format(round(difftime(Sys.time(), t_start, units = "mins"), 1))))

})
