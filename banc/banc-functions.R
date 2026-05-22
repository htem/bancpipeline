#' banc-functions — Shared helper functions for the BANC pipeline.
#'
#' Sourced by `banc/banc-startup.R`. Defines the per-script helpers used
#' across the pipeline: neuron filtering (`filter_valid_neurons`), SeaTable
#' I/O wrappers + retries, status-column append helpers, side-index and
#' laterality helpers, plot palettes, edgelist construction, and the
#' versioned-data write helpers (`write_connectome_data`). No top-level
#' side effects — pure function definitions.
#'
#' @section Notes:
#'   - Not invoked directly; always loaded via `banc/banc-startup.R`.

filter_valid_neurons <- function(df,
                                 only_proofread = TRUE,
                                 deduplicate = TRUE) {
  # Exclude non-neuronal super_classes
  if ("super_class" %in% colnames(df)) {
    df <- df %>%
      dplyr::filter(!grepl("glia|trachea|not_a_neuron|merge|debris",
                           super_class, ignore.case = TRUE))
  }
  # Exclude non-neuronal status values
  if ("status" %in% colnames(df)) {
    df <- df %>%
      dplyr::filter(!grepl("GLIA|TRACHEA|NOT_A_NEURON|DEBRIS|MERGE|DELETE",
                           status))
  }
  # Proofread filter
  if (only_proofread) {
    if ("proofread" %in% colnames(df) && "roughly_proofread" %in% colnames(df)) {
      df <- df %>%
        dplyr::filter(as.logical(proofread) %in% TRUE |
                        as.logical(roughly_proofread) %in% TRUE)
    } else if ("proofread" %in% colnames(df)) {
      df <- df %>%
        dplyr::filter(as.logical(proofread) %in% TRUE)
    }
  }
  # Deduplicate by root_id, keeping the row with the most non-NA values
  if (deduplicate) {
    id_col <- if ("root_id" %in% colnames(df)) "root_id"
    else if ("id" %in% colnames(df)) "id"
    else NULL
    if (!is.null(id_col)) {
      df <- df %>%
        dplyr::mutate(.n_nonna = rowSums(!is.na(dplyr::pick(dplyr::everything())))) %>%
        dplyr::arrange(dplyr::desc(.n_nonna)) %>%
        dplyr::distinct(!!rlang::sym(id_col), .keep_all = TRUE) %>%
        dplyr::select(-.n_nonna)
    }
  }
  df
}

# Write a data frame as feather or parquet
write_connectome_data <- function(data, path, format = c("feather", "parquet")) {
  format <- match.arg(format)
  if (format == "parquet") {
    arrow::write_parquet(data, path,
                         version = "2.6",
                         compression = "snappy",
                         chunk_size = 100000,
                         use_dictionary = TRUE)
  } else {
    arrow::write_feather(data, path)
  }
}

# function
calculate_influence_norms <- function(influence.df,
                                      const = -24){
  inf.threshold <- exp(const)
  if(!"target"%in%colnames(influence.df)){
    influence.df$target <- influence.df$id
    orig.target = FALSE
  }else{
    orig.target = TRUE
  }
  if(!"influence_syn_original"%in%colnames(influence.df)){
    influence.df$influence_syn_original <- 1
  }
  if(!"influence_original"%in%colnames(influence.df)){
    influence.df$influence_original <- influence.df$influence
  }
  if(!"influence_norm_original"%in%colnames(influence.df)){
    influence.df$influence_norm_original <- influence.df$influence_norm
  }
  influence.df <- influence.df %>%
    dplyr::ungroup() %>%
    dplyr::mutate(no_seeds = influence_original/influence_norm_original,
                  no_seeds = ifelse(is.na(no_seeds),1,no_seeds),
                  no_synapses = influence_original/influence_syn_norm) %>%
    dplyr::group_by(target) %>%
    dplyr::mutate(no_targets = length(unique(id))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      no_seeds = as.numeric(no_seeds),
      no_targets = as.numeric(no_targets)
    ) %>%
    dplyr::group_by(seed) %>%
    dplyr::mutate(influence_per_seed = influence_original*no_seeds,
                  influence_per_synapse = influence_original*no_synapses) %>%
    dplyr::group_by(target, seed) %>%
    dplyr::mutate(influence = sum(influence_original,na.rm = TRUE),
                  total_seeds = sum(no_seeds,na.rm=TRUE),
                  total_synapses = sum(no_synapses,na.rm=TRUE),
                  influence_norm = sum(influence_per_seed,na.rm = TRUE)/(total_seeds*no_targets),
                  influence = ifelse(influence<inf.threshold,inf.threshold,influence),
                  influence_norm = ifelse(influence_norm<inf.threshold,inf.threshold,influence_norm),
                  influence_norm = sum(influence_norm,na.rm = TRUE),
                  influence_syn_norm =  sum(influence_per_synapse,na.rm = TRUE)/(no_targets*total_synapses),
                  influence_syn_norm = ifelse(influence_syn_norm<inf.threshold,inf.threshold,influence_syn_norm),
                  influence_syn_norm = sum(influence_syn_norm,na.rm = TRUE),
                  influence_norm_log = log(influence_norm),
                  influence_log = log((influence/no_targets)),
                  influence_syn_norm_log = log(influence_syn_norm)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(influence_norm_log = influence_norm_log-const,
                  influence_log = influence_log-const,
                  influence_syn_norm_log = influence_syn_norm_log-const) %>%
    dplyr::group_by(seed) %>%
    dplyr::mutate(influence_log = ifelse(is.na(influence),0,influence_log)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(target, seed, .keep_all = TRUE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(influence = signif(influence,4),
                  influence_log = signif(influence_log,4),
                  influence_norm = signif(influence_norm,4),
                  influence_syn_norm = signif(influence_syn_norm,4),
                  influence_norm_log = signif(influence_norm_log,4),
                  influence_syn_norm_log = signif(influence_syn_norm_log,4)
    ) 
  if(!orig.target){
    influence.df$target <- NULL
  }
  influence.df
}

lighten_color <- function(color, factor = 1.4) {
  col_rgb <- col2rgb(color)
  col_rgb <- pmin(col_rgb * factor, 255)
  rgb(t(col_rgb), maxColorValue = 255)
}

extract_root_id <- function(string) {
  pattern <- "root_id_(\\d+)_"
  match <- str_match(string, pattern)
  if (!is.na(match[1,2])) {
    return(bit64::as.integer64(match[1,2]))
  } else {
    return(NA)
  }
}

# hidden
banc_cave_cell_types <- bancr:::banc_cave_cell_types

banc_nucelus_id_to_rootid <- function(nucleus_ids){
  nuclei <- banc_nuclei()
  nuclei$root_id <- as.character(nuclei$root_id)
  nuclei$nucleus_id <- as.character(nuclei$nucleus_id)
  ids <- nuclei$root_id[match(nucleus_ids, nuclei$nucleus_id)]
  ids[is.na(ids)] <- "0"
  ids
}

banc_updateids <- bancr:::banc_updateids

# Pull banc_meta from SeaTable with a guaranteed root_<version> column.
# Fallback chain when SeaTable is missing the column (e.g. immediately after
# a materialization bump, before the column is provisioned) or unreachable:
#   1. If `sql` provided AND the query errors, retry without the root_<ver> term.
#   2. If the result has no root_<ver> column, resolve it from supervoxel_id
#      via CAVE banc_rootid(version = <ver>) in-memory.
#   3. If the initial query fails AND `cache_feather` is provided, read that
#      feather instead and continue with step 2.
banc_meta_with_root_version <- function(version = banc.version,
                                        sql = NULL,
                                        cache_feather = NULL) {
  root_col <- paste0("root_", version)
  .strip_root_col <- function(s) {
    s <- gsub(sprintf(",\\s*%s\\b", root_col), "", s)
    gsub(sprintf("\\b%s\\s*,\\s*", root_col), "", s)
  }
  .run_query <- function(q) {
    if (is.null(q)) banctable_query() else banctable_query(q)
  }

  bc <- tryCatch(
    .run_query(sql),
    error = function(e) {
      msg <- conditionMessage(e)
      message(sprintf("banctable_query error: %s", msg))
      if (!is.null(sql) && grepl(root_col, msg)) {
        message(sprintf("  retrying without %s in SELECT", root_col))
        return(tryCatch(.run_query(.strip_root_col(sql)),
                        error = function(e2) {
                          message(sprintf("  retry also failed: %s",
                                          conditionMessage(e2)))
                          if (!is.null(cache_feather) && file.exists(cache_feather)) {
                            message(sprintf("  reading cache feather: %s",
                                            cache_feather))
                            return(arrow::read_feather(cache_feather))
                          }
                          stop(e2)
                        }))
      }
      if (!is.null(cache_feather) && file.exists(cache_feather)) {
        message(sprintf("  reading cache feather: %s", cache_feather))
        return(arrow::read_feather(cache_feather))
      }
      stop(e)
    }
  )

  if (!root_col %in% names(bc)) {
    if (!"supervoxel_id" %in% names(bc)) {
      stop(sprintf("Cannot backfill %s: supervoxel_id missing from response.",
                   root_col))
    }
    message(sprintf("Backfilling %s via CAVE banc_rootid(version=\"%s\")...",
                    root_col, version))
    svids <- as.character(bc$supervoxel_id)
    has_svid <- !is.na(svids) & nzchar(svids) & svids != "0"
    bc[[root_col]] <- NA_character_
    if (any(has_svid)) {
      bc[[root_col]][has_svid] <- as.character(
        bancr::banc_rootid(svids[has_svid], version = version))
    }
    n_ok <- sum(!is.na(bc[[root_col]]) & bc[[root_col]] != "0")
    message(sprintf("  resolved %s/%s rows",
                    format(n_ok, big.mark = ","),
                    format(nrow(bc), big.mark = ",")))
  }
  bc
}

# hidden
# Vectorized: works correctly when called from dplyr::mutate() on a column
append_status <- function(status, update){
  vapply(status, function(s) {
    s <- paste(c(s, update), collapse = ",")
    s <- paste(sort(unique(unlist(strsplit(s, split = ",|, ")))), collapse = ",")
    gsub("^,| ", "", s)
  }, character(1), USE.NAMES = FALSE)
}

# hidden
# Vectorized: works correctly when called from dplyr::mutate() on a column
subtract_status <- function(status, update, invert = FALSE){
  vapply(status, function(s) {
    statuses <- sort(unique(unlist(strsplit(s, split = ",|, "))))
    if (invert) {
      statuses <- sort(unique(intersect(statuses, update)))
    } else {
      statuses <- sort(unique(setdiff(statuses, update)))
    }
    result <- paste0(statuses, collapse = ",")
    gsub("^,| ", "", result)
  }, character(1), USE.NAMES = FALSE)
}

# Filter out non-neuron entries from a data frame.
# Handles both SeaTable status column and metadata super_class column.
banc_filter_neurons <- function(df) {
  if ("status" %in% names(df)) {
    df <- dplyr::filter(df, !grepl("DELETE|NOT_A_NEURON|DEBRIS|GLIA|TRACHEA|MERGE", status))
  }
  if ("super_class" %in% names(df)) {
    df <- dplyr::filter(df, !super_class %in% c("glia", "trachea", "not_a_neuron", "debris"))
  }
  df
}

# hidden
quiet_function <- function(expr) {
  #temp file
  f = file()
  
  #write output to that file
  sink(file = f)
  
  #evaluate expr in original environment
  y = eval(expr, envir = parent.frame())
  
  #close sink
  sink()
  
  #get rid of file
  close(f)
  y
}

# hidden
banc_nblast_cell_type_images <- function(nblast.files,
                                         dir.images,
                                         query.meta = banctable_query(),
                                         other.meta = query.meta,
                                         query.obj.save.path,
                                         other.obj.save.path = query.obj.save.path,
                                         query.id = "root_id",
                                         other.id = "root_id",
                                         numCores = 1,
                                         redo = TRUE,
                                         region = c("both","brain","vnc"),
                                         max.hits = 5,
                                         volume = NULL,
                                         update = FALSE,
                                         completed.images = NULL
){

  # Register cores
  region <- match.arg(region)
  cl <- setup_parallel()
  
  # Simplify meshes for ggplotting
  banc_neuropil <- Rvcg::vcgQEdecim(as.mesh3d(banc_neuropil.surf), percent = 0.05)
  banc_brain_neuropil <- Rvcg::vcgQEdecim(as.mesh3d(banc_brain_neuropil.surf), percent = 0.05)
  banc_vnc_neuropil <- Rvcg::vcgQEdecim(as.mesh3d(banc_vnc_neuropil.surf), percent = 0.05)
  
  # Create 2D screening files
  message("##### banc_nblast_images_cell_type: creating 2D screening files #####")
  if(!redo){
    if(is.null(completed.images)) completed.images <- basename(list.files(dir.images, recursive = TRUE))
    combos_vec <- unique(as.character(interaction(gsub(sprintf(".*_%s_|_%s_.*",query.id,other.id),"",basename(completed.images)),
                                              gsub("._nblast_score_.*","",basename(completed.images)),
                                              sep="_")))
    # Use hashed environment for O(1) lookup instead of %in% on large vector
    combos <- new.env(hash = TRUE, parent = emptyenv(), size = length(combos_vec))
    for(.cv in combos_vec) combos[[.cv]] <- TRUE
    rm(combos_vec, .cv)
    message("Computed images: ", length(completed.images))
  }else{
    combos <- NULL
  }
  
  # Read nblast matches
  nblast.df <- data.frame()
  if(length(nblast.files)){
    by.query <- foreach::foreach(mfile = nblast.files) %do% {
      id <- gsub(".*_root_id_|\\.csv","",basename(mfile))
      mdf <- readr::read_csv(mfile, col_types = banc.col.types, show_col_types = FALSE)
      numeric_columns <- sapply(mdf, is.numeric)
      mdf[numeric_columns] <- round(mdf[numeric_columns], digits = 2)
      mdf$query <- id
      mdf <- dplyr::arrange(mdf, dplyr::desc(nb))
      mdf$hit <- mdf[[other.id]]
      mdf[1:max.hits,]
    }
    by.query <- by.query[unlist(lapply(by.query,is.data.frame))]
    nb.df <- do.call(plyr::rbind.fill, by.query) %>%
      dplyr::rename(match=hit, root_id=query) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(nb = round(as.numeric(nb),2)) %>%
      dplyr::ungroup() %>%
      dplyr::distinct(root_id, match, nb)
    nblast.df <- rbind(nblast.df, nb.df)
  }
  nblast.df[,other.id] <-nblast.df$match
  nblast.df <- dplyr::left_join(nblast.df, other.meta[,c(other.id,"cell_type")], by = other.id) %>%
    dplyr::arrange(dplyr::desc(nb))
  if(update){
    nblast.df <- banc_updateids(nblast.df)
  }
  matched.ids <- unique(nblast.df$root_id)
  
  # Get cell type matches and pairs
  query.meta.ct <- query.meta %>%
    dplyr::mutate(banc_match = ifelse(is.na(banc_match),"",banc_match)) %>%
    dplyr::select(root_id, banc_match, side, cell_type) %>%
    # Filter away bad/no matches
    # dplyr::filter(!is.na(banc_match), !banc_match%in%c("NA","ISSUE!","ISSUE?")) %>%
    # Arrange by side, so right is first
    dplyr::arrange(dplyr::desc(side)) %>%
    # Create all possible combinations
    dplyr::bind_rows(
      select(., root_id, root_id_matched = banc_match, cell_type),
      select(., root_id = banc_match, root_id_matched = root_id, cell_type)
    ) %>%
    # Ensure each pair appears only once
    dplyr::distinct() %>%
    # Sort root_id and root_id_matched within each row
    dplyr::rowwise() %>%
    dplyr::mutate(
      sorted_id1 = min(root_id, root_id_matched),
      sorted_id2 = max(root_id, root_id_matched)
    ) %>%
    # Remove duplicates after sorting
    dplyr::distinct(sorted_id1, sorted_id2, .keep_all = TRUE) %>%
    # Rename columns to final format
    dplyr::transmute(
      root_id = sorted_id1,
      root_id_matched = sorted_id2,
      cell_type = cell_type
    ) %>%
    # Sort the result
    dplyr::arrange(cell_type, root_id) %>%
    dplyr::filter(root_id%in%matched.ids|root_id_matched%in%matched.ids)
  matched.ids <- unique(c(query.meta.ct$root_id,query.meta.ct$root_id_matched))
  matched.ids <- matched.ids[matched.ids!=""]
  
  # Parallel image process
  by.query <- foreach::foreach(i = 1:nrow(query.meta.ct),
                               .combine = 'c',
                               .errorhandling = 'pass') %do% {

                                 # Early exit check OUTSIDE try() so next works
                                 query.meta.pair <- query.meta.ct[i,]
                                 pair.ids <- c(query.meta.pair$root_id,query.meta.pair$root_id_matched)
                                 pair.ids <- pair.ids[pair.ids!=""]
                                 id <- paste(pair.ids, collapse="_")
                                 message("Working on: ", id)
                                 svid <- query.meta[match(pair.ids[1],query.meta$root_id),"supervoxel_id"][[1]]
                                 if(length(pair.ids)==2){
                                   svid2 <- query.meta[match(pair.ids[2],query.meta$root_id),"supervoxel_id"][[1]]
                                   svid <- paste(unique(c(svid, svid2)), collapse="_")
                                 }

                                 cts.hit <- nblast.df %>%
                                   dplyr::filter(root_id %in% pair.ids) %>%
                                   dplyr::group_by(cell_type) %>%
                                   dplyr::mutate(max_nb_ct = max(nb, na.rm = TRUE),
                                                 mean_nb_ct = mean(nb, na.rm = TRUE)) %>%
                                   dplyr::ungroup() %>%
                                   dplyr::arrange(dplyr::desc(max_nb_ct), dplyr::desc(nb))
                                 cts <- unique(cts.hit$cell_type)

                                 cell.class <- query.meta[match(pair.ids[1], query.meta[[query.id]]),"cell_class"][[1]]
                                 if(is.na(cell.class)){
                                   cell.class <- query.meta[match(pair.ids[1], query.meta[[query.id]]),"super_class"][[1]]
                                   if(is.na(cell.class)){
                                     cell.class = "unknown"
                                   }
                                 }
                                 cell.class <- gsub(" |,","_",cell.class)
                                 cell.class <- gsub("__","_",cell.class)

                                 skip_this <- FALSE
                                 if(!redo && !is.null(combos)){
                                   nclass_guess <- "mesh3d"
                                   iden_guess <- paste0(nclass_guess,"_",gsub("\\.csv","",id),"_supervoxel_id","_",svid)
                                   all_done <- TRUE
                                   for(pos_check in seq_along(cts)){
                                     ct_check <- cts[pos_check]
                                     hit_check <- subset(cts.hit, cts.hit$cell_type==ct_check)[,other.id][[1]]
                                     score_check <- subset(cts.hit, cts.hit$cell_type==ct_check)[,"max_nb_ct"][[1]]
                                     file_check <- sprintf("%s_cell_type_%s_nblast_score_%s_%s_%s_supervoxel_id_%s_hit_%s_%s.png",
                                                           pos_check, ct_check, score_check, query.id, id, svid, other.id, hit_check)
                                     combo_check <- as.character(interaction(gsub(sprintf(".*_%s_|_%s_.*",query.id,other.id),"",basename(file_check)),
                                                                             gsub("._nblast_score_.*","",basename(file_check)),
                                                                             sep="_"))
                                     if(!isTRUE(tryCatch(combos[[combo_check]], error = function(e) FALSE))){
                                       all_done <- FALSE
                                       break
                                     }
                                   }
                                   if(all_done){
                                     message("skipping, all ", length(cts), " cell-type images already exist for: ", id)
                                     skip_this <- TRUE
                                   }
                                 }
                                 if(skip_this) next

                                 try({

                                   # Get BANC meshes (only loaded after early-exit check)
                                   mesh.obj <- file.path(query.obj.save.path,paste0(pair.ids,".obj"))
                                   query.neurons <- nat::as.neuronlist(lapply(mesh.obj,function(x) readobj::read.obj(x, convert.rgl = TRUE)[[1]]))
                                   query.neurons <- nat::nlapply(query.neurons,
                                                                 Rvcg::vcgQEdecim,
                                                                 percent = 0.1)

                                   # Identifier
                                   nclass <- class(query.neurons[[1]])[1]
                                   iden <- paste0(nclass,"_",gsub("\\.csv","",id),"_supervoxel_id","_",svid)

                                   # Make savedir
                                   savedir.id <- file.path(dir.images,cell.class,iden)
                                   dir.create(savedir.id, showWarnings = FALSE, recursive = TRUE)

                                   # Make an image per cell type level hit
                                   for(position in 1:length(cts)){
                                     
                                     # Get the top hit, by which to name the file
                                     ct <- cts[position]
                                     hit <- subset(cts.hit, cts.hit$cell_type==ct)[,other.id][[1]]
                                     score <- subset(cts.hit, cts.hit$cell_type==ct)[,"max_nb_ct"][[1]]
                                     
                                     # Do we have this file?
                                     file <- sprintf("%s_cell_type_%s_nblast_score_%s_%s_%s_supervoxel_id_%s_hit_%s_%s.png", 
                                                     position, ct, score, query.id, id, svid, other.id, hit)
                                     filename <- file.path(savedir.id,file)[1]
                                     combo <- as.character(interaction(gsub(sprintf(".*_%s_|_%s_.*",query.id,other.id),"",basename(filename)),
                                                                       gsub("._nblast_score_.*","",basename(filename)), 
                                                                       sep="_"))
                                     if(!is.null(combos) && isTRUE(tryCatch(combos[[combo]], error = function(e) FALSE))){
                                       if(!redo){
                                         message("skippping already made: ", filename)
                                         next
                                       }
                                     }
                                     
                                     # Get all hits that had an NBLAST score
                                     hits <- subset(cts.hit, cts.hit$cell_type==ct)[[other.id]]
                                     if(length(hits)>3){
                                       hits <- hits[1:3]
                                     }
                                     
                                     # Get hit meshes
                                     mesh.obj <- file.path(other.obj.save.path,paste0(hits,".obj"))
                                     hit.neurons <- nat::as.neuronlist(lapply(mesh.obj,function(x) readobj::read.obj(x, convert.rgl = TRUE)[[1]]))
                                     hit.neurons <- nat::nlapply(hit.neurons,
                                                                 Rvcg::vcgQEdecim,
                                                                 percent = 0.1)
                                     
                                     # Prepare information
                                     cell_class <- query.meta[match(pair.ids[1], query.meta[[query.id]]),"cell_class"][[1]]
                                     cell_type <- query.meta[match(pair.ids[1], query.meta[[query.id]]),"cell_type"][[1]]
                                     hit_class <- other.meta[match(hit, other.meta[[other.id]]),"cell_class"][[1]]
                                     hit_type <- other.meta[match(hit, other.meta[[other.id]]),"cell_type"][[1]]
                                     
                                     # Get info
                                     neuron1.info <- sprintf("%s:%s\ncell_class:%s\ncell_type:%s", query.id, id, cell_class, cell_type)
                                     neuron2.info <- sprintf("%s:%s\ncell_class:%s\ncell_type:%s", other.id, hit, hit_class, hit_type)
                                     
                                     # Plot
                                     result <- banc_neuron_comparison_plot(neuron1 = query.neurons[1],
                                                                           neuron2 = hit.neurons,
                                                                           neuron3 = query.neurons[2],
                                                                           neuron1.info = neuron1.info,
                                                                           neuron2.info = neuron2.info,
                                                                           filename = filename,
                                                                           banc_neuropil = banc_neuropil,
                                                                           banc_brain_neuropil = banc_brain_neuropil,
                                                                           banc_vnc_neuropil = banc_vnc_neuropil,
                                                                           region = region,
                                                                           volume = volume)
                                   }
                                 })
                               }
}


# hidden
banc_nblast_images <- function(nblast.files,
                               dir.images,
                               query.meta,
                               other.meta = query.meta,
                               query.obj.save.path,
                               other.obj.save.path = query.obj.save.path,
                               query.id = "root_id",
                               other.id = "root_id",
                               numCores = 1,
                               redo = TRUE,
                               region = c("auto","both","brain","vnc"),
                               query.mirror = FALSE,
                               by_cell_type = FALSE,
                               max.hits = 5,
                               volume = NULL,
                               completed.images = NULL
){

  # Register cores
  region <- match.arg(region)
  cl <- setup_parallel()
  
  # Simplify meshes for ggplotting
  banc_neuropil <- Rvcg::vcgQEdecim(as.mesh3d(banc_neuropil.surf), percent = 0.05)
  banc_brain_neuropil <- Rvcg::vcgQEdecim(as.mesh3d(banc_brain_neuropil.surf), percent = 0.05)
  banc_vnc_neuropil <- Rvcg::vcgQEdecim(as.mesh3d(banc_vnc_neuropil.surf), percent = 0.05)
  
  # Create 2D screening files
  message("##### banc_nblast_images: creating 2D screening files #####")
  if(!redo){
    if(is.null(completed.images)){completed.images <- basename(list.files(dir.images, recursive = TRUE))}
    combos_vec <- unique(as.character(interaction(gsub(sprintf(".*_%s_|_supervoxel_id_.*|_nucleus_id_.*",query.id),"",basename(completed.images)),
                                         gsub(".*_hit_id_|_hit_supervoxel_id_.*|_hit_nucleus_id_.*","",basename(completed.images)),
                                         sep="_")))
    # Use hashed environment for O(1) lookup instead of %in% on large vector
    combos <- new.env(hash = TRUE, parent = emptyenv(), size = length(combos_vec))
    for(.cv in combos_vec) combos[[.cv]] <- TRUE
    rm(combos_vec, .cv)
    message("computed images: ", length(completed.images))
  }else{
    combos <- NULL
  }
  
  ## Does not seem to like doing this in parallel?
  by.query <- foreach::foreach(nfile = sample(nblast.files),
                               .combine = 'c',
                               .init = list(),
                               .errorhandling = 'pass',
                               .packages = c("readr","nat","bancr","Rvcg"),
                               .export = c("query.meta","other.meta","max.hits", "other.id", "query.id",
                                           "dir.images", "query.mirror", "region", "volume",
                                           "banc_neuropil","banc_brain_neuropil","banc_vnc_neuropil")) %do% {

                                             # Early exit check OUTSIDE try() so next works
                                             id <- gsub(sprintf(".*_%s_|\\.csv",query.id),"",basename(nfile))
                                             message("working on:", id)
                                             skip_this <- FALSE
                                             if(!redo){
                                               svid_check <- query.meta[match(id, query.meta[[query.id]]),"supervoxel_id"][[1]]
                                               nresult_check <- tryCatch(
                                                 readr::read_csv(nfile,
                                                                 col_types = readr::cols(.default = readr::col_character(),
                                                                                         nb = readr::col_number()),
                                                                 show_col_types = FALSE) %>%
                                                   dplyr::arrange(dplyr::desc(nb)),
                                                 error = function(e) NULL)
                                               if (!is.null(nresult_check) && !is.na(svid_check)) {
                                                 nresult_check <- nresult_check[1:max.hits,]
                                                 hits_check <- nresult_check[[other.id]]
                                                 all_done <- TRUE
                                                 for(h in hits_check){
                                                   h.clean <- gsub("^m|m$","",h)
                                                   sc <- round(nresult_check$nb[nresult_check[[other.id]]==h][1],2)
                                                   pos <- match(h, hits_check)
                                                   h_svid <- other.meta[match(h.clean, other.meta[[other.id]]),"supervoxel_id"][[1]]
                                                   ct <- query.meta[match(id, query.meta[[query.id]]),"cell_type"][[1]]
                                                   ht <- other.meta[match(h.clean, other.meta[[other.id]]),"cell_type"][[1]]
                                                   fg <- sprintf("%s_nblast_score_%s_%s_%s_supervoxel_id_%s_hit_id_%s_hit_supervoxel_id_%s_query_cell_type_%s_hit_cell_type_%s.png",
                                                                 pos, sc, query.id, id,
                                                                 nullToNA(svid_check), h, nullToNA(h_svid), ct, ht)
                                                   cc <- as.character(interaction(gsub(sprintf(".*_%s_|_supervoxel_id_.*|_nucleus_id_.*",query.id),"",basename(fg)),
                                                                                  gsub(".*_hit_id_|_hit_supervoxel_id_.*|_hit_nucleus_id_.*","",basename(fg)),
                                                                                  sep="_"))
                                                   if(!isTRUE(tryCatch(combos[[cc]], error = function(e) FALSE))){
                                                     all_done <- FALSE
                                                     break
                                                   }
                                                 }
                                                 if(all_done){
                                                   message("skipping, all ", length(hits_check), " hit images already exist for: ", id)
                                                   skip_this <- TRUE
                                                 }
                                               }
                                             }
                                             if(skip_this) next

                                             try({

                                               # Resolve region per-neuron from query metadata
                                               if(region == "auto" && "region" %in% colnames(query.meta)){
                                                 qreg <- query.meta[match(id, query.meta[[query.id]]), "region"][[1]]
                                                 if(is.na(qreg) || qreg %in% c("", "NA")){
                                                   neuron_region <- "both"
                                                 } else if(grepl("vnc|ventral_nerve_cord", qreg, ignore.case = TRUE)){
                                                   neuron_region <- "vnc"
                                                 } else if(grepl("brain", qreg, ignore.case = TRUE)){
                                                   neuron_region <- "brain"
                                                 } else {
                                                   neuron_region <- "both"
                                                 }
                                               } else if(region == "auto"){
                                                 neuron_region <- "both"
                                               } else {
                                                 neuron_region <- region
                                               }

                                               svid <- query.meta[match(id, query.meta[[query.id]]),"supervoxel_id"][[1]]
                                               if(is.na(svid)){
                                                 stop("ID not in query.meta, BANC ID is likely out of date: ", id)
                                               }

                                               # Get side
                                               side <- query.meta[match(id, query.meta[[query.id]]),"side"][[1]]
                                               if(is.na(side)|side%in%c("NA","", "")){
                                                 side <- "unknown"
                                               }

                                               # Get cell class
                                               cell.class <- query.meta[match(id, query.meta[[query.id]]),"cell_class"][[1]]
                                               if(is.na(cell.class)|cell.class%in%c("NA","", "")){
                                                 cell.class <- query.meta[match(id, query.meta[[query.id]]),"super_class"][[1]]
                                                 if(is.na(cell.class)|cell.class%in%c("NA","", "")){
                                                   cell.class = "unknown"
                                                 }
                                               }
                                               cell.class <- gsub(" |,","_",cell.class)
                                               cell.class <- gsub("__","_",cell.class)

                                               # Get NBLAST result
                                               nresult <- readr::read_csv(nfile,
                                                                          col_types = readr::cols(.default = readr::col_character(),
                                                                                                  nb = readr::col_number())
                                               ) %>%
                                                 dplyr::arrange(dplyr::desc(nb))
                                               nresult <- nresult[1:max.hits,]
                                               hits <- nresult[[other.id]]

                                               # Get query neuron mesh
                                               mesh.obj <- file.path(query.obj.save.path,paste0(id,".obj"))
                                               if(file.exists(mesh.obj)){
                                                 neuron1 <- tryCatch(
                                                   nat::as.neuronlist(readobj::read.obj(mesh.obj, convert.rgl = TRUE)),
                                                   error = function(e) {
                                                     warning("failed to read query mesh for ", id, ": ", e$message)
                                                     NULL
                                                   })
                                               }else{
                                                 warning("trying to read mesh directly from source, looking at BANC")
                                                 neuron1 <- tryCatch(
                                                   bancr::banc_read_neuron_meshes(id),
                                                   error = function(e) {
                                                     warning("mesh download failed for query ", id, ": ", e$message)
                                                     NULL
                                                   })
                                               }
                                               if (is.null(neuron1) || length(neuron1) == 0) {
                                                 warning("no mesh available for query ", id, ", skipping")
                                                 next
                                               }
                                               neuron1 <- Rvcg::vcgQEdecim(neuron1[[1]], percent = 0.1)
                                               
                                               # Identifier
                                               nclass <- class(neuron1)[1]
                                               iden <- paste0(nclass,"_",gsub("\\.csv","",basename(nfile)))
                                               
                                               # Make savedir
                                               savedir.id <- file.path(dir.images,iden) # cell.class,side
                                               dir.create(savedir.id, showWarnings = FALSE, recursive = TRUE)
                                               
                                               # Plot mesh and
                                               for(hit in hits){
                                                 message("----> Visualising on hit: ", hit)
                                                 try({
                                                   # Calculate by cell type?
                                                   if(by_cell_type){
                                                     stop("by_cell_type not yet implemented")
                                                   }else{
                                                     
                                                     # Prepare identifiers
                                                     hit.clean <- gsub("^m|m$","",hit)
                                                     score <- round(nresult$nb[nresult[[other.id]]==hit][1],2)
                                                     position <- match(hit, hits)
                                                     svid <- query.meta[match(id, query.meta[[query.id]]),"supervoxel_id"][[1]]
                                                     hit_svid <- other.meta[match(hit.clean, other.meta[[other.id]]),"supervoxel_id"][[1]]
                                                     
                                                     # Prepare information
                                                     cell_class <- query.meta[match(id, query.meta[[query.id]]),"cell_class"][[1]]
                                                     cell_type <- query.meta[match(id, query.meta[[query.id]]),"cell_type"][[1]]
                                                     hit_class <- other.meta[match(hit.clean, other.meta[[other.id]]),"cell_class"][[1]]
                                                     hit_type <- other.meta[match(hit.clean, other.meta[[other.id]]),"cell_type"][[1]]
                                                     
                                                     # Construct file name
                                                     file <- sprintf("%s_nblast_score_%s_%s_%s_supervoxel_id_%s_hit_id_%s_hit_supervoxel_id_%s_query_cell_type_%s_hit_cell_type_%s.png", 
                                                                     position, score, query.id, id, 
                                                                     nullToNA(svid), hit, nullToNA(hit_svid), cell_type, hit_type)
                                                     filename <- file.path(savedir.id,file)
                                                     combo <- as.character(interaction(gsub(sprintf(".*_%s_|_supervoxel_id_.*|_nucleus_id_.*",query.id),"",basename(filename)),
                                                                                  gsub(".*_hit_id_|_hit_supervoxel_id_.*|_hit_nucleus_id_.*","",basename(filename)),
                                                                                  sep="_"))
                                                     if(!is.null(combos) && isTRUE(tryCatch(combos[[combo]], error = function(e) FALSE)) && !redo){
                                                       message("skippping already made: ", filename)
                                                       next
                                                     }
                                                     
                                                     # Get info
                                                     neuron1.info <- sprintf("%s:%s\ncell_class:%s\ncell_type:%s", query.id, id, cell_class, cell_type)
                                                     neuron2.info <- sprintf("%s:%s\ncell_class:%s\ncell_type:%s", other.id, hit, hit_class, hit_type)
                                                     
                                                     # Get mesh?
                                                     hit.obj <- file.path(other.obj.save.path,paste0(hit.clean,".obj"))
                                                     if(file.exists(hit.obj)){
                                                       neuron2 <- tryCatch(
                                                         nat::as.neuronlist(readobj::read.obj(hit.obj, convert.rgl = TRUE)),
                                                         error = function(e) {
                                                           warning("failed to read hit mesh for ", hit.clean, ": ", e$message)
                                                           NULL
                                                         })
                                                     }else{
                                                       warning("No local mesh for hit ", hit.clean, ", skipping")
                                                       neuron2 <- NULL
                                                     }
                                                     if (is.null(neuron2) || length(neuron2) == 0) {
                                                       warning("no mesh available for ", hit.clean, ", skipping")
                                                       next
                                                     }

                                                     # Simplify  neuron meshes
                                                     neuron2 <- Rvcg::vcgQEdecim(neuron2[[1]], percent = 0.1)
                                                     
                                                     # Mirror
                                                     if(query.mirror){
                                                       neuron2 <- bancr::banc_mirror(neuron2)
                                                     }
                                                     
                                                     # mirror
                                                     if(grepl("^m|m$",hit)){
                                                       message("BANC neuron mirrored")
                                                       neuron3 <- bancr::banc_mirror(neuron2)
                                                       neuron3.info <- "mirrored"
                                                     }else{
                                                       neuron3 <- NULL
                                                       neuron3.info <- NULL
                                                     }
                                                     
                                                     # Plot
                                                     result <- banc_neuron_comparison_plot(neuron1 = neuron1,
                                                                                           neuron2 = neuron2,
                                                                                           neuron3 = neuron3,
                                                                                           neuron1.info = neuron1.info,
                                                                                           neuron2.info = neuron2.info,
                                                                                           filename = filename,
                                                                                           banc_neuropil = banc_neuropil,
                                                                                           banc_brain_neuropil = banc_brain_neuropil,
                                                                                           banc_vnc_neuropil = banc_vnc_neuropil,
                                                                                           region = neuron_region,
                                                                                           volume = volume)
                                                   }
                                                   
                                                   # # check
                                                   # library(png)
                                                   # img <- readPNG(filename)
                                                   # plot(1:2, type='n', xlab="", ylab="", asp=1)
                                                   # rasterImage(img, 1, 1, 2, 2)
                                                   
                                                 })
                                               }
                                             })
                                             NULL           
                                           }
  
  # Stop cores
  stop_parallel(cl)

  # Were there errors?
  # message("##### banc_nblast_images: displaying any errors from foreach loop #####")
  # for(i in 1:length(by.query)){
  #   if(!is.null(by.query[[i]])){
  #     message(by.query[[i]])
  #   }
  # }
  
  # Change permissions
  img.results <- file.path(dir.images,"images")
  imgs <- list.files(img.results, recursive = TRUE, full.names = TRUE)
  for(img in imgs){
    Sys.chmod(img, mode = "0777")
  }
  
  # Return
  invisible()
  
}

plot_png <- function(filename){
  img <- png::readPNG(filename)
  plot(1:2, type='n', xlab="", ylab="", asp=1)
  graphics::rasterImage(img, 1, 1, 2, 2)
  invisible()
}

move_data <- function(source_dir, dest_dir, search = "\\.obj$|\\.csv$|\\.swc$"){
  # Create the new directory if it doesn't exist
  if (!dir.exists(dest_dir)) {
    dir.create(dest_dir)
  }

  # List files matching the search pattern in the source directory
  matched_files <- list.files(source_dir, pattern = search, full.names = TRUE)

  # Move each file to the new directory
  for (file in matched_files) {
    new_path <- file.path(dest_dir, basename(file))
    file.rename(file, new_path)
  }
  invisible()
}


euclidean_distances <-function(A, B) {
  sqrt(rowSums((A - B)^2))
}

# re-root MANC neuron
manc_reroot <- function(x, id = NULL, mc.positions, ...) UseMethod("manc_reroot")
manc_reroot.neuron <- function(x, id = NULL, mc.positions, ...){
  if(is.null(id)){
    id <- x$id
  }
  df <- subset(mc.positions, mc.positions$bodyid==id)
  if(nrow(df)){
    if(!is.na(df$X)){
      soma <- nat::xyzmatrix(df)
      x <- nat::reroot(x = x, point = c(soma))
      x$tags$soma <- nat::rootpoints(x )
    } 
  }
  x
}
manc_reroot.neuronlist <- function(x, id = NULL, mc.positions, ...){
  if(is.null(id)){
    id <- names(x)
  }
  for(i in 1:length(id)){
    x[[i]]$id <- id[i] 
  }
  nat::nlapply(x, FUN = manc_reroot.neuron, mc.positions = mc.positions, id = NULL, ...)
}

# Write deformetrica files
create_dataset_xml <- function(data_folder, 
                               subject_files,
                               output_file, 
                               subject_ids = NULL, 
                               visit_ids = NULL,
                               dataset="") {
  # Create the root element
  root <- xml2::xml_new_root("data-set")
  
  # If subject_ids is not provided, use a default
  if (is.null(subject_ids)) {
    subject_ids <- list(neurons = subject_files)
  }
  
  # Iterate over subjects
  for (subject_id in names(subject_ids)) {
    subject <- xml2::xml_add_child(root, "subject", id = subject_id)
    
    # If visit_ids is not provided, use a default
    if (is.null(visit_ids)) {
      visit_ids <- list(experiment = subject_ids[[subject_id]])
    }
    
    # Iterate over visits
    for (visit_id in names(visit_ids)) {
      visit <- xml2::xml_add_child(subject, "visit", id = visit_id)
      
      # Add filename elements for each file in this visit
      for (file in visit_ids[[visit_id]]) {
        object_id <- tools::file_path_sans_ext(gsub(dataset,"",basename(file)))
        xml2::xml_add_child(visit, "filename", 
                            file.path(data_folder, file), 
                            object_id = object_id)
      }
    }
  }
  
  # Write the XML to file
  xml2::write_xml(root, output_file)
}

create_deformetrica_script <- function(data_dir = "", output_dir = "output/", script_path = "run_deformetrica.sh") {
  script_content <- sprintf('#!/bin/bash

# Run deformetrica
cd %s
source activate deformetrica
deformetrica estimate model.xml data_set.xml --p optimization_parameters.xml --output=\'%s\' -v DEBUG
conda deactivate
', data_dir, output_dir)
  
  # Write the script to a file
  writeLines(script_content, con = script_path)
  
  # Make the script executable (this works on Unix-like systems)
  if (.Platform$OS.type == "unix") {
    system(paste("chmod +x", script_path))
  }
  
  cat("Deformetrica script created at:", script_path, "\n")
}

create_model_xml <- function(sample_files, 
                             output_file,
                             deformable_object_type = "PolyLine",
                             kernel_type = "torch",
                             attachment_type = "Varifold",
                             noise_std = 1,
                             kernel_width = 11,
                             kernel_device = "cpu",
                             number_of_timepoints = 10,
                             dataset = "") {
  # Create the root element
  root <- xml2::xml_new_root("model")
  
  # Add model details
  xml2::xml_add_child(root, "model-type", "Registration")
  xml2::xml_add_child(root, "dimension", "3")
  
  # Add template element
  template <- xml2::xml_add_child(root, "template")
  
  # Ensure all parameter lists have the same length as sample_files
  n_files <- length(sample_files)
  deformable_object_type <- rep_len(deformable_object_type, n_files)
  kernel_type <- rep_len(kernel_type, n_files)
  attachment_type <- rep_len(attachment_type, n_files)
  noise_std <- rep_len(noise_std, n_files)
  kernel_width <- rep_len(kernel_width, n_files)
  kernel_device <- rep_len(kernel_device, n_files)
  
  # Add object elements for each sample file
  for (i in seq_along(sample_files)) {
    file <- sample_files[i]
    object_id <- tools::file_path_sans_ext(gsub(dataset,"",basename(file)))
    object <- xml2::xml_add_child(template, "object", id = object_id)
    xml2::xml_add_child(object, "deformable-object-type", deformable_object_type[i])
    xml2::xml_add_child(object, "attachment-type", attachment_type[i])
    xml2::xml_add_child(object, "noise-std", as.character(noise_std[i]))
    xml2::xml_add_child(object, "kernel-width", as.character(kernel_width[i]))
    xml2::xml_add_child(object, "kernel-type", kernel_type[i])
    xml2::xml_add_child(object, "kernel-device", kernel_device[i])
    xml2::xml_add_child(object, "filename", file)
  }
  
  # Add deformation parameters
  deform_params <- xml2::xml_add_child(root, "deformation-parameters")
  xml2::xml_add_child(deform_params, "kernel-width", as.character(max(kernel_width)))
  xml2::xml_add_child(deform_params, "kernel-type", kernel_type[1])
  xml2::xml_add_child(deform_params, "number-of-timepoints", as.character(number_of_timepoints))
  
  # Write the XML to file
  xml2::write_xml(root, output_file)
}

create_optimization_xml <- function(
    output_file = "optimization_parameters.xml",
    optimization_method_type = "ScipyLBFGS",
    max_iterations = 100,
    save_every_n_iters = 10,
    print_every_n_iters = 1,
    convergence_tolerance = 1e-4,
    initial_step_size = 1e-6,
    use_cuda = "Off",
    number_of_processes = 10,
    freeze_template = "On"
) {
  # Create the XML content
  xml_content <- paste0(
    '<?xml version="1.0"?>\n',
    '<optimization-parameters>\n\n',
    '    <optimization-method-type>', optimization_method_type, '</optimization-method-type>\n\n',
    '    <max-iterations>', max_iterations, '</max-iterations>\n',
    '    <save-every-n-iters>', save_every_n_iters, '</save-every-n-iters>\n',
    '    <print-every-n-iters>', print_every_n_iters, '</print-every-n-iters>\n\n',
    '    <convergence-tolerance>', format(convergence_tolerance, scientific = TRUE), '</convergence-tolerance>\n',
    '    <initial-step-size>', format(initial_step_size, scientific = TRUE), '</initial-step-size>\n\n',
    '    <use-cuda>', use_cuda, '</use-cuda>\n\n',
    '    <number-of-processes>', number_of_processes, '</number-of-processes>\n\n',
    '    <freeze-template>', freeze_template, '</freeze-template>\n\n',
    '</optimization-parameters>'
  )
  
  # Write the XML content to the file
  writeLines(xml_content, output_file)
  
  # Return the file path
  return(output_file)
}

# # For the dataset XML
# data_folder <- "data/"
# subject_files <- c("data/DNa01_left_banc.vtk", "data/DNa01_right_banc.vtk",
#                   "data/DNa03_left_banc.vtk", "data/DNa03_right_banc.vtk")
# 
# # Simple usage (same as before)
# create_dataset_xml(data_folder, subject_files, "dataset.xml")
# 
# # # Advanced usage with custom subject and visit IDs
# # subject_ids <- list(
# #   subject1 = c("DNa01_left_fafb.vtk", "DNa01_right_fafb.vtk"),
# #   subject2 = c("DNa03_left_fafb.vtk", "DNa03_right_fafb.vtk")
# # )
# # visit_ids <- list(
# #   visit1 = c("DNa01_left_fafb.vtk", "DNa03_left_fafb.vtk"),
# #   visit2 = c("DNa01_right_fafb.vtk", "DNa03_right_fafb.vtk")
# # )
# # create_dataset_xml(data_folder, subject_files, "dataset_advanced.xml",
# #                    subject_ids = subject_ids, visit_ids = visit_ids)
# 
# # For the model XML
# sample_files <- c("DNa01_left_fafb.vtk", "DNa01_right_fafb.vtk",
#                    "DNa03_left_fafb.vtk", "DNa03_right_fafb.vtk")
# 
# # Simple usage (same as before)
# create_model_xml(sample_files, "model.xml", deformable_object_type = "SurfaceMesh")
# 
# # # Advanced usage with custom parameters
# # create_model_xml(sample_files, "model_advanced.xml",
# #                  deformable_object_type = "SurfaceMesh",
# #                  kernel_type = "custom",
# #                  attachment_type = "Custom",
# #                  noise_std = 2,
# #                  kernel_width = 15,
# #                  kernel_device = "gpu",
# #                  number_of_timepoints = 20)

# hidden
write_neuron_to_vtk_paired <- function(neuron, file) {
  
  # Extract points from the neuron
  points <- nat::xyzmatrix(neuron)
  
  # Open the file for writing
  con <- file(file, "w")
  
  # Write VTK header
  writeLines("# vtk DataFile Version 3.0", con)
  writeLines("Neuron VTK file", con)
  writeLines("ASCII", con)
  writeLines("DATASET POLYDATA", con)
  
  # Write POINTS section
  writeLines(sprintf("POINTS %d float", nrow(points)), con)
  write.table(points, con, row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # Prepare LINES section
  line_data <- neuron$d %>%
    dplyr::mutate(Parent = Parent,
                  PointNo  = PointNo) %>%
    dplyr::mutate(from = (1:dplyr::n())-1,
                  to = from[match(Parent,PointNo)],
                  pair = 2) %>%
    dplyr::filter(!is.na(to)) %>%
    dplyr::select(pair, to, from)
  num_pairs <- nrow(line_data)
  
  # Write LINES section
  writeLines(sprintf("LINES %d %d", num_pairs, num_pairs * 3), con)
  write.table(line_data, con, row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # Close the file
  close(con)
  cat(sprintf("VTK file written to: %s\n", file))
}

list_files_age_sorted <- function(path = ".", ...) {
  files <- list.files(path, full.names = TRUE, ...)
  file_info <- file.info(files)
  sorted_files <- file_info[order(file_info$mtime), ]
  return(rownames(sorted_files))
}

# flywire read swc files
banc_read_swc <- function(swc, 
                          meta, 
                          ids = NULL, 
                          id = 'root_id', 
                          template = "BANC",
                          nams = NULL){
  if(!is.null(ids)){
    banc_skels <- nat::read.neurons(swc, 
                                    pattern = paste(ids,collapse="|"), 
                                    OmitFailures = TRUE)
  }else{
    banc_skels <- nat::read.neurons(swc,
                                    pattern = "\\.swc",
                                    OmitFailures = TRUE)
  }
  if(is.null(nams)){
    banc_skels <- banc_skels[!duplicated(names(banc_skels))]
    df <- data.frame(root_id=names(banc_skels))
  }else{
    new.nams <- nams[names(banc_skels)]
    new.nams[is.na(new.nams)] = names(banc_skels)[is.na(new.nams)]
    banc_skels <- banc_skels[!duplicated(new.nams)]
    df <- data.frame(root_id=new.nams)
  }
  colnames(df) <- id
  attr(banc_skels,"df") <- df
  banc_skels[,"dataset"] = template
  nat.templatebrains::regtemplate(banc_skels) = regtemplate
  meta <- meta[!duplicated(meta[[id]]),]
  banc_skels <- hemibrainr:::update_metdata(banc_skels, meta = meta, id = id)
  banc_skels
}

banc_read_metrics_csvs <- function(save.path, 
                                   delete.bad = FALSE, 
                                   delete.ids = NULL){
  files <- list.files(save.path, full.names = TRUE)
  message("Reading metrics: ")
  metrics.csvs <- pbapply::pblapply(files, function(file){
    root_id <- gsub("\\.csv","",basename(file))
    csv <- try(suppressMessages(readr::read_csv(file, 
                                                col_types = hemibrainr:::sql_col_types,
                                                progress = FALSE,
                                                show_col_types = FALSE)))
    if(delete.bad){
      csv.old <- nrow(csv)
      csv <- csv %>%
        dplyr::filter(!is.na(root_id),
                      #!is.na(total_inputs),
                      !is.na(cable_length)
        )
      if(!is.null(delete.ids)){
        csv <- csv %>%
          dplyr::filter(!root_id%in%delete.ids)
      }
      miss <- csv.old-nrow(csv)
      if(miss>=1){
        warning("deleting ", csv.old-nrow(csv)," entries from ", file)
      }
      readr::write_csv(x=csv,file=file)
    }
    csv
  })
  metrics.csvs <- metrics.csvs[sapply(metrics.csvs,is.data.frame)]
  do.call(plyr::rbind.fill, metrics.csvs)
}

banc_read_swc_split <-function(swc, 
                               synapses, 
                               meta, 
                               ids = NULL, 
                               regtemplate = "BANC"){
  message("reading skeletons: ")
  banc_skels <- banc_read_swc(swc = swc,
                              ids = ids,
                              meta = meta)
  pb <- progress::progress_bar$new(total = length(banc_skels))
  message("reading synapses: ")
  for(id in names(banc_skels)){
    pb$tick()
    file <- file.path(synapses,paste0(id,".csv"))
    csv <- tryCatch(suppressWarnings(readr::read_csv(file, col_types = hemibrainr:::sql_col_types)), 
                    error = function(e) NULL)
    if(is.null(csv)){
      csv <- utils::read.csv(file)
      idx <- names(hemibrainr:::sql_col_types$cols) %in% colnames(csv)
      sql_col_types <- hemibrainr:::sql_col_types
      sql_col_types$cols <- sql_col_types$cols[idx]
      csv <- readr::type_convert(csv, col_types = sql_col_types, guess_integer = TRUE) # should work, doesn't
    }
    if(!nrow(csv)){
      banc_skels[match(id,names(banc_skels))] <- NULL
      next
    }
    csv <- banc_process_synapses(csv)
    banc_skels[match(id,names(banc_skels))][[1]]$connectors <- csv
  }
  banc_skels <- nat::nlapply(banc_skels, function(neuron){
    class(neuron) <- c("catmaidneuron","neuprintneuron","neuron")
    neuron
  })
  banc_skels
}

assign_strahler <- function(x, ...) UseMethod("assign_strahler")
assign_strahler.neuronlist <- function(x, ...){
  nlapply(x, assign_strahler.neuron, ...)
}
assign_strahler.neuron <- function(x, ...){
  if(ifelse(!is.null(x$nTrees),x$nTrees!=1,FALSE)){
    warning("Neuron has multiple trees, calculating Strahler order for each subtree separately")
    x$d$strahler_order = 1
    for(tree in 1:x$nTrees){
      v = unique(unlist(x$SubTrees[tree]))
      if(length(v)<2){
        x$d[x$d$PointNo%in%v,]$strahler_order = 1
      }else{
        neuron = tryCatch(nat::prune_vertices(x,verticestoprune = v, invert = TRUE), error=function(e) NULL)
        if(sum(branchpoints(x)%in%v)==0){
          x$d[x$d$PointNo%in%v,]$strahler_order = 1
        }else if (!is.null(neuron)){
          s = nat::strahler_order(neuron)
          x$d[x$d$PointNo%in%v,]$strahler_order = s$points
        }
      }
    }
  }else{
    s = nat::strahler_order(x)
    x$d$strahler_order = s$points
  }
  if("synaptic"%in%class(x)){
    relevant.points = subset(x$d, PointNo%in%x$connectors$treenode_id)
    x$connectors$strahler_order = relevant.points[match(x$connectors$treenode_id,relevant.points$PointNo),]$strahler_order
  }
  x
}

# Rvcg
find_closest_region <- function(df, mesh_list, max.dist = 5000) {
  
  # Function to find distance from point to mesh
  point_to_mesh_distance <- function(point, mesh) {
    distances <- Rvcg::vcgClostKD(mesh=mesh, x=point)
    distances <- distances$quality
    point$distances <- distances
    return(point)
  }
  
  # Iterate through rows where region is NA
  distances <- pbapply::pblapply(mesh_list$RegionList, function(reg){
    p <- point_to_mesh_distance(nat::xyzmatrix(df), mesh = as.mesh3d(subset(mesh_list,reg)))
    p$distances
  })
  
  # Find nearest mesh
  distances.m <- do.call(cbind,distances)
  colnames(distances.m) <- mesh_list$RegionList
  chosen <- apply(abs(distances.m), 1, which.min)
  min.dists <- apply(abs(distances.m), 1, function(row) min(row)<max.dist)
  df$neuropil <- colnames(distances.m)[chosen]
  
  # Assign
  df <- df %>%
    dplyr::mutate(region = dplyr::case_when(
      grepl("vnc",neuropil) ~ "vnc",
      grepl("optic|^LO|^LOP|^ME|^AME",neuropil) ~ "optic_lobes",
      grepl("GNG|CAN|FLA|AMMC|SAD|PRW",neuropil) ~ "subesophageal_zone",
      TRUE ~ "central_brain",
    ))
  
  # Determine that some are outside
  df$neuropil[!min.dists] <- paste0("outside_",df$neuropil[!min.dists])
  df$region[!min.dists] <- paste0("outside_",df$region[!min.dists])
  
  # Return
  return(df)
}

# Find which neuropil surfaces synapses are nearest to
pointsnearby_banc <- function(x,id="id",only.missing = TRUE){
  
  # Get volume list
  volumes <- c(banc_brain_neuropils.surf,banc_vnc_neuropils.surf)
  
  # Neuropil missing
  if(only.missing){
    x.no.neuropil <- x %>%
      dplyr::filter((is.na(region)|is.na(neuropil)|grepl("^brain|outside",region))|grepl("outside",neuropil))
    x.neuropil <- x %>%
      dplyr::anti_join(x.no.neuropil, by=id)
    x.corrected <- find_closest_region(x.no.neuropil, volumes) 
  }else{
    x.corrected <- find_closest_region(x, volumes) 
  }
  
  # Re-combine and return
  rbind(x.neuropil,x.corrected)
  
}

# Find which neuropil synapses are inside of
pointsinside_banc <- function(x,
                              neuropils = list(banc_brain_neuropils.surf,
                                               banc_vnc_neuropils.surf),
                              volumes = list(neck = banc_neck_connective.surf,
                                             brain = banc_brain_neuropil.surf,
                                             optic_lobes = as.hxsurf(as.mesh3d(subset(banc_brain_neuropils.surf,"optic"))),
                                             sez = as.hxsurf(as.mesh3d(subset(banc_brain_neuropils.surf,"GNG|CAN|FLA|AMMC|SAD|PRW"))),
                                             central_brain = as.hxsurf(as.mesh3d(subset(banc_brain_neuropils.surf,"midbrain"))),
                                             vnc = banc_vnc_neuropil.surf),
                              alpha = 50000,
                              scaling = NULL){
  df = as.data.frame(x)
  df$neuropil = NA
  df$region = NA
  df$side = NA
  df$neuropil <- ""
  df$region <- ""
  df$side <- ""
  points = nat::xyzmatrix(df)
  if(!is.null(scaling)){
    points = points/scaling
  } 
  lrdiffs <- bancr:::banc_lr_position(points,units = "nm")
  sides <- ifelse(lrdiffs>0,"right","left")
  df$side <- sides
  for(vol in 1:length(volumes)){
    neuropil = volumes[[vol]]
    reg = names(volumes)[vol]
    if (!is.null(alpha)) {
      neuropil = alphashape3d::ashape3d(nat::xyzmatrix(neuropil), 
                                        alpha = alpha)
      a = alphashape3d::inashape3d(points = points, 
                                   as3d = neuropil, 
                                   indexAlpha = "ALL")
    }
    else {
      a = nat::pointsinside(x = points, surf = neuropil)
    }
    if(sum(a)) df$region[which(a == T)] = reg
  } 
  for(brain in neuropils){
    nps = sort(brain$RegionList)
    for (np in nps) {
      neuropil <- subset(brain, np)
      region <- NA
      if(is.na(region)) region = ifelse(np %in%banc_vnc_neuropils.surf$RegionList,"vnc",NA)
      if(is.na(region)) region = ifelse(grepl("^LO|^ME|^AME|^LOP",np),"optic_lobes",NA)
      if(is.na(region)) region = ifelse(grepl("^CAN|^GNG|^FLA|^AMMC|^SAD|^PRW",np),"suboesophageal_zone",NA)
      if(is.na(region)) region = ifelse(np %in%banc_brain_neuropils.surf$RegionList,"central_brain",NA)
      if (!is.null(alpha)) {
        neuropil = alphashape3d::ashape3d(nat::xyzmatrix(neuropil), 
                                          alpha = alpha)
        a = alphashape3d::inashape3d(points = points, 
                                     as3d = neuropil, 
                                     indexAlpha = "ALL")
      }
      else {
        a = nat::pointsinside(x = points, surf = neuropil)
      }
      if(sum(a)){
        df$neuropil[which(a)] = sapply(df$neuropil[which(a)], function(x) paste(unique(unlist(strsplit(paste(x,np,sep=","),split=","))),sep=",",collapse=","))
        unassigned <- which(a)[df$region[which(a)] == ""]
        if (length(unassigned)) df$region[unassigned] = region
      }
    } 
  }
  df <- df %>%
    dplyr::mutate(neuropil = ifelse(neuropil=="","outside",neuropil),
                  region = ifelse(region=="","outside",region)) %>%
    dplyr::mutate(neuropil = gsub("^,","",neuropil),
                  region = gsub("^,","",region)) 
  df
}

# round numbers
round_dataframe <- function(x, exclude=NULL, digits = 4, ...) {
  numcols <- names(x)[sapply(x, function(c) is.numeric(c) && !inherits(c, 'integer64'))]
  numcols <- setdiff(numcols, exclude)
  for(i in numcols) {
    col=x[[i]]
    # does it look like an int, if so, make it one
    intcol=try(checkmate::asInteger(col), silent = TRUE)
    if((sum(is.na(col))==length(col))){
      x[[i]]=col
    }else if(is.integer(intcol)){
      x[[i]]=intcol
    }
    else{
      x[[i]]= signif(col, digits)
    }
  }
  x
}

read_synapse_csvs <- function(save.path,
                              lookup = FALSE,
                              ids = NULL,
                              numCores = NULL){
  if(!is.null(numCores)){
    cl <- parallel::makeCluster(numCores)
  }else{
    cl = NULL
  }
  if(length(save.path)>1){
    files <- save.path
  }else{
    files <- list.files(save.path, full.names = TRUE)
  }
  message("Reading synapses: ")
  synapse.csvs <- pbapply::pblapply(files, function(file){
    root_id <- gsub("\\.csv","",basename(file))
    if(!is.null(ids)){
      good <- root_id %in% ids
    }else{
      good <- TRUE
    }
    if(good){
      csv <- try(suppressMessages(data.table::fread(file,
                                                    showProgress = FALSE)),
                  silent = TRUE)
      if(is.data.frame(csv)){
        if(nrow(csv)){
          csv$root_id <- root_id
          if(!is.null(ids)){
            csv <- csv %>%
              dplyr::filter(pre_id %in% ids, 
                            post_id %in% ids)
          }
          if(lookup){
            csv <- csv %>%
              dplyr::rename(label = Label) %>%
              dplyr::mutate(
                post_label = ifelse(prepost==1,label,NA)
              ) %>%
              dplyr::mutate(
                pre_label = ifelse(prepost==0,label,NA)
              ) %>%
              dplyr::distinct(connector_id, 
                              pre_id, post_id,
                              pre_label, post_label)
          }
          csv  
        }else{
          NULL
        }
      }else{
        NULL
      } 
    }
  }, 
  cl = cl)
  if(!is.null(cl)){
    parallel::stopCluster(cl)
  }
  synapse.csvs <- synapse.csvs[sapply(synapse.csvs,is.data.frame)]
  do.call(plyr::rbind.fill, synapse.csvs)
}

read_elist_csvs <- function(save.path,
                            split = FALSE,
                            lookup = NULL,
                            ids = NULL,
                            numCores = NULL){
  if(!is.null(numCores)){
    cl <- parallel::makeCluster(numCores)
  }else{
    cl = NULL
  }
  if(length(save.path)>1){
    files <- save.path
  }else{
    files <- list.files(save.path, full.names = TRUE)
  }
  message("Building edgelist: ")
  synapse.csvs <- pbapply::pblapply(files, function(file){
    root_id <- gsub("\\.csv","",basename(file))
    csv <- try(suppressWarnings(readr::read_csv(file, col_types = hemibrainr:::sql_col_types)))
    if(is.data.frame(csv)){
      if(!is.null(ids)){
        csv <- csv %>%
          dplyr::filter(pre_id %in% ids, post_id %in% ids)
      }
      if(nrow(csv)){
        csv$root_id <- root_id
        if(split){
          csv <- csv %>%
            dplyr::rename(label = Label) %>%
            dplyr::rename(post = post_id, pre = pre_id, post_label = label) %>% 
            dplyr::rowwise() %>%
            dplyr::left_join(lookup, by = "connector_id") %>%
            dplyr::ungroup() %>%
            dplyr::mutate(pre_count = sum(prepost==0)) %>%
            dplyr::mutate(post_count = sum(prepost==1)) %>%
            dplyr::mutate(pre_label = hemibrainr:::standard_compartments(pre_label)) %>%
            dplyr::mutate(post_label = hemibrainr:::standard_compartments(post_label)) %>%
            dplyr::filter(prepost==1) %>%
            dplyr::group_by(post, post_label) %>%
            dplyr::mutate(post_label_count = dplyr::n()) %>% 
            dplyr::group_by(pre, post, post_label, pre_label) %>%
            dplyr::mutate(count = dplyr::n()) %>% 
            dplyr::ungroup() %>%
            dplyr::mutate(norm = count/post_count,
                          norm_label = count/post_label_count) %>%
            dplyr::filter(count > 0) %>% 
            dplyr::distinct(post, 
                            pre,
                            post_svid,
                            pre_svid,
                            post_label, 
                            pre_label,
                            count, 
                            norm,
                            norm_label,
                            post_label_count,
                            post_count,
                            pre_count,
                            root_id,
                            .keep_all = FALSE) %>%
            as.data.frame()
          csv$norm = round(csv$norm, digits = 6)
          csv$norm_label = round(csv$norm_label, digits = 6)
          csv
        }else{
          csv <- csv %>%
            dplyr::rename(post = post_id, pre = pre_id) %>% 
            dplyr::mutate(pre_count = sum(prepost==0)) %>%
            dplyr::mutate(post_count = sum(prepost==1)) %>%
            dplyr::filter(prepost==1) %>%
            dplyr::group_by(pre, post) %>%
            dplyr::mutate(count = dplyr::n()) %>% 
            dplyr::ungroup() %>%
            dplyr::mutate(norm = count/post_count) %>%
            dplyr::filter(count > 0) %>% 
            dplyr::distinct(post, 
                            pre, 
                            post_svid,
                            pre_svid,
                            count, 
                            norm,
                            post_count,
                            pre_count,
                            root_id,
                            .keep_all = FALSE) %>%
            as.data.frame()
          csv$norm = round(csv$norm, digits = 6)
          csv
        }
      }else{
        NULL
      }
    }else{
      NULL
    }
  }, 
  cl = cl)
  if(!is.null(cl)){
    parallel::stopCluster(cl)
  }
  synapse.csvs <- synapse.csvs[sapply(synapse.csvs,is.data.frame)]
  do.call(plyr::rbind.fill, synapse.csvs)
}

read_metrics_csvs <- function(save.path){
  files <- list.files(save.path, full.names = TRUE)
  message("Reading metrics: ")
  metrics.csvs <- pbapply::pblapply(files, function(file){
    root_id <- gsub("\\.csv","",basename(file))
    csv <- try(suppressMessages(readr::read_csv(file, 
                                                col_types = hemibrainr:::sql_col_types,
                                                progress = FALSE,
                                                show_col_types = FALSE)))
  })
  metrics.csvs <- metrics.csvs[sapply(metrics.csvs,is.data.frame)]
  do.call(plyr::rbind.fill, metrics.csvs)
}

########################
### folder wrangling ###
########################

# Function to get all PNG files recursively, discard  = "/done$"
get_png_files <- function(base_dir, discard = NULL) {
  if(!is.null(discard)){
    dir_ls(base_dir, recurse = TRUE, glob = "*.png") %>%
      discard(~grepl(discard, .x))
  }else{
    dir_ls(base_dir, recurse = TRUE, glob = "*.png")
  }
}

# Updated function to move immediate subfolder containing the file
move_subfolder <- function(file, destination, base_images, verbose = FALSE) {
  # Get the immediate parent folder of the file
  parent_folder <- path_dir(file)
  
  # Calculate relative path, excluding 'todo' and 'done' from the path
  path_relative <- path_rel(parent_folder, base_images)
  path_relative <- gsub("/todo/|todo/|/done/|done/", "", path_relative)
  
  # Construct new path for the parent folder
  new_path <- path(destination, path_relative)
  
  # Create destination directory if it doesn't exist
  dir_create(new_path)
  
  # Move all contents of the parent folder to the new location
  files_to_move <- dir_ls(parent_folder)
  walk(files_to_move, function(f) {
    file_move(f, path(new_path, path_file(f)))
  })
  
  # Remove the original parent folder if it's now empty
  if (length(dir_ls(parent_folder)) == 0) {
    dir_delete(parent_folder)
  }
  if(verbose){
    cat("Moved folder:", parent_folder, "to", new_path, "\n")
  }
}

# Batch move files by grouping on parent directory
# Instead of calling move_subfolder per-file (which does dir_ls per call),
# this groups files by their parent folder and processes each unique parent once.
batch_move_subfolders <- function(files, destination, base_images, verbose = FALSE) {
  if (!length(files)) return(invisible())

  # Group files by parent directory
  parents <- path_dir(files)
  unique_parents <- unique(parents)

  n_moved <- 0
  for (parent in unique_parents) {
    # Skip if parent was already emptied by a previous iteration
    if (!dir_exists(parent)) next

    # Calculate relative path, stripping todo/done from path
    path_relative <- path_rel(parent, base_images)
    path_relative <- gsub("/todo/|todo/|/done/|done/", "", path_relative)

    # Construct and create destination
    new_path <- path(destination, path_relative)
    dir_create(new_path)

    # Move all contents of this parent folder at once
    sibling_files <- dir_ls(parent)
    for (f in sibling_files) {
      tryCatch(
        file_move(f, path(new_path, path_file(f))),
        error = function(e) NULL
      )
    }
    n_moved <- n_moved + length(sibling_files)

    # Remove parent if now empty
    remaining <- tryCatch(dir_ls(parent), error = function(e) character(0))
    if (length(remaining) == 0) {
      tryCatch(dir_delete(parent), error = function(e) NULL)
    }
  }

  if (verbose) {
    cat("Batch moved", n_moved, "files across", length(unique_parents),
        "directories to", destination, "\n")
  }
}

# Function to remove empty directories, excluding 'express' and 'done'
# Uses shell find command for speed instead of R dir_ls + walk
remove_empty_dirs <- function(base_dir, clean_thumbs = TRUE) {
  if (clean_thumbs) delete_thumbs_db(base_dir)
  # Use find to delete empty dirs, excluding express/ and done/
  # Run repeatedly until no more empty dirs (handles nested empties)
  for (i in 1:5) {
    result <- system2("find", c(
      base_dir,
      "-mindepth", "1",
      "-not", "-path", "*/express/*",
      "-not", "-path", "*/done/*",
      "-not", "-name", "express",
      "-not", "-name", "done",
      "-type", "d", "-empty", "-delete",
      "-print"
    ), stdout = TRUE, stderr = TRUE)
    if (length(result)) {
      cat("Removed empty directories:\n")
      cat(paste(result, collapse = "\n"), "\n")
    } else {
      break
    }
  }
}

extract_all_supervoxel_ids <- function(string) {
  # Split the string by underscores
  parts <- tryCatch(strsplit(string, "_")[[1]],error = function(e) NA)
  if(all(is.na(parts))){
    return(NA)
  }
  
  # Find all indices of "supervoxel" and "id"
  supervoxel_indices <- which(parts == "supervoxel")
  id_indices <- which(parts == "id")
  
  # Initialize a vector to store results
  results <- character()
  
  # Check each pair of "supervoxel" and "id"
  for (i in seq_along(supervoxel_indices)) {
    supervoxel_index <- supervoxel_indices[i]
    
    # Find the next "id" after this "supervoxel"
    next_id_index <- id_indices[id_indices > supervoxel_index][1]
    
    # Check if we found a matching "id"
    if (!is.na(next_id_index) && next_id_index == supervoxel_index + 1) {
      # Get the index of the ID
      id_position <- next_id_index + 1
      
      # Extract the ID if it exists
      if (id_position <= length(parts)) {
        id <- parts[id_position]
        # Check if the ID is numeric and add it to results
        if (grepl("^\\d+$", id)) {
          results <- c(results, id)
        }
      }
    }
  }
  
  # Return the results (or NA if no valid IDs found)
  if (length(results) > 0) {
    return(results)
  } else {
    return(NA)
  }
}

# Function to sort todo and done folders for each task
sort_todos <- function(base_images, base_correct,
                       use.latest = FALSE,
                       remove.duplicates = FALSE,
                       remove.old.roots = FALSE,
                       correct_queries = NULL){

  # Clean junk files once per directory (not again in remove_empty_dirs)
  cat("Deleting junk files  \n")
  delete_thumbs_db(base_correct)
  delete_thumbs_db(base_images)

  # Get all PNG files in both directories (single scan each)
  images_files <- get_png_files(base_images, discard  = "/express|/done")
  correct_files <- get_png_files(base_correct)

  # Remove duplicates
  if(remove.duplicates){
    images_dirs <-  list.dirs(base_images, recursive = TRUE)
    images_dirs <- images_dirs[!grepl("done|express|correct",images_dirs)]
    images_dirs <- images_dirs[grepl("root_id_",basename(images_dirs))]
    image_queries <- regmatches(images_dirs, regexpr("(?<=root_id_)\\d+", images_dirs, perl = TRUE))
    dupes <- images_dirs[duplicated(image_queries)]
    message("Deleted duplicated folders: ", length(dupes))
    unlink(dupes, recursive = TRUE)
  }
  if(remove.old.roots){
    # Re-use images_files from above instead of re-scanning
    image_queries <- regmatches(images_files, regexpr("(?<=root_id_)\\d+", images_files, perl = TRUE))
    olds <- images_files[!banc_islatest(image_queries)]
    message("Deleted outdated files: ", length(olds))
    file.remove(olds)
    # Re-scan after deletions
    images_files <- get_png_files(base_images, discard  = "/express|/done")
  }

  # Extract filenames without path, for correct
  images_filenames <- path_file(images_files)
  correct_filenames <- path_file(correct_files)
  if(is.null(correct_queries)){
    correct_queries <- regmatches(correct_filenames, regexpr("(?<=root_id_)\\d+", correct_filenames, perl = TRUE))
    if(use.latest){

      # Get files
      correct_files <- get_png_files(base_correct, discard  = "^[1-5]")
      correct_files <- as.character(basename(correct_files))
      correct_filenames <- path_file(correct_files)

      # Update by supervoxel IDs
      correct_svid_queries <- extract_all_supervoxel_ids(correct_filenames)
      correct_svid_queries <- setdiff(correct_svid_queries, "NA")

      # Update paired supervoxel IDs
      correct_filenames_svids <- grep("_supervoxel_id_",correct_filenames,value = TRUE)
      correct_sv_id_queries <- gsub(".*_supervoxel_id_|_hit_.*|_query_.*","",correct_filenames_svids)
      correct_sv_id_queries <- unique(unlist(strsplit(correct_sv_id_queries,split="_")))
      correct_sv_id_queries <- unique(c(correct_svid_queries,setdiff(correct_sv_id_queries, "NA")))
      correct_sv_id_queries <- banc_rootid(correct_sv_id_queries)
      correct_queries <- unique(c(correct_queries,correct_sv_id_queries))

      # Update by root IDs if not supervoxel IDs
      correct_filenames_roots <- grep("_nucleus_id_|_supervoxel_id_NA_hit_id",correct_filenames,value = TRUE)
      correct_root_id_queries <- gsub(".*_root_id_|_nucleus_id.*|_supervoxel_id_.*","",correct_filenames_roots)
      correct_root_id_queries <- unique(unlist(strsplit(correct_root_id_queries,split="_")))
      correct_root_id_queries <- pbapply::pbsapply(correct_root_id_queries, function(x) try(quiet_function(banc_updateids(x))))
      correct_queries <- unique(c(correct_queries,correct_root_id_queries))
      correct_queries <- correct_queries[correct_queries!="0"]
    }
    if(!length(correct_queries)){
      correct_queries <- "none"
    }
  }

  # --- Vectorized classification ---
  # Extract query root_ids from all image filenames at once
  cat("Sorting files in 'images': ",  length(images_files), " \n")

  if (!use.latest) {
    # Fast path: vectorized extraction + set membership
    # regmatches with regexpr drops non-matches, so use match positions to
    # keep alignment with images_files
    match_pos <- regexpr("(?<=root_id_)\\d+", images_filenames, perl = TRUE)
    has_match <- match_pos > 0
    queries <- rep(NA_character_, length(images_filenames))
    queries[has_match] <- regmatches(images_filenames, match_pos)

    # Classify: done (in correct_queries) vs todo (not in correct_queries) vs skip (no query)
    is_valid <- !is.na(queries) & queries != "0"
    is_done <- is_valid & queries %in% correct_queries
    is_todo <- is_valid & !is_done

    n_skipped <- sum(!is_valid)
    if (n_skipped > 0) {
      bad_files <- images_filenames[!is_valid]
      message("Skipping ", n_skipped, " files with no/invalid root_id")
    }

    cat("  -> done: ", sum(is_done), ", todo: ", sum(is_todo),
        ", skipped: ", n_skipped, "\n")

    # Batch move by parent directory
    batch_move_subfolders(images_files[is_done],
                          destination = path(base_images, "done"),
                          base_images = base_images, verbose = TRUE)
    batch_move_subfolders(images_files[is_todo],
                          destination = path(base_images, "todo"),
                          base_images = base_images)

  } else {
    # Slow path: use.latest requires per-file API calls to resolve current root_ids
    for(file in images_files){
      try({
        if(!file.exists(file)) next
        filename <- path_file(file)
        query.svid <- as.character(regmatches(filename, regexpr("(?<=supervoxel_id_)\\d+", filename, perl = TRUE)))
        query <- regmatches(filename, regexpr("(?<=root_id_)\\d+", filename, perl = TRUE))
        if(use.latest){
          if(length(query.svid)){
            query <- tryCatch(banc_rootid(query.svid), error = function(e){
              message("supervoxel_id could not be updated: ", query.svid, " using root_id")
              query <- regmatches(filename, regexpr("(?<=root_id_)\\d+", filename, perl = TRUE))
              query <- banc_updateids(query)
              query
            })
          }else{
            query <- regmatches(filename, regexpr("(?<=root_id_)\\d+", filename, perl = TRUE))
            query <- banc_updateids(query)
          }
        }
        # Fix: check length before is.na to avoid "argument is of length zero"
        if(length(query) == 0 || is.na(query) || query == "0"){
          message("Failed to find root ID for: ", basename(filename))
          next
        }
        if (query %in% correct_queries) {
          move_subfolder(file, dest=path(base_images, "done"), base_images=base_images, verbose = TRUE)
        } else {
          move_subfolder(file, dest=path(base_images, "todo"), base_images=base_images)
        }
      })
    }
  }

  # Remove empty directories (skip thumbs cleanup — already done above)
  cat("Removing empty directories from: ",  base_images, " \n")
  remove_empty_dirs(base_images, clean_thumbs = FALSE)

  # Return
  invisible()
}

# Function to sync files
sync_files <- function(src, dest, 
                       remote.source = FALSE,
                       move.old = FALSE, 
                       extensions = c("png","csv")) {
  
  if(remote.source){
    # Change permissions on source files if remote
    system(paste(c("ssh $(whoami)@transfer.rc.hms.harvard.edu find", src, "-type", "f", "-name", "*.png", "-o", "-name", "*.csv", "-exec", "chmod", "644", "{}", "+"),collapse=" "))
    system(paste0(c("ssh $(whoami)@transfer.rc.hms.harvard.edu find", src, "-type", "d", "-exec", "chmod", "755", "{}", "+"),collapse=" "))
  }else{
    # Change permissions on source files
    system2("find", c(src, "-type", "f", "-name", "*.png", "-o", "-name", "*.csv", "-exec", "chmod", "644", "{}", "+"))
    system2("find", c(src, "-type", "d", "-exec", "chmod", "755", "{}", "+")) 
  }
  
  # Construct the rsync command.
  # Arguments pass through TWO shells: local sh (from system()) then remote sh
  # (from ssh). ssh joins its args into a single string for the remote shell,
  # so we must double-quote: outer quotes for local shell, inner quotes survive
  # to remote shell. Use: "'arg with spaces'" — local shell strips outer quotes
  # and passes literal 'arg with spaces' to ssh; remote shell respects inner quotes.
  include_patterns <- paste(paste0("\"'--include=*.", extensions, "'\""), collapse = " ")
  rsync_cmd <- paste(
    "ssh $(whoami)@transfer.rc.hms.harvard.edu rsync",
    "-arlptv --no-g --recursive",
    "\"'--include=*/'\"",
    include_patterns,
    "\"'--exclude=*'\"",
    "--prune-empty-dirs",
    "--delete",
    "\"'--filter=P done/***'\"",
    "\"'--filter=P old***/***'\"",
    "\"'--filter=P express/***'\"",
    ifelse(move.old,"--backup",""),
    ifelse(move.old,sprintf("\"'--backup-dir=%s/%s%s/'\"",dest,"old_",Sys.Date()),""),
    shQuote(paste0(src,"/")),
    shQuote(paste0(dest,"/"))
  )
  cat("command: ",rsync_cmd)
  
  # Run rsync using system()
  result <- system(rsync_cmd, intern = TRUE)
  
  # # Remove empty directories in the source
  # remove_empty_dirs_cmd <- paste(
  #   "ssh $(whoami)@transfer.rc.hms.harvard.edu",
  #   "find", paste0(shQuote(dest),"/"),
  #   "-type d -empty -delete"
  # )
  # result2 <- system(remove_empty_dirs_cmd)
  
  # report
  cat(result)
  cat("Synced files from", src, "to", dest, "\n")
  cat(paste(result, collapse = "\n"), "\n")
}

# Organise neck connective task
move_files_preserve_structure <- function(source_dir, dest_dir, ids = NULL, pattern = "*.png") {
  
  # List all files in the source directory, including sub-directories
  files_to_move <- dir_ls(source_dir, recurse = TRUE, type = "file", glob = pattern)
  files_to_move.moved <- dir_ls(dest_dir, recurse = TRUE, type = "file", glob = pattern)
  files_to_move <- files_to_move[!basename(files_to_move)%in%basename(files_to_move.moved)]
  if(!is.null(ids)){
    queries <- regmatches(files_to_move, regexpr("(?<=root_id_)\\d+", files_to_move, perl = TRUE))
    correct.queries <- banc_updateids(queries)
    correct.queries <- correct.queries[correct.queries!="0"]
    files_to_move <- files_to_move[correct.queries%in%ids] 
  }
  files_to_move <- files_to_move[!is.na(files_to_move)]
  files_to_duplicates <- files_to_move[duplicated(basename(files_to_move))] 
  if(length(files_to_duplicates)){
    removed <- file.remove(files_to_duplicates)
    cat('removed duplicates:', table(removed))
  }
  
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

replace_broken_symlinks <- function(path_a, path_b) {
  # Get all symlinks in path_a
  symlinks <- dir_ls(path_a, recurse = TRUE, type = "symlink")
  
  # Filter for broken symlinks
  broken_symlinks <- symlinks[!file.exists(symlinks)]
  
  for (symlink in broken_symlinks) {
    # Get the filename of the broken symlink
    filename <- basename(symlink)
    
    # Search for the corresponding PNG file in path_b
    png_file <- dir_ls(path_b, recurse = TRUE, glob = paste0("*", filename))
    
    if (length(png_file) == 1) {
      # If a matching PNG file is found, replace the symlink
      file_delete(symlink)
      #link_delete(symlink)
      file_copy(path=png_file, new_path=symlink, overwrite = TRUE)
      cat("Replaced:", symlink, "with", png_file, "\n")
    } else if (length(png_file) > 1) {
      cat("Multiple matches found for:", symlink, "\n")
    } else {
      cat("No match found for:", symlink, "\n")
    }
  }
}

delete_thumbs_db <- function(root_dir) {
  # Find all Thumbs.db files
  thumbs_files <- dir_ls(root_dir, recurse = TRUE, glob = "*Thumbs.db")
  
  # Delete each Thumbs.db file
  for (file in thumbs_files) {
    file_delete(file)
    cat("Deleted:", file, "\n")
  }
  
  cat("Total files deleted:", length(thumbs_files), "\n")
}

banc_update_task <- function(task = c("meta","connectivity","matching","mirror","fafb","manc","hemibrain"), 
                             use.latest = FALSE){
  
  # Understand task
  task <- match.arg(task)
  message("Updating task: ", task)
  
  # O2 paths
  cat("Synchronising matching files between O2 and the fileserve \n")
  A <- '/n/data1/hms/neurobio/wilson/banc/'
  B <- '/n/files/Neurobio/wilsonlab/banc/'
  
  ### meta ###
  if(task=="meta"){
    sync_files(path(A, "meta/"), path(B, "meta/"), extensions = c("csv", "png"))
    remove_empty_dirs(path(A, "meta/"))
  }
  
  ### connectivity ###
  if(task=="connectivity"){
    sync_files(path(A, "connectivity/"), path(B, "connectivity/"), extensions = c("csv", "sqlite"))
    remove_empty_dirs(path(A, "connectivity/"))
  }
  
  ### hemibrain ###
  if(task=="mirror"){
    sync_files(path(B, "matching/mirror/correct/"), path(A, "matching/mirror/correct/"),  extensions = c("png"), move.old = FALSE)
    sort_todos(base_images=path(A, "matching/mirror/images/"),base_correct=path(A, "matching/mirror/correct/"), use.latest=use.latest)
    remove_empty_dirs(path(A, "matching/mirror/images/"))
    sync_files(path(A, "matching/mirror/images/"), path(B, "matching/mirror/images/"), extensions = c("png"), move.old = TRUE)
  }
  
  ### hemibrain ###
  if(task=="hemibrain"){
    sync_files(path(B, "matching/hemibrain/correct/"), path(A, "matching/hemibrain/correct/"), extensions = c("png"), move.old = FALSE)
    sort_todos(base_images=path(A, "matching/hemibrain/images/elastix_tpsreg_240721/"),base_correct=path(A, "matching/hemibrain/correct/"))
    remove_empty_dirs(path(A, "matching/hemibrain/images/elastix_tpsreg_240721/"))
    sync_files(path(A, "matching/hemibrain/images/elastix_tpsreg_240721/"), path(B, "matching/hemibrain/images/"), extensions = c("png"), move.old = TRUE)   
  }
  
  ### fafb ###
  if(task=="fafb"){
    sync_files(path(B, "matching/fafb/correct/"), path(A, "matching/fafb/correct/"), extensions = c("png"), move.old = FALSE)
    sort_todos(base_images=path(A, "matching/fafb/images/elastix_tpsreg_240721/"),base_correct=path(A, "matching/fafb/correct/"), use.latest=use.latest)
    remove_empty_dirs(path(A, "matching/fafb/images/elastix_tpsreg_240721/"))
    sync_files(path(A, "matching/fafb/images/elastix_tpsreg_240721/"), path(B, "matching/fafb/images/"), extensions = c("png"), move.old = TRUE)  
  }
  
  ### manc ###
  if(task=="manc"){
    sync_files(path(B, "matching/manc/correct/"), path(A, "matching/manc/correct/"), extensions = c("png"), move.old = FALSE)
    sort_todos(base_images=path(A, "matching/manc/images/elastix_tpsreg_240721/"),base_correct=path(A, "matching/manc/correct/"), use.latest=use.latest)
    remove_empty_dirs(path(A, "matching/manc/images/elastix_tpsreg_240721/"))
    sync_files(path(A, "matching/manc/images/elastix_tpsreg_240721/"), path(B, "matching/manc/images/"), extensions = c("png"), move.old = TRUE)
  }
  
  # Announce
  cat("File transfer and movement completed.\n")
  
  # Return
  invisible()
}

# Function to count files in a directory
count_files <- function(dir) {
  length(fs::dir_ls(dir, type = "file", recurse = TRUE))
}

# Function to recursively process directories
process_directories <- function(base_dir, recurse = 1) {
  # Get all subdirectories
  dirs <- fs::dir_ls(base_dir, type = "directory", recurse = recurse)
  # Create a data frame with directory paths and file counts
  result <- data.frame(
    path = c(base_dir, dirs),
    file_count = sapply(c(base_dir, dirs), count_files)
  )
  return(result)
}

# S3 method for prune_vertices on synapticneuron objects.
# hemibrainr's axonic_cable/dendritic_cable call nat::prune_vertices which
# dispatches on class. synapticneuron inherits from neuron but has a
# connectors field that must be preserved (same pattern as catnat's
# prune_vertices.catmaidneuron and hemibrainr's prune_vertices.neuprintneuron).
prune_vertices.synapticneuron <- function(x, verticestoprune, invert = FALSE, ...) {
  pruned <- nat::prune_vertices(structure(x, class = "neuron"),
                                 verticestoprune = verticestoprune,
                                 invert = invert, ...)
  if (!is.null(x$connectors) && nrow(x$connectors) > 0) {
    if ("treenode_id" %in% names(x$connectors)) {
      pruned$connectors <- x$connectors[x$connectors$treenode_id %in%
                                          pruned$d$PointNo, , drop = FALSE]
    } else {
      # Match connectors by nearest node position
      node_xyz <- as.matrix(pruned$d[, c("X", "Y", "Z")])
      syn_xyz <- as.matrix(x$connectors[, c("x", "y", "z")])
      nn <- RANN::nn2(node_xyz, syn_xyz, k = 1)
      keep <- nn$nn.dists[, 1] < 200
      pruned$connectors <- x$connectors[keep, , drop = FALSE]
    }
  }
  pruned$tags <- lapply(x$tags, function(t) t[t %in% pruned$d$PointNo])
  pruned$AD.segregation.index <- x$AD.segregation.index
  class(pruned) <- class(x)
  pruned
}
registerS3method("prune_vertices", "synapticneuron",
                  prune_vertices.synapticneuron, envir = asNamespace("nat"))

