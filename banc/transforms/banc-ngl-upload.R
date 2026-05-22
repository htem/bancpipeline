#' banc-ngl-upload — Upload Neuroglancer-compatible meshes + segment_properties to GCS.
#'
#' For each cross-dataset target (BANC, FAFB, MANC, hemibrain, FANC,
#' maleCNS, mirror) uploads precomputed-mesh fragments to the public GCS
#' neuroglancer layer and refreshes the per-dataset `segment_properties`
#' JSON. Per-dataset uploads are flagged at the top (`upload.banc <-` etc.).
#' If `BANC_NGL_SKIP_MESH_UPLOADS=1`, only the segment_properties refresh
#' runs.
#'
#' @section Reads:
#'   - per-dataset BANC-space meshes under `<banc.nblast.*.obj.save.path>/`
#'   - SeaTable `banc_meta`
#'   - `franken_meta()`
#'   - env var `BANC_NGL_SKIP_MESH_UPLOADS`
#'
#' @section Writes:
#'   - GCS `gs://lee-lab_..._/precomputed/<dataset>/mesh/` (precomputed mesh fragments)
#'   - GCS `gs://lee-lab_..._/precomputed/<dataset>/info` (+ `segment_properties`)
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_ngl.sh`,
#'   `o2/production/o2_banc_ngl_upload_array.sh`

#############################################
### Upload Neuroglancer Compatible meshes ###
#############################################
source("banc/banc-startup.R")
redo <- FALSE
version <- banc.nblast.version
# Env-var override (set BANC_NGL_SKIP_MESH_UPLOADS=1 to refresh segment_properties only).
.skip_uploads <- nzchar(Sys.getenv("BANC_NGL_SKIP_MESH_UPLOADS", ""))
upload.hemibrain <- FALSE
upload.fafb <- FALSE
upload.manc <- FALSE
upload.mirror <- FALSE
upload.fanc <- FALSE
upload.malecns <- FALSE
upload.banc <- !.skip_uploads
options(scipen = 999)

###############################################################################
### Helpers: tag-rich segment_properties JSON + GCS push                    ###
###                                                                         ###
### build_segment_properties(meta, id_col, label_col, tag_cols)             ###
###   → JSON string in neuroglancer_segment_properties format with both     ###
###     {id="label", type="label"} and {id="tags", type="tags"} entries.    ###
###     Tags are collected as the union of unique non-NA values across      ###
###     `tag_cols`; comma-separated cells are split. Output `values` for    ###
###     the tags property is a list-of-int-vectors indexing into `tags`.    ###
###                                                                         ###
### push_segment_properties(json_str, gs_path, local_path = NULL)           ###
###   Writes local_path (creating dirs) then `gsutil cp` to gs_path. Use    ###
###   gs_path WITHOUT the "gs://" scheme guard — caller supplies full URL.  ###
###                                                                         ###
### ensure_layer_info(gs_root, template_url = ...)                          ###
###   Idempotently posts a precomputed `info` file at <gs_root>/info if     ###
###   absent, copied from a sibling layer. Skips silently on success.       ###
###                                                                         ###
### read_compiled_meta(dataset)                                             ###
###   Reads gs://.../compiled_data/<dataset>/<dataset>_meta.feather; for    ###
###   "fanc" falls back to local /n/.../banc/meta/fanc_meta.csv built by    ###
###   fanc/fanc-meta.R.                                                     ###
###############################################################################

build_segment_properties <- function(meta, id_col, label_col,
                                     tag_cols = character(0),
                                     na_label = "unknown") {
  stopifnot(id_col %in% colnames(meta))
  meta <- as.data.frame(meta)
  meta <- meta[!is.na(meta[[id_col]]) & nzchar(as.character(meta[[id_col]])) &
                 as.character(meta[[id_col]]) != "0", , drop = FALSE]
  meta <- meta[!duplicated(meta[[id_col]]), , drop = FALSE]
  ids <- as.character(meta[[id_col]])

  # Label property (per-row single string)
  if (label_col %in% colnames(meta)) {
    label_vals <- as.character(meta[[label_col]])
  } else {
    label_vals <- rep(NA_character_, nrow(meta))
  }
  label_vals[is.na(label_vals) | !nzchar(label_vals)] <- na_label

  # Tags property: collect per-row tag sets (comma-split, trimmed, deduped)
  present_tag_cols <- intersect(tag_cols, colnames(meta))
  row_tag_lists <- lapply(seq_len(nrow(meta)), function(i) {
    vals <- unlist(lapply(present_tag_cols, function(cc) {
      x <- meta[[cc]][i]
      if (is.na(x)) return(character(0))
      x <- as.character(x)
      if (!nzchar(x)) return(character(0))
      strsplit(x, ",", fixed = TRUE)[[1]]
    }), use.names = FALSE)
    vals <- trimws(vals)
    vals <- vals[nzchar(vals) & !vals %in% c("NA", "na")]
    unique(vals)
  })

  global_tags <- sort(unique(unlist(row_tag_lists, use.names = FALSE)))
  # I() wraps array-valued fields so jsonlite::toJSON(auto_unbox=TRUE) does NOT
  # collapse length-1 arrays to bare scalars. Neuroglancer's segment_properties
  # parser requires each `tags.values[i]` to be a JSON array; without I(), rows
  # whose neuron carries a single tag get serialised as `371` (a bare integer)
  # instead of `[371]`, which trips the parser and yields "No label property".
  if (length(global_tags) == 0) {
    tag_values <- replicate(length(ids), I(integer(0)), simplify = FALSE)
  } else {
    tag_values <- lapply(row_tag_lists, function(tt)
      I(as.integer(sort(unique(match(tt, global_tags))) - 1L)))
  }

  properties <- list(
    list(id = "label", type = "label", values = I(as.character(label_vals)))
  )
  if (length(global_tags) > 0) {
    properties[[2]] <- list(id = "tags", type = "tags",
                            tags = I(as.character(global_tags)),
                            values = tag_values)
  }

  json_data <- list(
    `@type` = "neuroglancer_segment_properties",
    inline = list(ids = I(as.character(ids)), properties = properties)
  )
  jsonlite::toJSON(json_data, auto_unbox = TRUE, pretty = TRUE)
}

push_segment_properties <- function(json_str, gs_path, local_path = NULL) {
  stopifnot(grepl("^gs://", gs_path))
  if (!is.null(local_path)) {
    dir.create(dirname(local_path), recursive = TRUE, showWarnings = FALSE)
    writeLines(json_str, local_path)
    src <- local_path
  } else {
    src <- tempfile(fileext = ".info")
    writeLines(json_str, src)
  }
  rc <- system2("gsutil", c("cp", src, gs_path), stdout = "", stderr = "")
  if (rc != 0) warning("gsutil cp failed for ", gs_path, " (rc=", rc, ")")
  if (is.null(local_path)) try(file.remove(src), silent = TRUE)
  invisible(rc == 0)
}

ensure_layer_info <- function(gs_root,
                              template_url = "gs://lee-lab_brain-and-nerve-cord-fly-connectome/imported_meshes/fafb_783_meshes_elastix_tpsreg_240721/info") {
  gs_root <- sub("/$", "", gs_root)
  exists <- system2("gsutil", c("-q", "stat", paste0(gs_root, "/info")),
                    stdout = NULL, stderr = NULL) == 0
  if (exists) return(invisible(TRUE))
  message("Posting layer info to ", gs_root, "/info (from template)")
  rc <- system2("gsutil", c("cp", template_url, paste0(gs_root, "/info")),
                stdout = "", stderr = "")
  invisible(rc == 0)
}

read_compiled_meta <- function(dataset) {
  # dataset ∈ {"banc_888","malecns_09","fafb_783","manc_121","hemibrain_121","fanc"}
  if (dataset == "fanc") {
    return(suppressWarnings(readr::read_csv(
      file = file.path(banc.meta.save.path, "fanc_meta.csv"),
      show_col_types = FALSE)))
  }
  bkt <- arrow::gs_bucket("lee-lab_brain-and-nerve-cord-fly-connectome",
                          anonymous = FALSE)
  arrow::read_feather(bkt$path(
    sprintf("compiled_data/%s/%s_meta.feather", dataset, dataset)))
}

# Tag column list (applied across datasets; super_cluster only exists in BANC)
BANC_TAG_COLS  <- c("super_cluster", "super_class", "cell_function",
                    "body_part_sensory", "body_part_effector", "side", "region",
                    "cell_class", "hemilineage", "nerve", "cell_sub_class")
OTHER_TAG_COLS <- setdiff(BANC_TAG_COLS, "super_cluster")

# Make a version folder for the swc data
banc.nblast.hemibrain.obj.save.path.version <- file.path(banc.nblast.hemibrain.obj.save.path,version)
banc.nblast.fafb.obj.save.path.version <- file.path(banc.nblast.fafb.obj.save.path,version)
banc.nblast.manc.obj.save.path.version <- file.path(banc.nblast.manc.obj.save.path,version)
banc.nblast.fanc.obj.save.path.version <- file.path(banc.nblast.fanc.obj.save.path,version)
malecns.version <- banc.nblast.malecns.version
banc.nblast.malecns.obj.save.path.version <- file.path(banc.nblast.malecns.obj.save.path,malecns.version)

# Mesh-upload meta loads. The expensive ones (banctable_query, franken_meta)
# are skipped when no upload flag needs them — segment_properties below uses
# read_compiled_meta() and doesn't need these. Honour the SKIP_MESH_UPLOADS
# env var so a standalone segment_properties refresh doesn't pay SeaTable cost.
.need_franken <- any(c(upload.banc, upload.fafb, upload.manc, upload.malecns))
if (.need_franken) franken.meta <- franken_meta()

if (upload.banc) {
  banc.meta <- banctable_query()
  banc.meta <- banc.meta %>%
    banc_filter_neurons() %>%
    dplyr::filter(proofread=="TRUE" | roughly_proofread=="TRUE",
                  !is.na(root_888)) %>%
    dplyr::distinct(root_888, .keep_all = TRUE)
  banc.ids <- na.omit(sample(unique(banc.meta$root_888)))

  # SLURM array sharding: when run as `sbatch --array=0-(N-1)`, each task
  # picks a disjoint id slice by `as.integer(tail6digits) %% N == task_id`.
  # Disjoint slices means no two workers ever write the same GCS object
  # or the same per-id .obj scratch file — no flock needed.
  .shard_id <- suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "")))
  .shard_n  <- suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_COUNT", "")))
  if (!is.na(.shard_id) && !is.na(.shard_n) && .shard_n > 1L) {
    .ids_chr <- as.character(banc.ids)
    .key     <- as.integer(substr(.ids_chr, nchar(.ids_chr) - 5L, nchar(.ids_chr)))
    banc.ids <- banc.ids[!is.na(.key) & (.key %% .shard_n) == .shard_id]
    message(sprintf("SLURM array shard %d/%d -> %d ids",
                    .shard_id, .shard_n, length(banc.ids)))
  }
}

if (upload.fafb) {
  fw.meta <- franken.meta %>%
    dplyr::filter(!is.na(fafb_id))
  fw.ids <- sample(unique(fw.meta$fafb_id))
}

if (upload.hemibrain) {
  hb.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"hemibrain_meta.csv"),
                                              col_types = hemibrainr:::sql_col_types))
  hb.ids <- sample(unique(hb.meta$bodyid))
}

if (upload.manc) {
  mc.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"manc_meta.csv"),
                                              col_types = hemibrainr:::sql_col_types))
  mc.meta <- franken.meta %>%
    dplyr::filter(!is.na(manc_id))
  mc.ids <- sample(unique(mc.meta$manc_id))
}

if (upload.fanc) {
  fc.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"fanc_meta.csv"),
                                              col_types = hemibrainr:::sql_col_types))
  fc.ids <- sample(unique(fc.meta$root_id))
}

if (upload.malecns) {
  mcns.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"malecns_09_meta.csv"),
                                                col_types = hemibrainr:::sql_col_types))
  mcns.ids <- sample(unique(mcns.meta$malecns_09_id))
}

# Upload banc meshes
if(upload.banc){
  # In SLURM array mode, skip the upfront `gsutil ls -r` snapshot: it pulls
  # 270k+ entries and 20 workers all doing it would hammer GCS for redundant
  # info. The per-id `gsutil stat` pre-flight inside the loop suffices.
  .in_array <- nzchar(Sys.getenv("SLURM_ARRAY_TASK_ID", ""))
  if(!redo && !.in_array){
    command <- "gsutil ls -r gs://lee-lab_brain-and-nerve-cord-fly-connectome/neuron_meshes/meshes"
    output <- system(command, intern = TRUE)
    filenames <- unique(gsub("\\:.*","",basename(output)))
    banc.ids <- setdiff(banc.ids,filenames)
    message("Working on missing banc IDs: ", length(banc.ids))
  }
  for(banc.id in banc.ids){
    try({
      # Per-mesh GCS pre-flight: the upfront `gsutil ls -r` snapshot is one-shot
      # and goes stale during long runs / concurrent uploads. A fast per-id
      # stat (~50 ms) closes the concurrent-rework + partial-upload gaps.
      # Trigger fragment is `<id>:0`; a present `:0` means this mesh has at
      # least its primary LOD on GCS and we should skip unless `redo`.
      if (!redo) {
        gcs_target <- sprintf(
          "gs://lee-lab_brain-and-nerve-cord-fly-connectome/neuron_meshes/meshes/%s:0",
          banc.id)
        if (system2("gsutil", c("-q", "stat", gcs_target),
                    stdout = FALSE, stderr = FALSE) == 0L) {
          message(sprintf("skip: mesh already on GCS: %s", banc.id))
          next
        }
      }
      volume_name <- banc.meta[banc.meta$root_888==banc.id,]$cell_type[[1]]
      if(is.na(volume_name)){
        volume_name="unknown"
      }
      message("cell type: ", volume_name)
      obj.file <- file.path(banc.obj.save.path,sprintf("%s.obj",banc.id))
      if(!file.exists(obj.file)){
        mesh3d.banc <- banc_read_neuron_meshes(banc.id)
        bancr:::write_mesh3d_to_obj(mesh3d.banc[[1]], filename = obj.file)
        message("Written: ", obj.file)
      }
      message("file: ", obj.file)
      if(file.exists(obj.file)){
        message("Uploading to google storage: ", volume_name)
        mesh3d <- readobj::read.obj(obj.file, convert.rgl = TRUE)[[1]]
        mesh3d <- Rvcg::vcgQEdecim(mesh3d, percent = 0.1)
        banc_upload_mesh(
          mesh = obj.file,
          mesh_id = as.character(banc.meta$root_888[match(banc.id,banc.meta$root_888)]),
          vol = "precomputed://gs://lee-lab_brain-and-nerve-cord-fly-connectome/neuron_meshes/",
          compress = TRUE,
          overwrite = redo
        )
        banc.volume <- data.frame(volume_number=banc.id, volume_name = volume_name)
        rem <- file.remove(obj.file)
      }
    })
  }
}

# Segment properties: tag-rich JSON pushed straight to GCS.
#
# BANC keeps two parallel sets under the same mesh folder so users can pick
# between v626's (existing, at gs://.../neuron_meshes/segment_properties/info)
# and v888's annotations via a neuroglancer state's `segmentPropertiesUrl`.
# v888 lives at segment_properties_v888/info; segment_properties/info is
# also refreshed as the default for v888 (the "current" properties).
#
# In SLURM array mode: only shard 0 pushes — otherwise 20 workers would race
# to overwrite the same `info` file.
.do_segprops <- !nzchar(Sys.getenv("SLURM_ARRAY_TASK_ID", "")) ||
                identical(Sys.getenv("SLURM_ARRAY_TASK_ID"), "0")
if (.do_segprops) {
  ensure_layer_info("gs://lee-lab_brain-and-nerve-cord-fly-connectome/neuron_meshes/")
  banc_meta_for_segprop <- read_compiled_meta("banc_888")
  banc_json <- build_segment_properties(banc_meta_for_segprop,
                                        id_col   = "root_888",
                                        label_col = "cell_type",
                                        tag_cols  = BANC_TAG_COLS)
  push_segment_properties(banc_json,
    gs_path = "gs://lee-lab_brain-and-nerve-cord-fly-connectome/neuron_meshes/segment_properties_v888/info",
    local_path = "setup/neuron_meshes/segment_properties_v888/info")
  push_segment_properties(banc_json,
    gs_path = "gs://lee-lab_brain-and-nerve-cord-fly-connectome/neuron_meshes/segment_properties/info",
    local_path = "setup/neuron_meshes/segment_properties/info")
}

# Upload FANC meshes
if(upload.fanc){
  if(!redo){
    command <- sprintf("gsutil ls -r gs://lee-lab_brain-and-nerve-cord-fly-connectome/imported_meshes/fanc_1116_meshes_%s/meshes", version)
    output <- system(command, intern = TRUE)
    filenames <- unique(gsub("\\:.*","",basename(output)))
    fc.cellids <- as.character(fc.meta$cell_id[match(fc.ids,fc.meta$root_id)])
    fc.ids <- fc.ids[!fc.cellids%in%filenames]
    message("Working on missing FANC IDs: ", length(fc.ids))
  }
  fc.volumes.df <- data.frame()
  for(fc.id in fc.ids){
    try({
      volume_name <- fc.meta[fc.meta$root_id==fc.id,]$cell_type[[1]]
      if(is.na(volume_name)){
        volume_name="unknown"
      }
      obj.file <- file.path(banc.nblast.fanc.obj.save.path.version,sprintf("%s.obj",fc.id))
      if(file.exists(obj.file)){
        message("Uploading to google storage: ", volume_name)
        #message("Uploading: ", obj.file)
        mesh3d <- readobj::read.obj(obj.file, convert.rgl = TRUE)[[1]]
        mesh3d <- Rvcg::vcgQEdecim(mesh3d, percent = 0.1)
        banc_upload_mesh(
          mesh = obj.file, # obj.file
          mesh_id = as.character(fc.meta$cell_id[match(fc.id,fc.meta$root_id)]),
          vol = sprintf("precomputed://gs://lee-lab_brain-and-nerve-cord-fly-connectome/imported_meshes/fanc_1116_meshes_%s/", version),
          compress = TRUE,
          overwrite = redo
        )
        fc.volume <- data.frame(volume_number=fc.id, volume_name = volume_name)
        fc.volumes.df <- rbind(fc.volumes.df,fc.volume)
      }
    })
  }
}

# FANC segment_properties — local meta (compiled_data feather not yet published)
ensure_layer_info(sprintf("gs://lee-lab_brain-and-nerve-cord-fly-connectome/imported_meshes/fanc_1116_meshes_%s/", version))
fanc_meta_for_segprop <- read_compiled_meta("fanc")
fanc_json <- build_segment_properties(fanc_meta_for_segprop,
                                      id_col   = "cell_id",
                                      label_col = "cell_type",
                                      tag_cols  = OTHER_TAG_COLS)
push_segment_properties(fanc_json,
  gs_path = sprintf("gs://lee-lab_brain-and-nerve-cord-fly-connectome/imported_meshes/fanc_1116_meshes_%s/segment_properties/info", version),
  local_path = sprintf("setup/imported_meshes/fanc_1116_meshes_%s/segment_properties/info", version))

# Upload FAFB-FlyWire meshes
if(upload.fafb){
  if(!redo){
    command <- sprintf("gsutil ls -r gs://lee-lab_brain-and-nerve-cord-fly-connectome/imported_meshes/fafb_783_meshes_%s/meshes", version)
    output <- system(command, intern = TRUE)
    filenames <- unique(gsub("\\:.*","",basename(output)))
    fw.ids <- setdiff(fw.ids,filenames)
    message("Working on missing FAFB IDs: ", length(fw.ids))
  }
  fw.volumes.df <- data.frame()
  for(fw.id in fw.ids){
    try({
      volume_name <- fw.meta[fw.meta$fafb_id==fw.id,]$cell_type[[1]]
      if(is.na(volume_name)){
        volume_name="unknown"
      }
      message("Uploading to google storage: ", volume_name)
      obj.file <- file.path(banc.nblast.fafb.obj.save.path.version,sprintf("%s.obj",fw.id))
      if(file.exists(obj.file)){
        message("Uploading: ", obj.file)
        mesh3d <- readobj::read.obj(obj.file, convert.rgl = TRUE)[[1]]
        mesh3d <- Rvcg::vcgQEdecim(mesh3d, percent = 0.1)
        mesh_id <- as.character(fw.id)
        banc_upload_mesh(
          mesh = obj.file, # obj.file
          mesh_id = as.integer64(fw.id),
          vol = sprintf("precomputed://gs://lee-lab_brain-and-nerve-cord-fly-connectome/imported_meshes/fafb_783_meshes_%s/", version),
          compress = TRUE,
          overwrite = redo
        )
        fw.volume <- data.frame(volume_number=fw.id, volume_name = volume_name)
        fw.volumes.df <- rbind(fw.volumes.df,fw.volume)
      }
    })
  }
}

# FAFB-FlyWire segment_properties — compiled_data feather (id col: fafb_783_id)
ensure_layer_info(sprintf("gs://lee-lab_brain-and-nerve-cord-fly-connectome/imported_meshes/fafb_783_meshes_%s/", version))
fafb_meta_for_segprop <- read_compiled_meta("fafb_783")
fafb_json <- build_segment_properties(fafb_meta_for_segprop,
                                      id_col   = "fafb_783_id",
                                      label_col = "cell_type",
                                      tag_cols  = OTHER_TAG_COLS)
push_segment_properties(fafb_json,
  gs_path = sprintf("gs://lee-lab_brain-and-nerve-cord-fly-connectome/imported_meshes/fafb_783_meshes_%s/segment_properties/info", version),
  local_path = sprintf("setup/imported_meshes/fafb_783_meshes_%s/segment_properties/info", version))

# Upload MANC meshes
if(upload.manc){
  if(!redo){
    command <- sprintf("gsutil ls -r gs://lee-lab_brain-and-nerve-cord-fly-connectome/imported_meshes/manc_v1.2.1_meshes_%s/meshes", version)
    output <- system(command, intern = TRUE)
    filenames <- unique(gsub("\\:.*","",basename(output)))
    mc.ids <- setdiff(mc.ids,filenames)
    message("Working on missing MANC IDs: ", length(mc.ids))
  }
  mc.volumes.df <- data.frame()
  for(mc.id in mc.ids){
    try({
      volume_name <- mc.meta[mc.meta$manc_id==mc.id,"type"][[1]]
      if(is.na(volume_name)){
        volume_name="unknown"
      }
      message("Uploading to google storage: ", volume_name)
      obj.file <- file.path(banc.nblast.manc.obj.save.path.version,sprintf("%s.obj",mc.id))
      if(file.exists(obj.file)){
        message("Uploading: ", obj.file)
        banc_upload_mesh(
          mesh = obj.file,
          mesh_id = as.character(mc.id),
          vol = sprintf("precomputed://gs://lee-lab_brain-and-nerve-cord-fly-connectome/imported_meshes/manc_v1.2.1_meshes_%s/", version),
          compress = TRUE,
          overwrite = redo
        )
        mc.volume <- data.frame(volume_number=mc.id, volume_name = volume_name)
        mc.volumes.df <- rbind(mc.volumes.df,mc.volume)
      }
    })
  }

  # Convert to .json for neuroglancer indexing:
  json_data <- list(
    `@type` = "neuroglancer_segment_properties",
    inline = list(
      ids = as.character(mc.volumes.df$volume_number),
      properties = list(
        list(
          id = "mesh",
          type = "label",
          values = mc.volumes.df$volume_name
        )
      )
    )
  )
}

# MANC segment_properties — compiled_data feather (id col: manc_121_id)
ensure_layer_info(sprintf("gs://lee-lab_brain-and-nerve-cord-fly-connectome/imported_meshes/manc_v1.2.1_meshes_%s/", version))
manc_meta_for_segprop <- read_compiled_meta("manc_121")
manc_json <- build_segment_properties(manc_meta_for_segprop,
                                      id_col   = "manc_121_id",
                                      label_col = "cell_type",
                                      tag_cols  = OTHER_TAG_COLS)
push_segment_properties(manc_json,
  gs_path = sprintf("gs://lee-lab_brain-and-nerve-cord-fly-connectome/imported_meshes/manc_v1.2.1_meshes_%s/segment_properties/info", version),
  local_path = sprintf("setup/imported_meshes/manc_v1.2.1_meshes_%s/segment_properties/info", version))

# Upload hemibrain meshes
if(upload.hemibrain){
  if(!redo){
    command <- sprintf("gsutil ls -r gs://lee-lab_brain-and-nerve-cord-fly-connectome/imported_meshes/hemibrain_v1.2.1_meshes_%s/meshes", version)
    output <- system(command, intern = TRUE)
    filenames <- unique(gsub("\\:.*","",basename(output)))
    hb.ids <- setdiff(hb.ids,filenames)
    message("Working on missing hemibrain IDs: ", length(hb.ids))
  }
  hb.volumes.df <- data.frame()
  for(hb.id in hb.ids){
    volume_name <- hb.meta[hb.meta$bodyid==hb.id,"name"][[1]]
    if(is.na(volume_name)){
      volume_name="fragment"
    }
    message("Uploading to google storage: ", volume_name)
    obj.file <- file.path(banc.nblast.hemibrain.obj.save.path.version,sprintf("%s.obj",hb.id))
    if(file.exists(obj.file)){
      message("Uploading: ", obj.file)
      mesh3d <- readobj::read.obj(obj.file, convert.rgl = TRUE)[[1]]
      mesh3d <- Rvcg::vcgQEdecim(mesh3d, percent = 0.1)
      mesh_id <- as.character(hb.id)
      try({
        banc_upload_mesh(
          mesh = mesh3d, # obj.file,
          mesh_id = mesh_id,
          vol = sprintf("precomputed://gs://lee-lab_brain-and-nerve-cord-fly-connectome/imported_meshes/hemibrain_v1.2.1_meshes_%s/", version),
          compress = TRUE,
          overwrite = redo
        )
      })
      if(upload.mirror){
        message("Mirroring ", mesh_id)
        try({
          mesh3d.m <- banc_mirror(mesh3d)
          banc_upload_mesh(
            mesh = mesh3d.m,
            mesh_id = mesh_id,
            vol = sprintf("precomputed://gs://lee-lab_brain-and-nerve-cord-fly-connectome/imported_meshes/mirrored_hemibrain_v1.2.1_meshes_%s/", version),
            compress = TRUE,
            overwrite = redo
          )
        })
      }
      hb.volume <- data.frame(volume_number=hb.id, volume_name = volume_name)
      hb.volumes.df <- rbind(hb.volumes.df,hb.volume)
    }
  }
}

# Hemibrain segment_properties — compiled_data feather (id col: hemibrain_121_id)
ensure_layer_info(sprintf("gs://lee-lab_brain-and-nerve-cord-fly-connectome/imported_meshes/hemibrain_v1.2.1_meshes_%s/", version))
hemi_meta_for_segprop <- read_compiled_meta("hemibrain_121")
hemi_json <- build_segment_properties(hemi_meta_for_segprop,
                                      id_col   = "hemibrain_121_id",
                                      label_col = "cell_type",
                                      tag_cols  = OTHER_TAG_COLS)
push_segment_properties(hemi_json,
  gs_path = sprintf("gs://lee-lab_brain-and-nerve-cord-fly-connectome/imported_meshes/hemibrain_v1.2.1_meshes_%s/segment_properties/info", version),
  local_path = sprintf("setup/imported_meshes/hemibrain_v1.2.1_meshes_%s/segment_properties/info", version))

# Upload maleCNS meshes
if(upload.malecns){
  if(!redo){
    command <- sprintf("gsutil ls -r gs://lee-lab_brain-and-nerve-cord-fly-connectome/imported_meshes/malecns_v0.9_meshes_%s/meshes", malecns.version)
    output <- system(command, intern = TRUE)
    filenames <- unique(gsub("\\:.*","",basename(output)))
    mcns.ids <- setdiff(mcns.ids,filenames)
    message("Working on missing maleCNS IDs: ", length(mcns.ids))
  }
  mcns.volumes.df <- data.frame()
  for(mcns.id in mcns.ids){
    try({
      volume_name <- mcns.meta[mcns.meta$malecns_09_id==mcns.id,]$cell_type[[1]]
      if(is.na(volume_name)){
        volume_name="unknown"
      }
      message("Uploading to google storage: ", volume_name)
      obj.file <- file.path(banc.nblast.malecns.obj.save.path.version,sprintf("%s.obj",mcns.id))
      if(file.exists(obj.file)){
        # Guard: skip faceless OBJs. trimesh inside bikinibottom returns a
        # PointCloud for vertex-only OBJs and crashes with
        # "'PointCloud' object has no attribute 'faces'".
        mesh3d.check <- tryCatch(readobj::read.obj(obj.file, convert.rgl = TRUE)[[1]],
                                 error = function(e) NULL)
        if (is.null(mesh3d.check) || is.null(mesh3d.check$it) || ncol(mesh3d.check$it) == 0) {
          message("Skipping ", mcns.id, ": OBJ has no faces")
          next
        }
        message("Uploading: ", obj.file)
        banc_upload_mesh(
          mesh = obj.file,
          mesh_id = as.character(mcns.id),
          vol = sprintf("precomputed://gs://lee-lab_brain-and-nerve-cord-fly-connectome/imported_meshes/malecns_v0.9_meshes_%s/", malecns.version),
          compress = TRUE,
          overwrite = redo
        )
        mcns.volume <- data.frame(volume_number=mcns.id, volume_name = volume_name)
        mcns.volumes.df <- rbind(mcns.volumes.df,mcns.volume)
      }
    })
  }
}

# MaleCNS segment_properties — compiled_data feather (id col: malecns_09_id)
ensure_layer_info(sprintf("gs://lee-lab_brain-and-nerve-cord-fly-connectome/imported_meshes/malecns_v0.9_meshes_%s/", malecns.version))
mcns_meta_for_segprop <- read_compiled_meta("malecns_09")
mcns_json <- build_segment_properties(mcns_meta_for_segprop,
                                      id_col   = "malecns_09_id",
                                      label_col = "cell_type",
                                      tag_cols  = OTHER_TAG_COLS)
push_segment_properties(mcns_json,
  gs_path = sprintf("gs://lee-lab_brain-and-nerve-cord-fly-connectome/imported_meshes/malecns_v0.9_meshes_%s/segment_properties/info", malecns.version),
  local_path = sprintf("setup/imported_meshes/malecns_v0.9_meshes_%s/segment_properties/info", malecns.version))

