#' banc-ids — Discover BANC neurons to process and rebuild `banc_ids.csv`.
#'
#' Pulls SeaTable `banc_meta`, joins to CAVE proofreading-notes / status /
#' nuclei tables, and rebuilds `<banc.meta.save.path>/banc_ids.csv`. This
#' CSV is the canonical "which neurons exist for this pipeline run"
#' inclusion-set consumed by the downstream metric / NBLAST scripts.
#'
#' @section Reads:
#'   - SeaTable `banc_meta`
#'   - CAVE (`banc_proofreading_notes()`, `banc_nuclei()`, status tables)
#'   - existing `<banc.meta.save.path>/banc_ids.csv` (for delta detection)
#'
#' @section Writes:
#'   - `<banc.meta.save.path>/banc_ids.csv`
#'   - `logs/banc-ids_<date>.txt`
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_metrics.sh`,
#'   `o2/production/o2_banc_update.sh`,
#'   `o2/production/o2_banc_v888_rebuild.sh`,
#'   `o2/oneshots/o2_banc_v850_rebuild.sh`,
#'   `o2/oneshots/o2_banc_v890_rebuild.sh`

#################################
### Choose neurons to process ###
#################################
source("banc/banc-startup.R")

# Set up logging
log_dir <- "logs"
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
.log_file <- file.path(log_dir, sprintf("banc-ids_%s.txt", format(Sys.time(), "%Y-%m-%d")))
.log_con <- file(.log_file, open = "wt")
sink(.log_con, type = "output", split = TRUE)
sink(.log_con, type = "message")
message(sprintf("Log started: %s", Sys.time()))

# Direct us to the BANC dataset
bancr::choose_banc()

# Get old IDs, if applicable
dir.create(banc.meta.save.path, showWarnings = FALSE, recursive = TRUE)
banc.ids.file <- file.path(banc.meta.save.path,"banc_ids.csv")
if(file.exists(banc.ids.file)){
  banc.ids.orig <-  readr::read_csv(file=file.path(banc.meta.save.path,"banc_ids.csv"),
                                    col_types = banc.col.types,
                                    show_col_types = FALSE)
} else {
  banc.ids.orig <- NULL
}

# Read BANC meta seatable (retry on transient CAVE errors)
for (.try in 1:3) {
  tryCatch({
    bancr:::banctable_updateids()
    break
  }, error = function(e) {
    if (.try < 3) {
      message(sprintf("  banctable_updateids attempt %d failed: %s\n  Retrying in 30s...",
                       .try, conditionMessage(e)))
      Sys.sleep(30)
    } else {
      stop(e)
    }
  })
}

# Single broad query that serves all sections below (saves ~6 API round-trips).
# Each section filters/selects from bc.all in memory instead of re-querying.
# SeaTable SELECT — keep ALL historical root_<ver> columns plus the current
# banc.version column. Fall back gracefully if the current-version column
# doesn't exist yet (e.g. immediately after a materialization bump).
.bc_all_sql_full <- sprintf("SELECT _id, status, root_id, root_626, root_850, root_888, root_%s, supervoxel_id, position, side, proofread, roughly_proofread, nucleus_id, nucleus_supervoxel_id, nucleus_position, nucleus_position_nm, root_position_nm, root_position, root_region, region, super_class, cell_class, cell_sub_class, cell_type, fafb_cell_type, malecns_cell_type, hemibrain_cell_type, manc_cell_type, fanc_cell_type, banc_match, banc_match_supervoxel_id, banc_png_match, banc_png_match_supervoxel_id, banc_nblast_match, banc_nblast_match_supervoxel_id from banc_meta",
                            banc.version)
.bc_all_sql_legacy <- "SELECT _id, status, root_id, root_626, root_850, root_888, supervoxel_id, position, side, proofread, roughly_proofread, nucleus_id, nucleus_supervoxel_id, nucleus_position, nucleus_position_nm, root_position_nm, root_position, root_region, region, super_class, cell_class, cell_sub_class, cell_type, fafb_cell_type, malecns_cell_type, hemibrain_cell_type, manc_cell_type, fanc_cell_type, banc_match, banc_match_supervoxel_id, banc_png_match, banc_png_match_supervoxel_id, banc_nblast_match, banc_nblast_match_supervoxel_id from banc_meta"
bc.all <- tryCatch(
  banctable_query(.bc_all_sql_full),
  error = function(e) {
    msg <- conditionMessage(e)
    if (grepl(sprintf("root_%s", banc.version), msg, fixed = TRUE)) {
      message(sprintf("  banctable_query: root_%s column not in SeaTable yet, retrying without it. Run banc-ids%s.R first to add the column, or live without versioned root for this cycle.",
                      banc.version, banc.version))
      banctable_query(.bc_all_sql_legacy)
    } else stop(e)
  }
)
bc <- bc.all %>%
  dplyr::select(`_id`, status, root_id, supervoxel_id, position) %>%
  dplyr::filter(!is.na(root_id),!is.na(supervoxel_id))

# Get meta data
message("##### Building list of nucleus IDs and bona fide neurons #####")

# Get nuclei positions
banc.nulcei <- banc_nuclei(rawcoords = TRUE, table = "both")

# Get cell information
banc.cell.info <- banc_cell_info(rawcoords = TRUE,live=TRUE)
banc.peripheral.nerves <- banc_peripheral_nerves(live=TRUE)
banc.backbone.proofread <- banc_backbone_proofread(live=TRUE)
# banc_proofreading_notes() pipes the CAVE response through
# `mutate(across(ends_with("position"), safe_raw2nm_position))`. Both
# live=TRUE (initial bug, 2026-04-20) and live=FALSE (regression observed
# 2026-04-27) hit a vctrs cast failure when CAVE returns pt_position as a
# vctrs_list_of with mixed character/integer entries. Downstream we only
# use pt_position via paste0(., collapse=","), so we don't need the nm
# conversion at all. Try the bancr wrapper first; on failure, monkey-
# patch safe_raw2nm_position to identity (the caller flattens to string)
# and retry. The patch is reverted on exit so the rest of the script sees
# stock bancr.
banc.proofreading.notes <- local({
  attempt <- function(live_arg) {
    tryCatch(banc_proofreading_notes(live = live_arg),
             error = function(e) { message(sprintf("  banc_proofreading_notes(live=%s) failed: %s",
                                                    live_arg, conditionMessage(e))); NULL })
  }
  out <- attempt(FALSE)
  if (is.null(out)) out <- attempt(2)
  if (is.null(out)) {
    # Fallback: bypass safe_raw2nm_position. Replace with a no-op while we
    # call banc_proofreading_notes; pt_position survives as raw list/string
    # which paste0(., collapse=",") on line 112 handles either way.
    if (requireNamespace("bancr", quietly = TRUE) &&
        exists("safe_raw2nm_position", envir = asNamespace("bancr"))) {
      ns <- asNamespace("bancr")
      orig_fn <- get("safe_raw2nm_position", envir = ns)
      noop_fn <- function(x) x
      assignInNamespace("safe_raw2nm_position", noop_fn, ns = "bancr")
      on.exit(assignInNamespace("safe_raw2nm_position", orig_fn, ns = "bancr"),
              add = TRUE)
      out <- attempt(FALSE)
    }
  }
  if (is.null(out)) {
    warning("banc_proofreading_notes() unavailable; proceeding without notes-derived statuses.")
    out <- tibble::tibble(pt_root_id = character(),
                          pt_supervoxel_id = character(),
                          pt_position = list(),
                          tag = character())
  }
  out
})
banc.cell.ids <- banc_cell_ids()

# Get neck connective neurons
banc.neck.connective.neurons1 <- banc_neck_connective_neurons(table="neck_connective_y92500")
banc.neck.connective.neurons2 <- banc_neck_connective_neurons(table="neck_connective_y121000")
banc.neck.connective.neurons2 <- subset(banc.neck.connective.neurons2, 
                                        !banc.neck.connective.neurons2$pt_root_id%in%banc.neck.connective.neurons1$pt_root_id)
banc.neck.connective.neurons <- rbind.fill(banc.neck.connective.neurons1,banc.neck.connective.neurons2)

# Prepare all data frames
banc.nulcei.mod <- banc.nulcei %>%
  dplyr::mutate(nucleus_position_nm = nat::xyzmatrix2str(bancr::banc_raw2nm(nucleus_position),
                                                          format = "%.0f, %.0f, %.0f")) %>%
  dplyr::mutate(nucleus_id = as.character(nucleus_id)) %>%
  dplyr::ungroup() %>%
  dplyr::rename(nucleus_supervoxel_id = pt_supervoxel_id) %>%
  dplyr::select(nucleus_id, nucleus_supervoxel_id, root_id, nucleus_position, nucleus_position_nm) %>%
  dplyr::distinct(nucleus_id, .keep_all = TRUE) %>%
  dplyr::distinct(root_id, .keep_all = TRUE) %>%
  dplyr::select(root_id, nucleus_id, nucleus_supervoxel_id, nucleus_position, nucleus_position_nm)

############################################################
## SURGICAL DRIFT FIX (opt-in, gated on BANC_NUCLEUS_DRIFT_FIX)
##
## When BANC_NUCLEUS_DRIFT_FIX=TRUE, we (a) detect rows whose
## (nucleus_id, root_id) pair disagrees with CAVE's current
## banc_nuclei() table, (b) wipe nucleus_*, side from those rows
## directly in SeaTable so the existing CAVE-nuclei→banc.new join
## further down has a clean slate to re-attach the correct
## nucleus on each row's CURRENT root_id.
##
## Default OFF — first read the drift CSV report (written at the
## bottom of this script) to see the magnitude before flipping the
## flag.
############################################################
.banc_drift_fix_on <- isTRUE(as.logical(
  Sys.getenv("BANC_NUCLEUS_DRIFT_FIX", "FALSE")))
if (.banc_drift_fix_on) {
  message("\n##### Surgical drift fix: ENABLED #####")
  cave_nuc_for_fix <- banc.nulcei.mod %>%
    dplyr::transmute(nucleus_id   = as.character(nucleus_id),
                     cave_root_id = as.character(root_id))
  bc.live.fix <- bancr:::banctable_query(
    "SELECT _id, root_id, nucleus_id from banc_meta") %>%
    dplyr::filter(!is.na(nucleus_id), nucleus_id != "",
                  !is.na(root_id),    root_id    != "", root_id != "0") %>%
    dplyr::mutate(nucleus_id = as.character(nucleus_id),
                  root_id    = as.character(root_id))
  drift.fix <- bc.live.fix %>%
    dplyr::inner_join(cave_nuc_for_fix, by = "nucleus_id") %>%
    dplyr::filter(root_id != cave_root_id)
  if (nrow(drift.fix) > 0L) {
    message(sprintf("  Wiping nucleus_*, side from %d drifted SeaTable rows",
                    nrow(drift.fix)))
    wipe.df <- drift.fix %>%
      dplyr::transmute(`_id`                 = `_id`,
                       nucleus_id            = "",
                       nucleus_supervoxel_id = "",
                       nucleus_position      = "",
                       nucleus_position_nm   = "",
                       side                  = "")
    try(banctable_update_rows(base = 'banc_meta',
                              table = "banc_meta",
                              df = as.data.frame(wipe.df),
                              append_allowed = FALSE,
                              chunksize = 1000))
    message("  Wipe complete. Existing CAVE-nuclei join below will re-attach correctly.")
  } else {
    message("  No drift detected — nothing to wipe.")
  }
  rm(cave_nuc_for_fix, bc.live.fix, drift.fix)
} else {
  message("\nSurgical drift fix: DISABLED (set BANC_NUCLEUS_DRIFT_FIX=TRUE to enable)")
}

# Determine whether neuron is on the left or right of BANC
nuclei.points <- nat::xyzmatrix(banc.nulcei.mod$nucleus_position_nm)
lrdiffs <- bancr:::banc_lr_position(nuclei.points,units = "nm")
sides <- ifelse(lrdiffs>0,"right","left")
banc.nulcei.mod$side <- sides

# Process other cell info. NOTE: we no longer prefer the cell_info `side`
# tag over the geometric side computed at lines 180-183 above. cell_info
# is community-contributed and not curator-vetted (controlled vocabulary
# NOT enforced), so geometry alone is the load-bearing source for side.
banc.cell.info.mod <- banc_cave_cell_types()
banc.nulcei.mod <- banc.nulcei.mod %>%
  dplyr::distinct(root_id, nucleus_id, .keep_all = TRUE)
banc.cell.info.mod$side <- NULL

# Process proofreading notes
banc.proofreading.notes.mod <- banc.proofreading.notes %>%
  dplyr::filter(tag %in% c("thoroughly proofread","roughly proofread","arbor is damaged","soma is damaged", "large fragment"),
                !pt_root_id %in% !!banc.backbone.proofread$pt_root_id) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(pt_position = paste0(pt_position,collapse=",")) %>%
  dplyr::ungroup() %>%
  dplyr::rename(root_id = pt_root_id, supervoxel_id = pt_supervoxel_id, position = pt_position, status = tag) %>%
  dplyr::select(root_id, supervoxel_id, position, status) %>%
  dplyr::mutate(status = snakecase::to_screaming_snake_case(status)) %>%
  dplyr::arrange(status) %>%
  dplyr::distinct(root_id, supervoxel_id, .keep_all = TRUE)

# Neck connective neurons
banc.neck.connective.neurons.mod <- banc.neck.connective.neurons %>%
  dplyr::rowwise() %>%
  dplyr::mutate(pt_position = paste0(pt_position,collapse=",")) %>%
  dplyr::ungroup() %>%
  dplyr::rename(root_id = pt_root_id, supervoxel_id = pt_supervoxel_id, position = pt_position, neck_connective = tag) %>%
  dplyr::select(root_id, supervoxel_id, position, neck_connective) %>%
  dplyr::arrange(neck_connective) %>%
  dplyr::distinct(root_id, supervoxel_id, .keep_all = TRUE)

# Peripheral nerves
banc.peripheral.nerves.mod <- banc.peripheral.nerves %>%
  dplyr::rowwise() %>%
  dplyr::mutate(pt_position = paste0(pt_position,collapse=",")) %>%
  dplyr::ungroup() %>%
  dplyr::rename(root_id = pt_root_id, supervoxel_id = pt_supervoxel_id, position = pt_position, nerve = tag) %>%
  dplyr::select(root_id, supervoxel_id, position, nerve) %>%
  dplyr::arrange(nerve) %>%
  dplyr::distinct(root_id, supervoxel_id, .keep_all = TRUE)

# Mode the proofread table
banc.backbone.proofread.mod <- banc.backbone.proofread %>%
  dplyr::rowwise() %>%
  dplyr::mutate(pt_position = paste0(pt_position,collapse=",")) %>%
  dplyr::ungroup() %>%
  dplyr::rename(root_id = pt_root_id, supervoxel_id = pt_supervoxel_id, position = pt_position) %>%
  dplyr::distinct(root_id, supervoxel_id, position) %>%
  dplyr::mutate(proofread = TRUE)

# Combine
banc.cave.ids <- c(banc.nulcei.mod$root_id, 
                   banc.neck.connective.neurons.mod$root_id,
                   banc.proofreading.notes.mod$root_id,
                   banc.peripheral.nerves.mod$root_id,
                   banc.backbone.proofread.mod$root_id,
                   banc.cell.info.mod$root_id)
banc.ids <- banc.nulcei.mod %>%
  dplyr::full_join(banc.cell.info.mod, by = c("root_id")) %>% 
  dplyr::mutate(supervoxel_id = dplyr::case_when(
    is.na(supervoxel_id) ~ nucleus_supervoxel_id,
    TRUE ~ supervoxel_id
  )) %>%
  dplyr::mutate(position = dplyr::case_when(
    is.na(position) ~ nucleus_position,
    TRUE ~ position
  )) %>%
  dplyr::full_join(banc.neck.connective.neurons.mod, 
                   by = c("root_id","supervoxel_id","position"),
                   relationship = "many-to-many") %>%
  dplyr::full_join(banc.proofreading.notes.mod, 
                   by = c("root_id","supervoxel_id","position"),
                   relationship = "many-to-many") %>%
  dplyr::full_join(banc.peripheral.nerves.mod, 
                   by = c("root_id","supervoxel_id","position"),
                   relationship = "many-to-many") %>%
  dplyr::left_join(banc.backbone.proofread.mod[,c("root_id","proofread")], 
                   by = c("root_id"),
                   relationship = "many-to-many") %>%
  dplyr::rename(valid = proofread) %>%
  dplyr::full_join(banc.backbone.proofread.mod, 
                   by = c("root_id","supervoxel_id","position"),
                   relationship = "many-to-many") %>%
  dplyr::mutate(position_nm = nat::xyzmatrix2str(banc_raw2nm(position),
                                                  format = "%.0f, %.0f, %.0f"),
                proofread = dplyr::case_when(
                  valid|proofread ~ TRUE,
                  TRUE ~ FALSE
                )) %>%
  dplyr::select(-valid) %>%
  dplyr::mutate(cell_class = dplyr::case_when(
    is.na(cell_class) & root_id %in% banc.peripheral.nerves.mod$root_id ~ "sensory neuron", 
    TRUE ~ cell_class
  )) %>%
  dplyr::mutate(super_class = dplyr::case_when(
    is.na(super_class) & root_id %in% banc.peripheral.nerves.mod$root_id ~ "afferent",
    TRUE ~ super_class
  )) %>%
  # Backfill glia/trachea super_class from CAVE cell_info for neurons
  # whose seatable super_class is blank (new neurons not yet classified)
  dplyr::mutate(cave_super_class = super_class) %>%
  dplyr::arrange(super_class, cell_class, cell_type, side, proofread, nucleus_id, nucleus_position_nm, position) %>%
  dplyr::filter(proofread | root_id %in% banc.peripheral.nerves.mod$root_id | !is.na(cell_class)) %>%
  dplyr::distinct(root_id, supervoxel_id, .keep_all = TRUE) %>%
  dplyr::select(-super_class)

# Make root_id distinct, but prefer those present in banctable
result <- banc.ids %>%
  dplyr::filter(!is.na(root_id)) %>%
  dplyr::mutate(root_id = as.character(root_id),
                supervoxel_id = as.character(supervoxel_id)) %>%
  dplyr::left_join(bc[,c("root_id","supervoxel_id")] %>% 
                     dplyr::distinct(root_id, supervoxel_id), 
                   by = c("root_id")) %>%
  dplyr::group_by(root_id) %>%
  dplyr::slice(which.max(!is.na(supervoxel_id.y)) %||% 1) %>%
  dplyr::ungroup() %>%
  dplyr::rename(supervoxel_id = supervoxel_id.x) %>%
  dplyr::select(-supervoxel_id.y) %>%
  dplyr::distinct(root_id, .keep_all = TRUE) %>%
  dplyr::mutate(voxels = NA)

# Add anything else tracked in seatable
bc <- bc.all %>%
  dplyr::select(root_id, supervoxel_id, position, side)
bc <- dplyr::anti_join(bc, result, by="root_id")
if(nrow(bc)){
  bc$proofread <- FALSE
  result <- plyr::rbind.fill(result,bc)
}

# # Get old IDs for voxel count
# banc.ids.old <- readr::read_csv(file=file.path(banc.meta.save.path,"banc_ids.csv"),
#                             col_types = banc.col.types, 
#                             show_col_types = FALSE)
# if(nrow(banc.ids.old)){
#   result <- result %>%
#     dplyr::left_join(result, banc.ids.old[,c("root_id","voxels")], by = "root_id")
# }
# 
# # Update voxel counts and then re-write
# result <- result %>%
#   dplyr::filter(!is.na(root_id)) %>%
#   dplyr::rowwise() %>%
#   dplyr::mutate(voxels = dplyr::case_when(
#     is.na(voxels) ~ length(tryCatch(banc_leaves(root_id), error = function(e) NA)),
#     TRUE ~ voxels)) %>%
#   dplyr::ungroup()

# Neaten missing data
result <- as.data.frame(result)
result[is.na(result)] <- ""
result[result=="NA,NA,NA"] <- ""
result[result=="0"] <- ""

# Save as .csv file
readr::write_csv(result, file=file.path(banc.meta.save.path,"banc_ids.csv"))

# What is new, what is old
ids <- nrow(result)
new.ids <- result$root_id[!result$root_id %in% banc.ids.orig$root_id]

# Announce
message("##### BANCpipeline: banc root ids updated #####")
message(sprintf("##### new IDs: %s", length(new.ids)))

###################
## Fix version ####
###################

# Add root_id 626
bc <- bc.all %>%
  dplyr::select(`_id`, root_626, root_id, supervoxel_id, position) %>%
  dplyr::filter(is.na(root_626)|root_626==""|root_626=="0"|root_626==0|is.na(supervoxel_id)|supervoxel_id=="0"|supervoxel_id==0)
if(nrow(bc)){
  bc$supervoxel_id <- bancr::banc_xyz2id(bc$position,root=FALSE,rawcoords = TRUE)
  bc$root_626 <- bancr::banc_rootid(bc$supervoxel_id, version = "626")
  if(nrow(bc)){
    bc[is.na(bc)] <- ''
    banctable_update_rows(base='banc_meta',
                          table = "banc_meta",
                          df = bc[,c("_id","root_id", "root_626","supervoxel_id")],
                          append_allowed = FALSE,
                          chunksize = 1000)
  }
}

# Add root_id 850
bc850 <- bc.all %>%
  dplyr::select(`_id`, root_850, root_id, supervoxel_id, position) %>%
  dplyr::filter(is.na(root_850)|root_850==""|root_850=="0"|root_850==0|is.na(supervoxel_id)|supervoxel_id=="0"|supervoxel_id==0)
if(nrow(bc850)){
  bc850$supervoxel_id <- bancr::banc_xyz2id(bc850$position,root=FALSE,rawcoords = TRUE)
  bc850$root_850 <- bancr::banc_rootid(bc850$supervoxel_id, version = "850")
  if(nrow(bc850)){
    bc850[is.na(bc850)] <- ''
    banctable_update_rows(base='banc_meta',
                          table = "banc_meta",
                          df = bc850[,c("_id","root_id", "root_850","supervoxel_id")],
                          append_allowed = FALSE,
                          chunksize = 1000)
  }
}

# Add root_id 888
bc888 <- bc.all %>%
  dplyr::select(`_id`, root_888, root_id, supervoxel_id, position) %>%
  dplyr::filter(is.na(root_888)|root_888==""|root_888=="0"|root_888==0|is.na(supervoxel_id)|supervoxel_id=="0"|supervoxel_id==0)
if(nrow(bc888)){
  bc888$supervoxel_id <- bancr::banc_xyz2id(bc888$position,root=FALSE,rawcoords = TRUE)
  bc888$root_888 <- bancr::banc_rootid(bc888$supervoxel_id, version = "888")
  if(nrow(bc888)){
    bc888[is.na(bc888)] <- ''
    banctable_update_rows(base='banc_meta',
                          table = "banc_meta",
                          df = bc888[,c("_id","root_id", "root_888")],
                          append_allowed = FALSE,
                          chunksize = 1000)
  }
}

# Add root_id for the CURRENT banc.version (parallel pattern, parametrized).
# No-op when banc.version == "888" (covered above) or when the column wasn't
# fetched into bc.all (legacy SeaTable, before column existed).
.root.col <- paste0("root_", banc.version)
if (.root.col %in% names(bc.all) && banc.version != "888") {
  bc.ver <- bc.all %>%
    dplyr::select(`_id`, !!rlang::sym(.root.col), root_id, supervoxel_id, position) %>%
    dplyr::filter(
      is.na(.data[[.root.col]]) | .data[[.root.col]] == "" |
      .data[[.root.col]] == "0" | .data[[.root.col]] == 0 |
      is.na(supervoxel_id) | supervoxel_id == "0" | supervoxel_id == 0
    )
  message(sprintf("  rows needing %s backfill: %s",
                  .root.col, format(nrow(bc.ver), big.mark = ",")))
  if (nrow(bc.ver)) {
    bc.ver$supervoxel_id <- bancr::banc_xyz2id(bc.ver$position,
                                               root = FALSE, rawcoords = TRUE)
    bc.ver[[.root.col]] <- bancr::banc_rootid(bc.ver$supervoxel_id,
                                              version = banc.version)
    bc.ver[is.na(bc.ver)] <- ''
    if (nrow(bc.ver)) {
      banctable_update_rows(base = 'banc_meta',
                            table = "banc_meta",
                            df = bc.ver[, c("_id", "root_id", .root.col)],
                            append_allowed = FALSE,
                            chunksize = 1000)
    }
  }
} else if (!.root.col %in% names(bc.all) && banc.version != "888") {
  message(sprintf("  %s column not in SeaTable yet — run banc-ids%s.R after adding the column",
                  .root.col, banc.version))
}

##################################
## Fix supervoxel_id/position ####
##################################

# Update
bc <- bc.all %>%
  dplyr::select(`_id`, root_id, supervoxel_id, position)
bc.new <- bc %>%
  dplyr::filter(is.na(supervoxel_id)|is.na(position)) %>%
  dplyr::select(-supervoxel_id,-position) %>%
  dplyr::left_join(banc.ids %>%
                     dplyr::mutate(root_id = as.character(root_id),
                                   supervoxel_id = as.character(supervoxel_id)) %>%
                     dplyr::distinct(root_id, .keep_all = TRUE) %>%
                     dplyr::select(root_id, supervoxel_id, position),
                   by = "root_id") %>%
  dplyr::filter(!is.na(supervoxel_id),!is.na(position)) %>%
  as.data.frame() 

# Update
if(nrow(bc.new)){
  bc.new[is.na(bc.new)] <- ''
  banctable_update_rows(base='banc_meta', 
                        table = "banc_meta", 
                        df = bc.new, 
                        append_allowed = FALSE, 
                        chunksize = 1000)
}

########################
## Update proofread ####
########################

bc <- bc.all %>%
  dplyr::select(`_id`, root_id, proofread, roughly_proofread)
proofread.ids <- unique(as.character(
  subset(banc.backbone.proofread, valid == TRUE | valid == "t")$pt_root_id
))
roughly.proofread.ids <-  banc.proofreading.notes %>%
  dplyr::filter(tag %in% c("roughly proofread", "arbor is damaged", "large fragment")) %>%
  dplyr::mutate(pt_root_id = as.character(pt_root_id)) %>%
  dplyr::pull(pt_root_id) %>%
  unique()

# Coverage diagnostic: flag CAVE root_ids that have no SeaTable row to flag.
# These typically come from new proofreading work whose neuron isn't tracked
# in banc_meta yet (e.g., no nucleus, no L2 metrics, not yet in cell_info).
# Whether they should also be appended is decided downstream in
# banc-update-seatable.R.
.bc_roots <- as.character(bc$root_id)
.bp_missing <- setdiff(proofread.ids, .bc_roots)
.rp_missing <- setdiff(roughly.proofread.ids, .bc_roots)
message(sprintf(
  "  proofread          : CAVE has %d, %d match SeaTable rows, %d CAVE root_ids absent from SeaTable",
  length(proofread.ids),
  sum(.bc_roots %in% proofread.ids),
  length(.bp_missing)))
message(sprintf(
  "  roughly_proofread  : notes have %d, %d match SeaTable rows, %d notes root_ids absent from SeaTable",
  length(roughly.proofread.ids),
  sum(.bc_roots %in% roughly.proofread.ids),
  length(.rp_missing)))

# Update logic: only flip a value when the CAVE evidence is unambiguous.
#   proofread          - flip TRUE when root_id is in CAVE backbone_proofread.
#                        Otherwise keep whatever SeaTable already has (so a
#                        join-time miss in a previous run can't silently wipe
#                        a curator-asserted TRUE). Default brand-new blanks
#                        to "FALSE" so the column is well-formed.
#   roughly_proofread  - flip FALSE when proofread (the strong claim wins),
#                        flip TRUE when the proofreading_notes tag asserts it,
#                        otherwise preserve existing. Blanks default FALSE.
.is_blank <- function(x) is.na(x) | x == "" | toupper(as.character(x)) == "NA"
bc.new <- bc %>%
  dplyr::mutate(proofread = dplyr::case_when(
    root_id %in% proofread.ids ~ "TRUE",
    .is_blank(proofread) ~ "FALSE",
    TRUE ~ as.character(proofread)
  )) %>%
  dplyr::mutate(roughly_proofread = dplyr::case_when(
    root_id %in% proofread.ids ~ "FALSE",
    root_id %in% roughly.proofread.ids ~ "TRUE",
    .is_blank(roughly_proofread) ~ "FALSE",
    TRUE ~ as.character(roughly_proofread)
  ))

message(sprintf(
  "  proofread TRUE     : %d -> %d (delta %+d)",
  sum(bc$proofread %in% c("TRUE", TRUE), na.rm = TRUE),
  sum(bc.new$proofread == "TRUE", na.rm = TRUE),
  sum(bc.new$proofread == "TRUE", na.rm = TRUE) -
    sum(bc$proofread %in% c("TRUE", TRUE), na.rm = TRUE)))
message(sprintf(
  "  roughly_proofread T: %d -> %d (delta %+d)",
  sum(bc$roughly_proofread %in% c("TRUE", TRUE), na.rm = TRUE),
  sum(bc.new$roughly_proofread == "TRUE", na.rm = TRUE),
  sum(bc.new$roughly_proofread == "TRUE", na.rm = TRUE) -
    sum(bc$roughly_proofread %in% c("TRUE", TRUE), na.rm = TRUE)))

# Update
if(nrow(bc.new)){
  bc.new[is.na(bc.new)] <- ''
  banctable_update_rows(base='banc_meta', 
                        table = "banc_meta", 
                        df = bc.new, 
                        append_allowed = FALSE, 
                        chunksize = 1000)
}

#########################
## Suggest deletions ####
#########################

bc.ids <- unique(banc.ids$root_id)
info_cols <- c("super_class", "cell_class", "cell_sub_class", "cell_type",
               "fafb_cell_type", "malecns_cell_type", "hemibrain_cell_type",
               "manc_cell_type", "fanc_cell_type")
bc <- bc.all %>%
  dplyr::select(`_id`, status, root_id, supervoxel_id, position,
                dplyr::all_of(info_cols))
duplicates <- unique(bc$root_id[duplicated(bc$root_id)])
bc.new <- bc %>%
  # Score each row by how many info columns are non-empty
  dplyr::rowwise() %>%
  dplyr::mutate(info_score = sum(!is.na(dplyr::c_across(dplyr::all_of(info_cols))) &
                                   dplyr::c_across(dplyr::all_of(info_cols)) != "")) %>%
  dplyr::ungroup() %>%
  # Flag neurons missing from feeder CAVE tables
  dplyr::mutate(CAVE_MISSING = !root_id %in% banc.cave.ids | grepl("CAVE_MISSING", status)) %>%
  dplyr::mutate(DUPLICATED = root_id %in% duplicates) %>%
  dplyr::mutate(BAD = grepl("DELETE", status)) %>%
  dplyr::filter(DUPLICATED | CAVE_MISSING | BAD | grepl("DUPLICATED|CAVE_MISSING", status)) %>%
  # For duplicates, mark the row with least info as the one to flag
  dplyr::group_by(root_id) %>%
  dplyr::mutate(LEAST_INFO = DUPLICATED & info_score == min(info_score) &
                  dplyr::row_number(info_score) == 1) %>%
  dplyr::ungroup() %>%
  # Strip old flags then re-apply
  dplyr::rowwise() %>%
  dplyr::mutate(status = ifelse(grepl("DELETE", status),
                                subtract_status(status, "DELETE"),
                                status)) %>%
  dplyr::mutate(status = ifelse(grepl("DUPLICATED", status),
                                subtract_status(status, "DUPLICATED"),
                                status)) %>%
  dplyr::mutate(status = ifelse(grepl("CAVE_MISSING", status),
                                subtract_status(status, "CAVE_MISSING"),
                                status)) %>%
  dplyr::mutate(status = dplyr::case_when(
    !root_id %in% banc.cave.ids ~ append_status(status, "CAVE_MISSING"),
    LEAST_INFO ~ append_status(status, "DUPLICATED"),
    TRUE ~ status
  )) %>%
  dplyr::ungroup() %>%
  dplyr::select(-dplyr::all_of(info_cols), -info_score, -CAVE_MISSING, -DUPLICATED, -BAD, -LEAST_INFO) %>%
  dplyr::arrange(root_id, supervoxel_id, position, status) %>%
  as.data.frame()

# Update
if(nrow(bc.new)){
  try(banctable_update_rows(base='banc_meta',
                            table = "banc_meta",
                            df = bc.new[,c("_id","status")],
                            append_allowed = FALSE,
                            chunksize = 1000))
}

#####################
## Update nuclei ####
#####################

# Determine the neuropil region for the ID
bc <- bc.all %>%
  dplyr::select(`_id`, root_id, nucleus_id, nucleus_supervoxel_id,
                nucleus_position, nucleus_position_nm,
                root_position_nm, root_position, root_region, region, side)
banc.nulcei.mod$nucleus_supervoxel_id <- as.character(banc.nulcei.mod$nucleus_supervoxel_id)
bc.new <- bc %>%
  dplyr::left_join(banc.nulcei.mod %>%
                     dplyr::mutate(root_id=as.character(root_id)) %>%
                     dplyr::distinct(root_id, .keep_all = TRUE) %>%
                     dplyr::select(root_id,nucleus_id,nucleus_supervoxel_id,nucleus_position,nucleus_position_nm),
                   by = "root_id") %>%
  dplyr::mutate(nucleus_position = dplyr::coalesce(nucleus_position.y, nucleus_position.x)) %>%
  dplyr::mutate(nucleus_position_nm = dplyr::coalesce(nucleus_position_nm.y, nucleus_position_nm.x)) %>%
  dplyr::mutate(nucleus_id = dplyr::coalesce(nucleus_id.y, nucleus_id.x)) %>%
  dplyr::mutate(nucleus_supervoxel_id = dplyr::coalesce(nucleus_supervoxel_id.y, nucleus_supervoxel_id.x)) %>%
  dplyr::mutate(root_position = dplyr::case_when(
    (is.na(root_position) | root_position == "") & !is.na(nucleus_position) & nucleus_position != "" ~ nucleus_position,
    TRUE ~ root_position
  )) %>%
  dplyr::mutate(root_position_nm = dplyr::case_when(
    is.na(root_position_nm) ~ nucleus_position_nm,
    TRUE ~ root_position_nm
  )) %>%
  dplyr::mutate(nucleus_supervoxel_id = dplyr::case_when(
    is.na(nucleus_id) ~ NA,
    TRUE ~ nucleus_supervoxel_id
  )) %>%
  dplyr::select(`_id`, root_id, nucleus_id, nucleus_supervoxel_id,
                nucleus_position, nucleus_position_nm, root_position, root_position_nm, root_region, region) %>%
  as.data.frame()

# Guard: nucleus_supervoxel_id must be unique (each nucleus has one supervoxel)
nsv_vals <- bc.new$nucleus_supervoxel_id[!is.na(bc.new$nucleus_supervoxel_id) & bc.new$nucleus_supervoxel_id != ""]
dup_nsv <- unique(nsv_vals[duplicated(nsv_vals)])
if (length(dup_nsv) > 0) {
  # Keep only the row whose root_id matches CAVE nuclei data; clear the rest
  cave_owners <- banc.nulcei.mod %>%
    dplyr::mutate(root_id = as.character(root_id),
                  nucleus_supervoxel_id = as.character(nucleus_supervoxel_id)) %>%
    dplyr::filter(nucleus_supervoxel_id %in% dup_nsv) %>%
    dplyr::distinct(nucleus_supervoxel_id, .keep_all = TRUE) %>%
    dplyr::select(root_id, nucleus_supervoxel_id)
  is_dup <- bc.new$nucleus_supervoxel_id %in% dup_nsv & !is.na(bc.new$nucleus_supervoxel_id)
  is_owner <- paste(bc.new$root_id, bc.new$nucleus_supervoxel_id) %in%
    paste(cave_owners$root_id, cave_owners$nucleus_supervoxel_id)
  n_clear <- sum(is_dup & !is_owner)
  if (n_clear > 0) {
    message(sprintf("  Clearing %d rows with duplicated nucleus_supervoxel_id (%d unique values)",
                    n_clear, length(dup_nsv)))
    bc.new$nucleus_supervoxel_id[is_dup & !is_owner] <- NA
    bc.new$nucleus_id[is_dup & !is_owner] <- NA
    bc.new$nucleus_position[is_dup & !is_owner] <- NA
    bc.new$nucleus_position_nm[is_dup & !is_owner] <- NA
  }
}

# Push nucleus data changes to SeaTable
nucleus_cols <- c("_id", "nucleus_id", "nucleus_supervoxel_id", "nucleus_position", "nucleus_position_nm")
bc.nucleus.update <- bc.new[, nucleus_cols] %>%
  dplyr::anti_join(bc[, nucleus_cols], by = nucleus_cols)
if (nrow(bc.nucleus.update)) {
  bc.nucleus.update[is.na(bc.nucleus.update)] <- ""
  message(sprintf("Updating nucleus data for %d neurons in banctable", nrow(bc.nucleus.update)))
  try(banctable_update_rows(base = 'banc_meta',
                            table = "banc_meta",
                            df = bc.nucleus.update,
                            append_allowed = FALSE,
                            chunksize = 1000))
}

# Enforce agreement: always derive nm from raw. Use %.0f format so SeaTable
# never sees scientific notation (default %g flips to e+0N at ~1e6 nm — and
# BANC X/Y can exceed that).
has_root_pos <- !is.na(bc.new$root_position) & bc.new$root_position != ""
if (any(has_root_pos)) {
  bc.new$root_position_nm[has_root_pos] <- nat::xyzmatrix2str(
    banc_raw2nm(xyzmatrix(bc.new$root_position[has_root_pos])),
    format = "%.0f, %.0f, %.0f")
}

# Calculate
bc.root.pos <- as.data.frame(nat::xyzmatrix(bc.new$root_position_nm))
bc.root.pos$root_id <- bc.new$root_id
bc.root.pos$neuropil <- NA
bc.root.pos$region <- NA
# Drop rows with NA coordinates — pointsnearby_banc needs valid XYZ
bc.root.pos <- bc.root.pos[complete.cases(bc.root.pos[, c("X", "Y", "Z")]), ]

try({
  bc.regions <- pointsnearby_banc(bc.root.pos, id = "root_id")
  bc.regions <- bc.regions %>%
    dplyr::mutate(root_region = gsub("outside_","",neuropil)) %>%
    dplyr::select(root_id, root_region)

  # Add sides — use bc.regions$root_id to look up positions from bc.new
  # (bc.root.pos may have different row count after pointsnearby_banc)
  side_pos <- bc.new[match(bc.regions$root_id, bc.new$root_id), "root_position_nm"]
  side_xyz <- nat::xyzmatrix(side_pos)
  lrdiffs <- bancr:::banc_lr_position(side_xyz, units = "nm")
  bc.regions$side <- ifelse(lrdiffs > 0, "right", "left")

  # Carry root_position_nm (derived from root_position at lines 490-496)
  bc.regions$root_position_nm <- bc.new$root_position_nm[match(bc.regions$root_id, bc.new$root_id)]

  # EMPTY-ONLY SIDE FILL — never overwrite an existing non-blank side.
  # The geometric L/R cascade above is appropriate for previously-unset
  # rows, but curators may have set `midline`, `center`, `?` or a manually-
  # adjudicated L/R that arbor-level analysis disagrees with the soma
  # position on. Those decisions are load-bearing — preserve them by
  # carrying the stored side forward whenever it is non-empty. Rows with
  # a stored side then become side-no-ops in the anti_join below
  # (root_region / root_position_nm updates still propagate), while rows
  # with empty stored side pick up the freshly-computed L/R.
  bc.regions <- bc.regions %>%
    dplyr::left_join(bc %>% dplyr::select(root_id, .stored_side = side),
                     by = "root_id") %>%
    dplyr::mutate(side = dplyr::if_else(
      is.na(.stored_side) | .stored_side == "",
      side,
      .stored_side
    )) %>%
    dplyr::select(-.stored_side)

  # Update data frame
  bc.update <- bc.regions %>%
    dplyr::anti_join(bc, by = c("root_id",
                                "root_region",
                                "side",
                                "root_position_nm")) %>%
    dplyr::left_join(bc %>%
                       dplyr::select(`_id`,root_id),
                     by="root_id", relationship = "many-to-many") %>%
    dplyr::distinct(`_id`,.keep_all = TRUE)

  if(nrow(bc.update)){
    bc.update[is.na(bc.update)] <- ''
    banctable_update_rows(base='banc_meta',
                          table = "banc_meta",
                          df = bc.update,
                          append_allowed = FALSE,
                          chunksize = 1000)
  }
})

# New rows for nuclei?
# bc.nuclei.new.rows <- nuclei.regions %>%
#   dplyr::mutate(root_id=as.character(root_id)) %>%
#   dplyr::anti_join(banc.backbone.proofread.mod %>%
#                      dplyr::mutate(root_id=as.character(root_id)),
#                    by = c("root_id"))  %>%
#   dplyr::anti_join(bc, by = c("root_id"))

####################
## Update nerve ####
####################

# # Correct nerve
# banc.peripheral.nerves <- banc_peripheral_nerves()
# banc.peripheral.nerves$nerve <- banc.peripheral.nerves$tag
# banc.peripheral.nerves$pt_root_id <- as.character(banc.peripheral.nerves$pt_root_id)
# bc <- banctable_query("SELECT _id, root_id, status, super_class, nerve, side from banc_meta") %>%
#   dplyr::filter(!grepl("GLIA|NOT_A_NEURON|NOT_NEURON|DELETE",status))
# bc.new <- bc %>%
#   dplyr::left_join(banc.peripheral.nerves %>%
#                      dplyr::select(pt_root_id, nerve) %>%
#                      dplyr::distinct(pt_root_id, .keep_all = TRUE),
#                    by=c("root_id"="pt_root_id")) %>%
#   dplyr::mutate(nerve = dplyr::case_when(
#     grepl("REMOVE_NERVE_LABEL",status) ~ NA,
#     !grepl("sensory|efferent|endocrine|motor|afferent",super_class) ~ NA,
#     !is.na(nerve.y) ~ nerve.y,
#     TRUE ~ nerve.x
#   )) %>%
#   dplyr::mutate(side = dplyr::case_when(
#     side%in%c("left","right","midline") ~ side,
#     grepl("right",nerve) ~ "right",
#     grepl("left",nerve) ~ "left",
#     grepl("midline",nerve) ~ "midline",
#     TRUE ~ side
#   )) %>%
#   dplyr::ungroup() %>%
#   dplyr::select(-nerve.x, -nerve.y, -super_class) %>%
#   as.data.frame()
# bc.update <- bc.new %>%
#   dplyr::anti_join(bc,by=c("root_id","nerve","side"))
# 
# # Update
# if(nrow(bc.update)){
#   bc.update[is.na(bc.update)] <- ''
#   banctable_update_rows(base='banc_meta', 
#                         table = "banc_meta", 
#                         df = bc.update, 
#                         append_allowed = FALSE, 
#                         chunksize = 1000)
# }

# Convert valid to logical proofread column
banc.backbone.proofread <- banc.backbone.proofread %>%
  dplyr::mutate(proofread = valid == TRUE | valid == "t")

# Count daily
proofread_over_time <- banc.backbone.proofread %>%
  dplyr::distinct(pt_root_id, .keep_all = TRUE) %>%
  dplyr::filter(proofread) %>%
  dplyr::mutate(date = as.Date(created)) %>%
  dplyr::count(date) %>%
  dplyr::arrange(date) %>%
  dplyr::mutate(cumulative_n = cumsum(n))
start <- proofread_over_time[1,"n"]
proofread_over_time <- proofread_over_time[-1,]

# Long format for faceting
proofread_long <- tidyr::pivot_longer(
  proofread_over_time,
  cols = c(n, cumulative_n),
  names_to = "type",
  values_to = "count"
) %>%
  dplyr::mutate(type = dplyr::recode(
    type,
    n = "Daily count",
    cumulative_n = "Cumulative count"
  ))

# Faceted plot
g <- ggplot2::ggplot(proofread_long, ggplot2::aes(x = date, y = count)) +
  ggplot2::geom_line(size = 1, color = "#0072B2") +
  ggplot2::facet_wrap(~type, scales = "free_y", ncol = 1) +
  ggplot2::labs(title = "Proofread counts over time",
                x = "Date", y = "Count") +
  ggplot2::theme_minimal() +
  ggplot2::geom_vline(xintercept = as.numeric(Sys.Date()),
                      linetype = "dashed", color = "red", size = 1) +
  ggplot2::scale_x_date(
    date_breaks = "1 week",
    date_labels = "%b %d"
  ) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))

# Save
ggsave(plot = g,
       filename = "inst/images/banc_proofreading_over_time.png",
       width = 12,
       height = 8,
       dpi = 300,
       bg = "transparent")


################################
## SAVE CONNECTIVITY UPDATE ####
################################
try({
  bc <- bc.all %>%
    dplyr::select(`_id`, root_626, root_850, root_888, root_id) %>%
    dplyr::filter(!is.na(root_626), root_626!="", root_626!="0")
  bc.edgelist <- banc_edgelist()
  bc.edgelist <- bc.edgelist %>%
    dplyr::rename(count = n) %>%
    dplyr::mutate(pre_pt_root_id = as.character(pre_pt_root_id),
                  post_pt_root_id = as.character(post_pt_root_id) ) %>%
    dplyr::left_join(bc %>%
                       dplyr::distinct(root_id, .keep_all = TRUE) %>%
                       dplyr::distinct(pre_root_626=root_626, pre_root_850=root_850,
                                       pre_root_888=root_888, pre_pt_root_id=root_id),
                     by = "pre_pt_root_id") %>%
    dplyr::left_join(bc %>%
                       dplyr::distinct(root_id, .keep_all = TRUE) %>%
                       dplyr::distinct(post_root_626=root_626, post_root_850=root_850,
                                       post_root_888=root_888, post_pt_root_id=root_id),
                     by = "post_pt_root_id") %>%
    dplyr::group_by(post_pt_root_id) %>%
    dplyr::mutate(total_input = sum(count)) %>%
    dplyr::group_by(pre_pt_root_id, post_pt_root_id) %>%
    dplyr::mutate(norm = count/total_input) %>%
    dplyr::ungroup()
  readr::write_csv(x=bc.edgelist,
                   file = file.path(banc.connectivity.save.path,"banc_simple_edgelist.csv"))
})

#######################################################
## CHECK: nucleus_id ↔ root_id consistency in banc_meta
##
## Drift mode: a row in banc_meta has nucleus_id = N and root_id = R,
## but CAVE's nuclei table now reports nucleus N inside root_id R' (R != R').
## Causes: a split/merge moved the nucleus's supervoxel to a different
## root since the row was last updated, and our update logic missed it
## (e.g. update path triggered on root_id-changed but not nucleus-moved).
##
## REPORT-ONLY for now — write a CSV. We'll wire a fix in once we've
## seen the drift count and the patterns in the report.
#######################################################
try({
  message("\n##### Drift check: nucleus_id ↔ root_id #####")

  cave_nuc <- bancr::banc_nuclei(rawcoords = TRUE, table = "both") %>%
    dplyr::transmute(
      nucleus_id      = as.character(id),
      cave_root_id    = as.character(pt_root_id)
    ) %>%
    dplyr::filter(!is.na(nucleus_id), nucleus_id != "",
                  !is.na(cave_root_id), cave_root_id != "0")

  bc.live <- bancr:::banctable_query(
    "SELECT _id, root_id, nucleus_id, side, cell_class, cell_type from banc_meta") %>%
    dplyr::filter(!is.na(nucleus_id), nucleus_id != "",
                  !is.na(root_id), root_id != "", root_id != "0") %>%
    dplyr::mutate(nucleus_id = as.character(nucleus_id),
                  root_id    = as.character(root_id))

  drift <- bc.live %>%
    dplyr::inner_join(cave_nuc, by = "nucleus_id") %>%
    dplyr::filter(root_id != cave_root_id)

  message(sprintf("  banc_meta rows with both ids:    %s",
                  format(nrow(bc.live),  big.mark = ",")))
  message(sprintf("  CAVE nuclei rows w/ root_id:     %s",
                  format(nrow(cave_nuc), big.mark = ",")))
  message(sprintf("  DRIFT rows (root_id != CAVE):    %s",
                  format(nrow(drift),    big.mark = ",")))

  # Also: nuclei in banc_meta that no longer exist in CAVE (deleted/missing)
  missing_in_cave <- bc.live %>%
    dplyr::anti_join(cave_nuc, by = "nucleus_id")
  message(sprintf("  banc_meta nucleus_ids missing in CAVE: %s",
                  format(nrow(missing_in_cave), big.mark = ",")))

  # Also: many-to-one (multiple banc_meta rows pointing at same nucleus)
  dupes <- bc.live %>%
    dplyr::count(nucleus_id) %>%
    dplyr::filter(n > 1L)
  message(sprintf("  nucleus_ids assigned to >1 banc_meta row: %s",
                  format(nrow(dupes), big.mark = ",")))

  out_dir <- file.path(banc.meta.save.path, "drift_reports")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  stamp <- format(Sys.time(), "%Y-%m-%d")
  if (nrow(drift) > 0L) {
    drift_path <- file.path(out_dir,
                            sprintf("banc_meta_nucleus_root_drift_%s.csv", stamp))
    readr::write_csv(drift, drift_path)
    message(sprintf("  -> %s", drift_path))
  }
  if (nrow(missing_in_cave) > 0L) {
    miss_path <- file.path(out_dir,
                           sprintf("banc_meta_nucleus_missing_in_cave_%s.csv", stamp))
    readr::write_csv(missing_in_cave, miss_path)
    message(sprintf("  -> %s", miss_path))
  }
  if (nrow(dupes) > 0L) {
    dup_path <- file.path(out_dir,
                          sprintf("banc_meta_nucleus_dupes_%s.csv", stamp))
    readr::write_csv(dupes, dup_path)
    message(sprintf("  -> %s", dup_path))
  }
}, silent = FALSE)

# Close logging
message(sprintf("Log ended: %s", Sys.time()))
sink(type = "message")
sink(type = "output")
close(.log_con)
message(sprintf("Log saved: %s", .log_file))

