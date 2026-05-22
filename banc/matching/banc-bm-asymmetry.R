#' banc-bm-asymmetry — Diagnose BM sensory L/R asymmetry and propose retypes.
#'
#' For all BANC sensory neurons whose cell_type starts with `BM_`, report
#' L/R count asymmetry and propose retype corrections by scoring each
#' BANC neuron's downstream cell_type connectivity profile against FAFB BM
#' type centroids (cosine similarity). Suggests fafb_match via best NBLAST
#' hit of the new type. Generates per-candidate + summary neuroglancer URLs.
#'
#' @section Reads:
#'   - SeaTable `banc_meta` (BANC BM neurons + downstream partners)
#'   - `franken_meta()` (FAFB BM neurons for centroids)
#'   - `banc_fafb_783_nblast.feather` (best-match NBLAST scores)
#'
#' @section Writes:
#'   - `data/codex/bm_asymmetry.csv` — L/R count diagnostics
#'   - `data/codex/bm_retype_candidates.csv` — proposed retypes + scores

###########################################################
### BM sensory L/R asymmetry diagnostic + retype proposals
###
### For all BANC sensory neurons whose cell_type starts with
### "BM_", report L/R count asymmetry and propose retype
### corrections by scoring each BANC BM neuron's downstream
### connectivity profile against FAFB BM type centroids.
###
### Method:
###   1. Build per-neuron downstream cell_type profiles in
###      both BANC and FAFB by aggregating downstream partners
###      by their cell_type label (drop unlabelled partners),
###      then L1-normalising. BANC seatable and franken_meta
###      share a cell_type vocabulary, so partner-type vectors
###      are directly comparable across datasets and sides.
###   2. Build per-cell_type centroids from FAFB BM neurons
###      (mean of L1-normalised profiles). FAFB is the target.
###   3. Score every BANC BM neuron (both sides) against every
###      FAFB BM centroid by cosine similarity. Pick best,
###      runner-up. Margin = best - runner_up.
###   4. Where best != current AND margin >= confidence_min,
###      propose a retype. Otherwise keep the current type.
###   5. Re-enumerate L/R counts under the new assignment.
###   6. Pick the best NBLAST hit of the new type from
###      banc_fafb_783 NBLAST → suggested fafb_match.
###   7. Build neuroglancer review URLs:
###        - one per retype candidate (current vs proposed
###          type, each in its own segmentation_with_graph
###          layer with a fixed colour)
###        - two summary URLs covering all BANC BM neurons,
###          one layer per BM type with a unique colour:
###            (a) coloured by ORIGINAL cell_type
###            (b) coloured by RETYPED cell_type
###   8. Build a SeaTable push frame; hold back rows with
###      HAS_MANUAL_ANNOTATION.
###
### Inputs:
###   banc.meta (live SeaTable)
###   franken_meta()
###   banc_{version}_edgelist_simple.feather
###   fafb_783_simple_edgelist.feather (or flywire_783_edgelist.feather)
###   banc_fafb_783_nblast.feather
###
### Outputs:
###   data/codex/bm_asymmetry.csv             — L/R counts before+after
###   data/codex/bm_retype_candidates.csv     — all BANC BM scores + URLs
###   data/codex/bm_retype_push.csv           — eligible push frame
###   data/codex/bm_retype_manual_review.csv  — held-back rows
###   data/codex/bm_summary_ngl_urls.txt      — 2 summary NGL URLs
###########################################################
source("banc/banc-startup.R")

local({

message("### banc: BM sensory L/R asymmetry + retype candidates ###")

#######################
### Tunables        ###
#######################

confidence_min <- 0.1   # min margin (best - runner_up) to accept a retype
min_partners   <- 3     # min typed downstream partners to score a neuron
max_per_layer  <- 30    # max ref-pool neurons per per-neuron NGL layer
color_current  <- "#e64a4a"  # red — candidate's current type
color_proposed <- "#4ae64a"  # green — candidate's proposed type

###########################
### Read BANC + FAFB    ###
###########################

select_cols <- paste(c("_id", "root_888", "root_id", "cell_type", "super_class",
                       "side", "status", "fafb_match", "cell_type_source"),
                     collapse = ", ")
bc <- banctable_query(sprintf("SELECT %s FROM banc_meta", select_cols)) %>%
  dplyr::mutate(dplyr::across(dplyr::everything(), as.character)) %>%
  dplyr::filter(!is.na(root_888), root_888 != "") %>%
  dplyr::distinct(root_888, .keep_all = TRUE)
message(sprintf("  BANC seatable: %d rows", nrow(bc)))

message("  Loading franken_meta...")
fm <- franken_meta() %>%
  dplyr::mutate(dplyr::across(dplyr::everything(), as.character))

# BANC BM_* sensory neurons (both sides)
bm_bc <- bc %>%
  dplyr::filter(!is.na(super_class),
                grepl("sensory", super_class, ignore.case = TRUE),
                !is.na(cell_type), grepl("^BM_", cell_type))
message(sprintf("  BANC BM_* sensory neurons: %d", nrow(bm_bc)))

fm_bm <- fm %>%
  dplyr::filter(!is.na(fafb_id), fafb_id != "",
                !is.na(cell_type), grepl("^BM_", cell_type))
message(sprintf("  franken_meta BM_* rows:    %d", nrow(fm_bm)))

###########################
### Asymmetry (before)  ###
###########################

count_by_side <- function(df, dataset_name) {
  df %>%
    dplyr::mutate(side = ifelse(is.na(side) | side == "", "unknown", side)) %>%
    dplyr::count(cell_type, side, name = "n") %>%
    dplyr::mutate(dataset = dataset_name)
}

asym_old <- dplyr::bind_rows(count_by_side(bm_bc, "banc"),
                             count_by_side(fm_bm, "fafb")) %>%
  tidyr::pivot_wider(names_from = c(dataset, side), values_from = n, values_fill = 0L)
for (col in c("banc_left", "banc_right", "fafb_left", "fafb_right")) {
  if (!col %in% names(asym_old)) asym_old[[col]] <- 0L
}
asym_old <- asym_old %>%
  dplyr::mutate(
    banc_total   = banc_left + banc_right,
    fafb_total   = fafb_left + fafb_right,
    banc_lr_diff = banc_right - banc_left,
    fafb_lr_diff = fafb_right - fafb_left
  ) %>%
  dplyr::arrange(dplyr::desc(abs(banc_lr_diff)))

top_n <- min(25, nrow(asym_old))
message("\n  === BANC BM_* L/R imbalances (current cell_type) ===")
for (i in seq_len(top_n)) {
  r <- asym_old[i, ]
  message(sprintf("    %-30s  banc L/R = %3d/%-3d  fafb L/R = %3d/%-3d",
                  r$cell_type, r$banc_left, r$banc_right,
                  r$fafb_left, r$fafb_right))
}

###########################
### Edgelists           ###
###########################

on_o2 <- nzchar(Sys.getenv("SLURM_JOB_ID")) ||
  grepl("o2\\.rc\\.hms\\.harvard\\.edu|compute-", Sys.info()["nodename"])

cache_dir <- file.path(tempdir(), "bm_asymmetry_cache")
dir.create(cache_dir, showWarnings = FALSE)

resolve_edgelist <- function(local_basenames, gcs_uri) {
  candidates <- character(0)
  for (b in local_basenames) {
    candidates <- c(candidates,
      if (exists("banc.versioned.save.path")) file.path(banc.versioned.save.path, b),
      if (exists("banc.connectivity.save.path")) file.path(banc.connectivity.save.path, b),
      file.path("/n/data1/hms/neurobio/wilson/connectomes/banc",
                sprintf("banc_%s", banc.version), b),
      file.path("/n/data1/hms/neurobio/wilson/banc/connectivity", b))
  }
  for (p in candidates) if (!is.null(p) && file.exists(p)) return(p)
  cached <- file.path(cache_dir, basename(gcs_uri))
  if (!file.exists(cached)) {
    message(sprintf("    Downloading %s from GCS...", basename(gcs_uri)))
    system2("gsutil", c("cp", gcs_uri, cached), stdout = FALSE, stderr = FALSE)
  }
  if (!file.exists(cached)) stop("Could not locate any of: ",
                                  paste(local_basenames, collapse = ", "))
  cached
}

banc_elist_basename <- sprintf("banc_%s_edgelist_simple.feather", banc.version)
banc_elist_path <- resolve_edgelist(
  banc_elist_basename,
  sprintf("gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/banc_%s/%s",
          banc.version, banc_elist_basename))
message(sprintf("  Loading BANC edgelist: %s", banc_elist_path))
banc_el <- arrow::read_feather(banc_elist_path) %>%
  dplyr::mutate(pre = as.character(pre), post = as.character(post),
                count = as.numeric(count))
message(sprintf("    BANC edgelist: %d edges", nrow(banc_el)))

# FAFB edgelist — try the O2 name first, then the GCS name
fafb_elist_path <- resolve_edgelist(
  c("flywire_783_edgelist.feather", "fafb_783_simple_edgelist.feather"),
  "gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/fafb_783/fafb_783_simple_edgelist.feather")
message(sprintf("  Loading FAFB edgelist: %s", fafb_elist_path))
fafb_el <- arrow::read_feather(fafb_elist_path) %>%
  dplyr::mutate(pre = as.character(pre), post = as.character(post),
                count = as.numeric(count))
message(sprintf("    FAFB edgelist: %d edges", nrow(fafb_el)))

###########################
### Connectivity profiles
###########################

# Type lookups: post root -> partner cell_type. BANC and FAFB share vocab.
banc_type_lookup <- bc %>%
  dplyr::select(rid = root_888, partner_type = cell_type) %>%
  dplyr::filter(!is.na(partner_type), partner_type != "")

fafb_type_lookup <- fm %>%
  dplyr::filter(!is.na(fafb_id), fafb_id != "",
                !is.na(cell_type), cell_type != "") %>%
  dplyr::select(rid = fafb_id, partner_type = cell_type) %>%
  dplyr::distinct(rid, .keep_all = TRUE)

build_profiles <- function(elist, pre_ids, type_lookup) {
  elist %>%
    dplyr::filter(pre %in% pre_ids) %>%
    dplyr::inner_join(type_lookup, by = c("post" = "rid")) %>%
    dplyr::group_by(pre, partner_type) %>%
    dplyr::summarise(weight = sum(count), .groups = "drop") %>%
    dplyr::group_by(pre) %>%
    dplyr::mutate(weight = weight / sum(weight)) %>%
    dplyr::ungroup()
}

bm_pre_banc <- bm_bc$root_888
bm_pre_fafb <- fm_bm$fafb_id

message("  Building BANC BM downstream profiles...")
profiles_banc <- build_profiles(banc_el, bm_pre_banc, banc_type_lookup)
message("  Building FAFB BM downstream profiles...")
profiles_fafb <- build_profiles(fafb_el, bm_pre_fafb, fafb_type_lookup)

# Drop neurons with too few typed partners (noisy profiles)
np_banc <- profiles_banc %>% dplyr::count(pre, name = "n_partners")
np_fafb <- profiles_fafb %>% dplyr::count(pre, name = "n_partners")
keep_banc <- np_banc$pre[np_banc$n_partners >= min_partners]
keep_fafb <- np_fafb$pre[np_fafb$n_partners >= min_partners]
profiles_banc <- profiles_banc %>% dplyr::filter(pre %in% keep_banc)
profiles_fafb <- profiles_fafb %>% dplyr::filter(pre %in% keep_fafb)
message(sprintf("  BANC BM neurons with >= %d typed partners: %d / %d",
                min_partners, length(keep_banc), length(bm_pre_banc)))
message(sprintf("  FAFB BM neurons with >= %d typed partners: %d / %d",
                min_partners, length(keep_fafb), length(bm_pre_fafb)))

# Shared partner-type vocabulary (union of types observed in either dataset)
all_partner_types <- sort(unique(c(profiles_banc$partner_type,
                                   profiles_fafb$partner_type)))
all_pre_banc      <- sort(unique(profiles_banc$pre))
all_pre_fafb      <- sort(unique(profiles_fafb$pre))

M_banc <- Matrix::sparseMatrix(
  i = match(profiles_banc$pre, all_pre_banc),
  j = match(profiles_banc$partner_type, all_partner_types),
  x = profiles_banc$weight,
  dims = c(length(all_pre_banc), length(all_partner_types)),
  dimnames = list(all_pre_banc, all_partner_types))

M_fafb <- Matrix::sparseMatrix(
  i = match(profiles_fafb$pre, all_pre_fafb),
  j = match(profiles_fafb$partner_type, all_partner_types),
  x = profiles_fafb$weight,
  dims = c(length(all_pre_fafb), length(all_partner_types)),
  dimnames = list(all_pre_fafb, all_partner_types))

message(sprintf("  Profile matrices: BANC %d x %d, FAFB %d x %d (shared vocab of %d types)",
                nrow(M_banc), ncol(M_banc), nrow(M_fafb), ncol(M_fafb),
                length(all_partner_types)))

###########################
### FAFB centroids      ###
###########################

# Mean FAFB profile per BM cell_type
fm_bm_used <- fm_bm[match(all_pre_fafb, fm_bm$fafb_id), , drop = FALSE]
fafb_type_factor <- factor(fm_bm_used$cell_type)
ind_fafb <- Matrix::sparseMatrix(
  i = seq_len(nrow(M_fafb)),
  j = as.integer(fafb_type_factor),
  x = 1,
  dims = c(nrow(M_fafb), nlevels(fafb_type_factor)))
type_n     <- as.integer(table(fafb_type_factor))
type_sum   <- Matrix::t(ind_fafb) %*% M_fafb
centroid_M <- as.matrix(type_sum) / type_n
rownames(centroid_M) <- levels(fafb_type_factor)
message(sprintf("  FAFB BM centroids: %d types x %d partner types",
                nrow(centroid_M), ncol(centroid_M)))

# L2-normalise centroids and BANC profiles for cosine similarity
cent_norm <- sqrt(rowSums(centroid_M^2))
cent_norm[cent_norm == 0] <- 1
centroid_M_norm <- centroid_M / cent_norm

banc_norm <- sqrt(Matrix::rowSums(M_banc^2))
banc_norm[banc_norm == 0] <- 1
M_banc_norm <- as.matrix(M_banc / banc_norm)

###########################
### Score BANC neurons  ###
###########################

sim_mat <- M_banc_norm %*% t(centroid_M_norm)
rownames(sim_mat) <- all_pre_banc
colnames(sim_mat) <- rownames(centroid_M_norm)
message(sprintf("  Cosine similarity: %d BANC neurons x %d FAFB centroids",
                nrow(sim_mat), ncol(sim_mat)))

best_idx   <- max.col(sim_mat, ties.method = "first")
best_score <- sim_mat[cbind(seq_along(best_idx), best_idx)]
sim_mat2   <- sim_mat
sim_mat2[cbind(seq_along(best_idx), best_idx)] <- -Inf
ru_idx   <- max.col(sim_mat2, ties.method = "first")
ru_score <- sim_mat2[cbind(seq_along(ru_idx), ru_idx)]
ru_score[!is.finite(ru_score)] <- NA_real_

# Self-similarity (current type, if it's in the centroid set)
bm_meta <- bm_bc[match(all_pre_banc, bm_bc$root_888), , drop = FALSE]
self_idx <- match(bm_meta$cell_type, colnames(sim_mat))
self_score <- ifelse(is.na(self_idx), NA_real_,
                     sim_mat[cbind(seq_along(self_idx), self_idx)])

candidates <- data.frame(
  `_id`              = bm_meta$`_id`,
  root_888           = all_pre_banc,
  root_id            = bm_meta$root_id,
  side               = bm_meta$side,
  status             = bm_meta$status,
  current_fafb_match = bm_meta$fafb_match,
  cell_type_source   = bm_meta$cell_type_source,
  current_type       = bm_meta$cell_type,
  current_score      = self_score,
  best_type          = colnames(sim_mat)[best_idx],
  best_score         = best_score,
  runner_up_type     = colnames(sim_mat)[ru_idx],
  runner_up_score    = ru_score,
  stringsAsFactors   = FALSE,
  check.names        = FALSE)
candidates$confidence <- candidates$best_score -
  ifelse(is.na(candidates$runner_up_score), 0, candidates$runner_up_score)
candidates$delta_self <- candidates$best_score -
  ifelse(is.na(candidates$current_score), 0, candidates$current_score)

# Confidence rule:
#   best == current                                          → no change
#   best != current AND confidence >= confidence_min          → propose retype
#   otherwise (low confidence)                                → keep current
candidates$is_change <- !is.na(candidates$best_type) &
                        candidates$best_type != candidates$current_type &
                        !is.na(candidates$confidence) &
                        candidates$confidence >= confidence_min

# Final assigned type after rule application
candidates$new_type <- ifelse(candidates$is_change,
                              candidates$best_type,
                              candidates$current_type)

# BANC BM neurons that were dropped for too few partners stay at their current
# type (no change). Append them so the candidates frame covers ALL BM neurons.
unscored_ids <- setdiff(bm_bc$root_888, all_pre_banc)
if (length(unscored_ids) > 0) {
  unscored_meta <- bm_bc[match(unscored_ids, bm_bc$root_888), , drop = FALSE]
  unscored_df <- data.frame(
    `_id`              = unscored_meta$`_id`,
    root_888           = unscored_ids,
    root_id            = unscored_meta$root_id,
    side               = unscored_meta$side,
    status             = unscored_meta$status,
    current_fafb_match = unscored_meta$fafb_match,
    cell_type_source   = unscored_meta$cell_type_source,
    current_type       = unscored_meta$cell_type,
    current_score      = NA_real_,
    best_type          = NA_character_,
    best_score         = NA_real_,
    runner_up_type     = NA_character_,
    runner_up_score    = NA_real_,
    stringsAsFactors   = FALSE,
    check.names        = FALSE)
  unscored_df$confidence <- NA_real_
  unscored_df$delta_self <- NA_real_
  unscored_df$is_change  <- FALSE
  unscored_df$new_type   <- unscored_meta$cell_type
  candidates <- rbind(candidates, unscored_df)
  message(sprintf("  %d BANC BM neurons left at current type (too few typed partners)",
                  length(unscored_ids)))
}

n_changed <- sum(candidates$is_change)
n_low_conf <- sum(!is.na(candidates$best_type) &
                  candidates$best_type != candidates$current_type &
                  candidates$confidence < confidence_min, na.rm = TRUE)
message(sprintf("\n  BANC BM neurons scored:                %d", sum(!is.na(candidates$best_type))))
message(sprintf("  Proposed retypes (best != current):    %d", n_changed))
message(sprintf("  Rejected for low confidence (< %.2f):  %d", confidence_min, n_low_conf))

# Top source -> target retype flows
if (n_changed > 0) {
  flow_tab <- candidates %>%
    dplyr::filter(is_change) %>%
    dplyr::count(current_type, best_type, name = "n") %>%
    dplyr::arrange(dplyr::desc(n))
  message("\n  Top source -> target retype flows:")
  for (i in seq_len(min(20, nrow(flow_tab)))) {
    r <- flow_tab[i, ]
    message(sprintf("    %3d  %-30s -> %s", r$n, r$current_type, r$best_type))
  }
}

###########################
### Re-enumerate L/R    ###
###########################

count_by_side_simple <- function(types, sides) {
  side_norm <- ifelse(is.na(sides) | sides == "", "unknown", sides)
  df <- data.frame(cell_type = types, side = side_norm,
                   stringsAsFactors = FALSE)
  df %>% dplyr::count(cell_type, side, name = "n")
}

old_counts <- count_by_side_simple(candidates$current_type, candidates$side) %>%
  tidyr::pivot_wider(names_from = side, values_from = n, values_fill = 0L,
                     names_prefix = "old_")
new_counts <- count_by_side_simple(candidates$new_type, candidates$side) %>%
  tidyr::pivot_wider(names_from = side, values_from = n, values_fill = 0L,
                     names_prefix = "new_")

asym_new <- dplyr::full_join(old_counts, new_counts, by = "cell_type")
for (col in c("old_left", "old_right", "new_left", "new_right")) {
  if (!col %in% names(asym_new)) asym_new[[col]] <- 0L
  asym_new[[col]][is.na(asym_new[[col]])] <- 0L
}
asym_new <- asym_new %>%
  dplyr::mutate(
    old_total   = old_left + old_right,
    new_total   = new_left + new_right,
    old_lr_diff = old_right - old_left,
    new_lr_diff = new_right - new_left
  ) %>%
  dplyr::left_join(
    asym_old %>% dplyr::select(cell_type, fafb_left, fafb_right, fafb_total, fafb_lr_diff),
    by = "cell_type") %>%
  dplyr::mutate(dplyr::across(c(fafb_left, fafb_right, fafb_total, fafb_lr_diff),
                              ~ ifelse(is.na(.x), 0L, .x))) %>%
  dplyr::arrange(dplyr::desc(abs(new_lr_diff)))

asym_file <- "data/codex/bm_asymmetry.csv"
readr::write_csv(asym_new, asym_file)
message(sprintf("\n  Saved asymmetry table (before+after) to %s (%d types)",
                asym_file, nrow(asym_new)))

message("  Re-enumerated L/R (top 25 by |new_lr_diff|):")
for (i in seq_len(min(25, nrow(asym_new)))) {
  r <- asym_new[i, ]
  message(sprintf("    %-30s  banc OLD L/R = %3d/%-3d  NEW L/R = %3d/%-3d  fafb L/R = %3d/%-3d",
                  r$cell_type, r$old_left, r$old_right,
                  r$new_left, r$new_right, r$fafb_left, r$fafb_right))
}

###########################
### NBLAST -> fafb_match
###########################

message("\n  Loading FAFB NBLAST feather for fafb_match suggestion...")
nblast_cache <- file.path(tempdir(), "nblast_cache")
dir.create(nblast_cache, showWarnings = FALSE)
nblast_file <- file.path(nblast_cache, "banc_fafb_783_nblast.feather")
if (!file.exists(nblast_file)) {
  alt_file <- file.path("/tmp/nblast_cache", "banc_fafb_783_nblast.feather")
  if (file.exists(alt_file)) nblast_file <- alt_file
  else {
    message("    Downloading NBLAST feather from GCS...")
    system2("gsutil", c("cp",
      "gs://lee-lab_brain-and-nerve-cord-fly-connectome/nblast/banc_fafb_783_nblast.feather",
      nblast_file), stdout = FALSE, stderr = FALSE)
  }
}
fafb_nblast <- arrow::read_feather(nblast_file)
fafb_nblast_idx <- split(
  data.frame(match_id = as.character(fafb_nblast$match_id),
             score = fafb_nblast$score, stringsAsFactors = FALSE),
  as.character(fafb_nblast$root_888))
message(sprintf("  NBLAST: %d rows, indexed for %d BANC neurons",
                nrow(fafb_nblast), length(fafb_nblast_idx)))

# franken_meta cell_type -> fafb_ids index
fm_typed <- fm %>% dplyr::filter(!is.na(fafb_id), fafb_id != "",
                                 !is.na(cell_type), cell_type != "")
fm_type_ids <- split(as.character(fm_typed$fafb_id), fm_typed$cell_type)

candidates$suggested_fafb_match <- NA_character_
candidates$suggested_nblast     <- NA_real_

retype_mask <- candidates$is_change
for (i in which(retype_mask)) {
  rid  <- candidates$root_888[i]
  newt <- candidates$best_type[i]
  fids <- fm_type_ids[[newt]]
  if (is.null(fids) || length(fids) == 0) next
  hits <- fafb_nblast_idx[[rid]]
  if (is.null(hits) || nrow(hits) == 0) next
  mask <- hits$match_id %in% fids
  if (any(mask)) {
    th <- hits[mask, , drop = FALSE]
    j <- which.max(th$score)
    candidates$suggested_fafb_match[i] <- th$match_id[j]
    candidates$suggested_nblast[i]     <- th$score[j]
  }
}
message(sprintf("  Retype candidates with NBLAST-picked fafb_match: %d / %d",
                sum(!is.na(candidates$suggested_fafb_match)), sum(retype_mask)))

###############################################
### Neuroglancer URL helpers                ###
###############################################
### Decode the BANC base scene once and use  ###
### its `v850 neurons` segmentation_with_graph
### layer as a deep-copy template. New layers
### are appended via direct list assignment;
### segments and a per-segment colour map are
### written straight to each new layer's
### $segments / $segmentColors fields. We do
### NOT use fafbseg::ngl_add_colours, which
### writes segments to the FIRST seg layer
### regardless of the `layer=` argument.
###############################################

message("\n  Decoding base neuroglancer scene...")

ngl_url   <- "https://spelunker.cave-explorer.org/#!middleauth+https://global.daf-apis.com/nglstate/api/v1/4899366325190656"
ngl_url2  <- sub("#!middleauth+", "?", ngl_url, fixed = TRUE)
ngl_parts <- unlist(strsplit(ngl_url2, "?", fixed = TRUE))
ngl_json  <- fafbseg::flywire_fetch(ngl_parts[2], token = bancr:::banc_token(),
                                    return = "text", cache = TRUE)
sc_base   <- fafbseg::ngl_decode_scene(
  fafbseg::ngl_encode_url(ngl_json, baseurl = ngl_parts[1]))

# Find the BANC graphene segmentation layer. Match by name first, then by
# source URL as a fallback. Crucially, do NOT fall back to grep("banc")
# on the name — the BANC EM image layer is also named "BANC EM".
ngl_ls   <- fafbseg:::ngl_layer_summary(fafbseg::ngl_layers(sc_base))
banc_idx <- match("segmentation proofreading", ngl_ls$name)
if (is.na(banc_idx)) banc_idx <- match("v850 neurons", ngl_ls$name)
if (is.na(banc_idx)) {
  for (i in seq_along(sc_base$layers)) {
    src <- sc_base$layers[[i]]$source
    if (is.list(src)) src <- src$url
    if (!is.null(src) && is.character(src) &&
        grepl("graphene.*wclee_fly_cns", src)) {
      banc_idx <- i; break
    }
  }
}
if (is.na(banc_idx)) stop("Could not find BANC graphene segmentation layer in base scene")
banc_src <- sc_base$layers[[banc_idx]]$source
if (is.list(banc_src)) banc_src <- banc_src$url
banc_type <- sc_base$layers[[banc_idx]]$type
if (is.null(banc_type)) banc_type <- "<none>"
message(sprintf("  Base BANC layer: '%s' [%d] type=%s",
                ngl_ls$name[banc_idx], banc_idx, banc_type))
if (is.null(banc_src) || !grepl("graphene", banc_src)) {
  stop("Base BANC layer is not a graphene segmentation source: ",
       if (is.null(banc_src)) "<none>" else banc_src)
}

# Template = the base BANC layer with no segments / colours of its own.
banc_layer_template <- sc_base$layers[[banc_idx]]
banc_layer_template$segments         <- list()
banc_layer_template$hiddenSegments   <- NULL
banc_layer_template$segmentColors    <- NULL
banc_layer_template$segmentDefaultColor <- NULL

# Resolve a BANC neuron's NGL segment id (root_id preferred, fall back to root_888)
ngl_id <- function(root_id, root_888) {
  ifelse(!is.na(root_id) & root_id != "", root_id, root_888)
}

# Build a layer (deep copy of the template) holding `ids` all coloured `colour`.
# Both segmentDefaultColor (one colour for the whole layer) and segmentColors
# (per-segment dict) are set so the colour survives whichever neuroglancer
# build the user has open.
make_banc_layer <- function(name, ids, colour) {
  ids <- as.character(unique(ids))
  ids <- ids[!is.na(ids) & ids != ""]
  hex <- tolower(colour)
  L <- banc_layer_template
  L$name <- name
  L$segments <- as.list(ids)
  L$segmentDefaultColor <- hex
  if (length(ids) > 0) {
    L$segmentColors <- setNames(as.list(rep(hex, length(ids))), ids)
  } else {
    L$segmentColors <- NULL
  }
  L
}

# Wrap a list of new layers into a fresh scene; the original BANC layer is
# emptied so it doesn't draw anything on top of the new ones.
scene_with_layers <- function(layer_list) {
  sc <- sc_base
  sc$layers[[banc_idx]]$segments       <- list()
  sc$layers[[banc_idx]]$segmentColors  <- NULL
  sc$layers[[banc_idx]]$hiddenSegments <- NULL
  if (length(layer_list) > 0) {
    sc$layers <- c(sc$layers, layer_list)
  }
  sc
}

###############################################
### Per-neuron review URLs                  ###
###############################################

message("\n  Building per-neuron neuroglancer URLs for retype candidates...")

# Reference pool keyed by PRE-rework cell_type. The candidate sits in the
# current-type pool; the proposed-type pool shows what neurons currently
# labelled as the proposed type look like.
ref_pool <- bm_bc %>%
  dplyr::distinct(root_888, .keep_all = TRUE) %>%
  dplyr::mutate(ngl_id = ngl_id(root_id, root_888)) %>%
  dplyr::filter(!is.na(ngl_id), ngl_id != "")
ref_id_lookup <- split(ref_pool$ngl_id, ref_pool$cell_type)
message(sprintf("  Reference pool: %d BANC BM neurons across %d types",
                nrow(ref_pool), length(ref_id_lookup)))

build_per_neuron_url <- function(cand_ngl_id, cur_type, prop_type) {
  cur_pool  <- ref_id_lookup[[cur_type]]
  prop_pool <- ref_id_lookup[[prop_type]]
  if (is.null(cur_pool))  cur_pool  <- character(0)
  if (is.null(prop_pool)) prop_pool <- character(0)
  if (length(cur_pool)  > max_per_layer) cur_pool  <- sample(cur_pool,  max_per_layer)
  if (length(prop_pool) > max_per_layer) prop_pool <- sample(prop_pool, max_per_layer)

  cur_ids  <- unique(c(as.character(cand_ngl_id), as.character(cur_pool)))
  prop_ids <- as.character(prop_pool)

  layers <- list(
    make_banc_layer(paste0("BM_current:",  cur_type),  cur_ids,  color_current),
    make_banc_layer(paste0("BM_proposed:", prop_type), prop_ids, color_proposed))
  as.character(scene_with_layers(layers))
}

candidates$neuroglancer_url <- NA_character_
set.seed(42)
n_done <- 0L
first_err <- NULL
for (i in which(retype_mask)) {
  url <- tryCatch(
    build_per_neuron_url(
      ngl_id(candidates$root_id[i], candidates$root_888[i]),
      candidates$current_type[i],
      candidates$best_type[i]),
    error = function(e) {
      if (is.null(first_err)) first_err <<- conditionMessage(e)
      NA_character_
    })
  candidates$neuroglancer_url[i] <- url
  if (!is.na(url)) n_done <- n_done + 1L
}
if (!is.null(first_err)) message(sprintf("  First NGL error: %s", first_err))
message(sprintf("  Generated %d / %d per-neuron NGL URLs", n_done, sum(retype_mask)))

###############################################
### Two summary URLs                        ###
###############################################
### One layer per BM cell_type, with a unique
### fixed colour per layer. Built twice — once
### over the ORIGINAL labels, once over the
### RETYPED labels — covering ALL BANC BM
### neurons (both sides).

message("\n  Building summary neuroglancer URLs (one layer per BM type)...")

candidates$ngl_id <- ngl_id(candidates$root_id, candidates$root_888)

build_summary_url <- function(types_vec) {
  ok <- !is.na(candidates$ngl_id) & candidates$ngl_id != "" &
        !is.na(types_vec) & types_vec != ""
  type_to_ids <- split(candidates$ngl_id[ok], types_vec[ok])
  type_to_ids <- type_to_ids[order(names(type_to_ids))]
  n_types <- length(type_to_ids)
  if (n_types == 0) return(NA_character_)

  pal <- tryCatch(
    grDevices::palette.colors(n = max(n_types, 2), palette = "Polychrome 36"),
    error = function(e) grDevices::rainbow(max(n_types, 2)))
  if (length(pal) < n_types) pal <- rep_len(pal, n_types)
  cols <- as.character(pal[seq_len(n_types)])

  layers <- mapply(function(name, ids, col) {
    make_banc_layer(paste0("BM:", name), ids, col)
  }, names(type_to_ids), type_to_ids, cols, SIMPLIFY = FALSE)

  as.character(scene_with_layers(layers))
}

summary_url_original <- tryCatch(
  build_summary_url(candidates$current_type),
  error = function(e) {
    message("    summary (original) error: ", conditionMessage(e))
    NA_character_
  })
summary_url_retyped <- tryCatch(
  build_summary_url(candidates$new_type),
  error = function(e) {
    message("    summary (retyped) error: ", conditionMessage(e))
    NA_character_
  })

summary_file <- "data/codex/bm_summary_ngl_urls.txt"
writeLines(c(
  "# BANC BM_* summary neuroglancer URLs",
  sprintf("# Generated %s", Sys.time()),
  "",
  "## Coloured by ORIGINAL cell_type",
  if (is.na(summary_url_original)) "(failed to build)" else summary_url_original,
  "",
  "## Coloured by RETYPED cell_type (post FAFB-centroid scoring)",
  if (is.na(summary_url_retyped)) "(failed to build)" else summary_url_retyped),
  con = summary_file)
message(sprintf("  Saved summary URLs to %s", summary_file))

###########################
### Save / push         ###
###########################

cand_file <- "data/codex/bm_retype_candidates.csv"
readr::write_csv(candidates %>% dplyr::select(-ngl_id), cand_file)
message(sprintf("\n  Saved bm_retype_candidates.csv (%d rows)", nrow(candidates)))

has_manual <- grepl("HAS_MANUAL_ANNOTATION", candidates$status, ignore.case = TRUE)

push_df <- candidates[retype_mask & !has_manual, , drop = FALSE] %>%
  dplyr::transmute(
    `_id`,
    cell_type        = best_type,
    fafb_match       = ifelse(is.na(suggested_fafb_match), "", suggested_fafb_match),
    fafb_cell_type   = best_type,
    cell_type_source = append_status(cell_type_source, "bates"),
    neuroglancer_url) %>%
  as.data.frame()

push_file <- "data/codex/bm_retype_push.csv"
readr::write_csv(push_df, push_file)
message(sprintf("  SeaTable push frame: %d rows -> %s", nrow(push_df), push_file))
message(sprintf("    cell_type updated:    %d", sum(push_df$cell_type != "")))
message(sprintf("    fafb_match updated:   %d", sum(push_df$fafb_match != "")))

if (sum(retype_mask & has_manual) > 0) {
  manual_df <- candidates[retype_mask & has_manual, , drop = FALSE]
  manual_file <- "data/codex/bm_retype_manual_review.csv"
  readr::write_csv(manual_df, manual_file)
  message(sprintf("  Held back HAS_MANUAL_ANNOTATION: %d -> %s",
                  nrow(manual_df), manual_file))
}

# banctable_update_rows(base = 'banc_meta',
#                       table = "banc_meta",
#                       df = push_df,
#                       append_allowed = FALSE,
#                       chunksize = 200)
# message("  SeaTable update complete (BM retypes)")

message("\n### banc: BM asymmetry analysis complete ###")

})  # end local
