#' banc-sensory-classify — Classify proofread sensory neurons via cross-NBLAST against named GT.
#'
#' Classifies proofread sensory / afferent BANC neurons with placeholder /
#' missing cell_type by NBLAST against the BANC-native pool of named
#' sensory cell types. High-confidence GT matches → `NBLAST:<GT_type>`;
#' new clusters (intra-mean ≥ 0.20, n ≥ 3) → `NBLAST:SN<region>NN`.
#'
#' @section Reads:
#'   - SeaTable `banc_meta`
#'   - L2 SWCs under `<banc.l2swc.save.path>`
#'
#' @section Writes:
#'   - SeaTable `banc_meta`: col `cell_type` for matched / clustered neurons

###########################################################
### BANC sensory cell-type classification ###
###########################################################
# Classify proofread sensory/afferent BANC neurons that currently have
# placeholder or missing cell_type labels by NBLAST against the BANC-native
# pool of named sensory cell types.
#
# Inputs (pulled from banc_meta):
#   - Candidates: sensory/afferent neurons with placeholder cell_type
#       (NULL/blank, or contains "unknown" / "nerve" / "orphan", or ends "xx"), AND
#       (proofread = TRUE OR roughly_proofread = TRUE), AND
#       (nerve populated OR super_class contains "sensory" OR flow contains "afferent")
#   - Ground truth (GT): sensory/afferent neurons with a non-placeholder cell_type
#       (excludes the candidate patterns above + already NBLAST:* tagged), also proofread.
#
# Pipeline:
#   1. Pull both sets, fetch L2 skeletons
#   2. Mirror left-side neurons via banc_mirror(tpsreg) so L+R of a type co-cluster
#   3. nat::dotprops on all
#   4. Cross-NBLAST query=candidates vs target=GT(sampled per cell type)
#   5. Within-NBLAST among candidates (for clustering unmatched into new types)
#   6. Derive within-type / between-type score distributions from GT, choose threshold
#   7. Classify each candidate:
#         - high-confidence GT match (score >= TH_HI)  -> "NBLAST:<GT_type>"
#         - new cluster (intra mean >= 0.20, n >= 3)   -> "NBLAST:SN<region>NN"
#         - unassigned                                  -> kept as-is
#   8. Push to seatable with "NBLAST:" prefix
#         - only overwrite cell_type when it matches the candidate pattern (safe)
#         - new-cluster numbers continue from the highest already in seatable
#           per-prefix (cb/vnc/ol), with retry-on-rate-limit lookup
#
# Output: CSV with per-neuron classification, saved to
#   data/sensory_typing/banc_sensory_classification.csv
#
# To get the next-free SN<region>NN slot safely, use get_next_sn_numbers()
# below — it queries each prefix separately with retries to avoid the
# rate-limited combined-OR query.

source("banc/banc-startup.R")
suppressPackageStartupMessages({
  library(nat); library(nat.nblast); library(stringr)
})

out_dir <- "data/sensory_typing"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

###########################################################
### Helpers ###
###########################################################

# Returns the maximum trailing-digit number already used for a given
# SN<region> prefix (and its NBLAST:-prefixed variant). Robust to API
# rate-limits: retries the per-prefix query a few times.
get_next_sn_numbers <- function(prefixes = c("cb", "vnc", "ol"),
                                 retries = 4, sleep_seconds = 2) {
  out <- setNames(integer(length(prefixes)), prefixes)
  for (p in prefixes) {
    used_nums <- integer(0)
    for (variant in c(sprintf("SN%s%%", p), sprintf("NBLAST:SN%s%%", p))) {
      r <- NULL
      for (i in seq_len(retries)) {
        r <- tryCatch(
          banctable_query(sprintf("SELECT cell_type FROM banc_meta WHERE cell_type LIKE '%s'", variant)),
          error = function(e) NULL)
        if (!is.null(r) && nrow(r) > 0) break
        Sys.sleep(sleep_seconds * i)
      }
      if (!is.null(r) && nrow(r) > 0) {
        names_clean <- sub("^NBLAST:", "", unique(r$cell_type))
        hits <- as.integer(stringr::str_match(names_clean, sprintf("^SN%s([0-9]{2})$", p))[, 2])
        used_nums <- c(used_nums, hits[!is.na(hits)])
      }
    }
    out[p] <- if (!length(used_nums)) 0L else max(used_nums)
  }
  out
}

region_suffix <- function(r) {
  if (is.na(r) || r == "") return("xx")
  if (r == "central_brain")       return("cb")
  if (r == "ventral_nerve_cord")  return("vnc")
  if (r == "optic_lobe")          return("ol")
  "xx"
}

# Safety: only allow overwriting placeholder-pattern cell_types.
is_safe_to_overwrite <- function(x) {
  is.na(x) | x == "" |
    grepl("unknown", x, ignore.case = TRUE) |
    grepl("nerve",   x, ignore.case = TRUE) |
    grepl("orphan",  x, ignore.case = TRUE) |
    grepl("xx$",     x)
}

###################
### 1. Queries ###
###################

cand_sql <- paste(
  "SELECT _id, root_626, root_id, cell_type, super_class, flow, side, nerve, region,",
  "       proofread, roughly_proofread, status",
  "  FROM banc_meta WHERE",
  "(cell_type IS NULL OR cell_type = ''",
  "   OR cell_type LIKE '%unknown%'",
  "   OR cell_type LIKE '%nerve%'",
  "   OR cell_type LIKE '%orphan%'",
  "   OR cell_type LIKE '%xx')",
  "AND (proofread = 'TRUE' OR roughly_proofread = 'TRUE')",
  "AND (nerve IS NOT NULL OR super_class LIKE '%sensory%' OR flow LIKE '%afferent%')"
)
cand <- banctable_query(cand_sql) %>%
  dplyr::mutate(across(c(`_id`, root_626, root_id), as.character)) %>%
  dplyr::filter(!is.na(root_id), root_id != "", root_id != "0") %>%
  dplyr::distinct(root_id, .keep_all = TRUE)
message(sprintf("Candidates: %d", nrow(cand)))

gt_sql <- paste(
  "SELECT _id, root_626, root_id, cell_type, super_class, side, flow",
  "  FROM banc_meta WHERE",
  "(super_class LIKE '%sensory%' OR flow LIKE '%afferent%')",
  "AND cell_type IS NOT NULL AND cell_type != ''",
  "AND cell_type NOT LIKE '%unknown%' AND cell_type NOT LIKE '%nerve%'",
  "AND cell_type NOT LIKE '%orphan%' AND cell_type NOT LIKE '%xx'",
  "AND cell_type NOT LIKE 'NBLAST:%'",
  "AND (proofread = 'TRUE' OR roughly_proofread = 'TRUE')"
)
gt_all <- banctable_query(gt_sql) %>%
  dplyr::mutate(across(c(`_id`, root_626, root_id), as.character)) %>%
  dplyr::filter(!is.na(root_id), root_id != "", root_id != "0") %>%
  dplyr::distinct(root_id, .keep_all = TRUE)
message(sprintf("GT all: %d across %d cell types", nrow(gt_all), dplyr::n_distinct(gt_all$cell_type)))

# Subsample GT per cell type: keep up to N_PER_TYPE per type, balanced L/R
N_PER_TYPE <- 8
set.seed(42)
gt_sub <- gt_all %>%
  dplyr::group_by(cell_type) %>%
  dplyr::slice_sample(n = N_PER_TYPE) %>%
  dplyr::ungroup()
message(sprintf("GT subsampled: %d", nrow(gt_sub)))

# Combined ids
all_ids <- unique(c(cand$root_id, gt_sub$root_id))
message(sprintf("Total to NBLAST: %d (cand %d + GT %d)", length(all_ids), nrow(cand), nrow(gt_sub)))

#######################
### 2. Fetch skels ###
#######################

# Pull from cache where possible
get_cache <- function() {
  pool <- list()
  for (rds in c("/tmp/banc_mirror_skels.rds", "/tmp/banc_sn_nblast_skels.rds",
                "/tmp/banc_dcv_nblast_skels.rds", "/tmp/banc_ph_skels.rds")) {
    if (!file.exists(rds)) next
    nl <- readRDS(rds)
    for (n in names(nl)) if (is.null(pool[[n]])) pool[[n]] <- nl[[n]]
  }
  do.call(nat::as.neuronlist, list(pool))
}
cache <- get_cache()
have <- intersect(names(cache), all_ids)
need <- setdiff(all_ids, have)
message(sprintf("Cached: %d  Need to fetch: %d", length(have), length(need)))

if (length(need)) {
  fresh <- bancr::banc_read_l2skel(need, OmitFailures = TRUE)
  for (n in names(fresh)) cache[[n]] <- fresh[[n]]
}
skels <- cache[intersect(names(cache), all_ids)]
message(sprintf("Skels available: %d / %d", length(skels), length(all_ids)))
saveRDS(skels, "/tmp/banc_sensory_skels.rds")

#####################
### 3. Mirror L→R ###
#####################

side_lookup <- setNames(
  c(cand$side, gt_sub$side),
  c(cand$root_id, gt_sub$root_id)
)
left_ids <- intersect(names(skels), names(side_lookup)[!is.na(side_lookup) & side_lookup == "left"])
right_ids <- setdiff(names(skels), left_ids)
message(sprintf("Left: %d  Right/NA: %d", length(left_ids), length(right_ids)))

t_mir <- Sys.time()
left_m <- bancr::banc_mirror(skels[left_ids], banc.units = "nm", method = "tpsreg")
message(sprintf("Mirror done in %.1f min", as.numeric(difftime(Sys.time(), t_mir, units = "mins"))))

unified <- nat::as.neuronlist(c(as.list(left_m), as.list(skels[right_ids])))
names(unified) <- c(left_ids, right_ids)
saveRDS(unified, "/tmp/banc_sensory_unified.rds")

##########################################
### 4. Dotprops + NBLAST query vs target
##########################################

um <- nat::xform(unified, function(p, ...) p / 1e3)
dps <- nat::dotprops(um, k = 5, resample = 1, OmitFailures = TRUE)
message(sprintf("Dotprops: %d", length(dps)))
saveRDS(dps, "/tmp/banc_sensory_dps.rds")

cand_dps <- dps[intersect(names(dps), cand$root_id)]
gt_dps   <- dps[intersect(names(dps), gt_sub$root_id)]
message(sprintf("Cand dps: %d  GT dps: %d", length(cand_dps), length(gt_dps)))

# Cross NBLAST: candidates against GT (smaller of the two as query for speed)
t_nbl <- Sys.time()
cross <- nat.nblast::nblast(query = cand_dps, target = gt_dps, normalised = TRUE)
message(sprintf("Cross NBLAST: %dx%d  (%.1f min)", nrow(cross), ncol(cross),
                as.numeric(difftime(Sys.time(), t_nbl, units = "mins"))))

# Within-candidate NBLAST for clustering unmatched ones
t_nbl2 <- Sys.time()
within <- nat.nblast::nblast_allbyall(cand_dps, normalisation = "mean")
message(sprintf("Within-candidate NBLAST: %d x %d  (%.1f min)",
                nrow(within), ncol(within),
                as.numeric(difftime(Sys.time(), t_nbl2, units = "mins"))))

# GT-within: needed to compute threshold
t_nbl3 <- Sys.time()
gt_self <- nat.nblast::nblast_allbyall(gt_dps, normalisation = "mean")
message(sprintf("GT-self NBLAST: %dx%d  (%.1f min)", nrow(gt_self), ncol(gt_self),
                as.numeric(difftime(Sys.time(), t_nbl3, units = "mins"))))

saveRDS(list(cross = cross, within = within, gt_self = gt_self,
             cand = cand, gt_sub = gt_sub),
        "/tmp/banc_sensory_nblast.rds")

############################
### 5. Threshold from GT ###
############################

# Drop "xx-suffixed" cell types from GT (placeholder-y) before deriving the
# within-type / between-type distributions used for thresholding.
gt_clean <- gt_sub %>% dplyr::filter(!grepl("xx$", cell_type))
gt_groups <- split(gt_clean$root_id, gt_clean$cell_type)
gt_groups <- gt_groups[vapply(gt_groups, length, integer(1)) >= 2]
within_pairs <- c()
for (g in gt_groups) {
  sub <- gt_self[g, g]; diag(sub) <- NA
  within_pairs <- c(within_pairs, as.numeric(sub[upper.tri(sub)]))
}
within_pairs <- within_pairs[is.finite(within_pairs)]

named_ids <- unlist(gt_groups, use.names = FALSE)
name_of <- rep(names(gt_groups), vapply(gt_groups, length, integer(1)))
set.seed(7)
ii <- sample.int(length(named_ids), 20000, replace = TRUE)
jj <- sample.int(length(named_ids), 20000, replace = TRUE)
keep <- ii != jj & name_of[ii] != name_of[jj]
between_pairs <- gt_self[cbind(named_ids[ii][keep], named_ids[jj][keep])]
between_pairs <- between_pairs[is.finite(between_pairs)]

# Youden's J + high-precision
thr_grid <- quantile(c(within_pairs, between_pairs), seq(0.05, 0.99, 0.01))
best <- list(t = NA, j = -Inf)
for (t in thr_grid) {
  j <- mean(within_pairs >= t) - mean(between_pairs >= t)
  if (j > best$j) best <- list(t = t, j = j)
}
TH    <- best$t
TH_HI <- min(thr_grid[vapply(thr_grid, function(t) mean(between_pairs >= t) <= 0.05, logical(1))])
message(sprintf("Youden threshold: %.3f  High-precision (FPR<=0.05): %.3f", TH, TH_HI))

###################################
### 6. Classify each candidate ###
###################################

classify_one <- function(rid) {
  # nat.nblast::nblast(query, target) returns target-rows x query-cols, so a
  # candidate (query) lives in a COLUMN, and GT members live in ROWS.
  if (!rid %in% colnames(cross)) return(list(best="", score=NA_real_, runner="", runner_score=NA_real_))
  s <- cross[, rid]
  type_means <- vapply(gt_groups, function(g) mean(s[intersect(g, rownames(cross))], na.rm = TRUE),
                       numeric(1))
  ord <- order(type_means, decreasing = TRUE)
  list(best = names(type_means)[ord[1]],
       score = unname(type_means[ord[1]]),
       runner = if (length(ord) > 1) names(type_means)[ord[2]] else "",
       runner_score = if (length(ord) > 1) unname(type_means[ord[2]]) else NA_real_)
}
res <- lapply(cand$root_id, classify_one)
df_class <- data.frame(
  root_626 = cand$root_626,
  root_id = cand$root_id,
  cur_cell_type = cand$cell_type,
  super_class = cand$super_class,
  flow = cand$flow,
  side = cand$side,
  nerve = cand$nerve,
  best_match = vapply(res, function(x) x$best, character(1)),
  best_score = vapply(res, function(x) x$score, numeric(1)),
  runner = vapply(res, function(x) x$runner, character(1)),
  runner_score = vapply(res, function(x) x$runner_score, numeric(1)),
  stringsAsFactors = FALSE
) %>% dplyr::mutate(
  decision = dplyr::case_when(
    is.na(best_score)    ~ "no_data",
    best_score >= TH_HI  ~ "high",
    best_score >= TH     ~ "medium",
    TRUE                 ~ "no_match"
  )
)

message("Decisions:")
print(table(df_class$decision))
print(df_class %>% dplyr::filter(decision == "high") %>% dplyr::count(best_match, sort = TRUE) %>% head(15))

##########################################
### 7. Cluster unmatched -> new types ###
##########################################

unmatched_ids <- df_class$root_id[df_class$decision %in% c("no_match", "medium")]
new_clusters <- list()
if (length(unmatched_ids) >= 3) {
  sub <- within[unmatched_ids, unmatched_ids]
  hc <- stats::hclust(as.dist(1 - sub), method = "average")
  best <- list()
  for (k in seq(2, min(30, length(unmatched_ids) - 1))) {
    ct <- stats::cutree(hc, k = k)
    g_list <- split(names(ct), ct)
    for (g in g_list) {
      if (length(g) < 3) next
      m <- sub[g, g]; diag(m) <- NA
      intra <- mean(m, na.rm = TRUE)
      if (intra < TH) next
      key <- paste(sort(g), collapse=",")
      if (is.null(best[[key]]) || intra > best[[key]]$intra) {
        best[[key]] <- list(members = g, intra = intra)
      }
    }
  }
  best <- best[order(-vapply(best, function(x) x$intra, numeric(1)))]
  taken <- character(0)
  for (c in best) {
    if (any(c$members %in% taken)) next
    new_clusters[[length(new_clusters) + 1]] <- c
    taken <- c(taken, c$members)
  }
}
message(sprintf("New coherent clusters among unmatched: %d", length(new_clusters)))
for (i in seq_along(new_clusters)) {
  c1 <- new_clusters[[i]]
  message(sprintf("  cluster %d  n=%d  intra=%.3f", i, length(c1$members), c1$intra))
}

# Region-aware new-cluster naming: each cluster's region is the modal region
# of its members; the new name continues the SN<region>NN counter at the next
# free number per prefix (queried with retries, see get_next_sn_numbers()).
df_class$new_cluster <- NA_character_
df_class$region <- cand$region[match(df_class$root_id, cand$root_id)]

next_sn <- if (length(new_clusters)) {
  message("Querying next-free SN<region>NN numbers per prefix (cb/vnc/ol)...")
  get_next_sn_numbers(c("cb", "vnc", "ol"))
} else {
  c(cb = 0L, vnc = 0L, ol = 0L)
}
message("Highest existing SN<prefix>NN per region: ",
        paste(sprintf("%s=%d", names(next_sn), next_sn), collapse = ", "))

# Sort new clusters by intra (best first) so the "best" clusters get earlier IDs
new_clusters <- new_clusters[order(-vapply(new_clusters, function(x) x$intra, numeric(1)))]
counter <- next_sn  # per-prefix running counter (starts at max-used)
for (i in seq_along(new_clusters)) {
  members <- new_clusters[[i]]$members
  member_regions <- df_class$region[df_class$root_id %in% members]
  modal_region <- {
    tt <- table(member_regions[!is.na(member_regions) & member_regions != ""])
    if (length(tt)) names(sort(tt, decreasing = TRUE))[1] else NA_character_
  }
  suf <- region_suffix(modal_region)
  if (suf == "xx") {
    name_i <- sprintf("SNxx%02d", i)  # unrecognized region; user will inspect
  } else {
    counter[suf] <- counter[suf] + 1L
    name_i <- sprintf("SN%s%02d", suf, counter[suf])
  }
  for (rid in members) {
    df_class$new_cluster[df_class$root_id == rid] <- name_i
  }
}

#######################
### 8. Save results ###
#######################

readr::write_csv(df_class, file.path(out_dir, "banc_sensory_classification.csv"))
message(sprintf("Wrote: %s", file.path(out_dir, "banc_sensory_classification.csv")))

# Summary
message("\n=== Final per-class summary ===")
message(sprintf("  high-confidence match:  %d", sum(df_class$decision == "high")))
message(sprintf("  medium-confidence:      %d", sum(df_class$decision == "medium")))
message(sprintf("  no_match:               %d", sum(df_class$decision == "no_match")))
message(sprintf("  no_data:                %d", sum(df_class$decision == "no_data")))
message(sprintf("  in new clusters:        %d (across %d clusters)",
                sum(!is.na(df_class$new_cluster)), length(new_clusters)))

#######################
### 9. Seatable push ###
#######################
# Push proposed type as "NBLAST:<name>" with these rules:
#   - high-confidence GT match   -> NBLAST:<gt_type>
#   - new-cluster member          -> NBLAST:<SN<region>NN>
# Safety: only overwrite cell_type that matches the candidate pattern
# (blank / contains unknown|nerve|orphan / ends "xx"). All other current
# values are skipped to avoid clobbering manual annotations made in between.
#
# Set DO_PUSH=TRUE explicitly (e.g. via env) to actually write — default is dry.

DO_PUSH <- isTRUE(as.logical(Sys.getenv("SENSORY_NBLAST_PUSH", "FALSE")))

# Re-pull current cell_type values for the candidate set so the safety check
# uses fresh state (the user may have hand-fixed some between prep and push).
fresh <- banctable_query(
  sprintf("SELECT _id, root_id, cell_type FROM banc_meta WHERE root_id IN ('%s')",
          paste(unique(df_class$root_id), collapse = "','"))
) %>%
  dplyr::mutate(across(c(`_id`, root_id), as.character)) %>%
  dplyr::distinct(root_id, .keep_all = TRUE)

push_df <- df_class %>%
  dplyr::mutate(
    proposed = dplyr::case_when(
      decision == "high"      ~ paste0("NBLAST:", best_match),
      !is.na(new_cluster)     ~ paste0("NBLAST:", new_cluster),
      TRUE                    ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(proposed)) %>%
  dplyr::left_join(fresh, by = "root_id", suffix = c("", ".fresh")) %>%
  dplyr::filter(is_safe_to_overwrite(cell_type.fresh)) %>%
  dplyr::transmute(`_id` = `_id`, root_626 = root_626, cell_type = proposed)

message(sprintf("\nPush candidates after safety filter: %d", nrow(push_df)))
print(head(push_df, 10))

if (DO_PUSH && nrow(push_df) >= 1) {
  message("Pushing to seatable (3-col: _id + root_626 + cell_type)...")
  bancr::banctable_update_rows(
    table = "banc_meta",
    df = as.data.frame(push_df),
    append_allowed = FALSE,
    chunksize = 100L
  )
  message("Push complete.")
} else {
  message("Dry run (set SENSORY_NBLAST_PUSH=TRUE to push). No seatable writes performed.")
}

message("\nDone. CSV at: ", file.path(out_dir, "banc_sensory_classification.csv"))
