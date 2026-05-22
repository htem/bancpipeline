#' banc-class-fixes — Propagate classification columns for bates-sourced neurons.
#'
#' For BANC neurons whose `cell_type_source` contains "bates", pull
#' downstream classification columns from the matched FAFB neuron in
#' `franken_meta` (primary) or the matched maleCNS neuron (fallback). Writes
#' the propagation plan to CSV; SeaTable push is COMMENTED OUT pending review.
#'
#' @section Reads:
#'   - SeaTable `banc_meta` (live)
#'   - `franken_meta()` (FAFB target follower columns)
#'   - maleCNS SeaTable (fallback target)
#'
#' @section Writes:
#'   - `data/codex/class_fixes.csv` — proposed per-neuron column updates

###########################################################
### Propagate classification columns for bates-sourced neurons
###
### For BANC neurons where cell_type_source contains "bates",
### look up the fafb_match in franken_meta (primary) or
### malecns_match in malecns seatable (fallback) and pull
### follower columns from that specific matched neuron.
###
### Columns updated: flow, super_class, cell_class, cell_sub_class,
### cell_function, cell_function_detailed, peripheral_target_type,
### body_part_sensory, body_part_effector, hemilineage,
### neurotransmitter_verified, neuropeptide_verified,
### sexually_dimorphic
###
### Output: data/codex/class_fixes.csv
###########################################################
source("banc/banc-startup.R")

message("### banc: classification column propagation for bates-sourced neurons ###")

# Columns to propagate (seatable names)
target_cols <- c("flow", "super_class", "cell_class", "cell_sub_class",
                 "cell_function", "cell_function_detailed",
                 "peripheral_target_type", "body_part_sensory",
                 "body_part_effector", "hemilineage",
                 "neurotransmitter_verified", "neuropeptide_verified",
                 "sexually_dimorphic")

###########################
### Read BANC data      ###
###########################

select_cols <- paste(c("_id", "root_888", "root_id", "cell_type", "cell_type_source",
                        "fafb_match", "fafb_cell_type", "hemibrain_cell_type", "fanc_cell_type",
                        "malecns_match", "malecns_cell_type",
                        "manc_match", "manc_cell_type", "side", "status", target_cols), collapse = ", ")
bc <- banctable_query(sprintf("SELECT %s FROM banc_meta", select_cols)) %>%
  dplyr::mutate(dplyr::across(dplyr::everything(), as.character))

# Filter to bates-sourced neurons with a cell_type
bates <- bc %>%
  dplyr::filter(grepl("bates", cell_type_source, ignore.case = TRUE),
                !is.na(cell_type), cell_type != "")
message(sprintf("  bates-sourced neurons with cell_type: %d", nrow(bates)))

###########################
### Build lookup tables ###
###########################

# Primary: franken_meta — look up by fafb_match -> fafb_id
# Use banc_neurotransmitter_verified / banc_neuropeptide_verified
message("  Loading franken_meta...")
fm <- franken_meta() %>%
  dplyr::filter(!is.na(fafb_id), fafb_id != "") %>%
  dplyr::mutate(fafb_id = as.character(fafb_id))

# Rename BANC-specific columns to match seatable column names
if ("banc_neurotransmitter_verified" %in% names(fm)) {
  fm <- fm %>% dplyr::mutate(neurotransmitter_verified = banc_neurotransmitter_verified)
}
if ("banc_neuropeptide_verified" %in% names(fm)) {
  fm <- fm %>% dplyr::mutate(neuropeptide_verified = banc_neuropeptide_verified)
}

# Build lookup keyed by fafb_id (one row per fafb_id)
fm_cols <- intersect(target_cols, names(fm))
fm_lookup <- fm %>%
  dplyr::distinct(fafb_id, .keep_all = TRUE) %>%
  dplyr::select(fafb_id, fm_cell_type = cell_type, dplyr::all_of(fm_cols))
message(sprintf("  franken_meta lookup: %d fafb_ids, %d columns", nrow(fm_lookup), length(fm_cols)))

# Also build cell_type -> fafb_ids index for NBLAST-based suggestions
fm_type_ids <- split(fm_lookup$fafb_id, fm_lookup$fm_cell_type)

# Build lookup keyed by manc_id for ascending/VNC neurons
fm_manc <- fm %>%
  dplyr::filter(!is.na(manc_id), manc_id != "") %>%
  dplyr::mutate(manc_id = as.character(manc_id)) %>%
  dplyr::distinct(manc_id, .keep_all = TRUE) %>%
  dplyr::select(manc_id, fm_cell_type_manc = cell_type, dplyr::all_of(fm_cols))
message(sprintf("  franken_meta manc lookup: %d manc_ids", nrow(fm_manc)))

# Fallback: malecns seatable — look up by malecns_match -> malecns_09_id
message("  Loading malecns seatable...")
mcns <- tryCatch(
  banctable_query("SELECT * FROM malecns", base = "cns_meta") %>%
    dplyr::filter(!is.na(malecns_09_id), malecns_09_id != "") %>%
    dplyr::mutate(malecns_09_id = as.character(malecns_09_id)),
  error = function(e) {
    warning("malecns query failed: ", e$message)
    NULL
  }
)

mcns_lookup <- NULL
if (!is.null(mcns)) {
  mcns_cols <- intersect(target_cols, names(mcns))
  mcns_lookup <- mcns %>%
    dplyr::distinct(malecns_09_id, .keep_all = TRUE) %>%
    dplyr::select(malecns_09_id, mcns_cell_type = cell_type, dplyr::all_of(mcns_cols))
  message(sprintf("  malecns lookup: %d IDs, %d columns", nrow(mcns_lookup), length(mcns_cols)))
}

# MANC lookup. NB: bancr::franken_meta() does not currently populate the
# manc_id column (verified: 0 / 166k rows have a manc_id), so the fm_manc
# path that the bates VNC loop relies on is effectively dead. Load MANC
# metadata from the local CSV instead — it is the same source the rest of
# the pipeline uses (manc/manc-meta.R writes it from neuprint).
manc_lookup <- NULL
manc_csv <- file.path(banc.meta.save.path, "manc_meta.csv")
if (file.exists(manc_csv)) {
  message("  Loading MANC metadata from local CSV...")
  mc.meta <- suppressWarnings(readr::read_csv(
    manc_csv, col_types = hemibrainr:::sql_col_types, show_col_types = FALSE))
  manc_cols <- intersect(target_cols, names(mc.meta))
  manc_lookup <- mc.meta %>%
    dplyr::filter(!is.na(bodyid)) %>%
    dplyr::mutate(bodyid = as.character(bodyid)) %>%
    dplyr::distinct(bodyid, .keep_all = TRUE) %>%
    dplyr::select(manc_id = bodyid,
                  manc_cell_type = cell_type,
                  dplyr::all_of(manc_cols))
  message(sprintf("  manc lookup: %d IDs, %d follower columns",
                  nrow(manc_lookup), length(manc_cols)))
} else {
  message(sprintf("  manc CSV not found at %s — manc source path inactive",
                  manc_csv))
}

# Load FAFB NBLAST for suggesting better matches
message("  Loading FAFB NBLAST feather...")
nblast_cache <- file.path(tempdir(), "nblast_cache")
dir.create(nblast_cache, showWarnings = FALSE)
nblast_file <- file.path(nblast_cache, "banc_fafb_783_nblast.feather")
if (!file.exists(nblast_file)) {
  alt_file <- file.path("/tmp/nblast_cache", "banc_fafb_783_nblast.feather")
  if (file.exists(alt_file)) nblast_file <- alt_file
  else {
    message("    Downloading...")
    system2("gsutil", c("cp",
      "gs://lee-lab_brain-and-nerve-cord-fly-connectome/nblast/banc_fafb_783_nblast.feather",
      nblast_file), stdout = FALSE, stderr = FALSE)
  }
}
fafb_nblast <- arrow::read_feather(nblast_file)
message(sprintf("  FAFB NBLAST: %d rows", nrow(fafb_nblast)))

# Pre-index NBLAST by root_888
message("  Pre-indexing NBLAST...")
fafb_nblast_idx <- split(
  data.frame(match_id = as.character(fafb_nblast$match_id),
             score = fafb_nblast$score, stringsAsFactors = FALSE),
  as.character(fafb_nblast$root_888)
)
message(sprintf("  Index: %d neurons", length(fafb_nblast_idx)))

###########################
### Apply fixes         ###
###########################

message("  Matching...")

results <- bates
results$source <- NA_character_
results$fafb_match_cell_type <- NA_character_  # cell_type of the fafb_match in franken_meta
n_fm_fafb <- 0L
n_fm_manc <- 0L
n_mcns <- 0L
n_none <- 0L

for (i in seq_len(nrow(results))) {
  fafb_id <- results$fafb_match[i]
  mcns_id <- results$malecns_match[i]
  manc_id <- results$manc_match[i]
  sc <- results$super_class[i]
  st <- results$status[i]

  # Apply X_MATCH_PREFERRED: override cell_type from the preferred dataset's cell_type column
  if (!is.na(st) && grepl("_MATCH_PREFERRED", st)) {
    pref_map <- c(
      FAFB_MATCH_PREFERRED      = "fafb_cell_type",
      HEMIBRAIN_MATCH_PREFERRED  = "hemibrain_cell_type",
      MALECNS_MATCH_PREFERRED    = "malecns_cell_type",
      MANC_MATCH_PREFERRED       = "manc_cell_type",
      FANC_MATCH_PREFERRED       = "fanc_cell_type"
    )
    for (pref_status in names(pref_map)) {
      if (grepl(pref_status, st, fixed = TRUE)) {
        pref_col <- pref_map[[pref_status]]
        if (pref_col %in% names(results)) {
          pref_val <- results[[pref_col]][i]
          if (!is.na(pref_val) && pref_val != "") {
            results$cell_type[i] <- pref_val
          }
        }
        break
      }
    }
  }

  # For ascending / VNC neurons: use malecns_cell_type or manc_cell_type as cell_type,
  # and pull follower columns from manc_match -> manc_id in franken_meta
  is_vnc <- !is.na(sc) && grepl("ascending|ventral_nerve_cord", sc, ignore.case = TRUE)

  if (is_vnc) {
    # Determine effective cell_type: prefer malecns_cell_type, fallback to manc_cell_type
    vnc_ct <- results$malecns_cell_type[i]
    if (is.na(vnc_ct) || vnc_ct == "") vnc_ct <- results$manc_cell_type[i]
    if (!is.na(vnc_ct) && vnc_ct != "") results$cell_type[i] <- vnc_ct

    # Preserve sexually_dimorphic — BANC seatable is authoritative for VNC neurons
    orig_sd <- results$sexually_dimorphic[i]

    # Pull follower columns from manc_match -> franken_meta manc_id
    if (!is.na(manc_id) && manc_id != "") {
      fm_row <- fm_manc[fm_manc$manc_id == manc_id, ]
      if (nrow(fm_row) == 1) {
        for (col in fm_cols) {
          val <- fm_row[[col]][1]
          results[[col]][i] <- if (!is.na(val) && val != "") val else ""
        }
        results$sexually_dimorphic[i] <- orig_sd
        results$source[i] <- "franken_meta_manc"
        n_fm_manc <- n_fm_manc + 1L
        next
      }
    }

    # Fallback: malecns seatable via malecns_match
    if (!is.null(mcns_lookup) && !is.na(mcns_id) && mcns_id != "") {
      mcns_row <- mcns_lookup[mcns_lookup$malecns_09_id == mcns_id, ]
      if (nrow(mcns_row) == 1) {
        for (col in intersect(target_cols, names(mcns_row))) {
          val <- mcns_row[[col]][1]
          results[[col]][i] <- if (!is.na(val) && val != "") val else ""
        }
        results$sexually_dimorphic[i] <- orig_sd
        results$source[i] <- "malecns"
        n_mcns <- n_mcns + 1L
        next
      }
    }

    results$source[i] <- "no_match"
    n_none <- n_none + 1L
    next
  }

  # Non-VNC: try franken_meta via fafb_match
  if (!is.na(fafb_id) && fafb_id != "") {
    fm_row <- fm_lookup[fm_lookup$fafb_id == fafb_id, ]
    if (nrow(fm_row) == 1) {
      results$fafb_match_cell_type[i] <- fm_row$fm_cell_type[1]
      for (col in fm_cols) {
        val <- fm_row[[col]][1]
        # Overwrite unconditionally — blank out columns the source doesn't have
        results[[col]][i] <- if (!is.na(val) && val != "") val else ""
      }
      results$source[i] <- "franken_meta_fafb"
      n_fm_fafb <- n_fm_fafb + 1L
      next
    }
  }

  # Fallback: malecns via malecns_match
  if (!is.null(mcns_lookup) && !is.na(mcns_id) && mcns_id != "") {
    mcns_row <- mcns_lookup[mcns_lookup$malecns_09_id == mcns_id, ]
    if (nrow(mcns_row) == 1) {
      for (col in intersect(target_cols, names(mcns_row))) {
        val <- mcns_row[[col]][1]
        results[[col]][i] <- if (!is.na(val) && val != "") val else ""
      }
      results$source[i] <- "malecns"
      n_mcns <- n_mcns + 1L
      next
    }
  }

  results$source[i] <- "no_match"
  n_none <- n_none + 1L
}

###########################
### Report              ###
###########################

message(sprintf("\n  === Results ==="))
message(sprintf("  Total bates-sourced neurons: %d", nrow(results)))
message(sprintf("    franken_meta via fafb_match: %d", n_fm_fafb))
message(sprintf("    franken_meta via manc_match (ascending/VNC): %d", n_fm_manc))
message(sprintf("    malecns seatable: %d", n_mcns))
message(sprintf("    no match: %d", n_none))

# Count how many columns actually change per row
matched <- results %>% dplyr::filter(source != "no_match")
n_changes <- 0L
for (col in target_cols) {
  if (!col %in% names(bates)) next
  old_vals <- bates[[col]][results$source != "no_match"]
  new_vals <- matched[[col]]
  changed <- (!is.na(new_vals) & new_vals != "") &
             (is.na(old_vals) | old_vals == "" | old_vals != new_vals)
  n_col_changes <- sum(changed, na.rm = TRUE)
  if (n_col_changes > 0) {
    message(sprintf("    %s: %d changes", col, n_col_changes))
    n_changes <- n_changes + n_col_changes
  }
}
message(sprintf("  Total cell-level changes: %d", n_changes))

# Flag potential issues: super_class or cell_class changed from a non-blank value
check_cols <- c("super_class", "cell_class")
issue_rows <- logical(nrow(matched))
for (col in check_cols) {
  if (!col %in% names(bates)) next
  old_vals <- bates[[col]][results$source != "no_match"]
  new_vals <- matched[[col]]
  changed <- (!is.na(old_vals) & old_vals != "") &
             (!is.na(new_vals) & new_vals != "") &
             (old_vals != new_vals)
  issue_rows <- issue_rows | changed
}

potential_issues <- data.frame(
  root_888 = matched$root_888[issue_rows],
  root_id = matched$root_id[issue_rows],
  cell_type = matched$cell_type[issue_rows],
  fafb_match = matched$fafb_match[issue_rows],
  source = matched$source[issue_rows],
  old_super_class = bates$super_class[results$source != "no_match"][issue_rows],
  new_super_class = matched$super_class[issue_rows],
  old_cell_class = bates$cell_class[results$source != "no_match"][issue_rows],
  new_cell_class = matched$cell_class[issue_rows],
  stringsAsFactors = FALSE
)

message(sprintf("\n  === Potential issues (super_class or cell_class changed from non-blank): %d ===",
                nrow(potential_issues)))
if (nrow(potential_issues) > 0) {
  issues_file <- "data/codex/class_fixes_potential_issues.csv"
  readr::write_csv(potential_issues, issues_file)
  message(sprintf("  Saved to %s", issues_file))
  for (i in seq_len(min(20, nrow(potential_issues)))) {
    r <- potential_issues[i, ]
    message(sprintf("  %s  %s  super: %s -> %s  class: %s -> %s",
                    r$root_id, r$cell_type,
                    r$old_super_class, r$new_super_class,
                    r$old_cell_class, r$new_cell_class))
  }
}

###############################################
### fafb_match type mismatches             ###
### cell_type != fafb_match's cell_type    ###
### Suggest better fafb_match via NBLAST   ###
###############################################

fm_matched <- results %>%
  dplyr::filter(source == "franken_meta_fafb",
                !is.na(fafb_match_cell_type), fafb_match_cell_type != "",
                !is.na(cell_type), cell_type != "",
                cell_type != fafb_match_cell_type)

message(sprintf("\n  === fafb_match type mismatches (cell_type != fafb_match cell_type): %d ===",
                nrow(fm_matched)))

if (nrow(fm_matched) > 0) {
  # For each, find best NBLAST hit to a FAFB neuron of the correct cell_type
  fm_matched$suggested_fafb_match <- NA_character_
  fm_matched$suggested_score <- NA_real_

  pb <- txtProgressBar(min = 0, max = nrow(fm_matched), style = 3)
  for (i in seq_len(nrow(fm_matched))) {
    setTxtProgressBar(pb, i)
    rid <- fm_matched$root_888[i]
    target_type <- fm_matched$cell_type[i]

    # Get FAFB IDs of the correct type
    correct_ids <- fm_type_ids[[target_type]]
    if (is.null(correct_ids) || length(correct_ids) == 0) next

    # Get NBLAST hits for this BANC neuron
    hits <- fafb_nblast_idx[[rid]]
    if (is.null(hits) || nrow(hits) == 0) next

    # Find best hit of correct type
    type_mask <- hits$match_id %in% correct_ids
    if (any(type_mask)) {
      type_hits <- hits[type_mask, ]
      best_idx <- which.max(type_hits$score)
      fm_matched$suggested_fafb_match[i] <- type_hits$match_id[best_idx]
      fm_matched$suggested_score[i] <- type_hits$score[best_idx]
    }
  }
  close(pb)

  fafb_mismatch_df <- fm_matched %>%
    dplyr::select(root_888, root_id, cell_type, fafb_match, fafb_match_cell_type,
                  suggested_fafb_match, suggested_score)

  n_suggested <- sum(!is.na(fafb_mismatch_df$suggested_fafb_match))
  message(sprintf("  With NBLAST suggestion: %d / %d", n_suggested, nrow(fafb_mismatch_df)))

  mismatch_file <- "data/codex/class_fixes_fafb_match_mismatches.csv"
  readr::write_csv(fafb_mismatch_df, mismatch_file)
  message(sprintf("  Saved to %s", mismatch_file))

  for (i in seq_len(min(15, nrow(fafb_mismatch_df)))) {
    r <- fafb_mismatch_df[i, ]
    message(sprintf("  %s  banc=%s  fafb_match_type=%s  suggested=%s (score=%.2f)",
                    r$root_id, r$cell_type, r$fafb_match_cell_type,
                    ifelse(is.na(r$suggested_fafb_match), "none", r$suggested_fafb_match),
                    ifelse(is.na(r$suggested_score), 0, r$suggested_score)))
  }

  # SeaTable update: push corrected fafb_match + fafb_cell_type for neurons with suggestions
  suggested <- fm_matched %>%
    dplyr::filter(!is.na(suggested_fafb_match), suggested_fafb_match != "")

  if (nrow(suggested) > 0) {
    push_fafb_fix <- suggested %>%
      dplyr::transmute(
        `_id`,
        fafb_match = suggested_fafb_match,
        fafb_cell_type = cell_type,
        cell_type_source = append_status(cell_type_source, "bates")
      ) %>%
      as.data.frame()

    push_fafb_fix$fafb_match[is.na(push_fafb_fix$fafb_match)] <- ""
    push_fafb_fix$fafb_cell_type[is.na(push_fafb_fix$fafb_cell_type)] <- ""

    message(sprintf("\n  SeaTable push (fafb_match fixes): %d rows", nrow(push_fafb_fix)))
    message(sprintf("    fafb_match updated: %d", sum(push_fafb_fix$fafb_match != "")))
    message(sprintf("    fafb_cell_type updated: %d", sum(push_fafb_fix$fafb_cell_type != "")))

    # banctable_update_rows(base = 'banc_meta',
    #                       table = "banc_meta",
    #                       df = push_fafb_fix,
    #                       append_allowed = FALSE,
    #                       chunksize = 200)
    # message("  SeaTable update complete (fafb_match fixes)")
  }
}

###########################
### Save CSV            ###
###########################

out_file <- "data/codex/class_fixes.csv"
readr::write_csv(results, out_file)
message(sprintf("\n  Saved %d rows to %s", nrow(results), out_file))

actionable <- results %>% dplyr::filter(source != "no_match")
act_file <- "data/codex/class_fixes_actionable.csv"
readr::write_csv(actionable, act_file)
message(sprintf("  Saved %d actionable rows to %s", nrow(actionable), act_file))

###########################
### SeaTable update     ###
###########################

# Separate manual-annotation neurons for careful review
has_manual <- grepl("HAS_MANUAL_ANNOTATION", actionable$status, ignore.case = TRUE)

push_cols <- c("_id", target_cols)

# Main push: neurons without HAS_MANUAL_ANNOTATION
push_df <- actionable[!has_manual, ] %>%
  dplyr::select(dplyr::any_of(push_cols)) %>%
  as.data.frame()
for (col in target_cols) {
  if (col %in% names(push_df)) push_df[[col]][is.na(push_df[[col]])] <- ""
}

# Manual-annotation neurons: separate for review
push_manual <- actionable[has_manual, ] %>%
  dplyr::select(dplyr::any_of(c(push_cols, "root_888", "root_id", "cell_type", "status"))) %>%
  as.data.frame()
for (col in target_cols) {
  if (col %in% names(push_manual)) push_manual[[col]][is.na(push_manual[[col]])] <- ""
}

message(sprintf("\n  SeaTable push: %d rows, %d columns", nrow(push_df), ncol(push_df) - 1))
message(sprintf("  Held back (HAS_MANUAL_ANNOTATION): %d rows", nrow(push_manual)))

if (nrow(push_manual) > 0) {
  manual_file <- "data/codex/class_fixes_manual_review.csv"
  readr::write_csv(push_manual, manual_file)
  message(sprintf("  Saved to %s", manual_file))
}

# banctable_update_rows(base = 'banc_meta',
#                       table = "banc_meta",
#                       df = push_df,
#                       append_allowed = FALSE,
#                       chunksize = 200)
# message("  SeaTable update complete")

###############################################
### Neuromere update                        ###
### For ALL BANC neurons with a manc_match, ###
### pull neuromere from franken_meta.       ###
### Independent of bates cell_type_source.  ###
###############################################

message("\n### Neuromere update ###")

# Re-load unfiltered franken_meta — the earlier `fm` is filtered to fafb_id != "",
# which excludes VNC-only neurons that have manc_id but no fafb_id.
fm_full <- franken_meta()

if ("neuromere" %in% names(fm_full)) {
  # Build neuromere lookup keyed by manc_id (rename to avoid join collision)
  fm_neuromere <- fm_full %>%
    dplyr::filter(!is.na(manc_id), manc_id != "",
                  !is.na(neuromere), neuromere != "") %>%
    dplyr::mutate(manc_id = as.character(manc_id)) %>%
    dplyr::distinct(manc_id, .keep_all = TRUE) %>%
    dplyr::select(manc_id, neuromere_new = neuromere)
  message(sprintf("  franken_meta manc->neuromere: %d entries", nrow(fm_neuromere)))

  # All BANC neurons with a manc_match (re-query to ensure we have status + current neuromere)
  bc_all <- banctable_query(
    "SELECT _id, root_888, root_id, manc_match, neuromere, status FROM banc_meta"
  ) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), as.character)) %>%
    dplyr::filter(!is.na(manc_match), manc_match != "") %>%
    dplyr::rename(neuromere_old = neuromere)
  message(sprintf("  BANC neurons with manc_match: %d", nrow(bc_all)))

  # Some manc_match values may be comma-separated — use first ID
  bc_all$manc_match_key <- sapply(strsplit(bc_all$manc_match, ","), function(x) trimws(x[1]))

  # Join to franken_meta neuromere
  neuromere_join <- bc_all %>%
    dplyr::left_join(fm_neuromere, by = c("manc_match_key" = "manc_id")) %>%
    dplyr::filter(!is.na(neuromere_new), neuromere_new != "")
  message(sprintf("  Matched to franken_meta: %d", nrow(neuromere_join)))

  # Only push where neuromere is different (or currently blank)
  changed <- neuromere_join %>%
    dplyr::filter(is.na(neuromere_old) | neuromere_old == "" |
                  neuromere_old != neuromere_new)
  message(sprintf("  Neuromere changes: %d", nrow(changed)))

  # Separate HAS_MANUAL_ANNOTATION neurons for careful review
  has_manual_nm <- grepl("HAS_MANUAL_NEUROMERE", changed$status, ignore.case = TRUE)

  push_neuromere <- changed[!has_manual_nm, ] %>%
    dplyr::transmute(`_id`, root_id, neuromere = neuromere_new) %>%
    as.data.frame()
  push_neuromere$neuromere[is.na(push_neuromere$neuromere)] <- ""

  push_neuromere_manual <- changed[has_manual_nm, ] %>%
    dplyr::select(`_id`, root_888, root_id, manc_match,
                  neuromere_old, neuromere_new, status) %>%
    as.data.frame()

  message(sprintf("  SeaTable push (neuromere): %d rows", nrow(push_neuromere)))
  message(sprintf("  Held back (HAS_MANUAL_NEUROMERE): %d rows", nrow(push_neuromere_manual)))

  if (nrow(push_neuromere_manual) > 0) {
    nm_manual_file <- "data/codex/neuromere_manual_review.csv"
    readr::write_csv(push_neuromere_manual, nm_manual_file)
    message(sprintf("  Saved to %s", nm_manual_file))
  }

  # banctable_update_rows(base = 'banc_meta',
  #                       table = "banc_meta",
  #                       df = push_neuromere,
  #                       append_allowed = FALSE,
  #                       chunksize = 200)
  # message("  Neuromere update complete")
}

###############################################
### Sensory follower-column update          ###
### For all BANC neurons whose super_class  ###
### contains "sensory", pull follower       ###
### columns from the matched neuron via     ###
### ID-level joins:                         ###
###   region == ventral_nerve_cord OR       ###
###   super_class contains "ascending"      ###
###     -> prefer manc_match -> manc_id     ###
###   else                                  ###
###     -> prefer fafb_match -> fafb_id     ###
###   fallback (if both missing):           ###
###     malecns_match -> malecns_09_id      ###
### Neurons with no ID-level match are      ###
### offered separately as                   ###
### push_df_sensory_uncertain (filled with  ###
### franken_meta cell_type consensus).      ###
### HAS_MANUAL_ANNOTATION rows are held     ###
### back from both push frames.             ###
###############################################

message("\n### Sensory follower-column update ###")

# Re-query BANC for all sensory neurons. Include match-key columns and region
# in addition to target_cols.
sens_select_cols <- unique(c("_id", "root_888", "root_id", "cell_type",
                             "super_class", "region",
                             "fafb_match", "manc_match", "malecns_match",
                             "status", target_cols))
sens_select <- paste(sens_select_cols, collapse = ", ")
sens_bc <- banctable_query(sprintf("SELECT %s FROM banc_meta", sens_select)) %>%
  dplyr::mutate(dplyr::across(dplyr::everything(), as.character)) %>%
  dplyr::filter(!is.na(super_class),
                grepl("sensory", super_class, ignore.case = TRUE),
                !is.na(cell_type), cell_type != "")
message(sprintf("  BANC sensory neurons with cell_type: %d", nrow(sens_bc)))

# --- Determine source preference per row ---
prefer_manc <- (!is.na(sens_bc$region) & sens_bc$region == "ventral_nerve_cord") |
               (!is.na(sens_bc$super_class) &
                grepl("ascending", sens_bc$super_class, ignore.case = TRUE))
message(sprintf("  Prefer manc (VNC or ascending): %d", sum(prefer_manc)))
message(sprintf("  Prefer fafb (other):            %d", sum(!prefer_manc)))

# --- Vectorised lookups into the existing reference tables ---
fafb_idx <- match(sens_bc$fafb_match, fm_lookup$fafb_id)
manc_idx <- if (!is.null(manc_lookup))
  match(sens_bc$manc_match, manc_lookup$manc_id) else
  rep(NA_integer_, nrow(sens_bc))
mcns_idx <- if (!is.null(mcns_lookup))
  match(sens_bc$malecns_match, mcns_lookup$malecns_09_id) else
  rep(NA_integer_, nrow(sens_bc))

# Source code per row:  1=fafb, 2=manc, 3=malecns, 4=none
src_fafb_first <- ifelse(!is.na(fafb_idx), 1L,
                  ifelse(!is.na(manc_idx), 2L,
                  ifelse(!is.na(mcns_idx), 3L, 4L)))
src_manc_first <- ifelse(!is.na(manc_idx), 2L,
                  ifelse(!is.na(fafb_idx), 1L,
                  ifelse(!is.na(mcns_idx), 3L, 4L)))
src <- ifelse(prefer_manc, src_manc_first, src_fafb_first)

# --- Cell-type consistency check ---
# A row only counts as "certain" if the chosen source's own cell_type also
# agrees with BANC's cell_type. Otherwise the ID-level match is suspect
# (stale match, wrong join, or franken_meta carrying a different neuron's
# row for that ID) and we demote the row to the uncertain bucket so it
# gets the cell_type consensus fallback instead of a blind copy.
matched_ct <- rep(NA_character_, nrow(sens_bc))
m1 <- src == 1L
if (any(m1)) matched_ct[m1] <- fm_lookup$fm_cell_type[fafb_idx[m1]]
m2 <- src == 2L
if (any(m2) && !is.null(manc_lookup) && "manc_cell_type" %in% names(manc_lookup)) {
  matched_ct[m2] <- manc_lookup$manc_cell_type[manc_idx[m2]]
}
m3 <- src == 3L
if (any(m3) && !is.null(mcns_lookup) && "mcns_cell_type" %in% names(mcns_lookup)) {
  matched_ct[m3] <- mcns_lookup$mcns_cell_type[mcns_idx[m3]]
}

bc_ct <- sens_bc$cell_type
ct_mismatch <- src %in% c(1L, 2L, 3L) &
               (is.na(matched_ct) | matched_ct == "" |
                is.na(bc_ct) | bc_ct == "" |
                matched_ct != bc_ct)
n_demoted <- sum(ct_mismatch)
if (n_demoted > 0) {
  message(sprintf("  Demoted to uncertain (source cell_type != BANC cell_type): %d", n_demoted))
}
src[ct_mismatch] <- 4L
src_label <- c("franken_meta_fafb", "franken_meta_manc", "malecns", "none")[src]

message(sprintf("  Source breakdown (after cell_type consistency check):"))
message(sprintf("    franken_meta_fafb: %d", sum(src == 1L)))
message(sprintf("    franken_meta_manc: %d", sum(src == 2L)))
message(sprintf("    malecns:           %d", sum(src == 3L)))
message(sprintf("    no ID-level match: %d", sum(src == 4L)))

# --- Build per-column "new value" vectors from the chosen source ---
mcns_cols_avail <- if (!is.null(mcns_lookup)) {
  intersect(target_cols, names(mcns_lookup))
} else character(0)
manc_cols_avail <- if (!is.null(manc_lookup)) {
  intersect(target_cols, names(manc_lookup))
} else character(0)
all_new_cols <- unique(c(fm_cols, manc_cols_avail, mcns_cols_avail))

new_vals <- vector("list", length(all_new_cols))
names(new_vals) <- all_new_cols
for (col in all_new_cols) {
  v_fafb <- if (col %in% names(fm_lookup))
    fm_lookup[[col]][fafb_idx] else rep(NA_character_, nrow(sens_bc))
  v_manc <- if (!is.null(manc_lookup) && col %in% names(manc_lookup))
    manc_lookup[[col]][manc_idx] else rep(NA_character_, nrow(sens_bc))
  v_mcns <- if (!is.null(mcns_lookup) && col %in% names(mcns_lookup))
    mcns_lookup[[col]][mcns_idx] else rep(NA_character_, nrow(sens_bc))

  v <- rep(NA_character_, nrow(sens_bc))
  v[src == 1L] <- v_fafb[src == 1L]
  v[src == 2L] <- v_manc[src == 2L]
  v[src == 3L] <- v_mcns[src == 3L]
  new_vals[[col]] <- v
}

# Guard: never demote a sensory_* super_class. MANC (and sometimes other
# sources) drop the sensory_ prefix in their taxonomy — e.g. a BANC
# sensory_ascending neuron whose manc_match resolves to a franken_meta row
# with super_class = "ascending". Without this guard, the follower-col copy
# would silently rewrite sensory_ascending -> ascending.
if ("super_class" %in% all_new_cols) {
  old_sc <- sens_bc$super_class
  new_sc <- new_vals$super_class
  is_sens <- !is.na(old_sc) & grepl("sensory", old_sc, ignore.case = TRUE)
  bad_new <- !is.na(new_sc) & !grepl("sensory", new_sc, ignore.case = TRUE)
  block <- is_sens & bad_new
  if (any(block)) {
    message(sprintf("  Guard: blocked %d sensory->non-sensory super_class proposals",
                    sum(block)))
    new_vals$super_class[block] <- NA
  }
}

# --- Cell-type consensus fallback (for the uncertain bucket only) ---
mode_value <- function(x) {
  v <- x[!is.na(x) & x != ""]
  if (length(v) == 0) return(NA_character_)
  names(sort(table(v), decreasing = TRUE))[1]
}

fm_sensory <- fm %>%
  dplyr::filter(!is.na(super_class),
                grepl("sensory", super_class, ignore.case = TRUE),
                !is.na(cell_type), cell_type != "")
fm_sensory_consensus <- fm_sensory %>%
  dplyr::group_by(cell_type) %>%
  dplyr::summarise(dplyr::across(dplyr::all_of(fm_cols), mode_value),
                   n_fm = dplyr::n(), .groups = "drop")
message(sprintf("  franken_meta sensory rows: %d (in %d cell_types)",
                nrow(fm_sensory), nrow(fm_sensory_consensus)))

cons_idx <- match(sens_bc$cell_type, fm_sensory_consensus$cell_type)
uncertain <- src == 4L
n_uncertain_with_cons <- 0L
for (col in fm_cols) {
  v_cons <- fm_sensory_consensus[[col]][cons_idx]
  fill <- uncertain & !is.na(v_cons) & v_cons != ""
  new_vals[[col]][fill] <- v_cons[fill]
}
# How many uncertain rows got at least one consensus value?
if (sum(uncertain) > 0) {
  cons_any <- rep(FALSE, nrow(sens_bc))
  for (col in fm_cols) {
    cons_any <- cons_any | (uncertain & !is.na(new_vals[[col]]) & new_vals[[col]] != "")
  }
  n_uncertain_with_cons <- sum(cons_any)
  message(sprintf("  Uncertain rows with cell_type consensus available: %d / %d",
                  n_uncertain_with_cons, sum(uncertain)))
}

# --- Detect column-level changes (vectorised) ---
per_col_changes <- integer(length(all_new_cols))
names(per_col_changes) <- all_new_cols
needs_update <- logical(nrow(sens_bc))
for (col in all_new_cols) {
  if (!col %in% names(sens_bc)) next
  old_v <- sens_bc[[col]]
  new_v <- new_vals[[col]]
  diff <- !is.na(new_v) & new_v != "" &
          (is.na(old_v) | old_v == "" | old_v != new_v)
  per_col_changes[col] <- sum(diff, na.rm = TRUE)
  needs_update <- needs_update | diff
}

message(sprintf("  Sensory neurons with >=1 follower-col update: %d",
                sum(needs_update)))
for (col in all_new_cols) {
  if (per_col_changes[[col]] > 0) {
    message(sprintf("    %s: %d", col, per_col_changes[[col]]))
  }
}

# --- Assemble per-row metadata + new values ---
sens_aug <- sens_bc
sens_aug$source <- src_label
for (col in all_new_cols) {
  sens_aug[[paste0(col, "_new")]] <- new_vals[[col]]
}
sens_aug$needs_update <- needs_update

certain_mask   <- needs_update & src != 4L
uncertain_mask <- needs_update & src == 4L
has_manual <- grepl("HAS_MANUAL_ANNOTATION", sens_aug$status, ignore.case = TRUE)

message(sprintf("  Certain (ID-level match):     %d (manual held back: %d)",
                sum(certain_mask), sum(certain_mask & has_manual)))
message(sprintf("  Uncertain (no ID-level match): %d (manual held back: %d)",
                sum(uncertain_mask), sum(uncertain_mask & has_manual)))

# --- Helper: build a push frame from a row mask ---
build_push_df <- function(mask) {
  if (sum(mask) == 0) return(NULL)
  out <- data.frame(`_id` = sens_aug$`_id`[mask],
                    stringsAsFactors = FALSE, check.names = FALSE)
  for (col in all_new_cols) {
    val <- new_vals[[col]][mask]
    val[is.na(val)] <- ""
    out[[col]] <- val
  }
  out
}

# Certain push frame (eligible)
push_sens_df <- build_push_df(certain_mask & !has_manual)
if (!is.null(push_sens_df)) {
  push_file <- "data/codex/class_fixes_sensory.csv"
  readr::write_csv(push_sens_df, push_file)
  message(sprintf("  Saved push_sens_df to %s (%d rows, %d cols)",
                  push_file, nrow(push_sens_df), ncol(push_sens_df) - 1))
}

# Uncertain push frame (eligible) — cell_type consensus only
push_df_sensory_uncertain <- build_push_df(uncertain_mask & !has_manual)
if (!is.null(push_df_sensory_uncertain)) {
  uncertain_file <- "data/codex/class_fixes_sensory_uncertain.csv"
  readr::write_csv(push_df_sensory_uncertain, uncertain_file)
  message(sprintf("  Saved push_df_sensory_uncertain to %s (%d rows, %d cols)",
                  uncertain_file, nrow(push_df_sensory_uncertain),
                  ncol(push_df_sensory_uncertain) - 1))
}

# Manual-annotation review CSVs (one for certain, one for uncertain)
if (sum(certain_mask & has_manual) > 0) {
  certain_manual_file <- "data/codex/class_fixes_sensory_manual_review.csv"
  readr::write_csv(sens_aug[certain_mask & has_manual, ], certain_manual_file)
  message(sprintf("  Saved certain-manual-review rows to %s", certain_manual_file))
}
if (sum(uncertain_mask & has_manual) > 0) {
  uncertain_manual_file <- "data/codex/class_fixes_sensory_uncertain_manual_review.csv"
  readr::write_csv(sens_aug[uncertain_mask & has_manual, ], uncertain_manual_file)
  message(sprintf("  Saved uncertain-manual-review rows to %s", uncertain_manual_file))
}

# banctable_update_rows(base = 'banc_meta',
#                       table = "banc_meta",
#                       df = push_sens_df,
#                       append_allowed = FALSE,
#                       chunksize = 200)
# message("  Sensory follower-column update complete (certain)")

# # push_df_sensory_uncertain — review before pushing
# banctable_update_rows(base = 'banc_meta',
#                       table = "banc_meta",
#                       df = push_df_sensory_uncertain,
#                       append_allowed = FALSE,
#                       chunksize = 200)

