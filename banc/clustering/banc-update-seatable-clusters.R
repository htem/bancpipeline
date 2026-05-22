#' banc-update-seatable-clusters — Push spectral cns_cluster + cns_network to SeaTable.
#'
#' Reads the canonical v888 v2 spectral CSV, joins on `root_888` (NOT
#' `root_id` — root_id is live/dynamic, root_888 pins to the version the
#' spectral pipeline wrote), and overwrites both `cns_cluster` and
#' `cns_network` columns in `banc_meta`. Wipe + rewrite semantics: every
#' SeaTable row gets a new value (old labels never linger). Sensory / glia /
#' trachea / optic_lobe_intrinsic neurons are explicitly cleared.
#'
#' @section Reads:
#'   - `data/cns_network/spectral_clustering_..._v2.csv`
#'   - SeaTable `banc_meta`
#'
#' @section Writes:
#'   - SeaTable `banc_meta`: cols `cns_cluster` (CNS_NN form), `cns_network`
#'
#' @section CLI:
#'   --dry-run   don't push to SeaTable
#'
#' @section Paper:
#'   Methods §"Spectral clustering" — labels the columns the paper uses.

###########################################################
### Update cns_cluster and cns_network in BANC seatable
###
### Reads the canonical v888 v2 spectral CSV (cluster_count=13),
### joins on SeaTable's root_888 column (NOT root_id — root_id
### is live/dynamic; root_888 is version-pinned to match what
### the spectral pipeline wrote), and overwrites both columns
### in banc_meta.
###
### cns_network = the descriptive label in the spectral CSV's
### cns_network column (no hard-coded CNS_XX -> name mapping).
### cns_cluster = "CNS_NN" derived from spectral_cluster index.
###
### Neurons whose super_class matches
###   sensory | glia | trachea | optic_lobe_intrinsic
### are explicitly cleared (no label even if the spectral CSV
### picked them up).
###
### Wipe + rewrite semantics: every SeaTable row gets a new
### value. Rows not present in the spectral CSV (or matching
### the exclusion rule) get empty strings, so old labels never
### linger.
###
### Usage:
###   Rscript banc/clustering/banc-update-seatable-clusters.R --dry-run
###   Rscript banc/clustering/banc-update-seatable-clusters.R
###########################################################
source("banc/banc-startup.R")

local({
  message("### banc: update seatable cns_cluster + cns_network ###")
  t_start <- Sys.time()

  # ---------------------------------------------------------------
  # Configuration
  # ---------------------------------------------------------------
  dry_run <- "--dry-run" %in% commandArgs(trailingOnly = TRUE)
  if (dry_run) message("  *** DRY RUN — no seatable writes ***")

  # Retry-with-exponential-backoff wrapper. SeaTable cloud API has a 30s
  # per-chunk read timeout that fails sporadically under load; idempotent
  # update_rows tolerates re-execution so we just retry.
  .seatable_retry <- function(expr, label, max_tries = 6L, base_wait = 10) {
    .e <- substitute(expr)
    for (.try in seq_len(max_tries)) {
      ok <- tryCatch({ eval(.e, envir = parent.frame()); TRUE },
        error = function(e) {
          message(sprintf("  [%s] attempt %d/%d failed: %s",
                          label, .try, max_tries, conditionMessage(e)))
          FALSE
        })
      if (ok) {
        if (.try > 1L) message(sprintf("  [%s] succeeded on retry %d", label, .try))
        return(invisible(TRUE))
      }
      if (.try < max_tries) {
        wait <- min(300, base_wait * 2^(.try - 1))
        message(sprintf("  [%s] sleeping %ds before retry", label, wait))
        Sys.sleep(wait)
      }
    }
    stop(sprintf("[%s] exhausted %d retries", label, max_tries), call. = FALSE)
  }

  # Canonical spectral CSV location. Mirrors the path BANC-Project
  # consumes (R/startup/banc-startup.R .banc_spectral_csv) — both
  # repos point at this same file in the BANC-project tree.
  bancproject_root <- "/n/data1/hms/neurobio/wilson/banc/BANC-project"
  clustering_file <- file.path(
    bancproject_root, "data/cns_network",
    sprintf("spectral_clustering_min_connection_strength_1_banc_version_%s_cluster_count_13_cluster_seed_10_embedding_seed_3_v2.csv",
            banc.version)
  )
  if (!file.exists(clustering_file)) {
    # Fallback: bancpipeline working copy.
    clustering_file <- file.path(
      "data/cns_network",
      sprintf("spectral_clustering_min_connection_strength_1_banc_version_%s_cluster_count_13_cluster_seed_10_embedding_seed_3_v2.csv",
              banc.version)
    )
  }
  stopifnot(file.exists(clustering_file))
  message(sprintf("  Spectral CSV: %s", clustering_file))

  # Exclusion rule (super_class regex). Neurons matching this never
  # get a label, even if the spectral CSV included them.
  excl_re <- "sensory|glia|trachea|optic_lobe_intrinsic"

  # ---------------------------------------------------------------
  # 1. Load clustering CSV
  # ---------------------------------------------------------------
  message("Loading clustering CSV...")
  clustering <- readr::read_csv(clustering_file,
                                col_types = readr::cols(.default = "c")) %>%
    dplyr::mutate(
      spectral_cluster = as.integer(spectral_cluster),
      new_cns_cluster  = paste0("CNS_", stringr::str_pad(spectral_cluster, width = 2, pad = "0")),
      new_cns_network  = cns_network
    ) %>%
    dplyr::select(spectral_root_888 = root_id, new_cns_cluster, new_cns_network) %>%
    dplyr::distinct(spectral_root_888, .keep_all = TRUE)
  message(sprintf("  %s neurons in clustering CSV",
                  format(nrow(clustering), big.mark = ",")))

  # ---------------------------------------------------------------
  # 2. Query current SeaTable state (raw — drives a write, no cache)
  # ---------------------------------------------------------------
  message("Querying SeaTable banc_meta ...")
  bc.current <- banctable_query(
    "SELECT _id, root_id, root_888, supervoxel_id, super_class, region, cns_cluster, cns_network FROM banc_meta"
  )
  bc.current <- bc.current %>%
    dplyr::mutate(root_888 = as.character(root_888),
                  root_id  = as.character(root_id))
  message(sprintf("  %s SeaTable rows (non-empty root_888: %s)",
                  format(nrow(bc.current), big.mark = ","),
                  format(sum(!is.na(bc.current$root_888) & nzchar(bc.current$root_888)), big.mark = ",")))

  # ---------------------------------------------------------------
  # 3. Compose new values per row
  # ---------------------------------------------------------------
  message("Building update...")
  bc.update <- bc.current %>%
    dplyr::left_join(clustering, by = c("root_888" = "spectral_root_888")) %>%
    dplyr::mutate(
      old_cns_cluster = ifelse(is.na(cns_cluster) | cns_cluster == "", NA_character_, cns_cluster),
      old_cns_network = ifelse(is.na(cns_network) | cns_network == "", NA_character_, cns_network),
      excluded        = grepl(excl_re, super_class, ignore.case = TRUE),
      cns_cluster_new = dplyr::case_when(
        excluded ~ "",
        !is.na(new_cns_cluster) ~ new_cns_cluster,
        TRUE ~ ""
      ),
      cns_network_new = dplyr::case_when(
        excluded ~ "",
        !is.na(new_cns_network) ~ new_cns_network,
        TRUE ~ ""
      )
    )

  # ---------------------------------------------------------------
  # 3b. Diagnostics
  # ---------------------------------------------------------------
  identical_na_blank <- function(a, b) {
    a2 <- ifelse(is.na(a) | a == "", NA_character_, a)
    b2 <- ifelse(is.na(b) | b == "", NA_character_, b)
    (is.na(a2) & is.na(b2)) | (!is.na(a2) & !is.na(b2) & a2 == b2)
  }

  n_total       <- nrow(bc.update)
  n_excluded    <- sum(bc.update$excluded, na.rm = TRUE)
  n_label_now   <- sum(nzchar(bc.update$cns_network_new), na.rm = TRUE)
  n_label_was   <- sum(!is.na(bc.update$old_cns_network))
  n_clear       <- sum(bc.update$cns_network_new == "", na.rm = TRUE)
  n_in_spec     <- sum(!is.na(bc.update$new_cns_network))
  n_in_spec_excl <- sum(bc.update$excluded & !is.na(bc.update$new_cns_network), na.rm = TRUE)
  n_changed     <- sum(!identical_na_blank(bc.update$old_cns_network, bc.update$cns_network_new))

  message("\n--- Counts ---")
  message(sprintf("  SeaTable rows:                              %s", format(n_total, big.mark = ",")))
  message(sprintf("  Currently labelled (cns_network non-empty): %s", format(n_label_was, big.mark = ",")))
  message(sprintf("  Spectral CSV matched into SeaTable:         %s", format(n_in_spec, big.mark = ",")))
  message(sprintf("  Exclusion (%s) hits:        %s", excl_re, format(n_excluded, big.mark = ",")))
  message(sprintf("    of which were in spectral (dropped):      %s", format(n_in_spec_excl, big.mark = ",")))
  message(sprintf("  WILL be labelled after update:              %s", format(n_label_now, big.mark = ",")))
  message(sprintf("  WILL be cleared:                            %s", format(n_clear, big.mark = ",")))
  message(sprintf("  Rows whose cns_network value changes:       %s", format(n_changed, big.mark = ",")))

  message("\n--- Exclusion breakdown by super_class ---")
  ex_tab <- sort(table(bc.update$super_class[bc.update$excluded], useNA = "no"),
                 decreasing = TRUE)
  print(as.data.frame(ex_tab))

  message("\n--- New cns_network distribution ---")
  new_tab <- sort(table(bc.update$cns_network_new[nzchar(bc.update$cns_network_new)]),
                  decreasing = TRUE)
  print(as.data.frame(new_tab))

  message("\n--- Transition matrix (top 25 rows) ---")
  trans <- bc.update %>%
    dplyr::mutate(
      old = ifelse(is.na(old_cns_network), "(none)", old_cns_network),
      new = ifelse(nzchar(cns_network_new), cns_network_new, "(none)")
    ) %>%
    dplyr::count(old, new) %>%
    dplyr::arrange(old, dplyr::desc(n))
  print(as.data.frame(head(trans, 25)), row.names = FALSE)

  # ---------------------------------------------------------------
  # 4. Finalise update dataframe (keep only seatable columns)
  # ---------------------------------------------------------------
  bc.push <- bc.update %>%
    dplyr::transmute(
      `_id`        = `_id`,
      root_id      = root_id,
      cns_cluster  = cns_cluster_new,
      cns_network  = cns_network_new
    )
  bc.push <- as.data.frame(bc.push)
  bc.push[is.na(bc.push)] <- ""

  preview_file <- file.path("data/cns_network",
                            sprintf("seatable_update_preview_%s.csv", banc.version))
  if (!dir.exists(dirname(preview_file))) dir.create(dirname(preview_file), recursive = TRUE)
  readr::write_csv(bc.push, preview_file)
  message(sprintf("\nPreview CSV: %s (%s rows)",
                  preview_file, format(nrow(bc.push), big.mark = ",")))

  # ---------------------------------------------------------------
  # 5. Push to SeaTable
  # ---------------------------------------------------------------
  if (dry_run) {
    message("\n  DRY RUN complete — skipping SeaTable update")
  } else {
    # Restrict to rows that actually change to cut wall-time + collision
    # risk; SeaTable's update_rows is idempotent so this is safe either way.
    bc.push.changed <- bc.update %>%
      dplyr::filter(!identical_na_blank(old_cns_network, cns_network_new) |
                    !identical_na_blank(old_cns_cluster, cns_cluster_new)) %>%
      dplyr::transmute(`_id` = `_id`, root_id = root_id,
                       cns_cluster = cns_cluster_new, cns_network = cns_network_new) %>%
      as.data.frame()
    bc.push.changed[is.na(bc.push.changed)] <- ""
    message(sprintf("\nUpdating SeaTable (%s changed rows, chunksize 500) ...",
                    format(nrow(bc.push.changed), big.mark = ",")))
    .seatable_retry(
      banctable_update_rows(
        base = "banc_meta",
        table = "banc_meta",
        df = bc.push.changed,
        append_allowed = FALSE,
        chunksize = 500
      ),
      label = "banctable_update_rows(cns_network/cns_cluster)"
    )
    message("  SeaTable update complete")
  }

  message(sprintf("\n### banc: seatable cluster update complete [%s] ###",
                  format(round(difftime(Sys.time(), t_start, units = "mins"), 2))))
})
