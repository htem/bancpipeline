#' banc-kc-by-connectivity — Name Kenyon cells by connectivity to known DANs / MBONs.
#'
#' Picks the dominant DAN/MBON subtype group (αβ / α'β' / γ) per KC, assigns
#' ASCII-safe names (`KCab` / `KCa'b'` / `KCg` / `KC`). Existing more-specific
#' cell_types are preserved.
#'
#' @section Reads:
#'   - SeaTable `banc_meta`, CAVE KC-MBON / KC-DAN connectivity
#'
#' @section Writes:
#'   - SeaTable `banc_meta`: col `cell_type`

###############################################################################
### BANC annotation: Kenyon-cell subtype by connectivity to DANs/MBONs
###
### Counts KC ↔ DAN/MBON synapses, maps each DAN/MBON to a KC subtype group (αβ, α'β', γ), assigns ASCII-safe names (KCab / KCa'b' / KCg / KC), and pushes to SeaTable. Keeps any existing cell_type that is more specific than the new one.
###
### Originally part of banc/utilities/banc-tracing.R; split out 2026-05-21
### into per-task scripts under banc/annotations/ because the original was
### never meant to be run end-to-end — each section is a one-off curation
### task. Sections may share state; if so, run them in order.
###############################################################################
source("banc/banc-startup.R")

##############################################################################
### Name Kenyon cells by their connectivity to known DANs and MBONs       ###
###                                                                        ###
### Logic:                                                                 ###
###   1. For each KC, count synapses to/from named DANs and MBONs         ###
###   2. Map each DAN/MBON to a KC subtype group (αβ, α'β', γ)            ###
###   3. Aggregate per KC × group, pick the dominant group (with tol)      ###
###   4. Assign ASCII-safe names: KCab, KCa'b', KCg, or plain KC          ###
###   5. Keep existing cell_type if it is more specific than the new one   ###
###   6. Update seatable                                                   ###
##############################################################################

library(dplyr)
library(bancr)

####################
### Load lookups ###
####################

dan_map  <- read.csv("data/annotation/mb/dan_kc_subtype_preference.csv")
mbon_map <- read.csv("data/annotation/mb/mbon_kc_subtype_preference.csv")

# Unified lookup: DAN/MBON type → kc_group + specificity weight
dm_lookup <- bind_rows(
  dan_map  %>% select(type, kc_group, kc_fraction),
  mbon_map %>% select(type, kc_group, kc_fraction)
) %>% distinct()

# Greek → ASCII name mapping
kc_group_to_name <- c(
  "αβ-KCs"    = "KCab",
  "α'β'-KCs"  = "KCa'b'",
  "γ-KCs"     = "KCg"
)

####################
### Get neurons  ###
####################

bc <- banctable_query()

kcs    <- subset(bc, grepl("^KC|Kenyon", cell_type, ignore.case = TRUE))
kc_ids <- na.omit(unique(kcs$root_id))

dans_mbons <- subset(bc, cell_type %in% dm_lookup$type)
dm_ids     <- na.omit(unique(dans_mbons$root_id))

##################################
### Get KC ↔ DAN/MBON edges   ###
##################################

# KC outputs → MBONs
kc_out <- pbapply::pblapply(kc_ids, function(x)
  tryCatch(
    banc_partners(x, partners = "output", synapse_table = "synapses_v2"),
    error = function(e) NULL
  ))
kc_out <- do.call(plyr::rbind.fill, kc_out[!sapply(kc_out, is.null)])

# DAN inputs → KCs
kc_in <- pbapply::pblapply(kc_ids, function(x)
  tryCatch(
    banc_partners(x, partners = "input", synapse_table = "synapses_v2"),
    error = function(e) NULL
  ))
kc_in <- do.call(plyr::rbind.fill, kc_in[!sapply(kc_in, is.null)])

##################################
### Unify edges                ###
##################################

out_edges <- kc_out %>%
  mutate(
    pre_pt_root_id  = as.character(pre_pt_root_id),
    post_pt_root_id = as.character(post_pt_root_id)
  ) %>%
  filter(post_pt_root_id %in% dm_ids) %>%
  mutate(kc_id = pre_pt_root_id, partner_id = post_pt_root_id)

in_edges <- kc_in %>%
  mutate(
    pre_pt_root_id  = as.character(pre_pt_root_id),
    post_pt_root_id = as.character(post_pt_root_id)
  ) %>%
  filter(pre_pt_root_id %in% dm_ids) %>%
  mutate(kc_id = post_pt_root_id, partner_id = pre_pt_root_id)

all_edges <- bind_rows(out_edges, in_edges)

##################################
### Join partner type → kc_group ###
##################################

all_edges <- all_edges %>%
  left_join(
    bc %>%
      distinct(root_id, .keep_all = TRUE) %>%
      select(root_id, partner_cell_type = cell_type),
    by = c("partner_id" = "root_id")
  ) %>%
  left_join(dm_lookup, by = c("partner_cell_type" = "type"))

##################################
### Aggregate per KC × kc_group ###
##################################

kc_targets <- all_edges %>%
  dplyr::filter(!is.na(kc_group)) %>%
  dplyr::group_by(kc_id, kc_group) %>%
  dplyr::summarise(total_count = n(), .groups = "drop")

########################################
### Top-band with tolerance          ###
########################################

tol <- 0.3

kc_best <- kc_targets %>%
  dplyr::group_by(kc_id) %>%
  dplyr::mutate(
    max_count = max(total_count, na.rm = TRUE),
    thresh    = (1 - tol) * max_count,
    is_top    = total_count >= thresh
  ) %>%
  dplyr::filter(is_top) %>%
  dplyr::mutate(n_top = n()) %>%
  dplyr::ungroup()

#########################################
### Build connectivity-based names    ###
#########################################

kc_names_from_conn <- kc_best %>%
  dplyr::mutate(
    conn_cell_type = if_else(
      n_top == 1,
      unname(kc_group_to_name[kc_group]),
      "KC"
    )
  ) %>%
  dplyr::distinct(kc_id, conn_cell_type, n_top)

# KCs with no/insufficient signal → plain "KC"
kc_no_signal <- setdiff(kc_ids, kc_names_from_conn$kc_id)

kc_names_from_conn <- bind_rows(
  kc_names_from_conn,
  data.frame(
    kc_id          = kc_no_signal,
    conn_cell_type = "KC",
    n_top          = NA_integer_
  )
)

#########################################
### Compare with existing cell_type   ###
### Keep old if it is more specific   ###
#########################################

kc_compare <- kc_names_from_conn %>%
  dplyr::left_join(
    bc %>%
      distinct(root_id, .keep_all = TRUE) %>%
      select(kc_id = root_id, existing_type = cell_type),
    by = "kc_id"
  ) %>%
  dplyr::mutate(
    # The existing type is "more specific" if it starts with the new name
    # but has additional detail, e.g. old = "KCab-s", new = "KCab"
    old_is_more_specific = !is.na(existing_type) &
      existing_type != conn_cell_type &
      startsWith(existing_type, conn_cell_type) &
      nchar(existing_type) > nchar(conn_cell_type),
    # Final decision: keep old if more specific, otherwise use new
    final_cell_type = if_else(old_is_more_specific, existing_type, conn_cell_type),
    # Only flag rows that actually need updating
    needs_update = is.na(existing_type) | existing_type != final_cell_type
  ) %>%
  select(kc_id, existing_type, conn_cell_type, final_cell_type,
         old_is_more_specific, needs_update, n_top)

# Summary
message(sprintf(
  "KCs total: %d | Assigned subtype: %d | Kept more specific old: %d | Need update: %d",
  nrow(kc_compare),
  sum(kc_compare$conn_cell_type != "KC"),
  sum(kc_compare$old_is_more_specific),
  sum(kc_compare$needs_update)
))

#########################################
### Update seatable                   ###
#########################################

kc_to_update <- kc_compare %>%
  filter(needs_update)

bc.restored <- bc %>%
  distinct(`_id`, root_id) %>%
  left_join(
    kc_to_update %>%
      select(root_id = kc_id, cell_type = final_cell_type),
    by = "root_id"
  ) %>%
  filter(!is.na(cell_type))

bc.restored <- as.data.frame(bc.restored)
bc.restored[is.na(bc.restored)] <- ''

# banctable_update_rows(
#   base = 'banc_meta',
#   table = "banc_meta",
#   df = bc.restored,
#   append_allowed = FALSE,
#   chunksize = 1000
# )

# ## --- Setup ---------------------------------------------------------------
# pattern <- stringr::regex("gustatory|taste|labellum_bristle_neuron", ignore_case = TRUE)
# 
# keep_cols_banc <- c("root_id","cell_type","side","super_class","cell_class","cell_sub_class")
# keep_cols_fr   <- c("neuron_id","cell_type","side","super_class","cell_class","cell_sub_class")
# 
# # destination folder
# out_dir <- "~/Downloads/"
# if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
# 
# ## --- Subset & deduplicate ------------------------------------------------
# # If your franken.meta uses a different ID column, rename it to neuron_id first, e.g.:
# # franken.meta <- dplyr::rename(franken.meta, neuron_id = id)
# franken.meta.g <- franken.meta |>
#   dplyr::filter(
#     grepl("sensory",super_class),
#     !grepl("^AN|^IN|^SA_",cell_type),
#     dplyr::if_any(dplyr::all_of(c("cell_sub_class","cell_function")),
#                   ~ stringr::str_detect(stringr::str_to_lower(dplyr::coalesce(.x, "")), pattern))
#   ) |>
#   dplyr::distinct(neuron_id, .keep_all = TRUE) |>
#   dplyr::select(dplyr::all_of(keep_cols_fr))
# 
# 
# banc.meta.g <- banc.meta |>
#   dplyr::filter(
#     grepl("sensory",super_class),
#     !grepl("^AN|^IN",cell_type),
#     dplyr::if_any(dplyr::all_of(c("cell_sub_class","cell_function")),
#                   ~ stringr::str_detect(stringr::str_to_lower(dplyr::coalesce(.x, "")), pattern))|cell_type%in%franken.meta.g$cell_type
#   ) |>
#   dplyr::distinct(root_id, .keep_all = TRUE) |>
#   dplyr::select(dplyr::all_of(keep_cols_banc))
# 
# 
# ## --- Save the filtered tables -------------------------------------------
# readr::write_csv(banc.meta.g,     file.path(out_dir, "banc_meta_gustatory.csv"))
# readr::write_csv(franken.meta.g,  file.path(out_dir, "franken_meta_gustatory.csv"))
# 
# ## --- Build copy-paste summary by cell_type -------------------------------
# ct_banc <- banc.meta.g    |> dplyr::count(cell_type, name = "n_banc")
# ct_fr   <- franken.meta.g |> dplyr::count(cell_type, name = "n_franken")
# 
# ct_summary_tbl <- dplyr::full_join(ct_banc, ct_fr, by = "cell_type") |>
#   tidyr::replace_na(list(n_banc = 0L, n_franken = 0L)) |>
#   dplyr::mutate(
#     delta_n   = n_banc - n_franken,
#     delta_pct = dplyr::if_else(n_franken > 0, 100 * (n_banc - n_franken) / n_franken, NA_real_)
#   ) |>
#   dplyr::arrange(dplyr::desc(delta_n))
# 
# # Pretty (string) version for email copy/paste
# ct_summary_pretty <- ct_summary_tbl |>
#   dplyr::mutate(
#     delta_pct = dplyr::case_when(
#       is.na(delta_pct) ~ NA_character_,
#       TRUE             ~ scales::percent(delta_pct / 100, accuracy = 0.1)
#     )
#   )
# 
# # Save the table
# readr::write_csv(ct_summary_pretty, file.path(out_dir, "gustatory_celltype_counts_delta.csv"))
# 
# # (Optional) print a quick view
# print(as.data.frame(ct_summary_pretty))

# Motor


fw <- franken_meta("SELECT fafb_id, manc_id, cell_type, cell_function, cell_function_detailed, body_part_sensory, body_part_effector from franken_meta")
f <- fw %>%
  dplyr::distinct(fafb_783 = fafb_id, 
                  manc_id = manc_id,
                  cell_type, 
                  cell_function,
                  cell_function_detailed,
                  body_part_sensory,
                  body_part_effector) %>%
  dplyr::filter(!is.na(cell_function)|!is.na(cell_function_detailed)|!is.na(body_part_sensory)|!is.na(body_part_effector))
readr::write_csv(x = f, 
                 file = "~/Downloads/banc_project_cell_functions_for_fafb_manc_cell_types.csv")


franken.meta <- franken_meta()
banc.neck.motor <- banc.meta %>%
  dplyr::filter(cell_class=="neck_motor_neuron"|grepl("prosternal",body_part_sensory)) %>%
  dplyr::distinct(root_id, root_626, supervoxel_id, 
                  hemilineage, flow, super_class, cell_sub_class, cell_type, 
                  cell_function, cell_function_detailed, 
                  body_part_sensory, body_part_effector, peripheral_target_type) %>%
  dplyr::mutate(dataset = "BANC")

fafb.neck.motor <- franken.meta %>%
  dplyr::filter(cell_class=="neck_motor_neuron"|grepl("prosternal",body_part_sensory)&!is.na(fafb_id)) %>%
  dplyr::distinct(root_783 = fafb_id,
                  hemilineage, flow, super_class, cell_sub_class, cell_type, 
                  cell_function, cell_function_detailed, 
                  body_part_sensory, body_part_sensory_detailed, body_part_effector, peripheral_target_type) %>% 
  dplyr::mutate(dataset = "FAFB")

manc.neck.motor <- franken.meta %>%
  dplyr::filter(cell_class=="neck_motor_neuron"|grepl("prosternal",body_part_sensory)&!is.na(manc_id)) %>%
  dplyr::distinct(bodyid = manc_id,
                  hemilineage, flow, super_class, cell_sub_class, cell_type, 
                  cell_function, cell_function_detailed, 
                  body_part_sensory, body_part_sensory_detailed, body_part_effector, peripheral_target_type) %>% 
  dplyr::mutate(dataset = "MANC")


readr::write_csv(x = banc.neck.motor, 
                 file = "~/Downloads/banc_neck_motor_and_sensory.csv")
readr::write_csv(x = fafb.neck.motor, 
                 file = "~/Downloads/fafb_neck_motor_and_sensory.csv")
readr::write_csv(x = manc.neck.motor, 
                 file = "~/Downloads/manc_neck_motor_and_sensory.csv")

# Check

library(coconatfly)
library(bancr)
bancr::register_banc_coconat()
banc_meta_create_cache(use_seatable=TRUE)
cf_cosine_plot(cf_ids('/type:^LB1|^LB2|^LB3|^LB4', datasets = c("banc", "malecns"), expand = TRUE), partners = "inputs", interactive = TRUE)
cf_cosine_plot(cf_ids('/type:SNxx33|SNch09|SNch13|SNch11|SNch02|SNch07|SNch03|SNch12|SNch04|SNxx20|SNch16|SNch15|SNxx21', datasets = c("banc","malecns"), expand = TRUE), interactive = TRUE)

