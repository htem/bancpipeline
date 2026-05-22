###############################
### GET HEMIBRAIN META DATA ###
###############################
source("banc/banc-startup.R")
library(neuprintr)

# Get the meta data
all.neurons.meta <- neuprint_search("Traced",field="status",dataset="hemibrain:v1.2.1")
all.neurons.meta <- subset(all.neurons.meta, statusLabel!="Leaves")
all.bodyids <- all.neurons.meta$bodyid
all.bodyids <- unique(c(all.bodyids,
                       kc.ids,
                       apl.ids,
                       rn.ids,
                       orn.ids,
                       hrn.ids,
                       pn.ids,
                       upn.ids,
                       mpn.ids,
                       vppn.ids,
                       dan.ids,
                       mbon.ids,
                       alln.ids,
                       lhn.ids,
                       ton.ids,
                       lc.ids,
                       cent.ids,
                       dn.ids
))
all.bodyids <- unique(all.bodyids)
all.bodyids <- na.omit(all.bodyids)
hb.meta.orig <- neuprint_get_meta(all.bodyids, dataset="hemibrain:v1.2.1")

# Organise
colnames(hb.meta.orig) <- snakecase::to_snake_case(colnames(hb.meta.orig))
hb.meta <- hb.meta.orig %>%
  dplyr::mutate(cell_type = type) %>%
  dplyr::select(bodyid,
                cropped,
                post,
                pre,
                voxels,
                soma,
                cell_body_fiber,
                cell_type,
                type,
                name
  ) %>%
  dplyr::mutate(cell_class = case_when(
    grepl("^FB|^SA1|^SAF|^SA2|^SA3",cell_type) ~ "CX_tangential",
    grepl("^vDelta|^hDelta|^FC|^FS|^FR",cell_type) ~ "CX_columnar",
    grepl("^LNO|^LCNO|^GLNO",cell_type) ~ "CX_nodulus",
    grepl("^ER|^EL|^EPG|^Ex",cell_type) ~ "CX_ellipsoid_body",
    grepl("^SpsP|^P6|^P1-9|^IbSpsP|^Delta7",cell_type) ~ "CX_protocerebral_bridge",
    grepl("^PF|^PFN|^PFL|^PFG|^PFR",cell_type) ~ "CX_protocerebral_bridge",
    bodyid %in% hemibrainr::upn.ids ~ "uPN",
    bodyid %in% hemibrainr::vppn.ids ~ "THPN",
    bodyid %in% hemibrainr::mpn.ids ~ "mPN",
    bodyid %in% hemibrainr::pn.ids ~ "ALPN",
    bodyid %in% hemibrainr::mbon.ids ~ "MBON",
    bodyid %in% hemibrainr::kc.ids ~ "KC",
    bodyid %in% hemibrainr::dan.ids ~ "DAN",
    bodyid %in% hemibrainr::alln.ids ~ "ALLN",
    bodyid %in% hemibrainr::ton.ids ~ "TON",
    bodyid %in% hemibrainr::dn.ids ~ "DN",
    bodyid %in% hemibrainr::apl.ids ~ "MB",
    bodyid %in% hemibrainr::lc.ids ~ "VPN",
    bodyid %in% hemibrainr::hrn.ids ~ "HRN",
    bodyid %in% hemibrainr::orn.ids ~ "ORN",
    bodyid %in% hemibrainr::rn.ids ~ "RN",
    bodyid %in% hemibrainr::lhn.ids ~ "LHN",
    bodyid %in% hemibrainr::cent.ids ~ "CENT",
    TRUE ~ "unknown"),
    super_class = case_when(
      grepl("^FB|^vDelta|^hDelta|^FC|^SA1|^SAF|^SA2|^SA3|^PF|IbSpsP|^EPG|^PFN|^Delta|^LNO|^LCNO|^PFL|^PFG|^ER|^FS|P1-9|^EL$|^FR|^Ex|^SpsP|^PFR|^SpsP|^GLNO|^P6",cell_type) ~ "central_complex",
      bodyid %in% hemibrainr::upn.ids ~ "projection_neuron",
      bodyid %in% hemibrainr::vppn.ids ~ "projection_neuron",
      bodyid %in% hemibrainr::mpn.ids ~ "projection_neuron",
      bodyid %in% hemibrainr::pn.ids ~ "projection_neuron",
      bodyid %in% hemibrainr::mbon.ids ~ "mushroom_body",
      bodyid %in% hemibrainr::kc.ids ~ "mushroom_body",
      bodyid %in% hemibrainr::dan.ids ~ "mushroom_body",
      bodyid %in% hemibrainr::alln.ids ~ "antennal_lobe",
      bodyid %in% hemibrainr::lhn.ids ~ "lateral_horn",
      bodyid %in% hemibrainr::ton.ids ~ "third_order_chemosensory",
      bodyid %in% hemibrainr::dn.ids ~ "descending",
      bodyid %in% hemibrainr::apl.ids ~ "mushroom_body",
      bodyid %in% hemibrainr::lc.ids ~ "projection_neuron",
      bodyid %in% hemibrainr::hrn.ids ~ "sensory_neuron",
      bodyid %in% hemibrainr::orn.ids ~ "sensory_neuron",
      bodyid %in% hemibrainr::rn.ids ~ "sensory_neuron",
      bodyid %in% hemibrainr::cent.ids ~ "centrifugal",
      TRUE ~ "unknown"), 
    nucleus_id = "0") %>%
  dplyr::mutate(cell_function = dplyr::case_when(
    grepl("sensory",super_class) & cell_class=='olfactory' ~ "olfactory",
    grepl("sensory",super_class) & cell_class=='unknown_sensory' ~ "unknown_sensory",
    grepl("sensory",super_class) & cell_class=='visual' ~ "visual",
    grepl("sensory",super_class) & cell_class=='thermosensory' ~ "thermosensory",
    grepl("sensory",super_class) & cell_class=='mechanosensory' ~ "mechanosensory",
    grepl("sensory",super_class) & cell_class=='hygrosensory' ~ "hygrosensory",
    grepl("sensory",super_class) & cell_class=='enteric_gustatory' ~ "gustatory",
    grepl("sensory",super_class) & cell_class=='enteric gustatory' ~ "gustatory",
    grepl("sensory",super_class) & cell_class=='gustatory' ~ "gustatory",
    grepl("sensory",super_class) & cell_class=='ocellar' ~ "ocellar",
    grepl("sensory",super_class) ~ "unknown_sensory",
    grepl("endocrine",super_class) ~ "endocrine",
    grepl("motor",super_class) ~ "motor",
    grepl("motor",cell_class) ~ "motor",
    grepl("efferent",super_class) &cell_class=="pars_intercerebralis" ~ "endocrine",
    grepl("efferent", super_class) & cell_class=="pars_lateralis" ~ "endocrine",
    grepl("efferent", super_class) & cell_class=="ocellar" ~ "ocellar",
    TRUE ~ NA
  ))

# Get known nt data
ft.midbrain <- fafbseg::flytable_query("select cell_type, hemibrain_type, top_nt, known_nt, known_nt_source from info")
ft.optic <- fafbseg::flytable_query("select cell_type, hemibrain_type, top_nt, known_nt, known_nt_source from optic")
ft <- rbind(ft.midbrain, ft.optic)
ft.nt <- ft %>%
  dplyr::mutate(hemibrain_type = ifelse(is.na(hemibrain_type),cell_type,hemibrain_type)) %>%
  dplyr::filter(!is.na(hemibrain_type),!is.na(known_nt)) %>%
  dplyr::select(-cell_type) %>%
  dplyr::distinct(hemibrain_type, .keep_all = TRUE) %>%
  dplyr::arrange(known_nt, top_nt)
hb.meta <- dplyr::left_join(hb.meta, ft.nt, by=c("type"="hemibrain_type"))

# Save
readr::write_csv(hb.meta, file.path(banc.meta.save.path,"hemibrain_meta.csv"))

# Announce
message("##### BANCpipeline: manc meta updated #####")
message(sprintf("##### we have meta for: %s neurons", nrow(hb.meta)))

