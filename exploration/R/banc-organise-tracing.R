### -----------------------------------------------------------------
### Legacy exploratory analysis by a member of the Wilson Lab.
### Historical record only — not part of the bancpipeline release.
### -----------------------------------------------------------------

# load packages
library(bancr)

# Direct tools ar BANC dataset
choose_banc()

# Get our meta data
bc.meta <- banctable_query() %>%
  dplyr::filter(region=="neck_connective")

# Get save folder
tracing.folder <- "tracing"
dir.create(tracing.folder, showWarnings = FALSE)

#####################################################################
### compile list of neurons from neck connective that need review ###
#####################################################################

# Get folders of issuesome neurons on the lab fileserve
miror.match <- "/Volumes/neurobio/wilsonlab/banc/matching/mirror/correct"
red_wrong <- file.path(miror.match,"red_wrong")
blue_wrong <- file.path(miror.match,"blue_wrong")
investigate <- file.path(miror.match,"investigate")

# Get files
red_wrong.pngs <- list.files(red_wrong, pattern = "png$")
blue_wrong.pngs <- list.files(blue_wrong, pattern = "png$")
investigate.pngs <- list.files(investigate, pattern = "png$")

# Get the 'issue ID'
red.issues <-  gsub(".*_hit_id_|_hit_nucleus_id_.*","",red_wrong.pngs) # the red neuron is 'hit_id', i.e. the nblast hit
red.bluecomparison <-  gsub(".*_root_id_|_nucleus_id_.*","",red_wrong.pngs) # the blue neuron is 'root_id', i.e. the query
blue.issues <-  gsub(".*_root_id_|_nucleus_id_.*","",blue_wrong.pngs)
blue.redcomparison <-  gsub(".*_hit_id_|_hit_nucleus_id_.*","",blue_wrong.pngs)
investigate.ids <- gsub(".*_root_id_|_nucleus_id_.*","",investigate.pngs) # I guess people usually mean the blue neuron?
investigate.comparison <- gsub(".*_hit_id_|_hit_nucleus_id_.*","",investigate.pngs)

# Compile data frame
tracing.issues.pngs <- data.frame(
  reviewed_root_id = c(red.issues, blue.issues, investigate.ids),
  comparison_root_id = c(red.bluecomparison, blue.redcomparison, investigate.comparison),
  issue = c(rep("tracing:red neuron in file",length(red.issues)),
            rep("tracing: blue neuron in file",length(blue.issues)),
            rep("strange morphology needs follow-up",length(investigate.ids))
            )
)
tracing.issues.pngs$reviewed_root_id_latest <- banc_latestid(tracing.issues.pngs$reviewed_root_id)
tracing.issues.pngs$comparison_root_id_latest <- banc_latestid(tracing.issues.pngs$comparison_root_id)

# Cut the list down by ensuring we only have unique combinations
tracing.issues.pngs <- tracing.issues.pngs %>%
  dplyr::rowwise() %>%
  dplyr::mutate(comparison = paste(sort(c(reviewed_root_id,comparison_root_id)),collapse="_")) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(comparison, issue, .keep_all = TRUE) %>%
  dplyr::select(-comparison) %>%
  dplyr::mutate(notes="",
                see_png = TRUE) %>%
  dplyr::select(reviewed_root_id, comparison_root_id,
                reviewed_root_id_latest, comparison_root_id_latest,
                issue, notes)

##########################################
### combine with updates from seatable ###
##########################################
tracing.issues.st <- bc.meta %>%
  dplyr::filter(region=="neck_connective",
                grepl("TRACING|ISSUE|MERGE",status)) %>%
  dplyr::rename(reviewed_root_id=root_id,
                comparison_root_id=banc_match,
                issue=notes) %>%
  dplyr::distinct(reviewed_root_id,
                  comparison_root_id,
                  fafb_match,
                  manc_match,
                  hemibrain_match,
                  nucleus_id,
                  status,
                  issue) %>%  
  dplyr::rowwise() %>%
  dplyr::mutate(comparison = paste(sort(c(reviewed_root_id,comparison_root_id)),collapse="_")) %>%
  dplyr::ungroup()

# updates IDs
tracing.issues.st$reviewed_root_id_latest <- banc_updateids(tracing.issues.st$reviewed_root_id)
tracing.issues.st$comparison_root_id_latest[is.na(tracing.issues.st$comparison_root_id)] <- banc_updateids(tracing.issues.st$comparison_root_id[is.na(tracing.issues.st$comparison_root_id)])

# combine
tracing.issues <- plyr::rbind.fill(tracing.issues.pngs,tracing.issues.st) %>%
  dplyr::mutate(issue=ifelse(is.na(issue),"tracing issue",issue),
                notes=ifelse(is.na(issue),"",issue),
                see_png=reviewed_root_id%in%tracing.issues.pngs$reviewed_root_id) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(comparison = paste(sort(c(reviewed_root_id_latest,comparison_root_id_latest)),collapse="_")) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(comparison, .keep_all = TRUE) %>%
  dplyr::distinct(reviewed_root_id_latest, .keep_all = TRUE) %>%
  dplyr::mutate(comparison_root_id = ifelse(is.na(comparison_root_id),"0",comparison_root_id),
                comparison_root_id_latest = ifelse(is.na(comparison_root_id_latest),"0",comparison_root_id_latest),
                fafb_match = ifelse(is.na(fafb_match),"0",fafb_match),
                manc_match = ifelse(is.na(manc_match),"0",manc_match),
                hemibrain_match = ifelse(is.na(hemibrain_match),"0",hemibrain_match),
                nucleus_id = ifelse(is.na(nucleus_id),"0",nucleus_id),
                reviewed_root_id = ifelse(is.na(reviewed_root_id),"0",reviewed_root_id),
                reviewed_root_id_latest = ifelse(is.na(reviewed_root_id_latest),"0",reviewed_root_id_latest)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(ngl_link = bancsee(banc_ids = c(reviewed_root_id_latest,comparison_root_id_latest),
                                   fafb_ids = c(fafb_match),
                                   manc_ids = c(manc_match),
                                   nuclei_ids = c(nucleus_id))) %>%
  dplyr::ungroup() %>%
  dplyr::select(reviewed_root_id_latest, comparison_root_id_latest,
                reviewed_root_id, comparison_root_id,
                nucleus_id, fafb_match, manc_match,
                status, issue, see_png,
                ngl_link, notes) %>%
  dplyr::arrange(desc(comparison_root_id_latest),desc(see_png),status,issue)

# Save
current_datetime <- Sys.time()
datetime_string <- format(current_datetime, "%Y-%m-%d")
file <- file.path(tracing.folder,paste0(datetime_string,"_neck_connective_review.csv"))
readr::write_csv(tracing.issues, file = file)

#####################################################################
### compile list of neurons from neck connective that need nuclei ###
#####################################################################

# Get nuclei positions
banc.nulcei <- banc_nuclei(rawcoords = FALSE, table = "both")
banc.nulcei.dupes <- unique(banc.nulcei$pt_root_id[duplicated(banc.nulcei$pt_root_id)])
banc.nulcei.dupes <- as.character(banc.nulcei.dupes)
banc.nulcei.dupes <- banc.nulcei.dupes[banc.nulcei.dupes!="0"]

# Get lists of
nuclei.check <- bc.meta %>%
  dplyr::mutate(duplicated_nuclei = root_id %in% banc.nulcei.dupes) %>%
  dplyr::mutate(missing_nucleus = dplyr::case_when(
    grepl("afferent|sensory|glia",super_class) ~ FALSE,
    grepl("sensory|glia|innervates haltere",cell_class) ~ FALSE,
    grepl("sensory|glia",notes) ~ FALSE,
    grepl("sensory|glia",cell_type) ~ FALSE,
    !nucleus_id%in%c("0","NA",""," ") ~ FALSE,
    duplicated_nuclei ~ FALSE,
    TRUE ~ TRUE
  )) %>%
  dplyr::filter(duplicated_nuclei|missing_nucleus) %>%
  dplyr::distinct(root_id, cell_class, nucleus_id,
                  nucleus_position_nm, root_position_nm, 
                  duplicated_nuclei, missing_nucleus) %>%
  dplyr::arrange(duplicated_nuclei, missing_nucleus) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(ngl_link = bancsee(banc_ids = c(root_id),
                                   nuclei_ids = c(nucleus_id))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(notes="")
  
# Save
current_datetime <- Sys.time()
datetime_string <- format(current_datetime, "%Y-%m-%d")
file <- file.path(tracing.folder,paste0(datetime_string,"_neck_connective_nucleus_review.csv"))
readr::write_csv(nuclei.check, file = file)

#######################################################################
### compile list of neurons from neck connective that are too small ###
#######################################################################

# Get the tadpoles
tadpoles <-  bc.meta %>%
  dplyr::filter(l2_cable_length_um<1000,
                l2_cable_length_um!=0,
                !grepl("glia",cell_class)) %>%
  dplyr::rename(reviewed_root_id=root_id) %>%
  dplyr::mutate(reviewed_root_id_latest = banc_latestid(reviewed_root_id),
                l2_cable_length_um=l2_cable_length_um/100) %>%
  dplyr::distinct(reviewed_root_id, reviewed_root_id_latest, nucleus_id, cell_class, l2_cable_length_um) %>%
  dplyr::arrange(l2_cable_length_um) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(ngl_link = bancsee(banc_ids = c(reviewed_root_id_latest),
                                   nuclei_ids = c(nucleus_id))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(notes="")

# Save
current_datetime <- Sys.time()
datetime_string <- format(current_datetime, "%Y-%m-%d")
file <- file.path(tracing.folder,paste0(datetime_string,"_neck_connective_too_small.csv"))
readr::write_csv(tadpoles, file = file)