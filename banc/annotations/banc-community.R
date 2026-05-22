#' banc-community — Build the BANC ↔ MANC per-neuron community lookup.
#'
#' Joins NBLAST + manual matches for community-curated downstream analyses.
#'
#' @section Reads:
#'   - SeaTable `banc_meta`, `banc_manc_v1.2.1_nblast.feather`, `manc_neuprint_meta()`
#'
#' @section Writes:
#'   - `banc_manc_community.csv`

##########################
### ORGANISE COMMUNITY ###
##########################
source("banc/banc-startup.R")
library(malevnc)
library(neuprintr)

# Get our meta data
bc.meta <- banctable_query()
bc.manc <- bc.meta %>%
  dplyr::filter(!is.na(manc_match)) %>%
  dplyr::distinct(root_id,manc_match,manc_nblast) %>%
  dplyr::rename(banc_match = root_id,
                bodyid = manc_match,
                nblast_score = manc_nblast)

# Read NBLAST data
manc.meta <- manc_neuprint_meta()
manc.nblast.scores <- arrow::read_feather(file.path(banc.meta.save.path, "banc_manc_v1.2.1_nblast.feather")) %>%
  dplyr::rename(match_cell_type=cell_type) %>%
  dplyr::group_by(match_id) %>%
  dplyr::arrange(desc(score), .by_group = TRUE) %>%
  dplyr::slice_head(n = 3) %>%
  dplyr::mutate(banc_nblast_matches = paste0(pt_root_id,collapse=",")) %>%
  dplyr::mutate(banc_nblast_scores = paste0(score,collapse=",")) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(match_id, .keep_all = TRUE)

###############################
### adbominal neuron search ###
###############################

# Get neuprint data
anm.in <- neuprint_find_neurons(input_ROIs = "ANm",
                                dataset="manc:v1.2.1", 
                                all_segments = FALSE)
anm.out <- neuprint_find_neurons(input_ROIs = NULL,
                                 output_ROIs = "ANm",
                                 dataset="manc:v1.2.1", 
                                 all_segments = FALSE)
anm <- dplyr::full_join(anm.in,anm.out)
anm <- anm %>%
  dplyr::mutate(ANm.post=as.numeric(ANm.post),
                ANm.pre=as.numeric(ANm.pre)) %>%
  dplyr::filter(ANm.post>100|ANm.pre>100)
anm.meta <- manc_neuprint_meta(ids = anm$bodyid)
anm <- dplyr::left_join(anm.meta,bc.manc,by="bodyid") %>%
  dplyr::left_join(manc.nblast.scores,by=c("bodyid"="match_id")) %>%
  dplyr::select(bodyid, type, class, cell_function, entryNerve, exitNerve, hemilineage, predictedNt, transmission, post, pre, 
                banc_match, banc_nblast_matches, banc_nblast_scores) %>%
  dplyr::mutate(notes="",
                banc_match=ifelse(is.na(banc_match),"0",banc_match),
                banc_nblast_matches=ifelse(is.na(banc_nblast_matches),"0",banc_nblast_matches)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(ngl_link = bancsee(banc_ids = unique(c(banc_match,unlist(strsplit(banc_nblast_matches,split=",")))),
                                   manc_ids = c(bodyid))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(class = factor(class,levels=c("motor neuron", "efferent neuron", "efferent ascending", "sensory neuron", 
                                              "intrinsic neuron", "TBD", "sensory ascending", "ascending neuron", "descending neuron", 
                                              "glia"))) %>%
  dplyr::arrange(class)
colnames(anm) <- snakecase::to_snake_case(colnames(anm))

# Save
current_datetime <- Sys.time()
datetime_string <- format(current_datetime, "%Y-%m-%d")
file <- file.path(banc.meta.save.path,paste0(datetime_string,"_community_challenge_abdominal_neuromere.csv"))
readr::write_csv(anm, file = file)

###########################
### motor neuron search ###
###########################

# Get neuprint data
motor <- manc.meta %>%
  dplyr::filter(class%in%c("motor neuron","efferent neuron","efferent ascending")) %>%
  dplyr::left_join(bc.manc,by="bodyid") %>%
  dplyr::left_join(manc.nblast.scores,by=c("bodyid"="match_id")) %>%
  dplyr::select(bodyid, type, class, cell_function, entryNerve, exitNerve, hemilineage, predictedNt, transmission, post, pre, 
                banc_match, banc_nblast_matches, banc_nblast_scores) %>%
  dplyr::mutate(notes="",
                banc_match=ifelse(is.na(banc_match),"0",banc_match),
                banc_nblast_matches=ifelse(is.na(banc_nblast_matches),"0",banc_nblast_matches)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(ngl_link = bancsee(banc_ids = unique(c(banc_match,unlist(strsplit(banc_nblast_matches,split=",")))),
                                   manc_ids = c(bodyid))) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(type)
colnames(motor) <- snakecase::to_snake_case(colnames(motor))

# Save
current_datetime <- Sys.time()
datetime_string <- format(current_datetime, "%Y-%m-%d")
file <- file.path(banc.meta.save.path,paste0(datetime_string,"_manc_efferent_matches.csv"))
readr::write_csv(motor, file = file)

#############################
### sensory neuron search ###
#############################

# Get neuprint data
sensory <- manc.meta %>%
  dplyr::filter(class%in%c("sensory neuron","sensory ascending","Sensory TBD")) %>%
  dplyr::left_join(bc.manc,by="bodyid") %>%
  dplyr::left_join(manc.nblast.scores,by=c("bodyid"="match_id")) %>%
  dplyr::select(bodyid, type, class, cell_function, entryNerve, exitNerve, hemilineage, predictedNt, transmission, post, pre, 
                banc_match, banc_nblast_matches, banc_nblast_scores) %>%
  dplyr::mutate(notes="",
                banc_match=ifelse(is.na(banc_match),"0",banc_match),
                banc_nblast_matches=ifelse(is.na(banc_nblast_matches),"0",banc_nblast_matches)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(ngl_link = bancsee(banc_ids = unique(c(banc_match,unlist(strsplit(banc_nblast_matches,split=",")))),
                                   manc_ids = c(bodyid))) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(type)
colnames(sensory) <- snakecase::to_snake_case(colnames(sensory))

# Save
current_datetime <- Sys.time()
datetime_string <- format(current_datetime, "%Y-%m-%d")
file <- file.path(banc.meta.save.path,paste0(datetime_string,"_manc_sensory_matches.csv"))
readr::write_csv(sensory, file = file)

