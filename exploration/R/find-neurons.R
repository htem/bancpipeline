### -----------------------------------------------------------------
### Legacy exploratory analysis by a member of the Wilson Lab.
### Historical record only — not part of the bancpipeline release.
### -----------------------------------------------------------------

# Libraries
source("banc/banc-startup.R")

# Settings
version <- "elastix_tpsreg_240721"

# Get meta data
fw.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"flywire_meta.csv"), 
                                            col_types = hemibrainr:::sql_col_types))
fw.ids <- unique(fw.meta$root_783)
mc.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"manc_meta.csv"), 
                                            col_types = hemibrainr:::sql_col_types))
mc.ids <- unique(mc.meta$bodyid)

# flywire 783 IDs
fw.search <- c(
  # ANs upstream of FB8A/H
  "720575940635527092",# 720575941641744978
  
  # ANs upstream of FB8C
  "720575940629247682", # 720575941471999072 --> 12059 (AN27X009)
  "720575940633304601",
  "720575940623302988",# 720575941583214214
  
  # FB1C
  "720575940612718563", # 720575941544590697 --> (AN07B004)
  "720575940626462014", # 720575941652988377
  
  # FB1G
  "720575940618054981", # 720575941483610690 --> 12492 (AN27X017)
  "720575940612718563", # 720575941544590697
  
  # SA3/SAF
  "720575940614668347", # 720575941449241466
  "720575940623758377", # 
  
  # (DNp27->) SLP355 -> SA3/SAF
  "720575940621599741", # 720575941543458665 --> 13101 (ANXXX136)
  "720575940631295506", # 720575941468470135
  
  # SLP257 -> SLP355
  "720575940636263543", # 720575941533961498 --> 12460 (AN05B101)
  
  # SMP169 -> SLP355
  # SMP169 -> FB8C
  "720575940624810782", # 720575941488228809 --> 11789 (AN05B097)
  "720575940633764247", # 720575941459327888
  "720575940633292769", # 720575941579234425 --> 20129 (AN09B018)
  "720575940606762953", # 720575941512249667
  "720575940623302988", # 720575941583214214
  "720575940614595643", # 720575941547143028 --> 11358 (AN05B096)
  "720575940639521902", # 720575941503733925 --> 28408 (SAxx01)
  "720575940618831213", # 720575941471906602
  "720575940629247682", # 720575941471999072
  "720575940635527092", # 720575941641744978
  
  # Rachel's selection
  "720575940627490065", # 720575941512249667
  "720575940653877409", # 720575941512249667
  "720575940611003858", # 720575941512249667
  "720575940620686360", # 720575941520348919
  "720575940606762953" # 720575941512249667
)

# Get first order connectivity
# chosen <- c("^SA1|^SA2|^SA3|^SAF|^FB7B|FB8|^FB9A|FB1C|FB1D|FB1G|^CB.FBT|FB7E|^FB1E|FB7|FBTI|FB1E1|FB1I1|LCNOpm|PFNm|IbSpsP")
# ft <- fafbseg::flytable_cell_types()
# fw.meta.chosen <- ft %>%
#   dplyr::filter(grepl(chosen, cell_type)|grepl(chosen, hemibrain_type))
# fw.inputs = fafbseg::flywire_partner_summary(unique(fw.meta.chosen$root_783),
#                                                          partners = "inputs",
#                                                          threshold = 5,
#                                                          remove_autapses = TRUE,
#                                                          cleft.threshold = 0.5,
#                                                          details = TRUE,
#                                                          roots = TRUE,
#                                                          Verbose = FALSE)

# Compile FAFB NBLAST hits
fafb.nblast.folder <- file.path(banc.nblast.fafb.save.path,"results",version)
fafb.nblast.files <- list.files(fafb.nblast.folder, pattern = "\\.csv", full.names = TRUE, recursive = FALSE)
by.query <- foreach::foreach(mfile = fafb.nblast.files) %do% {
    id <- gsub(".*_root_id_|\\.csv","",basename(mfile))
    mdf <- readr::read_csv(mfile, col_types = banc.col.types, show_col_types = FALSE)
    numeric_columns <- sapply(mdf, is.numeric)
    mdf[numeric_columns] <- round(mdf[numeric_columns], digits = 2)
    mdf$query <- id
    mdf <- mdf[mdf$nb>0.3,]
    mdf <- dplyr::arrange(mdf, dplyr::desc(nb))
    mdf
}
by.query <- by.query[unlist(lapply(by.query,is.data.frame))]
banc.meta.fafb.nb <- do.call(plyr::rbind.fill, by.query) %>%
  dplyr::rename(fafb_nblast_match=root_783, root_id=query) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(fafb_nblast = round(as.numeric(nb),2)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(root_id, fafb_nblast_match, fafb_nblast) %>%
  dplyr::left_join(fw.meta[,c("root_783","ito_lee_hemilineage", 
                              "top_nt", "flow", "super_class", "cell_class", 
                              "cell_type", "hemibrain_type", "hemibrain_match")], 
                   by = c("fafb_nblast_match"="root_783")) 

# Find
fw.seach.in.banc <- banc.meta.fafb.nb %>%
  dplyr::filter(fafb_nblast_match%in%fw.search) %>%
  dplyr::arrange(desc(fafb_nblast)) %>%
  dplyr::distinct(fafb_nblast_match, .keep_all = TRUE)

# Compile manc NBLAST hits
manc.nblast.folder <- file.path(banc.nblast.manc.save.path,"results",version)
manc.nblast.files <- list.files(manc.nblast.folder, pattern = "\\.csv", full.names = TRUE, recursive = FALSE)
by.query <- foreach::foreach(mfile = manc.nblast.files) %do% {
  id <- gsub(".*_root_id_|\\.csv","",basename(mfile))
  mdf <- readr::read_csv(mfile, col_types = banc.col.types, show_col_types = FALSE)
  numeric_columns <- sapply(mdf, is.numeric)
  mdf[numeric_columns] <- round(mdf[numeric_columns], digits = 2)
  mdf$query <- id
  mdf <- mdf[mdf$nb>0.3,]
  mdf <- dplyr::arrange(mdf, dplyr::desc(nb))
  mdf
}
by.query <- by.query[unlist(lapply(by.query,is.data.frame))]
banc.meta.manc.nb <- do.call(plyr::rbind.fill, by.query) %>%
  dplyr::rename(manc_nblast_match=bodyid, root_id=query) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(manc_nblast = round(as.numeric(nb),2)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(root_id, manc_nblast_match, manc_nblast) %>%
  dplyr::left_join(mc.meta[,c("bodyid","hemilineage","conf_nt", "entry_nerve", "exit_nerve", "subclass", 
                              "type")], 
                   by = c("manc_nblast_match"="bodyid")) %>%
  dplyr::rename(cell_type = type,
                cell_class = subclass,
                top_nt = conf_nt)

# Find
fw.seach.in.manc <- banc.meta.manc.nb %>%
  dplyr::filter(root_id%in%fw.seach.in.banc$root_id) %>%
  dplyr::arrange(dplyr::desc(manc_nblast)) %>%
  dplyr::distinct(root_id, .keep_all = TRUE)