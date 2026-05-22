#' banc-hemilineages — Build the hemilineage mapping worksheet.
#'
#' One-off generator for the hemilineage mapping table used in the BANC
#' annotation cascade. Pulls `franken_meta()` and assembles per-hemilineage
#' summary tables (Ito-Lee + Hartenstein); most blocks are currently
#' commented out — uncomment to regenerate.
#'
#' @section Reads:
#'   - `franken_meta()`
#'
#' @section Writes:
#'   - hemilineage worksheet CSVs under `<banc.meta.save.path>/`

#########################################
### BANC CREATE HEMILINEAGE WORKSHEET ###
#########################################
source("banc/banc-startup.R")
redo <- FALSE
version <- banc.nblast.version

# # Create hemilineage table
# franken.meta <- franken_meta()
# mappings <- franken.meta %>%
#   dplyr::filter(!is.na(FAFB_hartenstein_hemilineage),
#                 !is.na(FAFB_ito_lee_hemilineage)) %>%
#   dplyr::distinct(hemilineage=FAFB_ito_lee_hemilineage, hartenstein_hemilineage=FAFB_hartenstein_hemilineage)
# midbrain.hl.df <- franken.meta %>%
#   dplyr::mutate( FAFB_ito_lee_hemilineage =  dplyr::case_when(
#     is.na(FAFB_ito_lee_hemilineage) ~ FAFB_hartenstein_hemilineage,
#     TRUE ~ FAFB_ito_lee_hemilineage
#   )
#   ) %>%
#   dplyr::filter(grepl("FAFB",dataset),
#                 !is.na(fafb_id),
#                 !is.na(FAFB_ito_lee_hemilineage)|!is.na(FAFB_hartenstein_hemilineage)) %>%
#   dplyr::distinct(hemilineage = FAFB_ito_lee_hemilineage,
#                   side = side,
#                   cell_type,
#                   fafb_id) %>%
#   dplyr::group_by(hemilineage, side) %>%
#   dplyr::summarise(
#     no_cell_types = n_distinct(cell_type),
#     no_neurons = n(),
#     .groups = "drop"
#   ) %>%
#   dplyr::ungroup() %>%
#   dplyr::left_join(mappings, by = 'hemilineage') %>%
#   dplyr::mutate(tract_complex = substr(sub("_.*", "", hartenstein_hemilineage), 1, 4),
#                 region = "central_brain")
# vnc.hl.df <- franken.meta %>%
#   dplyr::filter(grepl("MANC",dataset),
#                 !is.na(manc_id),
#                 !is.na(hemilineage)) %>%
#   dplyr::distinct(hemilineage = hemilineage,
#                   side = side,
#                   cell_type,
#                   manc_id) %>%
#   dplyr::group_by(hemilineage, side) %>%
#   dplyr::summarise(
#     no_cell_types = n_distinct(cell_type),
#     no_neurons = n(),
#     .groups = "drop"
#   ) %>%
#   dplyr::mutate(hartenstein_hemilineage = hemilineage,
#                 tract_complex = substr(sub("_.*", "", hemilineage), 1, 2),
#                 region = "vnc")
# hl.df <- rbind.fill(midbrain.hl.df,vnc.hl.df)
# banctable_append_rows(hl.df,
#                       table = "hemilineages", 
#                       base = "banc_meta"
# )  

# Read franken meta
franken.meta <- franken_meta() 
banc.hls <- banctable_query("SELECT _id, hemilineage, region, tract_complex, side, ngl_complex FROM hemilineages")
bc <- banctable_query("SELECT root_id, side, hemilineage, cell_type, fafb_cell_type, manc_cell_type, fafb_match, manc_match from banc_meta")

# Add complexes
franken.meta.mod <- franken.meta %>%
  dplyr::mutate(hemilineage =  dplyr::case_when(
    is.na(hemilineage) ~ FAFB_hartenstein_hemilineage,
    TRUE ~ hemilineage
  )) %>%
  dplyr::mutate(FAFB_hartenstein_hemilineage =  dplyr::case_when(
    is.na(FAFB_hartenstein_hemilineage) ~ FAFB_ito_lee_hemilineage,
    TRUE ~ FAFB_hartenstein_hemilineage
  )) %>%
  dplyr::mutate(tract_complex = dplyr::case_when(
    grepl("MANC",dataset) ~ substr(sub("_.*", "", hemilineage), 1, 2),
    grepl("FAFB",dataset) ~ substr(sub("_.*", "", FAFB_hartenstein_hemilineage), 1, 4),
    TRUE ~ NA
  )) %>%
  dplyr::filter(!is.na(tract_complex))

# Sides
hemilineage.complexes <- sort(na.omit(unique(franken.meta.mod$tract_complex)))
vnc.hls <- subset(banc.hls,region=="vnc")$hemilineage
sides <- c("left","right","center")
df <- data.frame()
for(sid in unique(sides)){

  # Iterate over complex
  for(hlc in unique(hemilineage.complexes)){
    
    # Get hemilinages
    franken.meta.hls <- franken.meta.mod %>%
      dplyr::filter(tract_complex==hlc)
    hls <- na.omit(sort(unique(franken.meta.hls$hemilineage)))
    if(!length(hls)){
      message("Skipping: ", hlc)
      next
    }
    cts <- sort(na.omit(unique(franken.meta.hls$cell_type)))
  
    # Construct ngl link
    hl.cols <- rainbow(length(hls))
    names(hl.cols) <- hls
    banc.scene <- data.frame()
    mdf <- data.frame()
    for(hl in hls){
      if(hl%in%vnc.hls){
        is.vnc <- TRUE
      }else{
        is.vnc <- FALSE
      }
      banc.col <- lighten_color(hl.cols[[hl]],1.75)
      other.col <- lighten_color(hl.cols[[hl]],0.5)
      if(is.vnc){
        url <- "https://spelunker.cave-explorer.org/#!middleauth+https://global.daf-apis.com/nglstate/api/v1/6473787793997824"
        rename.views <- c("segmentation proofreading", "manc v1.2.1 imported")
        match.ids <- subset(franken.meta.mod,hemilineage==hl&side==sid)$manc_id
        bc.hls <- bc %>%
          dplyr::mutate(cell_type = gsub("auto:","",cell_type),
                        fafb_cell_type = gsub("auto:","",fafb_cell_type),
                        manc_cell_type = gsub("auto:","",manc_cell_type)) %>%
          dplyr::filter(
            cell_type %in% cts | fafb_cell_type %in% cts | manc_cell_type %in% cts | manc_match %in% match.ids
          ) %>%
          dplyr::distinct(root_id, cell_type, side) %>%
          dplyr::left_join(franken.meta.hls %>%
                             dplyr::select(cell_type, hemilineage) %>%
                             dplyr::distinct(cell_type, .keep_all = TRUE),
                           by = "cell_type")
        banc.ids <- subset(bc.hls,hemilineage==hl&side==sid)$root_id
        hl.link <- bancsee(url = url,
                           banc_ids = c(banc.ids),
                           manc_ids = c(match.ids),
                           banc.cols = banc.col,
                           manc.cols = other.col)
      }else{
        url <- "https://spelunker.cave-explorer.org/#!middleauth+https://global.daf-apis.com/nglstate/api/v1/4512779479285760"
        rename.views <- c("segmentation proofreading", "fafb v783 imported")
        match.ids <- subset(franken.meta.mod,hemilineage==hl&side==sid)$fafb_id
        bc.hls <- bc %>%
          dplyr::mutate(cell_type = gsub("auto:","",cell_type),
                        fafb_cell_type = gsub("auto:","",fafb_cell_type),
                        manc_cell_type = gsub("auto:","",manc_cell_type)) %>%
          dplyr::filter(
            cell_type %in% cts | fafb_cell_type %in% cts | manc_cell_type %in% cts | fafb_match %in% match.ids
          ) %>%
          dplyr::distinct(root_id, cell_type, side) %>%
          dplyr::left_join(franken.meta.hls %>%
                             dplyr::select(cell_type, hemilineage) %>%
                             dplyr::distinct(cell_type, .keep_all = TRUE),
                           by = "cell_type")
        banc.ids <- subset(bc.hls,hemilineage==hl&side==sid)$root_id
        hl.link <- bancsee(url = url,
                           banc_ids = c(banc.ids),
                           fafb_ids = c(match.ids),
                           banc.cols = banc.col,
                           fafb.cols = other.col)
      }
      hl.link.decoded <- ngl_decode_scene(hl.link)
      chosen.views <- c()
      for(rv in rename.views){
        nv <- paste(hl, hl.link.decoded$layers[[rv]]$name,collapse=" ")
        hl.link.decoded$layers[[rv]]$name <- nv
        names(hl.link.decoded$layers)[names(hl.link.decoded$layers)==rv] <- nv
        chosen.views <- c(chosen.views,nv)
      }
      if(!length(banc.scene)){
        banc.scene <- hl.link.decoded
      }else{
        for(nv in chosen.views){
          banc.scene$layers[[nv]] <- hl.link.decoded$layers[[nv]]
        }
      }
      ldf <- data.frame(
        hemilineage = hl,
        ngl_hemilineage = hl.link,
        tract_complex = hlc,
        side = sid,
        no_banc_cell_types = length(unique( subset(bc.hls,hemilineage==hl&side==sid)$cell_type)),
        no_banc_neurons = length(banc.ids)
      )
      mdf <- rbind(mdf,ldf)
    }
    banc.url <- banc_shorturl(ngl_encode_url(banc.scene))
    mdf$ngl_complex <- banc.url
    df <- rbind(df,mdf)
  }
}
update <- dplyr::left_join(banc.hls %>%
                             dplyr::select(`_id`,hemilineage,side), 
                           df %>%
                             dplyr::distinct(hemilineage, side, .keep_all = TRUE), 
                          by = c("hemilineage","side")) %>%
  as.data.frame()
update$no_banc_cell_types <- as.character(update$no_banc_cell_types)
update$no_banc_neurons <- as.character(update$no_banc_neurons)
banctable_update_rows(update,
                      table = "hemilineages",
                      base = "banc_meta"
)




