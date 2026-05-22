#' banc-sensory-jo-orn — Find missing Johnston's-organ sensory neurons and fix ORN typing.
#'
#' Per-task curation: scans for missing JO sensory neurons via the tracing
#' sheet and applies ORN cell-type corrections via the proofreading-notes
#' table.
#'
#' @section Reads:
#'   - SeaTable `banc_meta`
#'   - `banc_proofreading_notes()`
#'   - Google Sheet referenced by `banc.keys` (tracing sheet)
#'
#' @section Writes:
#'   - SeaTable `banc_meta`: cols `cell_type`, `super_class`,
#'     `cell_class` (JO + ORN rows only)
#'
#' @section Notes:
#'   - Split out 2026-05-21 from `banc/utilities/banc-tracing.R`.

###############################################################################
### BANC annotation: Johnston's-organ sensory neurons + ORN corrections
###
### Finds missing JO sensory neurons via the tracing sheet and applies ORN cell-type corrections.
###
### Originally part of banc/utilities/banc-tracing.R; split out 2026-05-21
### into per-task scripts under banc/annotations/ because the original was
### never meant to be run end-to-end — each section is a one-off curation
### task. Sections may share state; if so, run them in order.
###############################################################################
source("banc/banc-startup.R")

#########################################
### Tracing sheet to find missing JOs ###
#########################################

## Search for JO neurons
library(bancr)
banc.proofreading.notes <- banc_proofreading_notes()
bc <- banctable_query()
jos <- subset(bc,grepl("sensory",super_class)&grepl("brain",region)&!grepl("^R",cell_type)&!is.na(cell_type)&cell_type!="sensory_fragment")
jos <- subset(bc,!grepl("^ORN|^HRN|^TRN|^THRN",cell_type)&grepl("sensory",super_class)&grepl("brain",region))
proof <- c(bc$root_id,as.character(banc.proofreading.notes$root_id)) # roughly proofread
jos.id <- bancr::banc_rootid(na.omit(unique(jos$supervoxel_id)))
jos.id <- setdiff(jos.id,"0")
conns <- pbapply::pblapply(jos.id, function(x) tryCatch({df<-banc_partners(x, partners = "output", synapse_table = "synapses_v2")
df$post_pt_root_id <- as.character(df$post_pt_root_id)
df$pre_pt_root_id <- as.character(df$pre_pt_root_id)
df
},
                                                        error=function(e)NULL))
conn <- do.call(plyr::rbind.fill,conns[!unlist(sapply(conns, is.null))])
conn.top <- conn %>%
  dplyr::group_by(post_pt_root_id) %>%
  dplyr::mutate(count = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(post_pt_root_id,count) %>%
  dplyr::arrange(desc(count)) %>%
  dplyr::filter(count>=100)
top.partners <- unique(conn.top$post_pt_root_id)
top.partners <- top.partners[top.partners%in%proof]
conns2 <- pbapply::pblapply(top.partners, function(x) tryCatch(banc_partners(x, partners = "input", synapse_table = "synapses_v2"),error= function(e)NULL))
conn2 <- do.call(plyr::rbind.fill,conns2[unlist(sapply(conns2, is.data.frame))])
conn.top2 <- conn2 %>%
  dplyr::group_by(pre_pt_root_id) %>%
  dplyr::mutate(count = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::select(pre_pt_root_id,pre_pt_position,pre_pt_supervoxel_id,count) %>%
  dplyr::distinct(pre_pt_root_id,.keep_all=TRUE) %>%
  dplyr::arrange(desc(count)) %>%
  dplyr::filter(count>=15,
                !as.character(pre_pt_root_id) %in% jos.id,
                !as.character(pre_pt_root_id) %in% proof)

# Write as .csv file
conn.top2$position = nat::xyzmatrix2str(conn.top2$pre_pt_position)
jo.search <- conn.top2 %>%
  as.data.frame() %>%
  dplyr::mutate(root_id = as.character(pre_pt_root_id),
                supervoxel_id = as.character(pre_pt_supervoxel_id),
                status = "",
                notes="",
                annotator="") %>%
  dplyr::distinct(root_id,supervoxel_id,position,count,status,notes,annotator)

# Get better root positions
jo.search$position <- NA
n <- nrow(jo.search)
pb <- txtProgressBar(min = 0, max = n, style = 3)
for (i in seq_len(n)) {
  id <- jo.search$root_id[i]
  res <- try({
    l2   <- bancr::banc_read_l2skel(id, rawcoords = TRUE)
    strh <- nat::strahler_order(l2[[1]])
    pos  <- nat::xyzmatrix(l2[[1]]$d)[which.max(strh$points), ]
    pos.raw <- nat:::xyzmatrix2str(bancr:::banc_nm2raw(pos))
    l2_nodes <- nrow(l2[[1]]$d)
    list(pos.raw=pos.raw,l2_nodes=l2_nodes)
  }, silent = TRUE)
  if (!inherits(res, "try-error")) {
    jo.search$position[i] <- res[["pos.raw"]]
    jo.search$l2_nodes <- res[["l2_nodes"]]
  } else {
    jo.search$position[i] <- NA
    jo.search$l2_nodes[i] <- NA
  }
  setTxtProgressBar(pb, i)
}
close(pb)

# Save
banc.tracing.save.path <- "tracing/tracing_issues/"
current_datetime <- Sys.time()
datetime_string <- format(current_datetime, "%Y-%m-%d")
readr::write_csv(jo.search, file.path(banc.tracing.save.path,paste0(datetime_string,"_banc_brain_sensory_search.csv")))

# Confident that these are mostly sensory neurons...
bancr::banctable_append_rows(base='banc_meta',
                             table = "banc_meta",
                             df = jo.search %>%
                               dplyr::filter(!root_id %in% bc$root_id) %>%
                               dplyr::select(root_id,supervoxel_id,position) %>%
                               dplyr::mutate(status="SENSORY_SEARCH_2"),
                             chunksize = 1000)  

# Search for missing VNC neurons
vncs <- subset(bc,grepl("sensory",super_class)&grepl("ventral|VNC",region))
vncs.id <- bancr::banc_rootid(na.omit(unique(vncs$supervoxel_id)))
vncs.id <- setdiff(vncs.id,"0")
conns <- pbapply::pblapply(vncs.id, function(x) tryCatch(banc_partners(x, partners = "output", synapse_table = "synapses_v2"),error=function(e)NULL))
conn <- do.call(plyr::rbind.fill,conns[!unlist(sapply(conns, is.null))])
conn.top <- conn %>%
  dplyr::group_by(post_pt_root_id) %>%
  dplyr::mutate(count = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(post_pt_root_id,count) %>%
  dplyr::arrange(desc(count)) %>%
  dplyr::filter(count>=100)
top.partners <- unique(conn.top$post_pt_root_id)
top.partners <- top.partners[top.partners%in%proof]
conns2 <- pbapply::pblapply(top.partners, function(x) tryCatch(banc_partners(x, partners = "input", synapse_table = "synapses_v2"),error= function(e)NULL))
conn2 <- do.call(plyr::rbind.fill,conns2[unlist(sapply(conns2, is.data.frame))])
conn.top2 <- conn2 %>%
  dplyr::group_by(pre_pt_root_id) %>%
  dplyr::mutate(count = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::select(pre_pt_root_id,pre_pt_position,pre_pt_supervoxel_id,count) %>%
  dplyr::distinct(pre_pt_root_id,.keep_all=TRUE) %>%
  dplyr::arrange(desc(count)) %>%
  dplyr::filter(count>=30,
                !as.character(pre_pt_root_id) %in% jo.search$root_id,
                !as.character(pre_pt_root_id) %in% vncs.id,
                !as.character(pre_pt_root_id) %in% proof)

# Write as .csv file
conn.top2$position = nat::xyzmatrix2str(conn.top2$pre_pt_position)
vnc.search <- conn.top2 %>%
  as.data.frame() %>%
  dplyr::mutate(root_id = as.character(pre_pt_root_id),
                supervoxel_id = as.character(pre_pt_supervoxel_id),
                status = "",
                notes="",
                annotator="") %>%
  dplyr::distinct(root_id,supervoxel_id,position,status,notes,annotator)
banc.tracing.save.path <- "tracing/tracing_issues/"
current_datetime <- Sys.time()
datetime_string <- format(current_datetime, "%Y-%m-%d")
readr::write_csv(vnc.search, file.path(banc.tracing.save.path,paste0(datetime_string,"_banc_vnc_sensory_search.csv")))

####################
### Correct ORNs ###
####################

# Search for ORN / TRN / HRN sensory neurons
bc <- banctable_query()

orns <- subset(
  bc,
  grepl("^ORN|^TRN|^HRN", cell_type)
)

# PNs: uniglomerular PNs + specified VP projection neurons
pns <- subset(
  bc,
  grepl("^uniglomerular_projection_neuron", cell_sub_class) |
    grepl(
      "VP4_vPN|VP1m\\+_lvPN|VP1d_il2PN|VP1l\\+_lvPN|VP3\\+_l2PN|VP5\\+_l2PN,VP5\\+VP2_l2PN|VP2\\+_adPN",
      cell_type
    )
)

pns.ids <- na.omit(unique(pns$root_id))
orns.id <- bancr::banc_rootid(na.omit(unique(orns$supervoxel_id)))
orns.id <- setdiff(orns.id, "0")

conns <- pbapply::pblapply(
  orns.id,
  function(x)
    tryCatch(
      banc_partners(x, partners = "output", synapse_table = "synapses_v2"),
      error = function(e) NULL
    )
)

conn <- do.call(plyr::rbind.fill, conns[!unlist(sapply(conns, is.null))])

conn.top <- conn %>%
  dplyr::mutate(
    pre_pt_root_id  = as.character(pre_pt_root_id),
    post_pt_root_id = as.character(post_pt_root_id)
  ) %>%
  dplyr::filter(post_pt_root_id %in% pns.ids) %>%
  dplyr::group_by(pre_pt_root_id, post_pt_root_id) %>%
  dplyr::mutate(count = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(pre_pt_root_id, post_pt_root_id, count) %>%
  dplyr::arrange(dplyr::desc(count)) %>%
  dplyr::filter(count > 1) %>%
  dplyr::left_join(
    bc %>%
      dplyr::distinct(root_id, .keep_all = TRUE) %>%
      dplyr::select(root_id, pre_cell_type = cell_type),
    by = c("pre_pt_root_id" = "root_id")
  ) %>%
  dplyr::left_join(
    bc %>%
      dplyr::distinct(root_id, .keep_all = TRUE) %>%
      dplyr::select(root_id, post_cell_type = cell_type),
    by = c("post_pt_root_id" = "root_id")
  )

