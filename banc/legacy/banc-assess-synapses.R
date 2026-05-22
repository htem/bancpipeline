###########################
### Assess synapse dump ###
###########################
source("banc/banc-startup.R")
library(vroom)
library(dtplyr)

########################################
### PROPORTION OF SYNAPSES PROOFREAD ###
########################################

# csv info
desired_columns <- c('id', 'size', 'pre_root_id', 'post_root_id', 'ctr_x', 'ctr_y', 'ctr_z')
col_types <- cols(
  id = col_character(),
  size = col_double(),
  pre_root_id = col_character(),
  post_root_id = col_character(),
  ctr_x = col_double(),
  ctr_y = col_double(),
  ctr_z = col_double(),
  .default = col_double()
)
column_names <- c('id', 'pre_x', 'pre_y', 'pre_z', 'post_x', 'post_y', 'post_z',
                  'ctr_x', 'ctr_y', 'ctr_z', 'size', 'pre_supervoxel_id',
                  'pre_root_id', 'post_supervoxel_id', 'post_root_id')

# # Get the synapse file remotely, you have it locally
version.path <- file.path(banc.save.path,"v821")
dir.create(version.path)
# system(sprintf("gsutil cp gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/v821/synapses_v2_human_readable.csv.gz %s ",
#                file.path(version.path,"synapses_v2_human_readable.csv.gz")))
# system(sprintf("gunzip %s",file.path(version.path,"synapses_v2_human_readable.csv.gz")))

# Get BANC meta
version.path <- file.path(banc.save.path,"v821")
banc.meta <- banctable_query() 
proof.svids <- banc.meta %>%
  dplyr::filter(proofread=="TRUE") %>%
  dplyr::pull(supervoxel_id)
roughly.svids <- banc.meta %>%
  dplyr::filter(roughly_proofread=="TRUE") %>%
  dplyr::pull(supervoxel_id)
bad.svids <- banc.meta %>%
  dplyr::filter(grepl("not_a_neuron|trachea|merge|glia",super_class)|grepl("NOT_A_NEURON|MERGE|TRACHEA|GLIA",status)) %>%
  dplyr::pull(supervoxel_id)
other.svids <- banc.meta %>%
  dplyr::filter(!supervoxel_id%in%bad.svids) %>%
  dplyr::pull(supervoxel_id)
svids <- unique(c(proof.svids,roughly.svids,other.svids))
identified.ids <- banc_rootid(svids,version="821")
bad.ids <-  banc_rootid(bad.svids,version="821")
neuron.ids <- c(bad.ids,identified.ids)
bancsynapses <- file.path(version.path,"synapses_v2_human_readable.csv")

# See what has been proofread
banc.syns.review <- vroom::vroom(file.path(version.path,"synapses_v2_human_readable.csv"),
                          col_names = column_names,
                          col_select = dplyr::all_of(desired_columns),
                          col_types = col_types,
                          skip = 1) %>%
  dplyr::rename(X=ctr_x,
                Y=ctr_y,
                Z=ctr_z) %>%
  dplyr::mutate(pre_status = case_when(
    pre_root_id %in% !!neuron.ids ~ "neuron",
    TRUE ~ "fragment"
  )) %>%
  dplyr::mutate(post_status = case_when(
    post_root_id %in% !!neuron.ids ~ "neuron",
    TRUE ~ "fragment"
  )) %>%
  tibble::as_tibble()

# Autapses
banc.syns.autapses <- banc.syns.review %>%
  dplyr::filter(pre_root_id==post_root_id) %>%
  nrow()
banc.syns.pre.good <- banc.syns.review %>%
  dplyr::filter(pre_status=="neuron") %>%
  nrow()
banc.syns.post.good <- banc.syns.review %>%
  dplyr::filter(post_status=="neuron") %>%
  nrow()

# Helper function for summary per root_id column
process_root_id <- function(df, root_col, status_col, label) {
  df2 <- df %>%
    dplyr::group_by(root_id = .data[[root_col]], pre_status = .data[[status_col]]) %>%
    dplyr::mutate(n_syn = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(pre_status, n_syn) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(n_syn)
  df_wide <- tidyr::pivot_wider(df2, names_from = pre_status, values_from = n, values_fill = 0)
  df_wide <- df_wide %>%
    dplyr::mutate(
      cum_neuron = if ("neuron" %in% names(.)) base::cumsum(neuron) else rep(0, dplyr::n()),
      cum_fragment = if ("fragment" %in% names(.)) base::cumsum(fragment) else rep(0, dplyr::n()),
      cum_total = cum_neuron + cum_fragment,
      pct_neuron = dplyr::if_else(cum_total > 0, cum_neuron / cum_total, NA_real_),
      root_id_col = label
    )
  return(df_wide)
}

# Calculate for both pre_root_id and post_root_id
dat_pre <- process_root_id(banc.syns.review, "pre_root_id", "pre_status", "pre_root_id")
dat_post <- process_root_id(banc.syns.review, "post_root_id", "post_status", "post_root_id")

# Combine for plotting
syn_comp <- dplyr::bind_rows(dat_pre, dat_post)

# Plot: two lines, one for each root_id type
g.ecdf <- ggplot2::ggplot(syn_comp, ggplot2::aes(x = n_syn, y = pct_neuron, color = root_id_col)) +
  ggplot2::geom_line(size = 1.2) +
  ggplot2::scale_x_log10() +
  ggplot2::annotation_logticks(sides = "b") +
  ggplot2::labs(
    x = "number of synapses per root_id (threshold)",
    y = "proportion of synapses in neurons (≤ threshold)",
    color = "",
    title = "cumulative share of synapses in neurons by root type"
  ) +
  scale_color_manual(values =  c(pre_root_id=paper.cols[["pre"]],
                                 post_root_id=paper.cols[["post"]])) +
  ggplot2::theme_minimal(base_size = 14) +
  theme(legend.position = "none")

# Save plot — both PNG (raster) and PDF (vector, for Illustrator linking).
# Same dir lives inside BANC-Project (figure_1/links/supplement) so a
# `git commit && push` from that repo ships the new figure.
print(g.ecdf)
banc.fig1.supp.path <- "/n/data1/hms/neurobio/wilson/banc/BANC-project/figures/figure_1/links/supplement"
dir.create(banc.fig1.supp.path, showWarnings = FALSE, recursive = TRUE)
for (.ext in c("png", "pdf")) {
  ggsave(plot = g.ecdf,
         filename = file.path(banc.fig1.supp.path,
                              sprintf("banc_synapse_proportion_on_cell.%s", .ext)),
         width = 6,
         height = 4,
         dpi = 300,
         bg = "transparent")
}

################################################
### Examine manually reviewed synapse sample ###
################################################

# Get regions
message("Examining false positive rate")
optic <- as.hxsurf(as.mesh3d(subset(banc_brain_neuropils.surf,"optic")))
midbrain <- as.hxsurf(as.mesh3d(subset(banc_brain_neuropils.surf,"midbrain")))
vnc <- banc_vnc_neuropil.surf
neck <- banc_neck_connective.surf

# Read results
data.new <- read_csv("data/tasks/2024-09-20_aelysia_synapse_sample_latest.csv")
data.orig <- read_csv("tracing/synapse_sample/2024-09-09_aelysia_synapse_sample_v2.csv")
data <- data.new %>%
  dplyr::select(id, `Coordinate 1`, Tags) %>%
  dplyr::left_join(data.orig %>%
                     dplyr::select(id, neuropil = Description, size, region), 
                   by = "id") %>%
  dplyr::filter(!is.na(Tags), !is.na(neuropil)) %>%
  dplyr::mutate(neuropil = gsub(" .*","",neuropil))  %>%
  dplyr::mutate(region = dplyr::case_when(
    neuropil %in% gsub(".*optic_","",banc_brain_neuropils.surf$RegionList[grepl("optic",banc_brain_neuropils.surf$RegionList)]) ~ "optic",
    neuropil %in% gsub(".*midbrain_","",banc_brain_neuropils.surf$RegionList[grepl("midbrain",banc_brain_neuropils.surf$RegionList)]) ~ "midbrain",
    neuropil %in% gsub(".*nerve_","",banc_vnc_nerves.surf$RegionList[grepl("nerve",banc_vnc_nerves.surf$RegionList)]) ~ "nerve",
    neuropil %in% gsub(".*vnc_","",banc_vnc_neuropils.surf$RegionList[grepl("vnc",banc_vnc_neuropils.surf$RegionList)]) ~ "vnc",
    TRUE ~ region
  )) %>%
  dplyr::mutate(region = factor(region, levels = c("optic","midbrain","vnc","nerve")))

# Assuming your data frame is named 'data'
plot_data1 <- data %>%
  # Create size bins with appropriate labels
  dplyr::mutate(size_bin = case_when(
    size <= 150 ~ cut(size, 
                      breaks = seq(0, 150, by = 5), 
                      labels = paste(seq(0, 145, by = 5), seq(4, 149, by = 5), sep="-"),
                      include.lowest = TRUE),
    size > 150 ~ "150+"
  )) %>%
  # Convert size_bin to factor to preserve order
  dplyr::mutate(size_bin = factor(size_bin, levels = c(paste(seq(0, 145, by = 5), seq(4, 149, by = 5), sep="-"), "150+"))) %>%
  # Convert Tags to a factor with levels True, False
  dplyr::mutate(Tags = factor(Tags, levels = c("True", "Ambiguous", "False"))) %>%
  # Group by size_bin and Tags
  dplyr::group_by(size_bin, Tags, region) %>%
  # Count occurrences
  dplyr::summarise(count = dplyr::n()) %>%
  # Calculate proportion within each size_bin
  dplyr::group_by(size_bin, region) %>%
  dplyr::mutate(proportion = count / sum(count))

# Create the plot
review.size <- ggplot(plot_data1, aes(x = size_bin, y = proportion, fill = Tags)) +
  facet_grid(region ~ ., scales = "free") +
  geom_col(position = "stack") +
  labs(x = "size bin", y = "proportion", title = "proportion of reviewed tags by size bin") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c(True = "#8FDA04", False = "#EE4244", Ambiguous  = "#1BB6AF"))

# Get the order of Description based on region
description_order <- data %>%
  dplyr::group_by(neuropil, region) %>%
  dplyr::summarise() %>%
  dplyr::arrange(region) %>%
  dplyr::pull(neuropil)
description_order <- unique(description_order)

# Assuming your data frame is named 'data'
plot_data2 <- data %>%
  # Convert neuropil to factor to preserve order
  dplyr::mutate(neuropil = factor(neuropil, description_order)) %>%
  # Convert Tags to a factor with levels True, False
  dplyr::mutate(Tags = factor(Tags, levels = c("True", "Ambiguous", "False"))) %>%
  # Group by neuropil and Tags
  dplyr::group_by(neuropil, Tags, region) %>%
  # Count occurrences
  dplyr::summarise(count = dplyr::n()) %>%
  # Calculate proportion within each neuropil
  dplyr::group_by(neuropil) %>%
  dplyr::mutate(proportion = count / sum(count))

# Create the plot
review.roi <- ggplot(plot_data2, aes(x = neuropil, y = count, fill = Tags)) +
  facet_grid(~ region, scales = "free") +
  geom_col(position = "stack") +
  labs(x = "neuropil", y = "proportion", title = "proportion of reviewed tags by neuropil") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c(True = "#8FDA04", False = "#EE4244", Ambiguous  = "#1BB6AF"))

# save
ggsave(plot = review.size,
       filename = file.path(banc.path,"images","synapses", "banc_reviewed_synapses_by_size.png"), 
       width = 24, height = 10, dpi = 300)

# save
ggsave(plot = review.roi,
       filename = file.path(banc.path,"images","synapses", "banc_reviewed_synapses_by_neuropil.png"), 
       width = 24, height = 10, dpi = 300)

# Write sample return for zetta
# readr::write_csv(x=data, file = file.path(banc.path,"data","synapses","2024-09-20_aelysia_synapse_sample_complete.csv"))

###########################
### Create sqlite table ###
###########################

# read ids
message("Making banc.ids")
banc.ids <- readr::read_csv(file=file.path(banc.meta.save.path,"banc_ids.csv"),
                            col_types = banc.col.types,
                            show_col_types = FALSE)
banc.root.ids <- unique(banc.ids$root_id)
banc.backbone.proofread <- banc_backbone_proofread()
proof.ids <- as.character(unique(banc.backbone.proofread$pt_root_id))

# lazy load the data
parallel <- FALSE
if(parallel){
  library(furrr)
  library(future)
  
  # For local parallelism:
  message("Available workers: ", parallel::detectCores() - 1)
  plan(multisession, workers = 10 - 1)
  
  # Read the data (this part stays sequential)
  message("parallel, reading: banc.syns")
  banc.syns <- vroom::vroom(bancsynapses,
                            col_names = column_names,
                            col_select = dplyr::all_of(desired_columns),
                            col_types = col_types,
                            skip = 1) %>%
    dplyr::rename(X=ctr_x,
                  Y=ctr_y,
                  Z=ctr_z) %>%
    dplyr::mutate(pre_status = case_when(
      pre_root_id %in% proof.ids ~ "proofread",
      pre_root_id %in% neuron.ids ~ "identified",
      TRUE ~ "fragment"
    )) %>%
    dplyr::mutate(post_status = case_when(
      post_root_id %in% proof.ids ~ "proofread",
      pre_root_id %in% neuron.ids ~ "identified",
      TRUE ~ "fragment"
    )) %>%
    tibble::as_tibble()
  
  # Define the chunk size
  message("making: banc.syns.np")
  chunk_size <- 100000 
  
  # Split the dataframe into chunks
  banc.syns <- banc.syns %>%
    dplyr::mutate(neuropil = NA,
                  region = NA,
                  side = NA,
                  chunk = ceiling(dplyr::row_number() / chunk_size)) %>%
    dplyr::group_by(chunk) %>%
    dplyr::group_split()
  
  # Process chunks in parallel - first function
  message("running: pointsinside_banc in parallel")
  banc.syns <- future_map(banc.syns, function(chunk) {
    message("Processing chunk for pointsinside_banc")
    pointsinside_banc(chunk)
  }, .progress = FALSE)
  
  # Process chunks in parallel - second function
  message("running: pointsnearby_banc in parallel")
  banc.syns <- future_map(banc.syns, function(chunk) {
    message("Processing chunk for pointsnearby_banc")
    pointsnearby_banc(chunk)
  }, .progress = FALSE)
  
  # Combine results
  banc.syns <- bind_rows(banc.syns) %>%
    dplyr::mutate(neuropil = gsub("ITO_optic_|ITO_midbrain_|COURT_vnc_","",neuropil)) %>%
    dplyr::arrange(region, neuropil)
  
  # Write to sql data base
  message("Saving synapse .sqlite file")
  dir.create(banc.connectivity.save.path, showWarnings = FALSE)
  con <- DBI::dbConnect(RSQLite::SQLite(),
                        file.path(banc.connectivity.save.path,"banc_data.sqlite"))
  DBI::dbWriteTable(con,
                    name = "synapses_250226_neuropils",
                    value = banc.syns,
                    overwrite = TRUE)
  dbDisconnect(con)
  message("Synapses saved")
  
  # Close parallel workers
  plan(sequential)
}else{
  message("non-parallel, reading: banc.syns")
  banc.syns <- vroom::vroom(bancsynapses,
                            col_names = column_names,
                            col_select = dplyr::all_of(desired_columns),
                            col_types = col_types,
                            skip = 1) %>%
    dplyr::rename(X=ctr_x,
                  Y=ctr_y,
                  Z=ctr_z) %>%
    dplyr::mutate(pre_status = case_when(
      pre_root_id %in% proof.ids ~ "identified",
      pre_root_id %in% neuron.ids ~ "identified",
      TRUE ~ "fragment"
    )) %>%
    dplyr::mutate(post_status = case_when(
      post_root_id %in% proof.ids ~ "identified",
      pre_root_id %in% neuron.ids ~ "identified",
      TRUE ~ "fragment"
    )) %>%
    tibble::as_tibble()
  
  # Get calculated BANC synapses
  if(file.exists(file.path(banc.connectivity.save.path,"banc_data.sqlite"))){
    con <- DBI::dbConnect(RSQLite::SQLite(),
                          file.path(banc.connectivity.save.path,"banc_data.sqlite"))
    banc.syns.nps.done <- dplyr::tbl(con, "synapses_250226_neuropils") %>%
      dplyr::select(id) %>%
      dplyr::collect()
    dbDisconnect(con) 
    done.ids <- unique(banc.syns.nps.done$id)
    message("calculated synapes: ", length(done.ids))
    banc.syns <- banc.syns %>%
      dplyr::filter(!id %in% done.ids)
    message("uncalculated synapes: ", nrow(banc.syns))
  }
  
  # Define the chunk size
  if(nrow(banc.syns)){
    message("making: banc.syns.np")
    chunk_size <- 1000000
    
    # Split the dataframe into a list of chunks
    banc.syns.np <- banc.syns %>%
      dplyr::mutate(neuropil = NA,
                    region = NA,
                    side = NA,
                    chunk = ceiling(dplyr::row_number() / chunk_size)) %>%
      dplyr::group_by(chunk) %>%
      dplyr::group_split()
    
    # Process each chunk
    banc.syns.np <- purrr::map(banc.syns.np, function(chunk) {
      message("running: pointsinside_banc")
      data <-pointsinside_banc(chunk)
      message("running: pointsnearby_banc")
      data <- pointsnearby_banc(data)
      data <- data %>%
        dplyr::mutate(neuropil = gsub("ITO_optic_|ITO_midbrain_|COURT_vnc_","",neuropil)) %>%
        dplyr::arrange(region, neuropil)
      con <- DBI::dbConnect(RSQLite::SQLite(),
                            file.path(banc.connectivity.save.path,"banc_data.sqlite"))
      DBI::dbWriteTable(con,
                        name = "synapses_250226_neuropils",
                        value = data,
                        overwrite = FALSE,
                        append = TRUE)
      dbDisconnect(con)
      message("synapse batch written")
      NULL
    }) 
  }
}

#########################################
### Send neuropil inclusion to bucket ###
#########################################

# Get BANC synapses
con <- DBI::dbConnect(RSQLite::SQLite(),
                      file.path(banc.connectivity.save.path,"banc_data.sqlite"))
banc.syns.nps <- dplyr::tbl(con, "synapses_250226") %>%
  dplyr::filter() %>%
  dplyr::collect()
dbDisconnect(con)
gc()

# # Large fragments
# large.post.ids <- banc.syns %>%
#   dplyr::group_by(post_root_id) %>%
#   dplyr::mutate(n = dplyr::n()) %>%
#   dplyr::ungroup() %>%
#   dplyr::filter(n > 1000) %>%
#   dplyr::pull(post_root_id)
# large.pre.ids <- banc.syns %>%
#   dplyr::group_by(pre_root_id) %>%
#   dplyr::mutate(n = dplyr::n()) %>%
#   dplyr::ungroup() %>%
#   dplyr::filter(n > 1000) %>%
#   dplyr::pull(pre_root_id)

# Get BANC synapses
banc.syns <- vroom::vroom(bancsynapses,
                          col_names = column_names,
                          col_select = dplyr::all_of(desired_columns),
                          col_types = col_types,
                          skip = 0) %>%
  dplyr::rename(X=ctr_x,
                Y=ctr_y,
                Z=ctr_z) %>%
  dplyr::filter(pre_root_id != post_root_id) %>%
  dplyr::distinct(id, .keep_all  = TRUE) %>%
  dplyr::mutate(pre_status = case_when(
    pre_root_id %in% identified.ids ~ "identified",
    TRUE ~ "fragment"
  )) %>%
  dplyr::mutate(post_status = case_when(
    post_root_id %in% identified.ids ~ "identified",
    TRUE ~ "fragment"
  )) %>%
  tibble::as_tibble()

# Clean up
banc.syns.nps.cleaned <- banc.syns.nps %>%
  dplyr::distinct(id, .keep_all= TRUE) %>%
  dplyr::filter(id %in% !!banc.syns$id) %>%
  dplyr::mutate(neuropil_full = neuropil) %>%
  dplyr::mutate(neuropil = gsub("MANC_.*","",neuropil),
                neuropil = gsub("outside_|vnc_","",neuropil),
                neuropil = gsub("\\,.*|_R$|_L$","",neuropil),
                neuropil = gsub("vncN","LN",neuropil),
                neuropil = gsub("outside|vnc","",neuropil),
                region = gsub("outside_|vnc_","",region)) %>%
  dplyr::mutate(neuropil = dplyr::case_when(
    neuropil == "H" ~ "LH",
    neuropil=="O" ~ "LO",
    neuropil=="OP" ~ "LOP",
    is.na(neuropil)|neuropil=='' ~ gsub("outside_|\\,.*|_R$|_L$|MANC_","",neuropil_full),
    TRUE ~ neuropil
  )) %>%
  dplyr::mutate(region = dplyr::case_when(
    grepl("outside",region_complex)|grepl("outside",neuropil_full) ~ "outside",
    neuropil %in% c("LO","LOP","AME","ME") ~ "optic_lobe",
    region =='sez' ~ 'central_brain',
    region == "optic_lobes" ~ "optic_lobe",
    region == "vnc" ~ "ventral_nerve_cord",
    TRUE ~ region
  ))

# Add missing synapse
missed <- data.frame(id = "205992435",
                     X = 602400, 
                     Y = 146192,
                     Z = 270720,
                     side = "left",
                     region = "central_brain",
                     neuropil = "SLP")
banc.syns.nps.cleaned <- plyr::rbind.fill(banc.syns.nps.cleaned,
                                          missed)

# Save
# con <- DBI::dbConnect(RSQLite::SQLite(),
#                       file.path(banc.connectivity.save.path,"banc_data.sqlite"))
# DBI::dbWriteTable(con,
#                   name = "synapses_250226",
#                   value = banc.syns.nps.cleaned,
#                   overwrite = TRUE,
#                   append = FALSE)
# dbDisconnect(con)

# Calculate some key numbers
summary <- banc.syns %>%
  dplyr::filter(!post_root_id %in% !!bad.ids,
                !pre_root_id %in% !!bad.ids) %>%
  dplyr::count(pre_status, post_status) %>%        
  dplyr::mutate(prop = n / sum(n),
                prop = round(prop,4))
print(summary)
write_csv(summary, file = file.path(banc.meta.save.path,"synapses_250226_gross_capture_rates.csv"))

# Not assigned
round(nrow(banc.syns.nps.cleaned)/nrow(banc.syns),6)

# By neuropil
banc.syns.np.full <- banc.syns %>%
  dplyr::left_join(banc.syns.nps.cleaned %>%
                     dplyr::select(id, side, region, neuropil), by = "id") %>%
  dplyr::mutate(side = ifelse(is.na(side),"uncalculated",side),
                region = ifelse(is.na(side),"uncalculated",region),
                neuropil = ifelse(is.na(side),"uncalculated",neuropil))

summary.no.optic <- banc.syns.np.full %>%
  dplyr::filter(!grepl("optic",region)) %>%
  dplyr::filter(!post_root_id %in% !!bad.ids,
                !pre_root_id %in% !!bad.ids) %>%
  dplyr::ungroup() %>%
  dplyr::count(pre_status, post_status) %>%        
  dplyr::mutate(prop = n / sum(n),
                prop = round(prop,4))
print(as.data.frame(summary.no.optic))

summary.inout <- banc.syns.np.full %>%
  dplyr::filter(!post_root_id %in% !!bad.ids,
                !pre_root_id %in% !!bad.ids) %>%
  dplyr::mutate(in_mesh = ifelse(region=='outside',"outside","inside")) %>%
  dplyr::group_by(in_mesh) %>%
  dplyr::count(pre_status, post_status, in_mesh) %>%        
  dplyr::mutate(prop = n / sum(n),
                prop = round(prop,4))
print(summary.inout)
write_csv(summary.inout, file = file.path(banc.meta.save.path,"synapses_250226_inout_of_mesh_capture_rates.csv"))

summary.reg.side <- banc.syns.np.full %>%
  dplyr::filter(!post_root_id %in% !!bad.ids,
                !pre_root_id %in% !!bad.ids) %>%
  dplyr::group_by(region, side) %>%
  dplyr::count(pre_status, post_status, side, region) %>%        
  dplyr::mutate(prop = n / sum(n),
                prop = round(prop,4))
print(as.data.frame(summary.reg.side))
write_csv(summary.reg.side, file = file.path(banc.meta.save.path,"synapses_250226_region_capture_rates.csv"))
summary.np <- banc.syns.np.full %>%
  dplyr::group_by(region, side, neuropil) %>%
  dplyr::count(pre_status, post_status, side, region, neuropil) %>%        
  dplyr::mutate(prop = n / sum(n),
                prop = round(prop,4))
print(summary.np)
write_csv(summary.np, file = file.path(banc.meta.save.path,"synapses_250226_neuropil_capture_rates.csv"))

# Send
readr::write_csv(banc.syns.nps.cleaned %>%
                   dplyr::distinct(id, .keep_all = TRUE) %>%
                   dplyr::select(id, X, Y, Z, side, region, neuropil), 
                 file =  file.path(banc.connectivity.save.path,"banc_synapses_to_neuropils_v2.csv"))
system("gsutil cp /n/data1/hms/neurobio/wilson/banc/connectivity/banc_synapses_to_neuropils_v2.csv gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/synapses")
message("Synapses sent to google bucket")

# Build large neural fragment list for proofreading
banc.syns.np.full.post.no.pre <- banc.syns.np.full %>%
  dplyr::filter(post_status == "identified",
                pre_status != "identified") %>%
  dplyr::group_by(pre_root_id) %>%
  dplyr::mutate(output_connections = dplyr::n()) %>%
  dplyr::distinct(pre_root_id, .keep_all = TRUE) %>%
  dplyr::select(root_id = pre_root_id,
                x=X,
                y=Y,
                z=Z,
                region,
                neuropil,
                output_connections) %>%
  dplyr::arrange(dplyr::desc(output_connections)) %>%
  dplyr::filter(output_connections>50) %>%
  dplyr::ungroup()
readr::write_csv(banc.syns.np.full.post.no.pre, 
                 file =  file.path(banc.connectivity.save.path,"banc_post_proofread_pre_not_proofread.csv"))
system("gsutil cp /n/data1/hms/neurobio/wilson/banc/connectivity/banc_post_proofread_pre_not_proofread.csv gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/synapses")
message("Synapses sent to google bucket")

# Build large neural fragment list for proofreading
banc.syns.np.full.pre.no.post <- banc.syns.np.full %>%
  dplyr::filter(pre_status == "identified",
                post_status != "identified") %>%
  dplyr::group_by(post_root_id) %>%
  dplyr::mutate(input_connections = dplyr::n()) %>%
  dplyr::distinct(post_root_id, .keep_all = TRUE) %>%
  dplyr::select(root_id = post_root_id,
                x=X,
                y=Y,
                z=Z,
                region,
                neuropil,
                input_connections) %>%
  dplyr::arrange(dplyr::desc(input_connections)) %>%
  dplyr::filter(input_connections>100) %>%
  dplyr::ungroup()
readr::write_csv(banc.syns.np.full.post.no.pre, 
                 file =  file.path(banc.connectivity.save.path,"banc_pre_proofread_post_not_proofread.csv"))
system("gsutil cp /n/data1/hms/neurobio/wilson/banc/connectivity/banc_post_proofread_pre_not_proofread.csv gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/synapses")
message("Synapses sent to google bucket")

####################
### Visual check ###
####################

np.dir <- "inst/images/synapse_neuropils"
dir.create(np.dir)
neuropils <- unique(banc.syns.nps.cleaned$neuropil)
for(np in neuropils){
  message("working on ", np)
  banc.syns.nps.single.np <- subset(banc.syns.nps.cleaned, neuropil==np)
  banc.syns.nps.single.np <- banc.syns.nps.single.np[sample(1:nrow(banc.syns.nps.single.np),size = min(10000,nrow(banc.syns.nps.single.np))),]
  cols <- unname(paper.cols[banc.syns.nps.single.np$side])
  reg <- names(rev(sort(table(banc.syns.nps.single.np$region))))
  if(all(grepl("brain|optic",reg))){
    view <- "front"
  }else if(all(grepl("vnc|ventral_nerve_cord",reg))){
    view <- "vnc"
  }else{
    view <- "main"
  }
  g <- banc_ggneuron(banc.syns.nps.single.np,banc_neuropil.surf,view = view, cols1=cols)
  ggsave(plot = g,
         filename = paste0("inst/images/synapse_neuropils/",gsub(" |,","",np),".png"), width = 8, height = 8, dpi = 300)
}

np.dir <- "inst/images/synapse_regions"
dir.create(np.dir)
regions <- unique(banc.syns.nps.cleaned$region)
for(np in regions){
  message("working on ", np)
  banc.syns.nps.single.np <- subset(banc.syns.nps.cleaned, region==np)
  banc.syns.nps.single.np <- banc.syns.nps.single.np[sample(1:nrow(banc.syns.nps.single.np),size = min(10000,nrow(banc.syns.nps.single.np))),]
  cols <- paper.cols[banc.syns.nps.single.np$side]
  reg <- unname(names(rev(sort(table(banc.syns.nps.single.np$region)))))
  if(all(grepl("brain|optic",reg))){
    view <- "front"
  }else if(all(grepl("vnc|ventral_nerve_cord",reg))){
    view <- "vnc"
  }else{
    view <- "main"
  }
  g <- banc_ggneuron(xyzmatrix(banc.syns.nps.single.np),banc_neuropil.surf,view = view, cols1=cols)
  ggsave(plot = g,
         filename = paste0("inst/images/synapse_regions/",gsub(" |,","",np),".png"), width = 8, height = 8, dpi = 300)
}

################################
### Plot synapse information ###
################################

# Get BANC synapses
con <- DBI::dbConnect(RSQLite::SQLite(),
                      file.path(banc.connectivity.save.path,"banc_data.sqlite"))
banc.syns <- dplyr::tbl(con, "synapses_250226") %>%
  dplyr::collect()
dbDisconnect(con)

# # summarise the data
# message("making: banc_synapses_on_proofread_neurons.png")
# summarized_data0 <- banc.syns %>%
#   dplyr::group_by(pre_status, post_status) %>% # post_status
#   dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
#   dplyr::mutate(percentage = count / sum(count) * 100,
#                 label = ifelse(percentage > 5, 
#                                paste0(round(percentage, 1), "% \n(", count, ")"), 
#                                "")) %>%
#   dplyr::ungroup()
# total_count <- sum(summarized_data0$count)
# 
# # create the plot using the summarized data
# g0 <- ggplot(summarized_data0, aes(x = pre_status, y = percentage, fill = post_status)) +
#   geom_bar(stat = "identity", position = "stack") +
#   geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 2, col = "black") +
#   ggplot2::scale_fill_manual(values = c("proofread" = hemibrainr:::hemibrain_bright_colors[["green"]], 
#                                         "fragment" = hemibrainr:::hemibrain_bright_colors[["cerise"]])) +
#   labs(title = paste0("distribution of synaptic status combinations (total: ", total_count, ")"),
#        x = "pre-synaptic status",
#        y = "percentage of total synapses",
#        fill = "post-synaptic status") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         legend.position = "right",
#         plot.title = element_text(hjust = 0.5)) +
#   coord_flip() +
#   scale_y_continuous(labels = scales::percent_format(scale = 1))
# 
# # save
# ggsave(plot = g0,
#        filename = "inst/images/banc_synapses_on_proofread_neurons.png", width = 8, height = 4, dpi = 300)

# Data processing with dplyr
message("making: banc_synapse_distribution_of_size_by_status.png")
processed_data1 <- banc.syns %>%
  dplyr::mutate(size_bin = cut(size, 
                               breaks = seq(0, 150, by = 5),
                               labels = as.character(seq(1, 150, by = 5)),
                               include.lowest = TRUE)) %>%
  dplyr::group_by(pre_status, post_status, size_bin) %>%
  dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(percentage = count / sum(count, na.rm = TRUE) * 100) %>%
  tidyr::complete(pre_status, post_status, size_bin, fill = list(count = 0,
                                                                 percentage = 0))
total_count <- sum(processed_data1$count)

# Create the plot
g1 <- ggplot(processed_data1, aes(x = size_bin, y = count, fill = size_bin)) +
  ggplot2::geom_bar(stat = "identity", position = "stack") +
  ggplot2::facet_grid(pre_status ~ post_status, scales = "free_y") +
  ggplot2::scale_fill_viridis_d() +
  ggplot2::labs(title = paste0("synapse size distribution by pre- and post-synaptic status (total: ", total_count, ")"),
                x = "synapse size", 
                y = "count",
                fill = "") +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1),
                 legend.position = "bottom",
                 strip.background = element_rect(fill = "lightgrey"),
                 strip.text = element_text(face = "bold")) +
  ggplot2::guides(fill = guide_legend(nrow = 3, byrow = TRUE))  +
  ggplot2::guides(fill = "none",
                  color = "none")

# save
ggsave(plot = g1,
       filename = "inst/images/banc_synapse_distribution_of_size_by_status.png", 
       width = 12, height = 12, dpi = 300)

# process the data
message("making: synapse_distribution_by_neuron_synaptic_count.png")
total_count <- nrow(banc.syns)
processed_data_output <- banc.syns %>%
  dplyr::ungroup() %>%
  dplyr::group_by(pre_root_id) %>%
  dplyr::mutate(pre_synapse_count = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(synapse_size_bin = cut(pre_synapse_count, 
                                       breaks = c(0, 5, 10, 50, 100, 500, 1000, 5000, Inf),
                                       labels = c("1-5", "6-10", "11-50", "51-100", "101-500", "501-1000", "1001-5000", "5000+"),
                                       include.lowest = TRUE)) %>%
  dplyr::select(pre_root_id, pre_synapse_count, synapse_size_bin) %>%
  dplyr::group_by(synapse_size_bin) %>%
  dplyr::summarise(connection_count = dplyr::n(),
                   proportion = round(connection_count / total_count,2),
                   .groups = "drop") %>%
  dplyr::mutate(prepost = "output")
processed_data_input <- banc.syns %>%
  dplyr::ungroup() %>%
  dplyr::group_by(post_root_id) %>%
  dplyr::mutate(post_synapse_count = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(synapse_size_bin = cut(post_synapse_count, 
                                       breaks = c(0, 5, 10, 50, 100, 500, 1000, 5000, Inf),
                                       labels = c("1-5", "6-10", "11-50", "51-100", "101-500", "501-1000", "1001-5000", "5000+"),
                                       include.lowest = TRUE)) %>%
  dplyr::select(post_root_id, post_synapse_count, synapse_size_bin) %>%
  dplyr::group_by(synapse_size_bin) %>%
  dplyr::summarise(connection_count = dplyr::n(),
                   proportion = round(connection_count / total_count,2),
                   .groups = "drop") %>%
  dplyr::mutate(prepost = "input") 
plot_data2 <- rbind(processed_data_output,processed_data_input)

# Step 3: Create the plot
g2 <- ggplot(plot_data2, aes(x = synapse_size_bin, y = proportion, fill = prepost)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(data = plot_data2 %>% filter(proportion >= 0.05),  # Only label bars > 5%
            aes(label = scales::percent(proportion, accuracy = 0.1)), 
            position = position_stack(vjust = 0.5)) +
  facet_wrap(~prepost,nrow=2) +
  labs(title = "proportion of connections by neuron type and synapse count",
       x = "number of synapses per neuron",
       y = "proportion of connections",
       fill = "direction") +
  ggplot2::scale_fill_manual(values = c("input" = hemibrainr:::hemibrain_bright_colors[["cyan"]], 
                                        "output" = hemibrainr:::hemibrain_bright_colors[["cerise"]])) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 12)) +
  scale_y_continuous(labels = scales::percent) +
  ggplot2::guides(fill = "none",
                  color = "none")

# Save the plot
ggsave(filename="inst/images/synapse_distribution_by_neuron_synaptic_count.png", 
       plot = g2, 
       width = 12, height = 8, dpi = 300)

########################
### Assign neuropils ###
########################

# data processing with dplyr
message("making: banc_synapses_sizes_by_neuropil.png")
processed_data3 <- banc.syns %>%
  dplyr::filter(!is.na(region), !is.na(neuropil)) %>%
  dplyr::mutate(neuropil = ifelse(is.na(neuropil),"outside",neuropil)) %>%
  dplyr::mutate(region = ifelse(neuropil=="outside",neuropil,region)) %>%
  dplyr::mutate(
    size_bin = cut(size, 
                   breaks = c(0, 10, 20, 30, 40, 50, 100, 200, Inf),
                   labels = c("1-10", "11-20", "21-30", "31-40", "41-50", "51-100","101-200", "201+"),
                   include.lowest = TRUE)) %>%
  dplyr::group_by(region, neuropil, side, size_bin) %>%
  dplyr::summarise(count = dplyr::n(), 
                   .groups = "drop") %>%
  tidyr::complete(region, neuropil, size_bin, fill = list(count = 0)) %>%
  dplyr::mutate(region = as.character(region),
                neuropil = as.character(neuropil)) %>%
  dplyr::filter(count > 0) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(neuropil = grep(pattern = "^ITO|^COURT|^outside", 
                                x = unlist(strsplit(neuropil,split=",")),value=T)[1]) %>%
  dplyr::mutate(neuropil = gsub("ITO_midbrain_|ITO_optic_|COURT_vnc_|_R$|_L$","", as.character(neuropil))) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(region, neuropil, side, size_bin) %>%
  dplyr::summarise(count = sum(count, na.rm = TRUE), 
                   .groups = "drop") %>%
  dplyr::group_by(region) %>%
  dplyr::mutate(neuropil = factor(neuropil, levels = sort(unique(neuropil))))

# create the plot
g3 <- ggplot(processed_data3, aes(x = neuropil, y = count, group = side, fill = size_bin)) +
  ggplot2::geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  ggplot2::scale_fill_viridis_d() +
  ggplot2::labs(title = "synapse size distribution by neuropil and region",
                x = "synapse size",
                y = "count",
                fill = "neuropil") +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1),
                 legend.position = "bottom",
                 strip.background = element_rect(fill = "lightgrey"),
                 strip.text = element_text(face = "bold")) +
  ggplot2::guides(fill = guide_legend(nrow = 3, byrow = TRUE))

# save
ggsave(plot = g3,
       filename = "inst/images/banc_synapses_sizes_by_neuropil.png", 
       width = 24, height = 10, dpi = 300)

# data processing with dplyr
message("making: banc_synapses_status_by_neuropil.png")
processed_data4 <- banc.syns %>%
  dplyr::filter(!is.na(region), !is.na(neuropil)) %>%
  dplyr::mutate(neuropil = ifelse(is.na(neuropil),"outside",neuropil)) %>%
  dplyr::mutate(region = ifelse(neuropil=="outside",neuropil,region)) %>%
  dplyr::group_by(region, neuropil, side, pre_status, post_status) %>%
  dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(count > 0) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(neuropil = grep(pattern = "^ITO|^COURT|^outside", 
                                x = unlist(strsplit(neuropil,split=",")),value=T)[1]) %>%
  dplyr::mutate(neuropil = gsub("ITO_midbrain_|ITO_optic_|COURT_vnc_|_R$|_L$","", as.character(neuropil))) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(region, neuropil, side, pre_status, post_status) %>%
  dplyr::summarise(count = sum(count, na.rm = TRUE), 
                   .groups = "drop") %>%
  dplyr::ungroup() %>%
  dplyr::mutate(percentage = count / sum(count, na.rm=TRUE) * 100) %>%
  dplyr::mutate(status_combination = paste(pre_status, post_status, sep = " - ")) %>%
  dplyr::group_by(region) %>%
  dplyr::mutate(neuropil = factor(neuropil, levels = sort(unique(neuropil))))

# create the plot
g4 <- ggplot(processed_data4, aes(x = neuropil, y = percentage, group = side, fill = status_combination)) +
  ggplot2::geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  #ggplot2::facet_wrap(~ region, scales = "free_x", ncol = 1) +
  ggplot2::scale_fill_manual(values = c("proofread - proofread" = hemibrainr:::hemibrain_bright_colors[["green"]], 
                                        "proofread - fragment" = hemibrainr:::hemibrain_bright_colors[["yellow"]],
                                        "fragment - proofread" = hemibrainr:::hemibrain_bright_colors[["paleorange"]],
                                        "fragment - fragment" = hemibrainr:::hemibrain_bright_colors[["cerise"]])) +
  ggplot2::labs(title = "synapse count by neuropil, region, and status combination",
                x = "neuropil",
                y = "count",
                fill = "status combination (pre -> post)") +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1),
                 legend.position = "bottom",
                 strip.background = element_rect(fill = "lightgrey"),
                 strip.text = element_text(face = "bold")) +
  ggplot2::guides(fill = guide_legend(nrow = 2, byrow = TRUE))  +
  scale_y_continuous(labels = scales::percent_format(scale = 1))

# save
ggsave(plot = g4,
       filename = "inst/images/banc_synapses_status_by_neuropil.png", 
       width = 24, height = 10, dpi = 300)

##########################
### Get synapse sample ###
##########################

# Sample synapses
message("writing: banc_synapse_sample_v4.csv")

# First, sample from groups
grouped_sample <- banc.syns.np.full %>%
  dplyr::filter(size <= 5) %>%
  dplyr::distinct(id, X, Y, .keep_all = TRUE) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(neuropil = grep(pattern = "^ITO|^COURT|^outside", 
                                x = unlist(strsplit(neuropil,split=",")),value=T)[1]) %>%
  dplyr::mutate(neuropil = gsub("ITO_midbrain_|ITO_optic_|COURT_vnc_|_R$|_L$","", as.character(neuropil))) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(neuropil, region, side, size) %>%
  dplyr::slice_sample(n = 25, replace = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!neuropil%in%c("NA","","outside"))

# Then, sample from "outside" category
outside_sample <- banc.syns %>%
  dplyr::filter(region == "outside", neuropil == "outside") %>%
  dplyr::distinct(id, X, Y, .keep_all = TRUE) %>%
  dplyr::anti_join(grouped_sample, by = "id") %>%  
  dplyr::slice_sample(n = 1000, replace = FALSE)

# Combine the samples
synapse.sample <- bind_rows(grouped_sample, outside_sample) %>%
  dplyr::collect()

# Check the total number of rows
total_rows <- nrow(synapse.sample)
print(paste("Total rows in the final sample:", total_rows))

# Check the number of "outside" rows
outside_rows <- sum(synapse.sample$region == "outside" & synapse.sample$neuropil == "outside")
print(paste("Number of 'outside' rows:", outside_rows))

# Make a neuroglancer annotation file
synapse.sample <- as.data.frame(synapse.sample, stringsAsFactors = FALSE)
synapse.sample$`Coordinate 1` = apply(nat::xyzmatrix(synapse.sample),1,function(x) hemibrainr:::paste_coords(x))
banc.scan <- data.frame(`Coordinate 1` = synapse.sample$`Coordinate 1`,
                        `Coordinate 2` = "",
                        `Ellipsoid Dimensions` = "",
                        tags = "",
                        Description = hemibrainr:::nullToNA(synapse.sample$neuropil),
                        `Segment IDs` = "",
                        `Parent ID` = "",
                        Type = "Point",
                        ID = "",
                        id = hemibrainr:::nullToNA(synapse.sample$id),
                        size = hemibrainr:::nullToNA(synapse.sample$size),
                        region = hemibrainr:::nullToNA(synapse.sample$region))
colnames(banc.scan) = gsub("\\."," ",colnames(banc.scan))
banc.scan$`Coordinate 1` = as.character(banc.scan$`Coordinate 1`)
datetime_string <- format(Sys.time(), "%Y-%m-%d")
file <- file.path("tracing",paste0(datetime_string,"_banc_synapse_sample_v4.csv"))
readr::write_excel_csv(banc.scan, file = file)

# Write to database
con <- DBI::dbConnect(RSQLite::SQLite(),
                      file.path(banc.connectivity.save.path,"banc_data.sqlite"))
DBI::dbWriteTable(con,
                  name = "synapse_sample_v4",
                  value = synapse.sample,
                  overwrite = TRUE)

# Close the database connection
dbDisconnect(con)

system(sprintf("gsutil cp %s gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/synapses",file))
message("Synapses sent to google bucket")

