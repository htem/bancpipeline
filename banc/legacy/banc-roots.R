########################################
### Calculate rootpoints for neurons ###
########################################
source("banc/banc-startup.R")
library(R.utils)
message("##### Working out BANC roots #####")
numCores <- 1

# Direct us to the BANC dataset
bancr::choose_banc()
bc <- banctable_query("SELECT _id, status, proofread, root_id, supervoxel_id, position, l2_nodes, region, root_region, root_position, root_position_nm, nucleus_id, nucleus_position, nucleus_position_nm from banc_meta") %>%
  dplyr::filter(!grepl("DELETE|NOT_A_NEURON|DEBRIS|GLIA|TRACHEA|MERGE",status))

# Split metrics folder
split.metrics.folder <- file.path(banc.l2split.save.path,"metrics")

# Read saved meta file
banc.metrics.old <- readr::read_csv(file=file.path(banc.save.path,"banc_root_positions.csv"),
                                    col_types = banc.col.types, 
                                    show_col_types = FALSE)

# Add in roots from seatable
banc.metrics.updated <- banc.metrics.old %>%
  dplyr::distinct(root_id, .keep_all = TRUE) %>%
  dplyr::left_join(bc %>% dplyr::distinct(root_id, 
                                          .keep_all = TRUE), 
                   by = "root_id",
                   relationship = "many-to-many") %>%
  dplyr::mutate(root_position = dplyr::case_when(
    is.na(root_position.x) ~ root_position.y,
    is.na(root_position.y) ~ root_position.x,
    !is.na(root_position.y) ~ root_position.y,
    !is.na(root_position.x) ~ root_position.x,
    TRUE ~ root_position.x
  ),
  root_position_nm = dplyr::case_when(
    is.na(root_position_nm.x) ~ root_position_nm.y,
    is.na(root_position_nm.y) ~ root_position_nm.x,
    !is.na(root_position_nm.y) ~ root_position_nm.y,
    !is.na(root_position_nm.x) ~ root_position_nm.x,
    TRUE ~ root_position_nm.x
  ),
  region = dplyr::case_when(
    is.na(region.x) ~ region.y,
    is.na(region.y) ~ region.x,
    !is.na(region.y) ~ region.y,
    !is.na(region.x) ~ region.x,
    TRUE ~ region.x
  ),
  root_region = dplyr::case_when(
    is.na(root_region.x) ~ root_region.y,
    is.na(root_region.y) ~ root_region.x,
    !is.na(root_region.y) ~ root_region.y,
    !is.na(root_region.x) ~ root_region.x,
    TRUE ~ root_region.x
  ),
  region = dplyr::case_when(
    region=='vnc' ~ "ventral_nerve_cord",
    region=='midbrain' ~ "central_brain",
    region=='optic' ~ "optic_lobe",
    TRUE ~ region
  )) %>%
  dplyr::mutate(region = dplyr::case_when(
    !is.na(region)|region!="" ~ region,
    grepl("midbrain",root_region) ~ "central_brain",
    grepl("optic",root_region) ~ "optic_lobe",
    grepl("vnc",root_region) ~ "ventral_nerve_cord",
    TRUE ~ region
  )) %>%
  dplyr::select(-ends_with(".x"), -ends_with(".y"))  %>%
  dplyr::filter(!is.na(region)&!is.na(root_id)) %>%
  dplyr::distinct(root_id, .keep_all = TRUE)

# Enforce agreement
banc.metrics.updated$root_position_nm <- apply(banc_raw2nm(xyzmatrix(banc.metrics.updated$root_position)), 1, bancr:::paste_coords)
banc.metrics.updated$root_position_nm <- gsub("\\(|\\)", "", banc.metrics.updated$root_position_nm)
banc.metrics.updated$root_position <- apply(banc_nm2raw(xyzmatrix(banc.metrics.updated$root_position_nm)), 1, bancr:::paste_coords)
banc.metrics.updated$root_position <- gsub("\\(|\\)", "", banc.metrics.updated$root_position)
banc.has.root.positions <- unique(subset(banc.metrics.updated, !is.na(banc.metrics.updated$root_position_nm))$root_id)

# Update what is new
too_many_digits <- function(x){sapply(strsplit(x, ","), function(parts) {
  any(nchar(trimws(parts)) > 7)
})}
banc.ids.todo.all <- bc  %>%
  dplyr::filter(is.na(region)|is.na(root_id)|too_many_digits(root_position_nm)) 
banc.root.ids <- unique(c(banc.ids.todo.all$root_id))
banc.metrics.updated <- banc.metrics.updated %>%
  dplyr::filter(!root_id %in% banc.root.ids)
message("Working on roots and regions for: ", length(banc.root.ids), " neurons")

# Batch for processing 
multiplier <- 100
todo <- 1:nrow(banc.ids.todo.all)
upper <- ifelse((numCores*multiplier)<length(todo),numCores*multiplier,length(todo))
batches <- split(todo, round(seq(from = 1, to = upper, length.out = length(todo))))
batches <- sample(batches)
for(batch in batches){
  try({
    banc.ids.todo <- banc.ids.todo.all[batch,]
    if(nrow(banc.ids.todo)){
      
      # Register cores
      cl <- setup_parallel()
      
      # Read L2 skeletons
      nfiles <- paste0(file.path(banc.l2swc.save.path,unique(banc.ids.todo$root_id)),".swc")
      nfiles <- nfiles[file.exists(nfiles)]
      nfiles.ids <- gsub(".swc","",basename(nfiles))
      
      # Read synapses
      sfiles <- paste0(file.path(banc.synapses.save.path,unique(banc.ids.todo$root_id)),".csv")
      sfiles <- sfiles[file.exists(sfiles)]
      sfiles.ids <- gsub(".csv","",basename(sfiles))
      
      # Get the soma/root points
      soma.positions <- banc.ids.todo %>%
        dplyr::mutate(root_position_nm = ifelse(is.na(root_position_nm),nucleus_position_nm,root_position_nm)) %>%
        dplyr::distinct(root_id, supervoxel_id, nucleus_id, root_position_nm) %>%
        dplyr::filter(!is.na(root_position_nm))
      if(nrow(soma.positions)){
        soma.positions[,c("X","Y","Z")] <- nat::xyzmatrix(soma.positions$root_position_nm)
        # soma.positions$root_position_nm <- soma.positions$nucleus_position_nm
        # soma.positions$root_position <- soma.positions$nucleus_position
        # soma.positions$root_position <- apply(banc_nm2raw(soma.positions$nucleus_position_nm), 1, bancr:::paste_coords)
        # soma.positions$root_position <- gsub("\\(|\\)", "", soma.positions$root_position)
        # Get roots
        soma.positions <- dplyr::distinct(soma.positions, supervoxel_id, root_id, .keep_all = TRUE)
      }else{
        soma.positions <- data.frame()
      }
      
      # First loop: Neurons with no nucleus, take leaf points
      banc.missing.nuclei <- subset(banc.ids.todo, (is.na(root_position_nm)))
      #banc.missing.nuclei <- subset(banc.missing.nuclei,!root_id%in%banc.has.root.positions)
      if(nrow(banc.missing.nuclei)){
        banc.missing.nuclei.ids <- banc.missing.nuclei$root_id
        not.root.files <- nfiles[nfiles.ids%in%banc.missing.nuclei.ids]
        message("calculating missing roots")
        with_progress({
          p <- progressor(steps = length(not.root.files))
          banc.missing.nuclei.ids.guess <- foreach(nrfile = not.root.files, 
                                                   .combine = rbind, 
                                                   .packages = c("nat", "bancr"), 
                                                   .errorhandling = 'pass') %do% {
                                                     p(sprintf("Processing ID: %s", nrfile))
                                                     tryCatch({
                                                       neuron <- read.neuron(nrfile)
                                                       id <- gsub(".swc","",basename(nrfile))
                                                       leaves <- nat::endpoints(neuron)
                                                       npoints1 <- nat::xyzmatrix(neuron)[leaves,]
                                                       if(nrow(npoints1) == 0) return(NULL)
                                                       npoints <- npoints1
                                                       pin <- nat::pointsinside(x = npoints, surf = bancr::banc_neuropil.surf)
                                                       npoints2 <- data.frame(npoints[!pin,])
                                                       if(nrow(npoints2) > 0 && sum(!pin) > 2) npoints <- npoints2
                                                       pin <- nat::pointsinside(x = npoints, surf = bancr::banc_neck_connective.surf)
                                                       npoints3 <- data.frame(npoints[!pin,])
                                                       if(nrow(npoints3) > 0 && sum(!pin) > 2) npoints <- npoints3
                                                       if(!nrow(npoints)){
                                                         next
                                                       }
                                                       if(is.null(nrow(npoints))){
                                                         npoints <- matrix(npoints, ncol = 3)
                                                       }
                                                       npoints <- as.data.frame(npoints)
                                                       npoints$nucleus_id <- 0
                                                       npoints$root_id <- id
                                                       npoints$root_position_nm <- apply(npoints, 1, function(x) paste(x[c("X","Y","Z")], collapse=","))
                                                       data.frame(npoints[1,])
                                                     },
                                                     error = function(e) {message("failed root estimation on: ",id); NULL})
                                                   } 
        })
        soma.positions <- plyr::rbind.fill(soma.positions, banc.missing.nuclei.ids.guess)
      }
      soma.positions$root_position <- apply(banc_nm2raw(xyzmatrix(soma.positions$root_position_nm)), 1, bancr:::paste_coords)
      soma.positions$root_position <- gsub("\\(|\\)", "", soma.positions$root_position)
      
      # Second loop: Work out what region neuron belongs in
      optic_lobes <- subset(banc_brain_neuropils.surf, "optic")
      central_brain <- subset(banc_brain_neuropils.surf, "midbrain|central_brain")
      message("calculating region of innervation")
      with_progress({
        p <- progressor(steps = length(nfiles))
        results2 <- foreach(nfile = nfiles, 
                            .combine = rbind, 
                            .packages = c("nat", "bancr"), 
                            .errorhandling = 'pass') %do% {
                              tryCatch({
                                p(sprintf("Processing ID: %s", nfile))
                                neuron <- read.neuron(nfile)
                                id <- gsub(".swc","",basename(nfile))
                                leaves <- nat::endpoints(neuron)
                                points <- nat::xyzmatrix(neuron)
                                in.vnc <- any(pointsinside(x = points, surf = banc_vnc_neuropil.surf))
                                in.brain <- any(pointsinside(x = points, surf = central_brain))
                                in.optic_lobes <- any(pointsinside(x = points, surf = optic_lobes))
                                in.neck <- any(pointsinside(x = points, surf = banc_neck_connective.surf))
                                in.greater.brain <- any(pointsinside(x = points, surf = banc_brain_neuropil.surf))
                                in.cortex <- any(pointsinside(x = points, surf = banc.surf))
                                if(in.brain && in.vnc) {
                                  region <- "neck_connective"
                                } else if(in.vnc) {
                                  region <- "ventral_nerve_cord"
                                } else if(in.optic_lobes) {
                                  region <- "optic_lobe"
                                } else if(in.brain) {
                                  region <- "central_brain"
                                } else if(in.neck) {
                                  region <- "neck_connective"
                                } else if(in.greater.brain){
                                  region <- "brain"
                                } else if(in.cortex){
                                  in.cortex <- "rind"
                                } else {
                                  region <- NA
                                }
                                return(data.frame(root_id = id, region = region))
                              }, 
                              error = function(e) data.frame(root_id = id, region = NA))
                            }
      })
      
      # Update soma.positions with the new region information
      soma.positions <- merge(soma.positions, results2, by = "root_id", all.x = TRUE)
      
      # Third loop: Work synapse side index
      message("calculating side index")
      with_progress({
        p <- progressor(steps = nrow(banc.ids.todo))
        # run
        results3 <- foreach(id = banc.ids.todo$root_id, 
                            .packages = c("nat", "bancr", "R.utils"), 
                            .errorhandling = 'pass') %do% {
                              p()
                              result <- tryCatch({
                                withTimeout({
                                  # -- begin the original tryCatch body here --
                                  syn.file <-  file.path(banc.synapses.save.path,paste0(id,".csv"))
                                  if(file.exists(syn.file)){
                                    syns <- hemibrainr:::suppress(read_csv(syn.file, col_types = banc.col.types, show_col_types = FALSE, progress = FALSE))
                                  }else{
                                    syns <- data.frame()
                                  }
                                  if(!nrow(syns)){
                                    in.syns <- bancr::banc_partners(id, partners = "input")  %>%
                                      dplyr::filter(valid=='t') %>%
                                      dplyr::select(id, 
                                                    pre_pt_supervoxel_id, pre_pt_root_id,
                                                    post_pt_supervoxel_id, post_pt_root_id,
                                                    position = pre_pt_position) %>%
                                      dplyr::rowwise() %>%
                                      dplyr::mutate(position = paste(position,collapse=", ")) %>%
                                      dplyr::ungroup()
                                    in.syns$prepost <- 1
                                    out.syns <- bancr::banc_partners(id, partners = "output") %>%
                                      dplyr::filter(valid=='t') %>%
                                      dplyr::select(id, 
                                                    pre_pt_supervoxel_id, pre_pt_root_id,
                                                    post_pt_supervoxel_id, post_pt_root_id,
                                                    position = pre_pt_position) %>%
                                      dplyr::rowwise() %>%
                                      dplyr::mutate(position = paste(position,collapse=", ")) %>%
                                      dplyr::ungroup()
                                    out.syns$prepost <- 0
                                    syns <- plyr::rbind.fill(in.syns,out.syns)
                                    write.csv(syns, file = file.path(banc.synapses.save.path,paste0(id,".csv")))
                                  }
                                  if(!nrow(syns)){
                                    return(NULL)
                                  }
                                  
                                  mitos <- tryCatch(banc_mitochondria(id), error = function(e) 0)
                                  mito.volume <- tryCatch(sum(mitos$volume, na.rm = TRUE), error = function(e) 0)
                                  
                                  pos <- nat::xyzmatrix(syns$position)
                                  lrdiffs <- bancr:::banc_lr_position(pos,units = "nm")
                                  sides <- ifelse(lrdiffs>0,"right","left")
                                  syns$side <- sides
                                  
                                  if(nrow(syns)){
                                    data.frame(
                                      root_id = id,
                                      mitochondria = nrow(mitos),
                                      mitochondria_volume = mito.volume,
                                      input_connections = sum(syns$prepost==1, na.rm = TRUE),
                                      output_connections = sum(syns$prepost==0, na.rm = TRUE),
                                      input_side_index = round((sum(syns$side=="right"&syns$prepost==1)-sum(syns$side=="left"&syns$prepost==1))/sum(syns$prepost==1),4),
                                      output_side_index = round((sum(syns$side=="right"&syns$prepost==0)-sum(syns$side=="left"&syns$prepost==0))/sum(syns$prepost==0),4)
                                    )
                                  }else{
                                    NULL
                                  }
                                  # -- end of original body --
                                }, timeout = 3600, onTimeout = "error")  # 3600 seconds = 1 hour
                              }, error = function(e) {
                                # On error (timeout, etc.) -- dummy row of NAs
                                data.frame(
                                  root_id = id,
                                  mitochondria = NA_real_,
                                  mitochondria_volume = NA_real_,
                                  input_connections = NA_real_,
                                  output_connections = NA_real_,
                                  input_side_index = NA_real_,
                                  output_side_index = NA_real_
                                )
                              })
                              result
                            }
      })
      
      # Update soma.positions with the new region information
      #results3[is.nan(results3)] <- NA
      soma.positions <- merge(soma.positions, as.data.frame(do.call(rbind,results3)), by = "root_id", all.x = TRUE) %>%
        mutate(across(everything(), ~ ifelse(is.nan(.), NA, .)))
      
      # Update with L2 skeleton stats
      banc.metrics.files <- file.path(banc.metrics.save.path, paste0(soma.positions$root_id,".csv"))
      banc.metrics.files <- banc.metrics.files[file.exists(banc.metrics.files)]
      
      # Read metrics files and combine
      if(length(banc.metrics.files)){
        by.query.metrics <- foreach::foreach(mfile = banc.metrics.files) %do% {
          mdf <- readr::read_csv(mfile, col_types = banc.col.types, show_col_types = FALSE)
          numeric_columns <- sapply(mdf, is.numeric)
          mdf[numeric_columns] <- round(mdf[numeric_columns], digits = 2)
          mdf[,c("l2_n_trees","l2_endpoints","l2_branchpoints","l2_segments","l2_root")] <- NULL
          mdf
        }
        by.query.metrics <- by.query.metrics[unlist(lapply(by.query.metrics,is.data.frame))]
        banc.metrics <- do.call(plyr::rbind.fill, by.query.metrics)
        banc.metrics <- banc.metrics %>%
          dplyr::mutate(l2_cable_length_um = round(l2_cable_length/1000,2)) %>%
          dplyr::distinct(root_id, l2_nodes, l2_cable_length_um) %>%
          dplyr::arrange(dplyr::desc(l2_nodes))
        soma.positions <- dplyr::left_join(soma.positions,banc.metrics,by="root_id")
      }
      
      # Read split metrics files and combine 
      banc.split.metrics.files <- file.path(split.metrics.folder, paste0(soma.positions$root_id,".csv"))
      banc.split.metrics.files <-banc.split.metrics.files[file.exists(banc.split.metrics.files)]
      if(length(banc.split.metrics.files)){
        by.query.split.metrics <- foreach::foreach(mfile = banc.split.metrics.files) %do% {
          mdf <- readr::read_csv(mfile, col_types = banc.col.types, show_col_types = FALSE)
          numeric_columns <- sapply(mdf, is.numeric)
          mdf[numeric_columns] <- round(mdf[numeric_columns], digits = 2)
          mdf[,c("n_trees","endpoints","branchpoints","segments","root")] <- NULL
          mdf
        }
        by.query.split.metrics <- by.query.split.metrics[unlist(lapply(by.query.split.metrics,is.data.frame))]
        by.query.split.metrics <- do.call(plyr::rbind.fill, by.query.split.metrics)
        by.query.split.metrics <- by.query.split.metrics %>%
          dplyr::mutate(cable_length_um = round(as.numeric(cable_length)/1000,2)) %>%
          dplyr::distinct(root_id, .keep_all = TRUE) %>%
          dplyr::arrange(dplyr::desc(nodes)) %>%
          dplyr::select(-nucleus_id, -nucleus_position_nm)
        soma.positions <- dplyr::left_join(soma.positions,
                                           by.query.split.metrics %>% dplyr::select(-root_position_nm),
                                           by="root_id")
      }
      
      # save
      message("Adding metrics for new neurons: ", nrow(soma.positions))
      soma.positions.nuclei <- nrow(subset(soma.positions, nucleus_id != "0"))
      soma.positions.no.nuclei <- subset(soma.positions, nucleus_id == "0")
      soma.positions.no.nuclei <- length(unique(soma.positions.no.nuclei$root_id))
      soma.positions <- plyr::rbind.fill(soma.positions, banc.metrics.updated) %>%
        dplyr::distinct(root_id, .keep_all = TRUE)
      message("Writing results")
      readr::write_csv(soma.positions, file=file.path(banc.save.path,"banc_root_positions.csv"))
      banc.metrics.updated <- soma.positions 
      
      # announce
      message("##### BANCpipeline: banc root positions updated #####")
      message(sprintf("##### soma positions for : %d neurons, proxy leaf nodes for: %d neurons", soma.positions.nuclei, soma.positions.no.nuclei))
      
      # stop parallel backend
      stopImplicitCluster()
    }
  })
}

#######################
### UPDATE SEATABLE ###
#######################

# Read saved meta file
banc.metrics <- readr::read_csv(file=file.path(banc.save.path,"banc_root_positions.csv"),
                                col_types = banc.col.types, 
                                show_col_types = FALSE) %>%
  dplyr::filter(!is.na(root_id))
banc.meta <- readr::read_csv(file=file.path(banc.meta.save.path,"banc_meta.csv"),
                             col_types = banc.col.types,
                             show_col_types = FALSE)  %>%
  dplyr::filter(!is.na(root_id))
banc.metrics2 <- dplyr::full_join(banc.meta %>%
                                        dplyr::distinct(root_id, .keep_all = TRUE) %>%
                                        dplyr::select(root_id,l2_cable_length_um, l2_nodes),
                  banc.metrics %>%
                     dplyr::distinct(root_id, .keep_all = TRUE) %>%
                     dplyr::select(root_id, region, #root_region,
                                   mitochondria, mitochondria_volume, 
                                   input_connections, output_connections, 
                                   input_side_index, output_side_index, 
                                   pd_width, segregation_index),
                   by="root_id")
bc <- banctable_query("SELECT _id, status, root_id, super_class, region, root_region, root_position, root_position_nm, l2_cable_length_um, l2_nodes, pd_width, input_connections, output_connections, input_side_index, output_side_index, segregation_index from banc_meta") 
metrics.update <- banc.metrics2 %>%
  dplyr::left_join(bc %>% 
                     dplyr::distinct(`_id`, .keep_all = TRUE) %>%
                     dplyr::select(`_id`, root_region, root_id, bc_region = region, super_class), 
                   by = "root_id") %>%
  dplyr::mutate(region = dplyr::case_when(
    !is.na(bc_region)&!bc_region%in%c("brain","outside","rind",""," ") ~ bc_region,
    !is.na(region)&!region%in%c("brain","rind","") ~ region,
    super_class %in% c("visual_centrifugal") ~ "central_brain",
    super_class %in% c("visual_projection") ~ "optic_lobe",
    super_class %in% c("ventral_nerve_cord_intrinsic") ~ "ventral_nerve_cord",
    super_class %in% c("central_brain_intrinsic") ~ "central_brain",
    super_class %in% c("central_brain_intrinsic") ~ "central_brain",
    grepl('ascending|descending', super_class) ~ "neck_connective",
    region %in% c("midbrain","central_brain") ~ "central_brain",
    region %in% c("vnc","ventral_nerve_cord") ~ "ventral_nerve_cord",
    region %in% c("optic") ~ "optic_lobe",
    grepl('neck_connective', region) ~ "neck_connective",
    grepl('CAN|FLA|GNG|AMMC|SAD|PRW',root_region) | region=="sez" ~ "central_brain",
    grepl('_midbrain_',root_region) ~ "central_brain",
    grepl('optic', root_region) ~ "optic_lobe",
    grepl('optic', region) ~ "optic_lobe",
    grepl("midbrain",root_region) ~ "central_brain",
    grepl("optic",root_region) ~ "optic_lobe",
    grepl("vnc",root_region) ~ "ventral_nerve_cord",
    TRUE ~ region
  )) %>%
  #dplyr::select(-root_region) %>%
  dplyr::distinct(`_id`, .keep_all = TRUE) %>%
  dplyr::filter(!is.na(`_id`)) %>%
  as.data.frame()

# Correct data types
metrics.update$segregation_index <- as.numeric(metrics.update$segregation_index)
metrics.update$region[is.na(metrics.update$region)] <- ""
metrics.update[is.na(metrics.update)] <- 0
metrics.update$mitochondria<-as.numeric(metrics.update$mitochondria)
metrics.update$mitochondria_volume<-as.numeric(metrics.update$mitochondria_volume)

# Run update
metrics.update$bc_region <- NULL
metrics.update$super_class <- NULL
if(nrow(metrics.update)){
  banctable_update_rows(base='banc_meta', 
                        table = 'banc_meta', 
                        df = metrics.update,
                        append_allowed = FALSE, 
                        chunksize = 1000)  
}

##############################
### PLOT METRICS SUMMARIES ###
##############################

# # Update sides
# roots <- nat::xyzmatrix(soma.positions$root_position)
# roots <- banc_raw2nm(roots)
# lrdiffs <- bancr:::banc_lr_position(roots,units = "nm")
# sides <- ifelse(lrdiffs>0,"right","left")
# soma.positions$side <- sides  
# 
# # Step 1: Process the data
# processed_data <- soma.positions %>%
#   dplyr::mutate(region = ifelse(is.na(region),"undetermined",region),
#                 side = ifelse(is.na(side),"undetermined",side)) %>%
#   dplyr::group_by(region, side) %>%
#   dplyr::summarise(count = dplyr::n()) %>%
#   dplyr::ungroup() %>%
#   dplyr::mutate(
#     group = interaction(region,side, sep="_"),
#     group = fct_inorder(group)) %>%
#   dplyr::mutate(
#     proportion = count / sum(count),
#     label = paste0(group,"\n",count, scales::percent(proportion, accuracy = 0.1), ")")
#   )
# 
# # Step 2: Define custom colors (replace with your desired hex codes)
# region_colors <- c(
#   "optic" = hemibrainr:::hemibrain_bright_colors[["marine"]],
#   "midbrain" = hemibrainr:::hemibrain_bright_colors[["cerise"]],
#   "vnc" = hemibrainr:::hemibrain_bright_colors[["orange"]],
#   "neck_connective" = hemibrainr:::hemibrain_bright_colors[["pink"]],
#   "undetermined" = "lightgrey"
# )
# 
# # Step 3: Create a function to lighten colors
# lighten_color <- function(color, factor = 1.4) {
#   col_rgb <- col2rgb(color)
#   col_rgb <- pmin(col_rgb * factor, 255)
#   rgb(t(col_rgb), maxColorValue = 255)
# }
# 
# # Step 4: Create the color palette
# color_palette <- c(
#   sapply(region_colors, function(color) color),
#   sapply(region_colors, lighten_color, factor = 1.2),
#   sapply(region_colors, lighten_color, factor = 1.4)
# )
# names(color_palette) <- c(
#   paste0(names(region_colors), "_left"),
#   paste0(names(region_colors), "_right"),
#   paste0(names(region_colors), "_undetermined")
# )
# 
# # Calculate the positions for the labels
# processed_data <- processed_data %>%
#   mutate(
#     fraction = count / sum(count),
#     ymax = cumsum(fraction),
#     ymin = c(0, head(ymax, n=-1)),
#     position = (ymax + ymin) / 2,
#     position_label = 1.05  # This will place the labels just outside the pie
#   )
# 
# # Step 5: Create the plot
# pie_chart <- ggplot(processed_data, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = group)) +
#   geom_rect() +
#   geom_text(aes(label = label, x=4, y=position), 
#             nudge_x = .5,
#             size = 2) +
#   geom_segment(aes(x = 4, y = position, xend = 4.05, yend = position), color = "black", size = 0.5) +
#   coord_polar(theta = "y") +
#   scale_fill_manual(values = color_palette) +  # Assuming you've defined color_palette as before
#   xlim(c(0, 4.5)) +  # Adjust this to leave more space for labels
#   theme_void() +
#   theme(legend.position = "none") +
#   labs(fill = "region and side", title = "distribution of BANC neurons by region and side")
# 
# # Display the plot
# print(pie_chart)
# 
# # Save
# ggsave(plot = pie_chart,
#        filename = "inst/images/metrics/banc_region_overview.png", width = 8, height = 8, dpi = 300)


