#############################################
### Calculate MANC flow centrality splits ###
#############################################
source("banc/banc-startup.R")
redo <- FALSE

# Choose MANC skeleton source
skeleton.folder <- file.path(banc.nblast.manc.swc.save.path, banc.nblast.version)
split.master.folder <- banc.nblast.manc.split.save.path

# Save folders
split.folder <- file.path(split.master.folder,"swc")
images.folder <- file.path(split.master.folder,"images")
synapses.folder <- file.path(split.master.folder,"synapses")
metrics.folder <- file.path(split.master.folder,"metrics")
dir.create(split.folder, recursive = TRUE, showWarnings = FALSE)
dir.create(images.folder, recursive = TRUE, showWarnings = FALSE)
dir.create(synapses.folder, recursive = TRUE, showWarnings = FALSE)
dir.create(metrics.folder, recursive = TRUE, showWarnings = FALSE)

# Get manc meta data
mc.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"manc_meta.csv"), 
                                            col_types = hemibrainr:::sql_col_types))
mc.ids <- unique(mc.meta$bodyid)
sensories <- subset(mc.meta, !is.na(receptor_type)|grepl("SN|SA|^DN",cell_type))$bodyid
motors <- subset(mc.meta, grepl("MN|^AN",cell_type)|grepl("ascending",cell_class))$bodyid

# For each  
swcs <- list.files(skeleton.folder, full.names = TRUE)
swc.ids <- gsub("\\.swc","", basename(swcs))
swcs <- swcs[swc.ids %in% mc.meta$bodyid]
swc.ids <- gsub("\\.swc","", basename(swcs))

# Do no redo what we have done
dones.skels <- list.files(split.folder)
dones.syns <- list.files(synapses.folder)
done.ids.syn <- gsub("\\.csv","",dones.syns)
done.ids.swc <- gsub("\\.swc","",dones.skels)
done.ids <- intersect(done.ids.syn,done.ids.swc)
swcs <- sample(swcs[!swc.ids%in%done.ids])
message("attempting to split ", length(swcs), " MANC neurons, we have done:", length(done.ids))

# Read root nodes
manc.roots <- mc.meta %>%
  dplyr::mutate(pt_position = dplyr::case_when(
    soma_location!="" ~ soma_location,
    tosoma_location!="" ~ tosoma_location,
    root_location!="" ~ root_location,
    TRUE ~ ""
  )) %>%
  dplyr::mutate(pt_position=gsub("list\\(|\\).*","",pt_position),
                nucleus_id = soma_location) %>%
  dplyr::distinct(bodyid, pt_position) %>%
  dplyr::filter(pt_position!="")
mc.positions.xyz <- nat::xyzmatrix(manc.roots$pt_position)*rep(8/1000, 3)
mc.positions.xyz.jrcvnc2018f <- nat.templatebrains::xform_brain(mc.positions.xyz, 
                                                                sample = "MANC", 
                                                                reference = "JRCVNC2018F")
mc.positions.xyz.banc <- bancr::banc_to_JRC2018F(mc.positions.xyz.jrcvnc2018f, 
                                                 region="vnc", 
                                                 method="tpsreg", 
                                                 banc.units = "nm", 
                                                 inverse = TRUE)
manc.roots$pt_position <- apply(mc.positions.xyz.banc,1,paste,collapse=",")
manc.roots$root_id <- manc.roots$bodyid

# Set up progress bar
iterations <- length(swc.ids)
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Batch the swc files into small batches
swcs <- sample(swcs)
multiplier <- 100
upper <- ifelse((numCores*multiplier)<length(swcs),numCores*multiplier,length(swcs))
batches <- split(swcs, round(seq(from = 1, to = upper, length.out = length(swcs))))

# register parallel
doMC.parallel <- FALSE

# Process splits achieved by other means
cl <- parallel::makeForkCluster(numCores)

# For dopar
doParallel::registerDoParallel(cl)

# Run
by.query <- foreach::foreach(batch = seq_along(batches), 
                             .verbose = numCores>1, 
                             .combine = 'c', 
                             .errorhandling='pass', 
                             .options.snow = opts) %do% {
                               
                               # Check if already processed
                               swc <- batches[[batch]]
                               target1 <- file.path(split.folder, basename(swc))
                               exists1 <- sapply(target1, file.exists)
                               target2 <- file.path(synapses.folder, gsub("\\.swc",".csv",basename(swc)))
                               exists2 <- sapply(target2, file.exists)
                               exists <- exists1 & exists2
                               swc <- swc[!exists]
                               if(!length(swc)){
                                 return(NULL)
                               }
                               ids <- gsub("\\.swc","", basename(swc))
                               message("running on: ", paste(ids, collapse = ", "))
                               
                               # Read neurons
                               neurons <- banc_read_swc(swc, 
                                                        ids = ids, 
                                                        meta = as.data.frame(mc.meta),
                                                        id = 'bodyid', 
                                                        template = "MANC")
                               
                               # resample neuron, .1 um
                               neurons <- nat:::resample.neuronlist(neurons, 100, OmitFailures = TRUE)
                               
                               # Re-root neurons
                               ids <- names(neurons)
                               neurons <- hemibrainr::add_field_seq(neurons, entries=ids, field="id")
                               mc.neurons.rerooted <- nlapply(neurons,
                                                              banc_reroot, 
                                                              id = NULL, 
                                                              roots = manc.roots,
                                                              estimate = FALSE)
                               
                               # Add synapses
                               ids <- names(mc.neurons.rerooted)
                               mc.neurons.rerooted <- hemibrainr::add_field_seq(mc.neurons.rerooted, entries=ids, field="id")
                               connectors <- neuprint_get_synapses(bodyids = ids, 
                                                                   dataset = 'manc:v1.2.1', 
                                                                   conn = manc_neuprint()) %>%
                                 dplyr::mutate(post_id=ifelse(prepost,bodyid,partner),
                                               pre_id=ifelse(prepost,partner,bodyid)) %>%
                                 dplyr::select(-bodyid,-partner) %>%
                                 dplyr::rename(X=x, Y=y, Z=z)
                               xyz.jrcvnc2018f <- nat.templatebrains::xform_brain(nat::xyzmatrix(connectors)*rep(8/1000, 3), 
                                                                                  sample = "MANC", 
                                                                                  reference = "JRCVNC2018F")
                               xyz.banc <- bancr::banc_to_JRC2018F(xyz.jrcvnc2018f, 
                                                                   region="vnc", 
                                                                   method="tpsreg", 
                                                                   banc.units = "nm", 
                                                                   inverse = TRUE)
                               nat::xyzmatrix(connectors) <- xyz.banc 
                               mc.neurons.syns <- nat::nlapply(mc.neurons.rerooted,
                                                               banc_add_synapses,
                                                               connectors = connectors,
                                                               update.id = FALSE,
                                                               id = NULL)
                               
                               # Drop those with no synapses
                               good.ids <- c()
                               for(i in 1:length(mc.neurons.syns)){
                                 if(!length(mc.neurons.syns[[i]]$connectors)){
                                   message("no synapses for: ", names(mc.neurons.syns)[i])
                                   next
                                 }else if(!nrow(mc.neurons.syns[[i]]$connectors)){
                                   message("no synapses for: ", names(mc.neurons.syns)[i])
                                   next
                                 }else{
                                   good.ids <- c(good.ids, i)
                                 }
                               }
                               mc.neurons.syns <- mc.neurons.syns[good.ids]
                               if(!length(mc.neurons.syns)){
                                 return(NULL)
                               }
                               
                               # Flag synapses so they do not affect flow_centrality calculation
                               mc.neurons.syns.rem <- mc.neurons.syns
                               somas <- sapply(mc.neurons.syns, function(x) !is.null(x$tags$soma))
                               for(i in 1:length(mc.neurons.syns)){
                                 mc.neuron.syns.rem <- remove_bad_synapses(mc.neurons.syns[i],
                                                                           meshes = banc_neuropil.surf,
                                                                           soma = somas[i],
                                                                           min.nodes.from.soma = 150, # 15 microns
                                                                           min.nodes.from.pnt = 10, # 0.1 microns
                                                                           primary.branchpoint = 0.25,
                                                                           OmitFailures = TRUE,
                                                                           .parallel = FALSE,
                                                                           wipe = TRUE)
                                 mc.neurons.syns.rem[[i]] <- mc.neuron.syns.rem[[1]]
                               }
                               
                               # Split skeletons
                               mc.neurons.flow <- mc.neurons.syns.rem
                               ids.sensories <- intersect(ids,sensories)
                               if (length(ids.sensories)){
                                 for(is in ids.sensories){
                                   mc.neurons.flow[[is]] <- hemibrainr:::add_Label(mc.neurons.flow[[is]], Label = 2)
                                 }
                               }
                               ids.motors <- intersect(ids,motors)
                               if (length(ids.motors)){
                                 for(is in ids.motors){
                                   mc.neurons.flow[[is]] <- hemibrainr:::add_Label(mc.neurons.flow[[is]], Label = 3)
                                 }
                               }
                               others <- setdiff(ids,c(motors,sensories))
                               if(length(others)){
                                 for(is in names(mc.neurons.syns.rem)){
                                   mc.neurons.syns.rem[[is]] <- hemibrainr:::add_Label(mc.neurons.syns.rem[[is]], Label = 0)
                                 }
                                 mc.neurons.flow[others] <- hemibrainr::flow_centrality(mc.neurons.syns.rem[others], 
                                                                                        mode = mode, 
                                                                                        polypre = polypre, 
                                                                                        split = split, 
                                                                                        .parallel = FALSE, 
                                                                                        OmitFailures = FALSE)
                                 fails <- unlist(sapply(mc.neurons.flow,function(x) class(x)=="try-error"))
                                 if(sum(fails)){
                                   mc.neurons.flow[fails] <- mc.neurons.syns.rem[fails]
                                 }
                               }
                               
                               # Remove bad synapses now knowing the flow_centrality calc
                               mc.neurons.flow.rem <- mc.neurons.flow
                               for(i in 1:length(mc.neurons.flow)){
                                 mc.neuron.flow.rem <- remove_bad_synapses(mc.neurons.flow[i],
                                                                           meshes = NULL,#?
                                                                           soma = TRUE,
                                                                           min.nodes.from.soma = 50, # 5 microns
                                                                           min.nodes.from.pnt = 5,# 0.5 microns
                                                                           primary.branchpoint = 0.25,
                                                                           OmitFailures = TRUE,
                                                                           .parallel = FALSE,
                                                                           wipe = TRUE)
                                 mc.neurons.flow.rem[[i]] <- mc.neuron.flow.rem[[1]]
                               }
                               
                               # Assign Strahler
                               mc.neurons.split <- assign_strahler(mc.neurons.flow, 
                                                                   .parallel = FALSE, 
                                                                   OmitFailures = TRUE)
                               
                               # Make sure connectors also have labels
                               mc.neurons.split <- nat::nlapply(mc.neurons.split, 
                                                                hemibrainr:::carryover_labels, 
                                                                .parallel = FALSE, 
                                                                OmitFailures = TRUE)
                               mc.neurons.split[,"bodyid"] <- names(mc.neurons.split)
                               
                               # Save synapses in another folder
                               message("writing data: ")
                               write.it <- capture.output(
                                 write.neurons(nl = mc.neurons.split,
                                               dir = split.folder,
                                               files = names(mc.neurons.split),
                                               Force = TRUE,
                                               format = "swc", 
                                               include.data.frame = FALSE) 
                               )
                               for(n in names(mc.neurons.split)){
                                 csv <- mc.neurons.split[n][[1]]$connectors
                                 write.csv(x = csv, 
                                           file = file.path(synapses.folder,paste0(n,".csv")), 
                                           row.names = FALSE)    
                               }
                               
                               # segregation index
                               mc.neurons.batch <- hemibrainr:::segregation_index.neuronlist(mc.neurons.split)
                               
                               # Scale
                               mc.neurons.split.microns <- hemibrainr::scale_neurons(mc.neurons.batch, 
                                                                                     scaling = 1/1000, 
                                                                                     .parallel = FALSE, 
                                                                                     OmitFailures = TRUE) # convert to microns
                               mets <- hemibrain_compartment_metrics(mc.neurons.split.microns, 
                                                                     OmitFailures = TRUE, 
                                                                     .parallel = FALSE, 
                                                                     delta = 5, 
                                                                     resample = NULL, 
                                                                     locality = FALSE)
                               colnames(mets) <- snakecase::to_snake_case(colnames(mets))
                               mc.metrics <- mets
                               mc.metrics$bodyid <- mc.metrics$id
                               mc.metrics$id <- NULL
                               mc.metrics$projection_score <- NULL
                               mc.metrics <- mc.metrics[!duplicated(mc.metrics$bodyid),]
                               
                               # Projection scores
                               proj.scores <- hemibrainr:::projection_score.neuronlist(mc.neurons.split.microns,
                                                                                       .parallel = FALSE,
                                                                                       OmitFailures = TRUE)
                               colnames(proj.scores) <- c("bodyid", "projection_score")
                               mc.metrics <- dplyr::left_join(mc.metrics, proj.scores, by = "bodyid")
                               
                               # summary of skeleton structure
                               summaries <- summary(mc.neurons.split.microns)[,c("nodes","cable.length")]
                               colnames(summaries) <- snakecase::to_snake_case(colnames(summaries))
                               summaries$bodyid <- rownames(summaries)
                               summaries$nsoma <- NULL
                               keep <- c("bodyid", setdiff(colnames(mc.metrics),colnames(summaries)))
                               mc.metrics <- dplyr::full_join(mc.metrics[,keep], summaries, by ="bodyid")
                               colnames(mc.metrics) <- snakecase::to_snake_case(colnames(mc.metrics))
                               
                               # Simplify                    
                               mc.metrics.sq <- round_dataframe(mc.metrics, digits = 4)
                               
                               # Write
                               message("metrics calculated for neurons: ", nrow(mc.metrics.sq))
                               nmets <- length(list.files(metrics.folder))+1
                               write <- readr::write_csv(mc.metrics.sq, 
                                                         file = file.path(metrics.folder,sprintf("manc_metrics_%d.csv",nmets)),
                                                         col_names = TRUE,
                                                         append = FALSE)
                               
                               # Take an image of each neuron
                               for(id in names(mc.neurons.split)){
                                 try({
                                   img <- banc_neuron_comparison_plot(filename = file.path(images.folder,paste0(id,".png")),
                                                                      neuron1 = mc.neurons.split[id],
                                                                      neuron1.info = id)
                                 })
                               }
                               
                               # Return
                               if(length(mc.neurons.split)==1){
                                 message("completed: ", names(mc.neurons.split), " : ", which(swcs==swc), "/",length(swcs)) 
                               }else{
                                 message("completed batch: ", batch, ": neuron count: ", length(mc.neurons.split)) 
                               }
                             }

# Stop cores
parallel::stopCluster(cl)

# Were there errors?
for(i in 1:length(by.query)){
  if(!is.null(by.query[[i]])){
    message(by.query[[i]])
  }
}

