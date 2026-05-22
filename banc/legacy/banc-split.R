#############################################
### Calculate BANC flow centrality splits ###
#############################################
source("banc/banc-startup.R")

# Choose BANC skeleton source
skeleton.folder <- banc.l2swc.save.path
split.master.folder <- banc.l2split.save.path
# png.split.folder <- paste0(banc.l2split.save.path,"_png")
# dir.create(png.split.folder)
#skeleton.folder <- banc.swc.save.path
#split.master.folder <- banc.split.save.path

# Save folders
split.folder <- file.path(split.master.folder,"swc")
images.folder <- file.path(split.master.folder,"images")
synapses.folder <- file.path(split.master.folder,"synapses")
metrics.folder <- file.path(split.master.folder,"metrics")
dir.create(split.folder, recursive = TRUE, showWarnings = FALSE)
dir.create(images.folder, recursive = TRUE, showWarnings = FALSE)
dir.create(synapses.folder, recursive = TRUE, showWarnings = FALSE)
dir.create(metrics.folder, recursive = TRUE, showWarnings = FALSE)

# Read IDs
# banc.meta <- readr::read_csv(file=file.path(banc.meta.save.path,"banc_meta.csv"),
#                              col_types = banc.col.types, 
#                              show_col_types = FALSE) 
bc <- banctable_query()
banc.meta <- bc %>%
  dplyr::filter(!grepl("DELETE",status),!grepl("NOT_A_NEURON",status),!grepl("GLIA|TRACHEA",status))
banc.root.ids <- unique(banc.meta$root_id)
sensories <- subset(banc.meta, grepl("sensory|afferent",super_class)|grepl("sensory|afferent",cell_class))$root_id
motors <- subset(banc.meta, grepl("motor|efferent|endocrine|visceral",super_class)|grepl("motor|efferent|endocrine|visceral",cell_class))$root_id

# For each  
swcs <- list.files(skeleton.folder, full.names = TRUE)
swc.ids <- gsub("\\.swc","", basename(swcs))
swcs <- swcs[swc.ids %in% banc.meta$root_id]
swc.ids <- gsub("\\.swc","", basename(swcs))

# Do no redo what we have done
dones.skels <- list.files(split.folder)
dones.syns <- list.files(synapses.folder)
done.ids.syn <- gsub("\\.csv","",dones.syns)
done.ids.swc <- gsub("\\.swc","",dones.skels)
done.ids <- intersect(done.ids.syn,done.ids.swc)
swcs <- sample(swcs[!swc.ids%in%done.ids])
message("attempting to split ", length(swcs), " BANC neurons, we have done: ", length(done.ids))

# Read root nodes
banc.roots <- readr::read_csv(file=file.path(banc.save.path,"banc_root_positions.csv"), 
                              col_types = banc.col.types, 
                              show_col_types = FALSE)
banc.roots <- as.data.frame(banc.roots)
banc.roots <- banc.roots %>% 
  dplyr::distinct(root_id, 
                  root_position,
                  root_position_nm)
banc.confirmed.roots <- bc %>% 
  dplyr::distinct(root_id, 
                 root_position,
                 root_position_nm) %>%
  rbind(banc.roots) %>%
  dplyr::distinct(root_id, 
                  root_position,
                  root_position_nm) %>%
  dplyr::filter(!is.na(root_position), 
                !is.na(root_position_nm),
                !grepl("NA",root_position_nm),
                !grepl("NA",root_position))
  
# Set up progress bar
iterations <- length(swc.ids)
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Batch the swc files into small batches
swcs <- sample(swcs)
multiplier <- 10
upper <- ifelse((numCores*multiplier)<length(swcs),numCores*multiplier,length(swcs))
batches <- split(swcs, round(seq(from = 1, to = upper, length.out = length(swcs))))

# register parallel
doMC.parallel <- FALSE

# Process splits achieved by other means
numCores <- 10
cl <- setup_parallel()

# Run
by.query <- foreach::foreach(batch = seq_along(batches), 
                             .verbose = numCores>1, 
                             .combine = 'c', 
                             .errorhandling='pass', 
                             .options.snow = opts) %dopar% {
                               
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
                                                        meta = as.data.frame(bc))
                               
                               # resample neuron, .1 um
                               neurons <- nat:::resample.neuronlist(neurons, 100, OmitFailures = TRUE)
                               
                               # Re-root neurons
                               id.updated <- banc_latestid(names(neurons))
                               neurons <- hemibrainr::add_field_seq(neurons, entries=id.updated, field="id")
                               bc.neurons.rerooted <- nlapply(neurons,
                                                              banc_reroot, 
                                                              id = NULL, 
                                                              roots = banc.confirmed.roots,
                                                              estimate = FALSE)

                               # Add synapses
                               id.updated <- banc_latestid(names(bc.neurons.rerooted))
                               bc.neurons.rerooted <- hemibrainr::add_field_seq(bc.neurons.rerooted, entries=id.updated, field="id")
                               bc.neurons.syns <- nlapply(bc.neurons.rerooted,
                                                          banc_add_synapses, 
                                                          OmitFailures = TRUE,
                                                          id = NULL)
                               
                               # Drop those with no synapses
                               good.ids <- c()
                               for(i in 1:length(bc.neurons.syns)){
                                 if(!nrow(bc.neurons.syns[[i]]$connectors)){
                                   message("no synapses for: ", names(bc.neurons.syns)[i])
                                   next
                                 }else if(!length(bc.neurons.syns[[i]]$connectors)){
                                   message("no synapses for: ", names(bc.neurons.syns)[i])
                                   next
                                 }else{
                                   good.ids <- c(good.ids, i)
                                 }
                               }
                               bc.neurons.syns <- bc.neurons.syns[good.ids]
                               if(!length(bc.neurons.syns)){
                                 return(NULL)
                               }
                               
                               # Flag synapses so they do not affect flow_centrality calculation
                               bc.neurons.syns.rem <- bc.neurons.syns
                               somas <- sapply(bc.neurons.syns, function(x) !is.null(x$tags$soma))
                               for(i in 1:length(bc.neurons.syns)){
                                 rem <- remove_bad_synapses(bc.neurons.syns[i],
                                                            meshes = NULL,
                                                            soma = somas[i],
                                                            min.nodes.from.soma = 150, # 15 microns
                                                            min.nodes.from.pnt = 10, # 0.1 microns
                                                            primary.branchpoint = 0.25,
                                                            OmitFailures = TRUE,
                                                            .parallel = FALSE,
                                                            wipe = TRUE)
                                 bc.neurons.syns.rem[[i]] <- rem[[1]]
                               }
                               ## check: 720575941559458675
                               
                               # Split skeletons
                               bc.neurons.flow <- hemibrainr::flow_centrality(bc.neurons.syns.rem, 
                                                                              mode = mode, 
                                                                              polypre = polypre, 
                                                                              split = split, 
                                                                              .parallel = FALSE, 
                                                                              OmitFailures = FALSE)
                               
                               # Remove bad synapses now knowing the flow_centrality calc
                               bc.neurons.flow.rem <- bc.neurons.flow
                               for(i in 1:length(bc.neurons.flow)){
                                 bc.neuron.flow.rem <- remove_bad_synapses(bc.neurons.flow[i],
                                                                           meshes = NULL,
                                                                           soma = TRUE,
                                                                           min.nodes.from.soma = 50, # 5 microns
                                                                           min.nodes.from.pnt = 5,# 0.5 microns
                                                                           primary.branchpoint = 0.25,
                                                                           OmitFailures = TRUE,
                                                                           .parallel = FALSE,
                                                                           wipe = TRUE)
                                 tryCatch({
                                   bc.neurons.flow.rem[[i]] <- bc.neuron.flow.rem[[1]]
                                 }, error = function(e){
                                   bc.neurons.flow.rem <- bc.neurons.flow[[1]]
                                 })
                               }
                               bc.neurons.flow <- bc.neurons.flow.rem
                               bc.neurons.flow <- bc.neurons.flow[sapply(bc.neurons.flow, function(x) is.neuron(x))]
                               
                               # Assert label
                               ids.sensories <- intersect(names(bc.neurons.flow),sensories)
                               if (length(ids.sensories)){
                                 for(is in ids.sensories){
                                   bc.neurons.flow[[is]] <- hemibrainr:::add_Label(bc.neurons.flow[[is]], Label = 2)
                                 }
                               }
                               ids.motors <- intersect(names(bc.neurons.flow),motors)
                               if (length(ids.motors)){
                                 for(is in ids.motors){
                                   bc.neurons.flow[[is]] <- hemibrainr:::add_Label(bc.neurons.flow[[is]], Label = 3)
                                 }
                               }
                               
                               # remove nulls
                               for(id in names(bc.neurons.flow)){
                                 # Get data
                                 neuron <- bc.neurons.flow[[id]]
                                 syns.df <- neuron$connectors
                                 
                                 # Fix
                                 d <- neuron$d
                                 g <- nat::as.ngraph(neuron)
                                 labeled_nodes = which(d$Label != 0)
                                 unlabeled_nodes = which(d$Label == 0)
                                 if(length(unlabeled_nodes)){
                                   if(length(unlabeled_nodes) > 0 && length(labeled_nodes) > 0) {
                                     dm = igraph::distances(g, v=unlabeled_nodes, to=labeled_nodes)
                                     closest_labeled_idx = apply(dm, 1, which.min)
                                     assigned_label = d$Label[labeled_nodes[closest_labeled_idx]]
                                     d$Label[unlabeled_nodes] = assigned_label
                                     neuron$d = d
                                   }
                                   if(nrow(syns.df)){
                                     node_coords <- as.matrix(neuron$d[,c("X","Y","Z")])
                                     syn_coords <- as.matrix(syns.df[,c("x","y","z")])
                                     nearest_node_idx <- RANN::nn2(node_coords, syn_coords, k = 1)$nn.idx[,1]
                                     syns.df$Label <- neuron$d$Label[nearest_node_idx]
                                     neuron$connectors <- syns.df
                                   } 
                                 }
                               }
                               
                               # Assign Strahler
                               bc.neurons.split <- assign_strahler(bc.neurons.flow, 
                                                                   .parallel = FALSE, 
                                                                   OmitFailures = TRUE)
                               
                               # Make sure connectors also have labels
                               bc.neurons.split <- nat::nlapply(bc.neurons.split, 
                                                                hemibrainr:::carryover_labels, 
                                                                .parallel = FALSE, 
                                                                OmitFailures = TRUE)
                               bc.neurons.split[,"root_id"] <- names(bc.neurons.split)
                               
                               # Save synapses in another folder
                               message("writing data: ")
                               write.it <- capture.output(
                                 write.neurons(nl = bc.neurons.split,
                                               dir = split.folder,
                                               files = names(bc.neurons.split),
                                               Force = TRUE,
                                               format = "swc", 
                                               include.data.frame = FALSE) 
                               )
                               for(n in names(bc.neurons.split)){
                                 csv <- bc.neurons.split[n][[1]]$connectors
                                 write.csv(x = csv, 
                                           file = file.path(synapses.folder,paste0(n,".csv")), 
                                           row.names = FALSE)    
                               }
                               
                               # segregation index
                               bc.neurons.batch <- hemibrainr:::segregation_index.neuronlist(bc.neurons.split)
                               
                               # Scale
                               bc.neurons.split.microns <- hemibrainr::scale_neurons(bc.neurons.batch, 
                                                                                     scaling = 1/1000, 
                                                                                     .parallel = FALSE, 
                                                                                     OmitFailures = TRUE) # convert to microns
                               bc.neurons.split.microns[,] <- bc.neurons.split.microns[,c("root_id","nucleus_id","root_position_nm")]
                               mets <- hemibrain_compartment_metrics(bc.neurons.split.microns, 
                                                                     OmitFailures = TRUE, 
                                                                     .parallel = FALSE, 
                                                                     delta = 5, 
                                                                     resample = NULL, 
                                                                     locality = FALSE)
                               colnames(mets) <- snakecase::to_snake_case(colnames(mets))
                               bc.metrics <- mets
                               bc.metrics$root_id <- bc.metrics$id
                               bc.metrics$id <- NULL
                               bc.metrics$projection_score <- NULL
                               bc.metrics <- bc.metrics[!duplicated(bc.metrics$root_id),]
                               
                               # projection scores
                               proj.scores <- hemibrainr:::projection_score.neuronlist(bc.neurons.split.microns,
                                                                                       .parallel = FALSE,
                                                                                       OmitFailures = TRUE)
                               colnames(proj.scores) <- c("root_id", "projection_score")
                               bc.metrics <- dplyr::left_join(bc.metrics, proj.scores, by = "root_id")
                               
                               # summary of skeleton structure
                               summaries <- summary(bc.neurons.split.microns)[,c("nodes","cable.length")]
                               colnames(summaries) <- snakecase::to_snake_case(colnames(summaries))
                               summaries$root_id <- rownames(summaries)
                               summaries$nsoma <- NULL
                               keep <- c("root_id", setdiff(colnames(bc.metrics),colnames(summaries)))
                               bc.metrics <- dplyr::full_join(bc.metrics[,keep], summaries, by ="root_id")
                               colnames(bc.metrics) <- snakecase::to_snake_case(colnames(bc.metrics))
                               
                               # simplify                    
                               bc.metrics.sq <- round_dataframe(bc.metrics, digits = 4)
                               
                               # write
                               message("metrics calculated for neurons: ", nrow(bc.metrics.sq))
                               nmets <- length(list.files(metrics.folder))+1
                               write <- readr::write_csv(bc.metrics.sq, 
                                                         file = file.path(metrics.folder,sprintf("banc_metrics_%d.csv",nmets)),
                                                         col_names = TRUE,
                                                         append = FALSE)
                               
                               # Take an image of each neuron
                               for(id in names(bc.neurons.split)){
                                 try({
                                   #neuron3 <- banc_read_neuron_meshes(id)
                                   img <- banc_neuron_comparison_plot(filename = file.path(images.folder,paste0(id,".png")),
                                                                      neuron1 = bc.neurons.split[id],
                                                                      #neuron3 = neuron3,
                                                                      neuron1.info = id)
                                 })
                               }
                               
                               # Return
                               if(length(bc.neurons.split)==1){
                                 message("completed: ", names(bc.neurons.split), " : ", which(swcs==swc), "/",length(swcs)) 
                               }else{
                                 message("completed batch: ", batch, ": neuron count: ", length(bc.neurons.split)) 
                               }
                             }

# Stop cores
stop_parallel(cl)

# Were there errors?
for(i in 1:length(by.query)){
  if(!is.null(by.query[[i]])){
    message(by.query[[i]])
  }
}

