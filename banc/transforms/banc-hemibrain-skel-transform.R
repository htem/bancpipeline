#' banc-hemibrain-skel-transform — Transform hemibrain skeletons into BANC space.
#'
#' For each hemibrain bodyid, fetches the skeleton (neuprint) and applies
#' the JRCFIB2018F → JRC2018F → BANC registration chain. Saves the
#' transformed SWC.
#'
#' @section Reads:
#'   - `<banc.meta.save.path>/hemibrain_meta.csv`
#'   - hemibrain SWCs from `<hemibrain.save.path>`
#'
#' @section Writes:
#'   - `<banc.nblast.hemibrain.swc.save.path>/<version>/<bodyid>.swc`

#############################################
### TRANSFORM hemibrain SKELETONS TO BANC ###
#############################################
source("banc/banc-startup.R")
nat.jrcbrains::register_saalfeldlab_registrations()
version <- banc.nblast.version
redo <- FALSE

# Create save folder
banc.nblast.hemibrain.swc.save.path.version <- file.path(banc.nblast.hemibrain.swc.save.path,version)
dir.create(banc.nblast.hemibrain.swc.save.path.version, recursive = TRUE, showWarnings = FALSE)

# Get hemibrain meta data
meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"hemibrain_meta.csv"), 
                                            col_types = hemibrainr:::sql_col_types))
ids <- unique(meta$bodyid)
# positions <- meta %>%
#   dplyr::mutate(root_pos = dplyr::case_when(
#     !is.na(soma_location)&soma_location!="" ~ soma_location,
#     !is.na(tosoma_location)&tosoma_location!="" ~ tosoma_location,
#     !is.na(root_location)&root_location!="" ~ root_location,
#     TRUE ~ ""
#   )) %>%
#   dplyr::mutate(root_pos=gsub("list\\(|\\).*","",root_pos)) %>%
#   dplyr::distinct(bodyid, root_pos, .keep_all = TRUE)

# Stable skeletons are just fetched from neuprint, so we only save the transforms
swc.full.files <- unique(unlist(lapply(ids, function(id) file.path(banc.nblast.hemibrain.swc.save.path.version,paste0(id,".swc")))))
swc.full.files <- swc.full.files[file.exists(swc.full.files)]
swc.ids <- gsub(pattern = "\\.swc$", "", basename(swc.full.files))

# What SWC files have we not downloaded?
if(!redo){
  todo.ids <- setdiff(ids, swc.ids)
}else{
  todo.ids <- ids
}
message("There are ", length(todo.ids), " todo, we have downloaded and transformed ", length(swc.ids), " hemibrain .swc files")

# Get neurons in hemibrain space
hb.split.JRCFIB2018Fraw <- read.neuronlistfh(file.path(hemibrain.save.path,"JRCFIB2018Fraw/hemibrain_all_neurons_flow_JRCFIB2018Fraw.rds"))

# Register cores
cl <- setup_parallel()

# Transforming skeletons
message("##### Transforming skeletons #####")
by.query <- foreach(id = sample(todo.ids),
                    .combine = 'c',
                    .packages = c('nat.jrcbrains','Morpho','rJava'),
                    .errorhandling = 'pass') %dopar% {
                      
                      # File
                      new.file <- file.path(banc.nblast.hemibrain.swc.save.path.version, paste0(id,".swc"))
                      
                      # Get neuron. Note: Microns will fit it into hemibrain.surf
                      neuron <- hb.split.JRCFIB2018Fraw[id]

                      # pos <- subset(positions, positions$bodyid==id)[1,]
                      # soma <- nat::xyzmatrix(pos$root_pos) # in raw coordinates
                      # if(nrow(soma)&!is.na(soma[1])){
                      #   soma <- soma*rep(8, 3)
                      #   neuron <- nat::reroot(x = neuron, point = c(soma))
                      #   neuron$tags$soma <- nat::rootpoints(neuron) 
                      # }
                      
                      # Change to nm
                      # nat::xyzmatrix(neuron) <- nat::xyzmatrix(neuron)*rep(8, 4)
                      
                      # Transform to JRCVNC2018F
                      neuron <- scale_neurons(neuron)
                      neuron.jrc2018f <- nat.templatebrains::xform_brain(neuron, sample = "JRCFIB2018F", reference = "JRC2018F")
                      
                      # Transform to BANC
                      neuron.banc <- bancr::banc_to_JRC2018F(neuron.jrc2018f, region="brain", method="tpsreg", banc.units = "nm", inverse = TRUE)
                      
                      # Save as a .swc file
                      nat::write.neuron(neuron.banc[[1]], file = new.file, Force = TRUE)
                      
                      # Return nothing
                      NULL
                    }                    

# Stop cores
stop_parallel(cl)

# Were there errors?
message("##### Displaying any errors from foreach loop #####")
for(i in 1:length(by.query)){
  if(!is.null(by.query[[i]])){
    message(by.query[[i]])
  }
}

# Announce
done <- list.files(banc.nblast.hemibrain.obj.save.path)
message("##### BANCpipeline: hemibrain-to-banc skeleton transforms calculated #####")
message(sprintf("##### we new transformed skeletons: %s of %s neurons", nrow(swc.full.files)), length(done))