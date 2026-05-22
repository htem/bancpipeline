#' banc-fafb-skel-transform — Transform FAFB (FlyWire v783) skeletons into BANC space.
#'
#' For each FAFB-v783 root_id, picks the detailed (split + rerooted)
#' skeleton if available else the L2 skeleton, and applies the FlyWire-
#' JRC2018F-BANC registration chain. Saves the transformed SWC.
#'
#' @section Reads:
#'   - `<banc.meta.save.path>/flywire_meta.csv`
#'   - SWCs from `<flywire.split.save.path>` (preferred) or `<flywire.swc.save.path>`
#'
#' @section Writes:
#'   - `<banc.nblast.fafb.swc.save.path>/<version>/<root_783>.swc`

####################################################
### TRANSFORM FAFB-FLYWIRE 783 SKELETONS TO BANC ###
####################################################
source("banc/banc-startup.R")
choose_segmentation("flywire31")
redo <- FALSE
version <- banc.nblast.version

# Make a version folder for the swc data
banc.nblast.fafb.swc.save.path.version <- file.path(banc.nblast.fafb.swc.save.path,version)
dir.create(banc.nblast.fafb.swc.save.path.version, recursive = TRUE, showWarnings = FALSE)

# Get fafb meta data
fw.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"flywire_meta.csv"), 
                                            col_types = hemibrainr:::sql_col_types))
fw.ids <- unique(fw.meta$root_783)

# Use the more detailed skeletons where we can
## Split and re-rooted: "/n/data1/hms/neurobio/wilson/fafbz//hemibrainr_data/flywire_neurons/783/split" 
# fw.swc.full.files <- list.files(flywire.split.save.path, pattern = ".swc$", full.names = TRUE)
fw.swc.full.files <- unique(unlist(lapply(fw.ids, function(id) file.path(flywire.split.save.path,paste0(id,".swc")))))
fw.swc.full.files <- fw.swc.full.files[file.exists(fw.swc.full.files)]
fw.swc.ids <- gsub(pattern = "\\.swc$", "", basename(fw.swc.full.files))

# What SWC files have we not downloaded?
if(redo){
  fw.swc.todo.files <- fw.ids
}else{
  fw.swc.todo.files <- setdiff(fw.ids, fw.swc.ids)
}
message("There are ", length(fw.swc.todo.files), " todo, we have downloaded ", length(fw.swc.ids), " FAFB-FlyWire .swc files")

# Get the leftover from the l2 files
fw.l2.full.files <- unique(unlist(lapply(fw.swc.todo.files, function(id) file.path(flywire.l2.save.path,paste0(id,".swc")))))
fw.l2.full.files <- fw.l2.full.files[file.exists(fw.l2.full.files)]
fw.swc.full.files <- c(fw.swc.full.files,fw.l2.full.files)

# Sort save directory
dir.create(banc.nblast.fafb.swc.save.path.version, recursive = TRUE, showWarnings = FALSE)
if(!redo){
  done <- list.files(banc.nblast.fafb.swc.save.path.version)
  fw.swc.full.files <- fw.swc.full.files[!basename(fw.swc.full.files)%in%done]
}
message("There are ", length(fw.swc.full.files), " SWC to transform to BANC space")

# Register cores
cl <- setup_parallel()

# Transforming skeletons
## There is some issue using parallelising the xform here?
message("##### Transforming skeletons #####")
by.query <- foreach(swc.file = sample(fw.swc.full.files),
                    .combine = 'c',
                    .packages = c('nat.jrcbrains','Morpho','rJava'),
                    .errorhandling = 'pass') %do% {
                      
                      # File
                      new.file <- file.path(banc.nblast.fafb.swc.save.path.version, basename(swc.file))
                      if(file.exists(new.file)&!redo){
                        next
                      }
                      
                      # Get neuron
                      neuron <- nat::read.neuron(swc.file)
                      
                      # If an L2 neuron, needs to be re-rooted
                      if(swc.file%in%fw.l2.full.files){
                        id <- gsub("\\.swc","",basename(swc.file))
                        pos <- subset(fw.meta, fw.meta$root_783==id)[1,]
                        if(!is.na(pos$soma_x)&pos$soma_x!=""){
                          soma <- nat::xyzmatrix(pos[,c("soma_x","soma_y","soma_z")])
                        }else{
                          soma <- nat::xyzmatrix(pos[,c("pos_x","pos_y","pos_z")])
                        }
                        if(nrow(soma)&!is.na(soma[1])){
                          neuron <- nat::reroot(x = neuron, point = c(soma))
                          neuron$tags$soma <- nat::rootpoints(neuron) 
                        }
                      }
                      
                      # Transform to JRC2018F
                      neuron.jrc2018f <- nat.templatebrains::xform_brain(neuron, sample = "FAFB14", reference = "JRC2018F")
                      
                      # Mirror, because of the FAFB-flip
                      neuron.jrc2018f <- nat.templatebrains::mirror_brain(neuron.jrc2018f, 
                                                                          brain = nat.flybrains::JRC2018F, 
                                                                          transform = "flip")
                      
                      # Transform to BANC
                      neuron.banc <- bancr::banc_to_JRC2018F(neuron.jrc2018f, method = "tpsreg", inverse = TRUE)
                      
                      # Save as a .swc file
                      nat::write.neuron(neuron.banc, file = new.file, Force = redo)
                      
                      # Return nothing
                      message("Transformed: ", new.file)
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
done <- list.files(banc.nblast.fafb.swc.save.path.version)
message("##### BANCpipeline: fafb-to-banc skeleton transforms calculated #####")
message(sprintf("##### we new transformed skeletons: %s of %s neurons", nrow(fw.swc.full.files)), length(done))
