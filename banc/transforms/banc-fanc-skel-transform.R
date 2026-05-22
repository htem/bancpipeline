#' banc-fanc-skel-transform — Transform FANC neuron skeletons into BANC space.
#'
#' For each FANC v1.116 root_id, fetches the skeleton via `fancr`/neuprint
#' and applies the FANC → JRCVNC2018F → BANC registration chain. Saves
#' the transformed SWC to the per-version BANC NBLAST skel tree.
#'
#' @section Reads:
#'   - `<banc.meta.save.path>/fanc_meta.csv`
#'   - FANC skeletons via `fancr`
#'
#' @section Writes:
#'   - `<banc.nblast.fanc.swc.save.path>/<version>/<root_id>.swc`
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_nblast_dataset.sh`

########################################
### TRANSFORM FANC SKELETONS TO BANC ###
########################################
source("banc/banc-startup.R")
library(fancr)
nat.jrcbrains::register_saalfeldlab_registrations()
version <- banc.nblast.version
redo <- FALSE

# Create save folder
banc.nblast.fanc.swc.save.path.version <- file.path(banc.nblast.fanc.swc.save.path,version)
dir.create(banc.nblast.fanc.swc.save.path.version, recursive = TRUE, showWarnings = FALSE)

# Get fanc meta data
fc.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"fanc_meta.csv"), 
                                            col_types = hemibrainr:::sql_col_types))
fc.ids <- unique(as.character(fc.meta$root_id))

# Stable skeletons are just fetched from neuprint, so we only save the transforms
fc.swc.full.files <- unique(unlist(lapply(fc.ids, function(id) file.path(banc.nblast.fanc.swc.save.path.version,paste0(id,".swc")))))
fc.swc.full.files <- fc.swc.full.files[file.exists(fc.swc.full.files)]
fc.swc.ids <- gsub(pattern = "\\.swc$", "", basename(fc.swc.full.files))

# What SWC files have we not downloaded?
if(!redo){
  fc.swc.todo.files <- setdiff(fc.ids, fc.swc.ids)
}else{
  fc.swc.todo.files <- fc.ids
}
message("There are ", length(fc.swc.todo.files), 
        " todo, we have downloaded and transformed ",
        length(fc.swc.ids), " fanc .swc files")

# Register cores
cl <- setup_parallel()

# Transforming skeletons
message("##### Transforming skeletons #####")
fanc.voxdims <- fancr::fanc_voxdims()
by.query <- foreach(swc.file = fc.swc.todo.files,
                    .combine = 'c',
                    .packages = c('nat.jrcbrains','Morpho','rJava','fancr'),
                    .errorhandling = 'pass') %do% {
                      
                      try({
                        # File
                        new.file <- file.path(banc.nblast.fanc.swc.save.path.version, paste0(basename(swc.file),".swc"))
                        id <- gsub("\\.swc","",basename(swc.file))
                        
                        # Get neuron. Note: Microns will fit it into MANC.surf
                        neuron <- fanc_read_l2skel(id)[[1]]
                        
                        # May need to reroot, based on 'raw' points.
                        pos <- subset(fc.meta, fc.meta$root_id==id)[1,]
                        soma <- fanc_raw2nm(nat::xyzmatrix(pos$root_position)) 
                        if(nrow(soma)&!is.na(soma[1])){
                          neuron <- nat::reroot(x = neuron, point = c(soma))
                          neuron$tags$soma <- nat::rootpoints(neuron) 
                        }
                        
                        # Change to MANC space
                        xyz.nm <- nat::xyzmatrix(neuron)
                        xyz.manc <- transform_fanc2manc(xyz.nm)/1000
                        nat::xyzmatrix(neuron) <- xyz.manc
                        
                        # Transform to JRCVNC2018F
                        neuron.jrcvnc2018f <- nat.templatebrains::xform_brain(neuron, 
                                                                              sample = "MANC",
                                                                              reference = "JRCVNC2018F")
                        
                        # Transform to BANC
                        neuron.banc <- bancr::banc_to_JRC2018F(neuron.jrcvnc2018f, 
                                                               region="vnc", 
                                                               method="tpsreg", 
                                                               banc.units = "nm", 
                                                               inverse = TRUE)
                        
                        # Save as a .swc file
                        nat::write.neuron(neuron.banc, file = new.file, Force = TRUE)
                        
                        # Return nothing
                        NULL
                      })
                      
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
done <- list.files(banc.nblast.fanc.obj.save.path)
message("##### BANCpipeline: fanc-to-banc skeleton transforms calculated #####")
message(sprintf("##### we new transformed skeletons: %s of %s neurons", nrow(fc.swc.full.files)), length(done))