#' banc-fanc-mesh-transform — Transform FANC neuron meshes into BANC space.
#'
#' For each FANC v1.116 root_id, fetches the OBJ mesh via `fancr` and
#' applies the FANC → JRCVNC2018F → BANC registration chain. Saves to the
#' per-version BANC NBLAST mesh tree.
#'
#' @section Reads:
#'   - `<banc.meta.save.path>/fanc_meta.csv`
#'   - FANC meshes via `fancr`
#'
#' @section Writes:
#'   - `<banc.nblast.fanc.obj.save.path>/<version>/<root_id>.obj`
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_nblast_dataset.sh`

######################################
### TRANSFORM FANC NEURONS TO BANC ###
######################################
source("banc/banc-startup.R")
library(fancr)
nat.jrcbrains::register_saalfeldlab_registrations()
malevnc:::choose_malevnc_dataset('MANC')
redo <- FALSE
version <- banc.nblast.version

# Make a version folder for the swc data
banc.nblast.fanc.obj.save.path.version <- file.path(banc.nblast.fanc.obj.save.path,version)
dir.create(banc.nblast.fanc.obj.save.path.version, recursive = TRUE, showWarnings = FALSE)

# Get meta data
fc.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"fanc_meta.csv"), 
                                            col_types = hemibrainr:::sql_col_types))
fc.ids <- unique(fc.meta$root_id)

# Count the number of 
fc.obj.full.files <- list.files(banc.nblast.fanc.obj.save.path.version, pattern = ".obj$", full.names = TRUE)
fc.obj.ids <- gsub(pattern = "\\.obj$", "", basename(fc.obj.full.files))

# What obj files have we not downloaded?
if(redo){
  fc.obj.todo.files <- fc.ids
}else{
  fc.obj.todo.files <- setdiff(fc.ids, fc.obj.ids)
}
message("There are ", length(fc.obj.todo.files), 
        " todo, we have downloaded ", length(fc.obj.ids), " FANC .obj files")

# Transforming meshes
## There is some issue using parallelising the xform here?
message("##### Transforming meshes #####")
by.query <- foreach(id = sample(fc.obj.todo.files),
                    .combine = 'c',
                    .packages = c('nat.jrcbrains','Morpho','rJava','fancr'),
                    .errorhandling = 'pass') %do% {
                      
                      # File
                      new.file <- file.path(banc.nblast.fanc.obj.save.path.version, paste0(id,".obj"))
                      if(file.exists(new.file)&!redo){
                        message("Skipping written file: ", new.file)
                        NULL
                      }else{
                        
                        # Get mesh3d
                        fanc.mesh <- read_fanc_meshes(id)
                        
                        # Change to MANC space
                        xyz.nm <- nat::xyzmatrix(fanc.mesh)
                        xyz.manc <- transform_fanc2manc(xyz.nm)/1000
                        nat::xyzmatrix(fanc.mesh) <- xyz.manc
                        
                        # Transform into JRCVNC2918F
                        fanc.mesh.jrcvnc2018f=xform_brain(fanc.mesh, 
                                                          reference = "JRCVNC2018F", 
                                                          sample="MANC")
                        
                        # Transform into the BANC
                        mesh3d.banc <- banc_to_JRC2018F(fanc.mesh.jrcvnc2018f, region="vnc", method="tpsreg", 
                                                        banc.units = "nm", inverse = TRUE)
                        
                        # Simplify meshes
                        # mesh3d.banc.simp <- Rvcg::vcgQEdecim(mesh3d.banc[[1]], percent = 0.2)
                        
                        # Save as a .obj file
                        bancr:::write_mesh3d_to_obj(mesh3d.banc[[1]], filename = new.file)
                        message("Written: ", new.file)
                        
                        # Return nothing
                        NULL
                      }
                    }                    

# Were there errors?
message("##### Displaying any errors from foreach loop #####")
for(i in 1:length(by.query)){
  if(!is.null(by.query[[i]])){
    message(by.query[[i]])
  }
}

# Announce
done <- list.files(banc.nblast.fanc.obj.save.path.version)
message("##### BANCpipeline: fanc meshes in BANC space, updated #####")
message(sprintf("##### fanc meshes transformed .obj files: %s", length(done)))
message(sprintf("##### we have new transformed meshes: %s of %s, neurons", nrow(fc.obj.full.files)), length(done))

