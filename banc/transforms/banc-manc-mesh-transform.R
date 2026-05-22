#' banc-manc-mesh-transform — Transform MANC meshes into BANC space.
#'
#' For each MANC bodyid, fetches the OBJ mesh and applies the
#' MANC → JRCVNC2018F → BANC registration chain. Saves the transformed OBJ.
#'
#' @section Reads:
#'   - `<banc.meta.save.path>/manc_meta.csv`
#'   - MANC meshes via `malevnc`
#'
#' @section Writes:
#'   - `<banc.nblast.manc.obj.save.path>/<version>/<bodyid>.obj`
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_ngl.sh`

######################################
### TRANSFORM MANC NEURONS TO BANC ###
######################################
source("banc/banc-startup.R")
nat.jrcbrains::register_saalfeldlab_registrations()
malevnc:::choose_malevnc_dataset('MANC')
redo <- FALSE
version <- banc.nblast.version

# Make a version folder for the swc data
banc.nblast.manc.obj.save.path.version <- file.path(banc.nblast.manc.obj.save.path,version)
dir.create(banc.nblast.manc.obj.save.path.version, recursive = TRUE, showWarnings = FALSE)

# Get meta data
mc.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"manc_meta.csv"), 
                                            col_types = hemibrainr:::sql_col_types))
mc.ids <- unique(mc.meta$bodyid)

# Count the number of 
mc.obj.full.files <- list.files(banc.nblast.manc.obj.save.path.version, pattern = ".obj$", full.names = TRUE)
mc.obj.ids <- gsub(pattern = "\\.obj$", "", basename(mc.obj.full.files))

# What obj files have we not downloaded?
if(redo){
  mc.obj.todo.files <- mc.ids
}else{
  mc.obj.todo.files <- setdiff(mc.ids, mc.obj.ids)
}
message("There are ", length(mc.obj.todo.files), " todo, we have downloaded ", length(mc.obj.ids), " MANC .obj files")

# Transforming meshes
## There is some issue using parallelising the xform here?
message("##### Transforming meshes #####")
by.query <- foreach(id = sample(mc.obj.todo.files),
                    .combine = 'c',
                    .packages = c('nat.jrcbrains','Morpho','rJava'),
                    .errorhandling = 'pass') %dopar% {
                      
                      # File
                      new.file <- file.path(banc.nblast.manc.obj.save.path.version, paste0(id,".obj"))
                      if(file.exists(new.file)&!redo){
                        message("Skipping written file: ", new.file)
                        NULL
                      }else{
                        
                        # Get mesh3d
                        manc.mesh <- read_manc_meshes(id)
                        
                        # Transform into JRCVNC2918F
                        manc.mesh.jrcvnc2018f=xform_brain(manc.mesh/1e3, reference = "JRCVNC2018F", sample="MANC")
                        
                        # Transform into the BANC
                        mesh3d.banc <- banc_to_JRC2018F(manc.mesh.jrcvnc2018f, region="vnc", method="tpsreg", banc.units = "nm", inverse = TRUE)
                        
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
done <- list.files(banc.nblast.manc.obj.save.path.version)
message("##### BANCpipeline: manc meshes in BANC space, updated #####")
message(sprintf("##### manc meshes transformed .obj files: %s", length(done)))
message(sprintf("##### we have new transformed meshes: %s of %s, neurons", nrow(mc.obj.full.files)), length(done))

