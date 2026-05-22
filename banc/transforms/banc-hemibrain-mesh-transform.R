#' banc-hemibrain-mesh-transform — Transform hemibrain meshes into BANC space.
#'
#' For each hemibrain bodyid, fetches the OBJ mesh and applies the
#' JRCFIB2018F → JRC2018F → BANC registration chain. Saves the
#' transformed OBJ.
#'
#' @section Reads:
#'   - `<banc.meta.save.path>/hemibrain_meta.csv`
#'   - hemibrain meshes from `<hemibrain.save.path>/obj`
#'
#' @section Writes:
#'   - `<banc.nblast.hemibrain.obj.save.path>/<version>/<bodyid>.obj`
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_ngl.sh`

###########################################
### TRANSFORM HEMIBRAIN NEURONS TO BANC ###
###########################################
source("banc/banc-startup.R")
choose_segmentation("flywire31")
redo <- FALSE
version <- banc.nblast.version

# Make a version folder for the swc data
banc.nblast.hemibrain.obj.save.path.version <- file.path(banc.nblast.hemibrain.obj.save.path,version)
dir.create(banc.nblast.hemibrain.obj.save.path.version, recursive = TRUE, showWarnings = FALSE)

# Get meta data
hb.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"hemibrain_meta.csv"), 
                                            col_types = hemibrainr:::sql_col_types))
hb.ids <- unique(hb.meta$bodyid)

# Count the number of 
hb.obj.full.files <- list.files(banc.nblast.hemibrain.obj.save.path.version, pattern = ".obj$", full.names = TRUE)
hb.obj.ids <- gsub(pattern = "\\.obj$", "", basename(hb.obj.full.files))

# What obj files have we not downloaded?
if(redo){
  hb.obj.todo.files <- hb.ids
}else{
  hb.obj.todo.files <- setdiff(hb.ids, hb.obj.ids)
}
message("There are ", length(hb.obj.todo.files), " todo, we have downloaded ", length(hb.obj.ids), " BANC .obj files")

# Transforming meshes
## There is some issue using paralleling the xform here?
message("##### Transforming meshes #####")
by.query <- foreach(id = sample(hb.obj.todo.files),
                    .combine = 'c',
                    .packages = c('nat.jrcbrains','Morpho','rJava'),
                    .errorhandling = 'pass') %dopar% {
                      
                      # File
                      new.file <- file.path(banc.nblast.hemibrain.obj.save.path.version, paste0(id,".obj"))
                      if(file.exists(new.file)&!redo){
                        message("Skipping written file: ", new.file)
                        NULL
                      }else{
                        # Get mesh3d
                        mesh3d <- hemibrainr::hemibrain_neuron_meshes(id, 
                                                                      dataset="hemibrain:v1.2.1",
                                                                      cloudvolume.url="precomputed://gs://neuroglancer-janelia-flyem-hemibrain/v1.2/segmentation") 
                        
                        # Transforms to JRC2018F, divide by 1000 to reach microns for JRCFIB2018F
                        mesh3d.jrc2018f <- xform_brain(mesh3d/1000, sample = "JRCFIB2018F", reference = "JRC2018F")
                        
                        # Transform to BANC
                        mesh3d.banc <- bancr::banc_to_JRC2018F(mesh3d.jrc2018f, method = "tpsreg", inverse = TRUE)
                        
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
done <- list.files(banc.nblast.hemibrain.obj.save.path.version)
message("##### BANCpipeline: Hemibrain meshes in BANC space, updated #####")
message(sprintf("##### hemibrain meshes transformed .obj files: %s", length(done)))
message(sprintf("##### we have new transformed meshes: %s of %s, neurons", nrow(hb.obj.full.files)), length(done))
