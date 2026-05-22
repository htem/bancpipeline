#' banc-fafb-mesh-transform — Transform FAFB (FlyWire v783) meshes into BANC space.
#'
#' For each FAFB-v783 root_id, fetches the OBJ mesh and applies the
#' FlyWire-JRC2018F-BANC registration chain. Saves the transformed OBJ to
#' the per-version BANC NBLAST mesh tree. Resumable.
#'
#' @section Reads:
#'   - `<banc.meta.save.path>/flywire_meta.csv`
#'   - FAFB OBJ meshes from `<flywire.obj.save.path>`
#'
#' @section Writes:
#'   - `<banc.nblast.fafb.obj.save.path>/<version>/<root_783>.obj`
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_ngl.sh`

##################################################
### TRANSFORM FAFB-FLYWIRE 783 NEURONS TO BANC ###
##################################################
source("banc/banc-startup.R")
choose_segmentation("flywire31")
redo <- FALSE
version <- banc.nblast.version

# Make a version folder for the swc data
banc.nblast.fafb.obj.save.path.version <- file.path(banc.nblast.fafb.obj.save.path,version)
dir.create(banc.nblast.fafb.obj.save.path.version, recursive = TRUE, showWarnings = FALSE)

# Get meta data
fw.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"flywire_meta.csv"), 
                                            col_types = hemibrainr:::sql_col_types))
fw.ids <- unique(fw.meta$root_783)

# # Count the number of 
# ## Ideally: "/n/data1/hms/neurobio/wilson/fafbz//hemibrainr_data/flywire_neurons/783/obj"
# fw.obj.full.files <- list.files(flywire.obj.save.path, pattern = ".obj$", full.names = TRUE)
# fw.obj.ids <- gsub(pattern = "\\.obj$", "", basename(fw.obj.full.files))
# 
# # What obj files have we not downloaded?
# fw.obj.todo.files <- setdiff(fw.ids, fw.obj.ids)
# message("There are ", length(fw.obj.todo.files), " todo, we have downloaded ", length(fw.obj.ids), "  flywie 783 .obj files")
# 
# # Download any meshes from release 783 that we lack
# if(length(fw.obj.todo.files)){
#   message("downloading ", length(fw.obj.todo.files), " .obj files")
#   hemibrainr:::download_neuron_obj_batch(ids = sample(fw.obj.todo.files),
#                                          numCores = numCores,
#                                          save.obj = flywire.obj.save.path)
# }

# Recount the number of meshes
fw.obj.full.files <- unique(unlist(lapply(fw.ids, function(id) file.path(flywire.obj.save.path,paste0(id,".obj")))))
fw.obj.full.files <- fw.obj.full.files[file.exists(fw.obj.full.files)]
if(!redo){
  done <- list.files(banc.nblast.fafb.obj.save.path.version)
  fw.obj.full.files <- fw.obj.full.files[!basename(fw.obj.full.files)%in%done]
}
message("Atempting to transform ", length(fw.obj.full.files), " FAFB-FlyWire-783 meshes")

# Register cores
cl <- setup_parallel()

# Transforming meshes
message("##### Transforming meshes #####")
by.query <- foreach(obj.file = sample(fw.obj.full.files),
                      .combine = 'c',
                      .packages = c('nat.jrcbrains','Morpho','rJava'),
                      .errorhandling = 'pass') %do% {
                        
    # File
    new.file <- file.path(banc.nblast.fafb.obj.save.path.version,basename(obj.file))
    if(file.exists(new.file)&!redo){
      message("Skipping written file: ", new.file)
      NULL
    }else{
      
      # Get mesh3d
      mesh3d <- readobj::read.obj(obj.file, convert.rgl=TRUE)[[1]] 
      
      # Transform to JRC2018F
      mesh3d.jrc2018f <- nat.templatebrains::xform_brain(mesh3d, sample = "FAFB14", reference = "JRC2018F")
      
      # Mirror, because of the FAFB flip
      mesh3d.jrc2018f <- nat.templatebrains::mirror_brain(mesh3d.jrc2018f, brain = nat.flybrains::JRC2018F, transform = "flip")
      
      # Transform to BANC
      mesh3d.banc <- bancr::banc_to_JRC2018F(mesh3d.jrc2018f, method = "tpsreg", inverse = TRUE)
      
      # Simplify meshes
      # mesh3d.banc.simp <- Rvcg::vcgQEdecim(mesh3d.banc, percent = 0.2)
      
      # Save as a .obj file
      bancr:::write_mesh3d_to_obj(mesh3d.banc, filename = new.file)
      message("Written: ", new.file)
      
      # Return nothing
      NULL 
    }
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
done <- list.files(banc.nblast.fafb.obj.save.path.version)
message("##### BANCpipeline: FAFB meshes in BANC space, updated #####")
message(sprintf("##### FAFB-FlyWire transformed .obj files: %s", length(done)))
message(sprintf("##### we new transformed meshes: %s of %s, neurons", nrow(fw.obj.full.files)), length(done))


