#' banc-obj — Download per-neuron OBJ meshes for proofread BANC neurons.
#'
#' @section Reads:
#'   - SeaTable `banc_meta`
#'
#' @section Writes:
#'   - `<banc.obj.save.path>/<root_888>.obj`
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_meshes.sh`, `o2/production/o2_banc_split.sh`.

##################################
### Save OBJ files for neurons ###
##################################
source("banc/banc-startup.R")
redo <- FALSE

# Direct us to the BANC dataset
bancr::choose_banc()
numCores <- 10

# Read IDs
banc.meta <- banctable_query() %>%
  banc_filter_neurons() %>%
  dplyr::filter(proofread=="TRUE"|roughly_proofread=="TRUE")
banc.root.ids <- na.omit(unique(banc.meta$root_888))

# Get IDs for downloaded meshes not in SWC folder
dir.create(banc.obj.save.path, showWarnings = FALSE, recursive = TRUE)
obj.todo.files <- file.path(banc.obj.save.path, paste0(banc.root.ids,".obj"))
obj.full.files <- list.files(banc.obj.save.path, pattern = ".obj$", full.names = TRUE)
obj.ids <- gsub(pattern = "\\.obj$", "", basename(obj.full.files))
obj.todo.files <- gsub(pattern = "\\.obj$", "", basename(obj.todo.files))

# What obj files have we not downloaded?
if(redo){
  obj.todo.files <- obj.todo.files
}else{
  obj.todo.files <- setdiff(obj.todo.files,obj.ids)
}  
message("There are ", length(obj.todo.files), " todo, we have downloaded ", length(obj.ids), "  BANC .obj files")

# Download new obj files
message("##### Downloading missing BANC .obj files #####")
by.query <- foreach(obj.file = sample(obj.todo.files),
                    .combine = 'c',
                    .errorhandling = 'pass') %do% {
                      
    # File
    new.file <- file.path(banc.obj.save.path, paste0(basename(obj.file),".obj"))
    if(file.exists(new.file)&!redo){
      message("Skipping written file: ", new.file)
      NULL
    }else{
                          
      # Get mesh3d
      mesh3d.banc <- banc_read_neuron_meshes(obj.file)
      
      # Simplify meshes
      # mesh3d.banc.simp <- Rvcg::vcgQEdecim(mesh3d.banc, percent = 0.2)
      
      # Save as a .obj file
      bancr:::write_mesh3d_to_obj(mesh3d.banc[[1]], filename = new.file)
      message("Written: ", new.file)
      # Return nothing
      invisible() 
  }
}

# Calculate what is new
obj.full.files.new <- list.files(banc.obj.save.path, pattern = ".obj$", full.names = TRUE)
obj.full.files.new <- length(obj.full.files.new)-length(obj.full.files)

# Announce
message("##### BANCpipeline: new banc .obj files downloaded #####")
message(sprintf("##### we have metrics for : %s neurons", obj.full.files.new))


