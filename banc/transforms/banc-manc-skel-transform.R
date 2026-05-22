#' banc-manc-skel-transform — Transform MANC neuron skeletons into BANC space.
#'
#' For each MANC bodyid, fetches the skeleton via `malevnc::read_manc_skeletons`
#' and applies the MANC → JRCVNC2018F → BANC registration chain. Saves
#' the transformed SWC to the per-version BANC NBLAST skel tree.
#'
#' @section Reads:
#'   - `<banc.meta.save.path>/manc_meta.csv`
#'   - MANC skeletons via `malevnc`
#'
#' @section Writes:
#'   - `<banc.nblast.manc.swc.save.path>/<version>/<bodyid>.swc`

########################################
### TRANSFORM MANC SKELETONS TO BANC ###
########################################
source("banc/banc-startup.R")
#malevnc::download_manc_registrations()
nat.jrcbrains::register_saalfeldlab_registrations()
malevnc:::choose_malevnc_dataset('MANC')
version <- banc.nblast.version
redo <- FALSE

# Create save folder
banc.nblast.manc.swc.save.path.version <- file.path(banc.nblast.manc.swc.save.path,version)
dir.create(banc.nblast.manc.swc.save.path.version, recursive = TRUE, showWarnings = FALSE)

# Get manc meta data
mc.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"manc_meta.csv"), 
                                            col_types = hemibrainr:::sql_col_types))
mc.ids <- unique(mc.meta$bodyid)
mc.positions <- mc.meta %>%
  dplyr::mutate(root_pos = dplyr::case_when(
    !is.na(soma_location)&soma_location!="" ~ soma_location,
    !is.na(tosoma_location)&tosoma_location!="" ~ tosoma_location,
    !is.na(root_location)&root_location!="" ~ root_location,
    TRUE ~ ""
  )) %>%
  dplyr::mutate(root_pos=gsub("list\\(|\\).*","",root_pos)) %>%
  dplyr::distinct(bodyid, root_pos, .keep_all = TRUE)

# Stable skeletons are just fetched from neuprint, so we only save the transforms
mc.swc.full.files <- unique(unlist(lapply(mc.ids, function(id) file.path(banc.nblast.manc.swc.save.path.version,paste0(id,".swc")))))
mc.swc.full.files <- mc.swc.full.files[file.exists(mc.swc.full.files)]
mc.swc.ids <- gsub(pattern = "\\.swc$", "", basename(mc.swc.full.files))

# What SWC files have we not downloaded?
if(!redo){
  mc.swc.todo.files <- setdiff(mc.ids, mc.swc.ids)
}else{
  mc.swc.todo.files <- mc.ids
}
message("There are ", length(mc.swc.todo.files), " todo, we have downloaded and transformed ", length(mc.swc.ids), " MANC .swc files")

# Register cores
cl <- setup_parallel()

# Transforming skeletons
message("##### Transforming skeletons #####")
by.query <- foreach(swc.file = mc.swc.todo.files,
                    .combine = 'c',
                    .packages = c('nat.jrcbrains','Morpho','rJava'),
                    .errorhandling = 'pass') %dopar% {
                      
                      # File
                      new.file <- file.path(banc.nblast.manc.swc.save.path.version, paste0(basename(swc.file),".swc"))
                      
                      # Get neuron. Note: Microns will fit it into MANC.surf
                      neuron <- manc_read_neurons(ids = swc.file, units="nm")[[1]]
                      
                      # May need to reroot, based on 'raw' points.
                      id <- gsub("\\.swc","",basename(swc.file))
                      pos <- subset(mc.positions, mc.positions$bodyid==id)[1,]
                      soma <- nat::xyzmatrix(pos$root_pos) # in raw coordinates
                      if(nrow(soma)&!is.na(soma[1])){
                        soma <- soma*rep(8, 3)
                        neuron <- nat::reroot(x = neuron, point = c(soma))
                        neuron$tags$soma <- nat::rootpoints(neuron) 
                      }
                      
                      # Change to nm
                      # nat::xyzmatrix(neuron) <- nat::xyzmatrix(neuron)*rep(8, 4)
  
                      # Transform to JRCVNC2018F
                      neuron.jrcvnc2018f <- nat.templatebrains::xform_brain(neuron/1000, sample = "MANC", reference = "JRCVNC2018F")
                      
                      # Transform to BANC
                      neuron.banc <- bancr::banc_to_JRC2018F(neuron.jrcvnc2018f, region="vnc", method="tpsreg", banc.units = "nm", inverse = TRUE)
                      
                      # Save as a .swc file
                      nat::write.neuron(neuron.banc, file = new.file, Force = TRUE)
                      
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
done <- list.files(banc.nblast.manc.obj.save.path)
message("##### BANCpipeline: manc-to-banc skeleton transforms calculated #####")
message(sprintf("##### we new transformed skeletons: %s of %s neurons", nrow(mc.swc.full.files)), length(done))