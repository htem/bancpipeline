#########################################
### Calculate accurate BANC skeletons ###
#########################################
source("banc/banc-startup.R")

# Direct us to the BANC dataset
choose_segmentation("flywire31") # must do this due to bug in skeletor R code?

# Read IDs
banc.meta <- banctable_query() %>%
  dplyr::filter(!grepl("GLIA|TRACHEA|NOT_A_NEURON|DELETE|DEBRIS",status)) %>%
  dplyr::filter(proofread=="TRUE"|roughly_proofread=="TRUE")
banc.root.ids <- unique(banc.meta$root_id)

# Get IDs for downloaded meshes not in SWC folder
obj.full.files <- list.files(banc.obj.save.path, pattern = ".obj$", full.names = TRUE)
obj.ids <- gsub(pattern = "\\.obj$", "", basename(obj.full.files))
obj.full.files <- obj.full.files[obj.ids%in%banc.root.ids]
obj.ids <- intersect(obj.ids,banc.root.ids)

# How many obj files not converted to accurate skeletons?
swc.files = gsub(pattern = "datahmneurobiwilsobanob|\\.swc$", "", list.files(banc.swc.save.path, pattern = ".swc$", full.names = FALSE))
todo.files = sample(obj.full.files[!obj.ids%in%swc.files])
message("There are ", length(todo.files), " of ", length(obj.full.files), " .obj files not converted to accurate skeletons")

# Have these been tried and failed before?
failed <- tryCatch(read_csv(file.path(banc.save.path,"banc_failed.csv"), col_types=cols("c","c")), error = function(e) NULL)
if(!is.null(failed)){
  failed = subset(failed, file.missing == "swc")
  fails = as.character(failed$root_id)
  todo.ids = gsub("\\.obj","",basename(todo.files))
  todo.good.files = todo.files[!todo.ids %in% fails]
}else{
  todo.good.files = todo.files
}
obj.not.failed.before <- length(todo.good.files)
message("There are ", obj.not.failed.before, " .obj that have not failed skeletonisation before. Attempting now...")
if(weekdays(Sys.time()) == "Saturday"){ # On Saturdays,re-try failures
  message("But on Saturdays, we also retry some of our past failures... ", length(fails)," to be exact")
  todo.good.files = todo.files
}
message("There are ", length(todo.good.files), " of ", length(todo.files), " .obj files have not been attempted before")

# Skeletonise
if(length(todo.good.files)){
  message("atttempting to skeletonise ", length(todo.good.files)," .obj files...")
  multiplier <- 1000
  upper <- ifelse((numCores*multiplier)<length(todo.good.files),numCores*multiplier,length(todo.good.files))
  batches <- split(sample(todo.good.files), round(seq(from = 1, to = upper, length.out = length(todo.good.files))))
  
  # Register cores
  cl <- setup_parallel()
  
  # Set up progress bar
  iterations <- length(batches)
  pb <- utils::txtProgressBar(max = iterations, style = 3)
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # Run foreach loop
  by.query <- foreach::foreach (batch = seq_along(batches),
                                .combine = 'c',
                                .errorhandling='pass',
                                .options.snow = opts) %dopar% {
    neuron.ids <- batches[[batch]]
    for(neuron.id in neuron.ids){
      try({
        skels <- fafbseg::skeletor(obj = neuron.id, 
                                   #cloudvolume.url = bancr:::banc_cloudvolume_url(),
                                   clean = TRUE,
                                   method = "wavefront",
                                   save.obj = NULL, 
                                   mesh3d = FALSE,
                                   waves = 1,
                                   k.soma.search = 100,
                                   radius.soma.search = 10000,
                                   heal = TRUE,
                                   reroot = TRUE,
                                   heal.threshold = 1000000,
                                   heal.k = 10L,
                                   reroot_method = "density",
                                   brain = bancr::banc_neuropil.surf,
                                   elapsed = 10000,
                                   resample = 1000)
        skels[,"id"] <- names(skels) <- basename(gsub("\\.obj|\\/datahmneurobiwilsobanob","",names(skels)))
        nat::write.neurons(skels, dir=banc.swc.save.path, format='swc', Force = TRUE)
        # message("completed: ", length(skels), " skeletonisations")
        # length(skels) 
      })
    }
    NULL
  }
  
  # Stop cores
  stop_parallel(cl)
  
  # Were there errors?
  for(i in 1:length(by.query)){
    if(!is.null(by.query[[i]])){
      message(by.query[[i]])
    }
  }
}

# Work out which neurons we just wrote
swc.files.new <- gsub(pattern = "\\.swc$", "", list.files(banc.swc.save.path, pattern = ".swc$", full.names = FALSE))
swc.files.new <- setdiff(swc.files.new, swc.files)
new.swc <- length(swc.files.new)

# Record which IDs failed, do not re-attempt
successful <- gsub("\\.swc","",basename(swc.files))
new.ids <- gsub("\\.swc","",basename(swc.files.new))
failed <- data.frame(stringsAsFactors = FALSE)
fw.ids.tried <- gsub(".obj","",basename(todo.good.files))
fail.obj <- setdiff(fw.ids.tried,obj.ids)
if(length(fail.obj)){
  failed.obj <- data.frame(root_id = fail.obj,
                           file.missing = "obj")
  failed <- rbind(failed,failed.obj)
}
fail.swc <- setdiff(fw.ids.tried,successful)
if(length(fail.swc)){
  failed.swc <- data.frame(root_id = fail.swc,
                           file.missing = "swc")
  failed <- rbind(failed,failed.swc)
}
if(nrow(failed)){
  failed$root_id <- as.character(failed$root_id)
  write.csv(failed,file <- file.path(banc.save.path,"banc_failed.csv"), row.names = FALSE)
}

# Announce
message("##### BANCpipeline: new banc skeleton files calculated #####")
message("succesfully skeletonised ", new.swc," .obj files out of ", length(todo.good.files))

# Rename files (unsure why I need to)
files_to_rename = list.files(banc.swc.save.path, pattern = ".swc$", full.names = TRUE)
new_names = gsub("datahmneurobiwilsobanob","",files_to_rename)
file.rename(from = files_to_rename, to = new_names)








