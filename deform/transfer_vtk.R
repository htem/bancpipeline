##############################################
### ORGANISE DEFORMETRICA RESULTS MATCHING ###
##############################################
old <- getwd()
setwd("/home/ab714/bancpipeline/")
source("banc/banc-startup.R")
library(fs)
library(purrr)
library(dplyr)
library(processx)

#################
### Functions ###
#################

# Function to sync files
sync_files <- function(src, dest, 
                       move.old = TRUE, 
                       extensions = c("png","csv")) {
  
  # Change permissions on source files
  system2("find", c(src, "-type", "f", "-name", "*.png", "-o", "-name", "*.csv", "-exec", "chmod", "644", "{}", "+"))
  system2("find", c(src, "-type", "d", "-exec", "chmod", "755", "{}", "+"))
  
  # Construct the rsync command
  include_patterns <- paste(paste0("--include='*.", extensions, "'"), collapse = " ")
  rsync_cmd <- paste(
    "ssh $(whoami)@transfer.rc.hms.harvard.edu rsync",
    "-rlptv --recursive",
    "--copy-links",
    "--keep-dirlinks",  # Add this option to preserve relative paths for symlinks
    "--include='*/'",
    include_patterns,
    "--exclude='*/done/*'",
    "--exclude='*/old/*'",
    "--exclude='old/*'",
    "--exclude='*/express/*'",
    "--exclude='*'",
    "--prune-empty-dirs",
    ifelse(move.old,"--backup",""),
    ifelse(move.old,sprintf("--backup-dir='%s/%s%s/'",dest,"old_",Sys.Date()),""),
    "--delete",
    shQuote(paste0(src,"/")),
    shQuote(paste0(dest,"/"))
  )
  cat("command: ",rsync_cmd)
  
  # Run rsync using system()
  result <- system(rsync_cmd, intern = TRUE)
  
  # # Remove empty directories in the source
  # remove_empty_dirs_cmd <- paste(
  #   "ssh $(whoami)@transfer.rc.hms.harvard.edu",
  #   "find", paste0(shQuote(dest),"/"),
  #   "-type d -empty -delete"
  # )
  # result2 <- system(remove_empty_dirs_cmd)
  
  # report
  cat(result)
  cat("Synced files from", src, "to", dest, "\n")
  cat(paste(result, collapse = "\n"), "\n")
}

############
### Main ###
############

# Remote paths
cat("Synchronising matching files between O2 and the fileserve")
A <- '/n/data1/hms/neurobio/wilson/banc/deformetrica/'
B <- '/n/files/Neurobio/wilsonlab/banc/deformetrica/'

# Remote synchronization of meta data
sync_files(path(A, "fafb/"), path(B, "fafb/"), extensions = c("vtk", "xml"))
setwd(old)