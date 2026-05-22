#!/usr/bin/env Rscript
source("banc/banc-startup.R")

# Download SWC files
out_dir <- "/n/data1/hms/neurobio/wilson/banc/matching/malecns/JRC2018U"
base::dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
gcs_path <- "gs://flyem-male-cns/v0.9/segmentation/skeletons-unisex-template/skeletons-swc/"

# Build and run gsutil command from R
cmd <- base::sprintf(
  "gsutil -m cp -r %s %s",
  base::shQuote(gcs_path),
  base::shQuote(out_dir)
)
base::message("Running: ", cmd)
status <- base::system(cmd)
if (status != 0L) {
  base::stop("gsutil command failed with status ", status)
} else {
  base::message("Download complete.")
}

# Download SWC files
out_dir2 <- "/n/data1/hms/neurobio/wilson/banc/matching/malecns/JRCFIB2022M"
base::dir.create(out_dir2, recursive = TRUE, showWarnings = FALSE)
gcs_path <- "gs://flyem-male-cns/v0.9/segmentation/skeletons-malecns/skeletons-swc/"

# Build and run gsutil command from R
cmd <- base::sprintf(
  "gsutil -m cp -r %s %s",
  base::shQuote(gcs_path),
  base::shQuote(out_dir2)
)
base::message("Running: ", cmd)
status <- base::system(cmd)
if (status != 0L) {
  base::stop("gsutil command failed with status ", status)
} else {
  base::message("Download complete.")
}

# Convert to BANC
banc_dir <- "/n/data1/hms/neurobio/wilson/banc/matching/malecns/banc_space_split/swc"
done.swc <- list.files(banc_dir)
swcs <- list.files(file.path(out_dir, "skeletons-swc"), full.names = TRUE)
swcs <- sample(swcs[!gsub("\\.swc", "", basename(swcs)) %in%
                      gsub("\\.swc", "", basename(done.swc))])
batch_size <- 1000L
idx <- seq_along(swcs)
batches <- split(idx, ceiling(idx / batch_size))
pb <- utils::txtProgressBar(min = 0, max = length(swcs), style = 3)
processed <- 0L
for (batch in batches) {
  for (i in batch) {
    swc <- swcs[i]
    try(
      suppressMessages({
        neuron <- nat::read.neuron(swc)
        in.brain <- any(
          nat::pointsinside(
            x    = nat::xyzmatrix(neuron),
            surf = JRC2018U.surf
          )
        )
        
        if (in.brain) {
          neuronf <- xform_brain(
            neuron,
            sample    = "JRC2018U",
            reference = "JRC2018F"
          )
          neuronb <- bancr::banc_to_JRC2018F(
            neuronf,
            region  = "brain",
            inverse = TRUE
          )
        } else {
          neuronf <- xform_brain(
            neuron,
            sample    = "JRCVNC2018U",
            reference = "JRCVNC2018F"
          )
          neuronb <- bancr::banc_to_JRC2018F(
            neuronf,
            region  = "vnc",
            inverse = TRUE
          )
        }
        
        nat::write.neuron(neuronb, dir = banc_dir, format = "swc", Force = TRUE)
      }),
      silent = TRUE
    )
    
    processed <- processed + 1L
    utils::setTxtProgressBar(pb, processed)
  }
  gc()
}
close(pb)

########################################################
### TRANSFORM maleCNS NATIVE SKELETONS (JRCFIB2022M) ###
### TO BANC (nm) USING PYTHON navis / flybrains      ###
########################################################
#
# The JRCFIB2022M-space SWC files are in 8 nm voxel coordinates.
# We multiply by 8 to get nm, then use the Python TPS transform
# (navis.xform_brain, source="JRCFIB2022M", target="BANC") since
# no R-native JRCFIB2022M → BANC registration exists.
#
# reticulate/Python calls are not fork-safe, so we run sequentially.

library(reticulate)
navis     <- import("navis", convert = FALSE)
flybrains <- import("flybrains")
np        <- import("numpy", convert = FALSE)
message("Python transform backend ready (navis + flybrains)")

version <- banc.nblast.malecns.version
banc.nblast.malecns.swc.save.path.version <- file.path(banc.nblast.malecns.swc.save.path, version)
dir.create(banc.nblast.malecns.swc.save.path.version, recursive = TRUE, showWarnings = FALSE)

# Find native-space SWC files (downloaded above to out_dir2)
mcns.swcs <- list.files(file.path(out_dir2, "skeletons-swc"), pattern = "\\.swc$", full.names = TRUE)

# Skip already-transformed files
done.banc.swc <- list.files(banc.nblast.malecns.swc.save.path.version, pattern = "\\.swc$")
mcns.swcs <- mcns.swcs[!gsub("\\.swc$", "", basename(mcns.swcs)) %in%
                          gsub("\\.swc$", "", done.banc.swc)]
mcns.swcs <- sample(mcns.swcs)
message("maleCNS native → BANC: ", length(mcns.swcs), " skeletons to transform (",
        length(done.banc.swc), " already done)")

# Transform sequentially (reticulate is not fork-safe)
pb2 <- utils::txtProgressBar(min = 0, max = length(mcns.swcs), style = 3)
processed2 <- 0L
batch_size2 <- 1000L
idx2 <- seq_along(mcns.swcs)
batches2 <- split(idx2, ceiling(idx2 / batch_size2))
for (batch in batches2) {
  for (i in batch) {
    swc <- mcns.swcs[i]
    try(
      suppressMessages({
        neuron <- nat::read.neuron(swc)

        # Convert from 8 nm voxel coordinates to nm
        nat::xyzmatrix(neuron) <- nat::xyzmatrix(neuron) * 8

        # Transform JRCFIB2022M (nm) → BANC (nm) via Python TPS
        coords.nm <- nat::xyzmatrix(neuron)
        coords.py <- np$array(r_to_py(coords.nm))
        coords.banc.py <- navis$xform_brain(coords.py, source = "JRCFIB2022M", target = "BANC")
        coords.banc <- py_to_r(coords.banc.py)
        nat::xyzmatrix(neuron) <- coords.banc

        # Save transformed skeleton
        new.file <- file.path(banc.nblast.malecns.swc.save.path.version, basename(swc))
        nat::write.neuron(neuron, file = new.file, Force = TRUE)
      }),
      silent = TRUE
    )
    processed2 <- processed2 + 1L
    utils::setTxtProgressBar(pb2, processed2)
  }
  gc()
}
close(pb2)

done2 <- list.files(banc.nblast.malecns.swc.save.path.version, pattern = "\\.swc$")
message("##### maleCNS native → BANC skeleton transforms complete #####")
message(sprintf("##### Transformed skeletons: %s of %s", length(done2), length(mcns.swcs) + length(done.banc.swc)))

