#' banc-malecns-mesh-transform — Transform maleCNS meshes into BANC space (navis TPS).
#'
#' For each maleCNS bodyid, fetches the OBJ mesh and applies the
#' JRCFIB2022M → BANC TPS transform via Python `navis` + `flybrains`
#' (no R-native registration exists yet for this chain).
#'
#' @section Reads:
#'   - `<banc.meta.save.path>/malecns_meta.csv`
#'   - maleCNS meshes via `malecns` package
#'
#' @section Writes:
#'   - `<banc.nblast.malecns.obj.save.path>/<version>/<bodyid>.obj`
#'
#' @section Invoked by:
#'   `o2/production/o2_banc_ngl.sh`
#'
#' @section Notes:
#'   - Requires `reticulate` Python with `navis`, `flybrains`, `numpy`.

##########################################
### TRANSFORM maleCNS NEURONS TO BANC  ###
##########################################
source("banc/banc-startup.R")
library(malecns)
library(reticulate)
redo <- FALSE
version <- banc.nblast.malecns.version

# Make a version folder for the obj data
banc.nblast.malecns.obj.save.path.version <- file.path(banc.nblast.malecns.obj.save.path, version)
dir.create(banc.nblast.malecns.obj.save.path.version, recursive = TRUE, showWarnings = FALSE)

# Setup Python transform backend (navis/flybrains)
# No R-native registration exists for JRCFIB2022M → BANC,
# so we use the Python TPS via reticulate.
navis     <- import("navis", convert = FALSE)
flybrains <- import("flybrains")
np        <- import("numpy", convert = FALSE)
message("Python transform backend ready (navis + flybrains)")

# Get meta data
mcns.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path, "malecns_09_meta.csv"),
                                              col_types = hemibrainr:::sql_col_types))
mcns.ids <- unique(mcns.meta$malecns_09_id)

# Count existing obj files
mcns.obj.full.files <- list.files(banc.nblast.malecns.obj.save.path.version, pattern = ".obj$", full.names = TRUE)
mcns.obj.ids <- gsub(pattern = "\\.obj$", "", basename(mcns.obj.full.files))

# What obj files have we not yet transformed?
if (redo) {
  mcns.obj.todo.files <- mcns.ids
} else {
  mcns.obj.todo.files <- setdiff(mcns.ids, mcns.obj.ids)
}
message("There are ", length(mcns.obj.todo.files), " todo, we have downloaded ", length(mcns.obj.ids), " maleCNS .obj files")

# Transforming meshes
# Note: reticulate/Python calls are not fork-safe, so we use %do% (sequential).
# The Python TPS transform is fast (~0.3-2s per mesh), so the bottleneck is
# mesh reading from DVID, not the transform itself.
message("##### Transforming meshes #####")
mcns.obj.todo.files <- sample(mcns.obj.todo.files)
by.query <- foreach(id = sample(mcns.obj.todo.files),
                    .combine = 'c',
                    .errorhandling = 'pass') %do% {

                      # File
                      new.file <- file.path(banc.nblast.malecns.obj.save.path.version, paste0(id, ".obj"))
                      if (file.exists(new.file) & !redo) {
                        message("Skipping written file: ", new.file)
                        NULL
                      } else {

                        # Get mesh3d in nm (read_mcns_meshes handles the x8 voxel conversion)
                        mcns.mesh <- read_mcns_meshes(id, units = "nm")

                        # Transform vertex coordinates: JRCFIB2022M (nm) → BANC (nm)
                        verts.nm <- t(mcns.mesh[[1]]$vb[1:3, ])
                        verts.py <- np$array(r_to_py(verts.nm))
                        verts.banc.py <- navis$xform_brain(verts.py, source = "JRCFIB2022M", target = "BANC")
                        verts.banc <- py_to_r(verts.banc.py)

                        # Replace vertices in mesh
                        mesh3d.banc <- mcns.mesh[[1]]
                        mesh3d.banc$vb[1, ] <- verts.banc[, 1]
                        mesh3d.banc$vb[2, ] <- verts.banc[, 2]
                        mesh3d.banc$vb[3, ] <- verts.banc[, 3]

                        # Simplify meshes
                        # mesh3d.banc <- Rvcg::vcgQEdecim(mesh3d.banc, percent = 0.2)

                        # Guard: skip vertex-only meshes (no faces). write_mesh3d_to_obj
                        # would otherwise silently produce a faceless .obj that trimesh
                        # later loads as a PointCloud and the uploader rejects.
                        if (is.null(mesh3d.banc$it) || ncol(mesh3d.banc$it) == 0) {
                          message("Skipping ", id, ": mesh has no faces (vertex-only mesh from DVID)")
                          NULL
                        } else {
                          # Save as a .obj file
                          bancr:::write_mesh3d_to_obj(mesh3d.banc, filename = new.file)
                          message("Written: ", new.file)

                          # Return nothing
                          NULL
                        }
                      }
                    }

# Were there errors?
message("##### Displaying any errors from foreach loop #####")
for (i in seq_along(by.query)) {
  if (!is.null(by.query[[i]])) {
    message(by.query[[i]])
  }
}

# Announce
done <- list.files(banc.nblast.malecns.obj.save.path.version, pattern = ".obj$")
message("##### BANCpipeline: malecns meshes in BANC space, updated #####")
message(sprintf("##### malecns meshes transformed .obj files: %s", length(done)))
message(sprintf("##### total neurons in metadata: %s", length(mcns.ids)))
