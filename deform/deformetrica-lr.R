library(bancr)

# Deformetrica folder
deform.folder <- "/Users/abates/projects/flyconnectome/bancpipeline/deformetrica/mirror"
deform.folder.points <- file.path(deform.folder,"Points")
deform.folder.meshes <- file.path(deform.folder,"NeuronMesh")
deform.folder.skeletons <- file.path(deform.folder,"Neuron")
deform.folder.surface <- file.path(deform.folder,"Surface")
deform.folder.imgs <- file.path(deform.folder,"images")
dir.create(deform.folder.imgs, showWarnings = FALSE)
dir.create(deform.folder.points, showWarnings = FALSE)
dir.create(deform.folder.meshes, showWarnings = FALSE)
dir.create(deform.folder.skeletons, showWarnings = FALSE)
dir.create(deform.folder.surface, showWarnings = FALSE)
deform.dirs <- c(deform.folder.surface, deform.folder.points, deform.folder.skeletons, deform.folder.meshes)

################################
### Matched nuclei positions ###
################################

# Add some matched nuclei
bc.meta <- banctable_query()
bc.meta.matches <- bc.meta %>%
  dplyr::filter(!is.na(cell_type), cell_type!="NA", side %in% c("right","left")) %>%
  dplyr::group_by(cell_type) %>%
  dplyr::filter(all(c("right","left")%in%side)) %>%
  dplyr::mutate(nct = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(nct==2) %>%
  dplyr::arrange(cell_type,side)

# Get perfectly matched nuclei coordinates
matched.nuclei <- data.frame()
cts <- unique(bc.meta.matches$cell_type)
for(ct in cts){
  dfr <- subset(bc.meta.matches, cell_type==ct & side == "right")
  dfl <- subset(bc.meta.matches, cell_type==ct & side == "left")
  if(nrow(dfr)&nrow(dfl)){
    df <- rbind(dfr, dfl)[,c("side","nucleus_position_nm")]
    matched.nuclei <- rbind(matched.nuclei,df)
  }
}
matched.nuclei <- matched.nuclei[!is.na(matched.nuclei$X),]
matched.nuclei.right <- nat::xyzmatrix(matched.nuclei$nucleus_position_nm[matched.nuclei$side=="right"])
matched.nuclei.left <- nat::xyzmatrix(matched.nuclei$nucleus_position_nm[matched.nuclei$side=="left"])

# Save as .vtk for deformetrica
matched.nuclei.right.xyz.m <- banc_mirror(matched.nuclei.right)
matched.nuclei.left.xyz.m <- banc_mirror(matched.nuclei.left)

# Remove things that are too far from each other
distances <- euclidean_distances(matched.nuclei.right, matched.nuclei.left.xyz.m)
keep <- distances<=100000
matched.nuclei.right <- matched.nuclei.right[keep,]
matched.nuclei.left.xyz.m <- matched.nuclei.left.xyz.m[keep,]
distances <- euclidean_distances(matched.nuclei.left, matched.nuclei.right.xyz.m)
keep <- distances<=100000
matched.nuclei.left <- matched.nuclei.left[keep,]
matched.nuclei.right.xyz.m <- matched.nuclei.right.xyz.m[keep,]

# Save as .vtk for deformetrica
matched.nuclei.xyz <- rbind(matched.nuclei.right, matched.nuclei.left)
matched.nuclei.xyz.m <- rbind(matched.nuclei.left.xyz.m,matched.nuclei.right.xyz.m)
deformetricar::write.vtk(matched.nuclei.xyz.m/1000, filename = file.path(deform.folder.points,"matched_nuclei_mirrored.vtk"))
deformetricar::write.vtk(matched.nuclei.xyz/1000, filename = file.path(deform.folder.points,"matched_nuclei_original.vtk"))

#################################
### Matched neuropil surfaces ###
#################################

banc_brain_neuropil.surf.m <- banc_mirror(banc_brain_neuropil.surf)
banc_vnc_neuropil.surf.m <- banc_mirror(banc_brain_neuropil.surf)
write_mesh3d_to_vtk(mesh = as.mesh3d(banc_vnc_neuropil.surf)/1000, filename = file.path(deform.folder.surface,"banc_vnc_neuropil_original.vtk"))
write_mesh3d_to_vtk(mesh = as.mesh3d(banc_brain_neuropil.surf)/1000, filename = file.path(deform.folder.surface,"banc_brain_neuropil_original.vtk"))
write_mesh3d_to_vtk(mesh = as.mesh3d(banc_vnc_neuropil.surf.m)/1000, filename = file.path(deform.folder.surface,"banc_vnc_neuropil_mirrored.vtk"))
write_mesh3d_to_vtk(mesh = as.mesh3d(banc_brain_neuropil.surf.m)/1000, filename = file.path(deform.folder.surface,"banc_brain_neuropil_mirrored.vtk"))

#######################
### Matched neurons ###
#######################

# Organise 1:1 hits
bc.meta.matches.unique <- bc.meta %>%
  dplyr::filter(!is.na(banc_match),!banc_match%in%c("NA",""," ","ISSUE!","ISSUE?")) %>%
  dplyr::filter(side=="right") %>%
  dplyr::arrange(dplyr::desc(as.numeric(banc_nblast))) %>%
  dplyr::distinct(root_id, banc_match) %>%
  dplyr::distinct(root_id, .keep_all = TRUE) %>%
  dplyr::distinct(banc_match, .keep_all = TRUE)

# Read data
choose_banc()
cts <- sort(unique(bc.meta.matches$cell_type))
for(ct in cts){

 try({

  # Get IDs
  if(ct=="giant_fiber"){
   next
  }
  bc.meta.matches.ct <- subset(bc.meta.matches, bc.meta.matches$cell_type==ct)
  right.id <- bc.meta.matches.ct[1,]$root_id
  left.id <- bc.meta.matches.ct[1,]$banc_match
  if(is.na(ct)|grepl("ISSUE|NA",ct)|ct%in%c(""," ","NA")){
    next
  }
  ct <- gsub(" ","_",ct)
  ct.left <- paste0(ct,"_left")
  ct.right <- paste0(ct,"_right")
  message("Working on: ", ct)

  # Read meshes
  right.mesh <- banc_read_neuron_meshes(right.id)
  left.mesh <- banc_read_neuron_meshes(left.id)

  # Simplify meshes
  right.mesh.simp <- Rvcg::vcgQEdecim(right.mesh[[1]], percent = 0.1)
  left.mesh.simp <- Rvcg::vcgQEdecim(left.mesh[[1]], percent = 0.1)

  # Cut our the parts within the neuropil
  right.mesh.simp.p <- banc_in_neuropil.mesh3d(right.mesh.simp, surf = banc_neuropil.surf, invert = TRUE)
  left.mesh.simp.p <- banc_in_neuropil.mesh3d(left.mesh.simp, surf = banc_neuropil.surf, invert = TRUE)

  # Mirror
  right.mesh.simp.m <- banc_mirror(right.mesh.simp)
  left.mesh.simp.m <- banc_mirror(left.mesh.simp)
  right.mesh.simp.m.p <- banc_mirror(right.mesh.simp.p)
  left.mesh.simp.m.p <- banc_mirror(left.mesh.simp.p)

  # Save meshes
  write_mesh3d_to_vtk(mesh = right.mesh.simp.p/1000, filename = file.path(deform.folder.meshes, paste0(ct.right,"_original",".vtk")), simplify = TRUE, percent = 1)
  write_mesh3d_to_vtk(mesh = left.mesh.simp.p/1000, filename = file.path(deform.folder.meshes, paste0(ct.left,"_original",".vtk")), simplify = TRUE, percent = 1)
  write_mesh3d_to_vtk(mesh = right.mesh.simp.m.p/1000, filename = file.path(deform.folder.meshes, paste0(ct.right,"_mirrored",".vtk")), simplify = TRUE, percent = 1)
  write_mesh3d_to_vtk(mesh = left.mesh.simp.m.p/1000, filename = file.path(deform.folder.meshes, paste0(ct.left,"_mirrored",".vtk")), simplify = TRUE, percent = 1)

  # Read skeletons
  right.nucleus <- nat::xyzmatrix(subset(bc.meta, root_id==right.id)$nucleus_position_nm)
  left.nucleus <- nat::xyzmatrix(subset(bc.meta, root_id==left.id)$nucleus_position_nm)
  if(is.matrix(left.nucleus)&&is.matrix(right.nucleus)){
    if(ncol(left.nucleus)==3&&ncol(right.nucleus)==3){

      # Read skeletons
      right.skel <- banc_read_l2skel(right.id)
      left.skel <- banc_read_l2skel(left.id)

      # Re-root
      right.skel <- banc_reroot(right.skel[[1]], id = right.id, banc_nuclei = bc.meta)
      left.skel <- banc_reroot(left.skel[[1]], id = left.id, banc_nuclei = bc.meta)

      # Simplify skeletons
      right.skel.simp <- nat::simplify_neuron(right.skel, n=4)
      left.skel.simp <- nat::simplify_neuron(left.skel, n=4)

      # Cut our the parts within the neuropil
      right.skel.simp.p <- hemibrainr:::primary_neurite(right.skel.simp)
      left.skel.simp.p <- hemibrainr:::primary_neurite(left.skel.simp)

      # Mirror
      right.skel.simp.m <- banc_mirror(right.skel.simp)
      left.skel.simp.m <- banc_mirror(left.skel.simp)
      right.skel.simp.m.p <- banc_mirror(right.skel.simp.p)
      left.skel.simp.m.p <- banc_mirror(left.skel.simp.p)

      # Write neurons as .vtk
      nat:::write.vtk.neuron(right.skel.simp.p/1000, file = file.path(deform.folder.skeletons, paste0(ct.right,"_original",".vtk")), format = "vtk")
      nat:::write.vtk.neuron(left.skel.simp.p/1000, file = file.path(deform.folder.skeletons, paste0(ct.left,"_original",".vtk")), format = "vtk")
      nat:::write.vtk.neuron(right.skel.simp.m.p/1000, file = file.path(deform.folder.skeletons, paste0(ct.right,"_mirrored",".vtk")), format = "vtk")
      nat:::write.vtk.neuron(left.skel.simp.m.p/1000, file = file.path(deform.folder.skeletons, paste0(ct.left,"_mirrored",".vtk")), format = "vtk")

    }
  }

  # Save 2D comparison plot\g
  rotation_matrix <- banc_rotation_matrices[["main"]]
  g <- ggplot2::ggplot() +
    #geom_neuron(x = banc_neuropil.surf, rotation_matrix = rotation_matrix, alpha = 0.05, cols = c("grey90", "grey50")) +
    geom_neuron(x=right.mesh.simp.p, rotation_matrix = rotation_matrix, cols = "blue", alpha = 0.5, linewidth = 0.3) +
    geom_neuron(x=left.mesh.simp.p, rotation_matrix = rotation_matrix, cols = "blue", alpha = 0.5, linewidth = 0.3) +
    geom_neuron(x=right.mesh.simp, rotation_matrix = rotation_matrix, cols = "blue", alpha = 0.1, linewidth = 0.3) +
    geom_neuron(x=left.mesh.simp, rotation_matrix = rotation_matrix, cols = "blue", alpha = 0.1, linewidth = 0.3) +
    geom_neuron(x=right.skel.simp.p, rotation_matrix = rotation_matrix, cols = "cyan", alpha = 0.5, linewidth = 0.3) +
    geom_neuron(x=left.skel.simp.p, rotation_matrix = rotation_matrix, cols = "cyan", alpha = 0.5, linewidth = 0.3) +
    geom_neuron(x=right.mesh.simp.m.p, rotation_matrix = rotation_matrix, cols = "red", alpha = 0.5, linewidth = 0.3) +
    geom_neuron(x=left.mesh.simp.m.p, rotation_matrix = rotation_matrix, cols = "red", alpha = 0.5, linewidth = 0.3) +
    geom_neuron(x=right.mesh.simp.m, rotation_matrix = rotation_matrix, cols = "red", alpha = 0.1, linewidth = 0.3) +
    geom_neuron(x=left.mesh.simp.m, rotation_matrix = rotation_matrix, cols = "red", alpha = 0.1, linewidth = 0.3) +
    geom_neuron(x=right.skel.simp.m.p, rotation_matrix = rotation_matrix, cols = "orange", alpha = 0.5, linewidth = 0.3) +
    geom_neuron(x=left.skel.simp.m.p, rotation_matrix = rotation_matrix, cols = "orange", alpha = 0.1, linewidth = 0.3) +
    ggplot2::coord_fixed() +
    ggplot2::theme_void() +
    ggplot2::guides(fill="none",color="none") +
    ggplot2::theme(legend.position = "none",
                   plot.title = ggplot2::element_text(hjust = 0, size = 8, face = "bold", colour = "black"),
                   axis.title.x=ggplot2::element_blank(),
                   axis.text.x=ggplot2::element_blank(),
                   axis.ticks.x=ggplot2::element_blank(),
                   axis.title.y=ggplot2::element_blank(),
                   axis.text.y=ggplot2::element_blank(),
                   axis.ticks.y=ggplot2::element_blank(),
                   axis.line = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(0, 0, 0, 0),
                   panel.spacing = ggplot2::unit(0, "cm"),
                   panel.border = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   plot.background = ggplot2::element_blank())
  ggplot2::ggsave(g, filename=file.path(deform.folder.imgs,paste0(ct,".png")))

 })
}

#######################
### registration.sh ###
#######################

# Create registration.sh
original.filelist.all <- c()
flipped.filelist.all <- c()
params <- c()
for(dd in deform.dirs){

  # Get paired files
  original.filelist <- list.files(dd, pattern = "original", full.names = TRUE)
  flipped.filelist <- gsub("original","mirrored",original.filelist)
  fexist <- file.exists(flipped.filelist)
  original.filelist <- basename(original.filelist[fexist])
  flipped.filelist <- basename(flipped.filelist[fexist])
  original.filelist.all <- c(original.filelist.all,paste0(basename(dd),"/",original.filelist))
  flipped.filelist.all <- c(flipped.filelist.all,paste0(basename(dd),"/",flipped.filelist))
  if(!length(original.filelist)){
    next
  }

  # Create registration file
  outfile <- file.path(deform.folder,tolower(paste0(basename(dd),"_registration.sh")))
  p <- sprintf("param%s.xml",basename(dd))
  params < c(params, rep(p,length(original.filelist)))
  lines <- c("#!/bin/bash\n", "/Applications/deformetrica-2.1/deformetrica/bin/sparseMatching3 paramDiffeos.xml")
  for (target in original.filelist){
    target <- paste0(basename(dd),"/",target)
    if(!basename(dd)%in%c("Surface","Points")){
      if(grepl("right",target)){
        match <- gsub("right","left",target)
      }else if(grepl("left",target)){
        match <- gsub("left","right",target)
      }else{
        next
      }
    }
    match <- gsub("original","mirrored",match)
    line <- c(p, match, target)
    lines <- c(lines, paste(line, collapse = " "))
  }
  cat(lines, file = outfile)
}

# Write master regisration.sh
overall.outfile <- file.path(deform.folder,"registration.sh")
p <- "paramNeuronMesh.xml"
lines <- c("#!/bin/bash\n", "/Applications/deformetrica-2.1/deformetrica/bin/sparseMatching3 paramDiffeos.xml")
for (i in 1:length(original.filelist.all)){
  target <- original.filelist.all[i]
  p <- params[i]
  if(grepl("right",target)){
    match <- gsub("right","left",target)
  }else if(grepl("left",target)){
    match <- gsub("left","right",target)
  }else{
    next
  }
  match <- gsub("original","mirrored",match)
  line <- c(p, match, target)
  lines <- c(lines, paste(line, collapse = " "))
}
cat(lines, file = overall.outfile)

#########################
### USE DEFORMETRICAR ###
#########################

# orig.points <- data.frame()
# mirrored.filelist <- list.files(deform.dir, pattern = "mirrored", full.names = TRUE)
# for(mf in mirrored.filelist){
#   dat <- read.vrk(mf)
#   xyz <- nat::xyzmatrix(dat)
#   xyz$PointNo <- 1:nrow(xyz)
#   xyz$file <- basename(mf)
#   orig.points <- rbind(orig.points,xyz)
# }
# orig.points <- orig.points %>%
#   dplyr::arrange(file, PointNo)
#
# deform.points <- data.frame()
# mirrored.filelist <- list.files(deform.dir, pattern = "__t_9.vtk", full.names = TRUE)
# for(mf in mirrored.filelist){
#   dat <- read.vrk(mf)
#   xyz <- nat::xyzmatrix(dat)
#   xyz$PointNo <- 1:nrow(xyz)
#   xyz$file <- gsub("__t_9","",basename(mf))
#   deform.points <- rbind(deform.points,xyz)
# }
# deform.points <- deform.points %>%
#   dplyr::arrange(file, PointNo)
#
# chosen.fies <- unique(deform.points$file)
# orig.points <- subset(orig.points, orig.points$file %in% chosen.fies)
#
#
#
