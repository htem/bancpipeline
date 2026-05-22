library(bancr)
library(malevnc)
library(nat.jrcbrains)
nat.jrcbrains::register_saalfeldlab_registrations()

# Deformetrica folder
deform.folder <- "/Users/abates/projects/flyconnectome/bancpipeline/deformetrica/manc"
deform.folder.points <- file.path(deform.folder,"Points")
deform.folder.meshes <- file.path(deform.folder,"NeuronMesh")
deform.folder.skeletons <- file.path(deform.folder,"Neuron")
deform.folder.surface <- file.path(deform.folder,"Surface")
deform.folder.imgs <- file.path(deform.folder,"images")
dir.create(deform.folder.points, showWarnings = FALSE, recursive = TRUE)
dir.create(deform.folder.meshes, showWarnings = FALSE, recursive = TRUE)
dir.create(deform.folder.skeletons, showWarnings = FALSE, recursive = TRUE)
dir.create(deform.folder.imgs, showWarnings = FALSE, recursive = TRUE)
dir.create(deform.folder.surface, showWarnings = FALSE, recursive = TRUE)
deform.dirs <- c(deform.folder.surface, deform.folder.points,
                 deform.folder.skeletons, deform.folder.meshes)

################################
### Matched nuclei positions ###
################################

# Add some matched nuclei
bc.meta <- banctable_query()
bc.meta.matches <- bc.meta %>%
  dplyr::filter(!is.na(cell_type),
                cell_type!="NA",
                !is.na(manc_match),
                ! manc_match %in% c("NA",""," "),
                side %in% c("right","left")) %>%
  dplyr::group_by(cell_type) %>%
  dplyr::filter(all(c("right","left")%in%side), proofread==TRUE) %>%
  dplyr::mutate(nct = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(nct==2) %>%
  dplyr::arrange(cell_type,side) %>%
  dplyr::distinct(root_id, .keep_all = TRUE) %>%
  dplyr::distinct(root_id, manc_match, side, cell_type)

################################
### Matched nuclei positions ###
################################

# NOT IMPLEMENTED

#################################
### Matched neuropil surfaces ###
#################################

manc.jrcvnc2018f=xform_brain(MANC.surf, reference = "JRCVNC2018F", sample="MANC")
manc.b <- banc_to_JRC2018F(manc.jrcvnc2018f, region="vnc", method="tpsreg", banc.units = "nm", inverse = TRUE)
write_mesh3d_to_vtk(mesh = as.mesh3d(banc_vnc_neuropil.surf)/1000, filename = file.path(deform.folder.surface,"vnc_neuropil_banc.vtk"))
write_mesh3d_to_vtk(mesh = as.mesh3d(manc.b)/1000, filename = file.path(deform.folder.surface,"vnc_neuropil_manc.vtk"))

#######################
### Matched neurons ###
#######################

# Read data
for(i in 1:nrow(bc.meta.matches)){

  try({

  # Get IDs
  banc.id <- bc.meta.matches[i,]$root_id
  manc.id <- bc.meta.matches[i,]$manc_match
  ct <- bc.meta.matches[i,]$cell_type
  side <- bc.meta.matches[i,]$side
  if(is.na(ct)|grepl("ISSUE|NA",ct)|ct%in%c(""," ","NA")){
    next
  }
  ct <- gsub(" ","_",ct)
  ct.manc <- paste0(ct,"_",side,"_manc")
  ct.banc <- paste0(ct,"_",side,"_banc")
  message("Working on: ", ct,"_",side)

  ### BANC ###

  # BANC meshes
  banc.mesh <- banc_read_neuron_meshes(banc.id)
  banc.mesh <- banc_decapitate(banc.mesh, invert = FALSE)
  banc.mesh.simp <- Rvcg::vcgQEdecim(banc.mesh[[1]], percent = 0.1)
  banc.mesh.simp.p <- banc_in_neuropil.mesh3d(banc.mesh.simp, surf = banc_neuropil.surf, invert = TRUE)

  # BANC skeletons
  banc.skel <- banc_read_l2skel(banc.id)
  soma.exists <- nrow(subset(bc.meta, root_id==banc.id & nucleus_id !='0'))
  banc.skel <- banc_reroot(banc.skel, id = banc.id, banc_nuclei = bc.meta)
  banc.skel <- banc_decapitate(banc.skel, invert = FALSE)
  banc.skel <- nat::resample(banc.skel, 1000)
  banc.skel.p <- hemibrainr:::primary_neurite(banc.skel)

  # Write neurons as .vtk
  if(soma.exists){
    nat:::write.vtk.neuron(banc.skel.p[[1]]/1000, file = file.path(deform.folder.skeletons, paste0(ct.banc,".vtk")), format = "vtk")
  }
  bancr:::write_mesh3d_to_vtk(mesh = banc.mesh.simp.p/1000, filename = file.path(deform.folder.meshes, paste0(ct.banc,".vtk")), simplify = TRUE, percent = 1)

  ### MANC ###

  #  MANC meshes
  manc.mesh <- read_manc_meshes(manc.id)
  manc.neuron <- manc_read_neurons(manc.id, unit = "nm")

  # Transform into JRCVNC2918F
  manc.mesh.jrcvnc2018f=xform_brain(manc.mesh/1e3, reference = "JRCVNC2018F", sample="MANC")
  manc.neuron.jrcvnc2018f=xform_brain(manc.neuron/1e3, reference = "JRCVNC2018F", sample="MANC")

  # Transform into the BANC
  manc.mesh.b <- banc_to_JRC2018F(manc.mesh.jrcvnc2018f, region="vnc", method="tpsreg", banc.units = "nm", inverse = TRUE)
  manc.neuron.b <- banc_to_JRC2018F(manc.neuron.jrcvnc2018f, region="vnc", method="tpsreg", banc.units = "nm", inverse = TRUE)
  manc.mesh.b.simp <- Rvcg::vcgQEdecim(manc.mesh.b[[1]], percent = 0.1)
  manc.neuron.b <- nat::resample(manc.neuron.b, 1000)

  # Write neurons as .vtk
  manc.neuron.b.p <- hemibrainr:::primary_neurite(manc.neuron.b)
  manc.mesh.b.simp.p <- banc_in_neuropil.mesh3d(manc.mesh.b.simp, surf = banc_neuropil.surf, invert = TRUE)
  if(soma.exists){
    nat:::write.vtk.neuron(manc.neuron.b.p[[1]]/1000, file = file.path(deform.folder.skeletons, paste0(ct.manc,".vtk")), format = "vtk")
  }
  write_mesh3d_to_vtk(mesh = manc.mesh.b.simp.p/1000, filename = file.path(deform.folder.meshes, paste0(ct.manc,".vtk")), simplify = TRUE, percent = 1)

  # Save 2D comparison plot\g
  g <- ggplot2::ggplot() +
    geom_neuron(x=banc.mesh.simp.p, rotation_matrix = rotation_matrix, cols = "blue", alpha = 0.5, linewidth = 0.3) +
    geom_neuron(x=banc.mesh.simp, rotation_matrix = rotation_matrix, cols = "blue", alpha = 0.1, linewidth = 0.3) +
    geom_neuron(x=banc.skel.p, rotation_matrix = rotation_matrix, cols = "cyan", alpha = 0.5, linewidth = 0.3) +
    geom_neuron(x=manc.mesh.b.simp.p, rotation_matrix = rotation_matrix, cols = "red", alpha = 0.5, linewidth = 0.3) +
    geom_neuron(x=manc.mesh.b.simp, rotation_matrix = rotation_matrix, cols = "red", alpha = 0.1, linewidth = 0.3) +
    geom_neuron(x=manc.neuron.b.p, rotation_matrix = rotation_matrix, cols = "orange", alpha = 0.1, linewidth = 0.3) +
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
  ggplot2::ggsave(g, filename=file.path(deform.folder.imgs,paste0(ct,"_",side,".png")))
  ggsave(g, filename=file.path(deform.folder.imgs,paste0(ct,"_",side,".png")))

  })
}

#######################
### registration.sh ###
#######################

# Create registration.sh
banc.filelist.all <- c()
flipped.filelist.all <- c()
params <- c()
for(dd in deform.dirs){

  # Get paired files
  banc.filelist <- list.files(dd, pattern = ".vtk", full.names = TRUE)
  flipped.filelist <- gsub("banc","manc",banc.filelist)
  fexist <- file.exists(flipped.filelist)
  banc.filelist <- basename(banc.filelist[fexist])
  flipped.filelist <- basename(flipped.filelist[fexist])
  banc.filelist.all <- c(banc.filelist.all,paste0(basename(dd),"/",banc.filelist))
  flipped.filelist.all <- c(flipped.filelist.all,paste0(basename(dd),"/",flipped.filelist))
  if(!length(banc.filelist)){
    next
  }

  # Create registration file
  outfile <- file.path(deform.folder,tolower(paste0(basename(dd),"_registration.sh")))
  p <- sprintf("param%s.xml",basename(dd))
  params < c(params, rep(p,length(banc.filelist)))
  lines <- c("#!/bin/bash\n", "/Applications/deformetrica-2.1/deformetrica/bin/sparseMatching3 paramDiffeos.xml")
  for (target in banc.filelist){
    target <- paste0(basename(dd),"/",target)
    if(grepl("right",target)){
      match <- gsub("right","left",target)
    }else if(grepl("left",target)){
      match <- gsub("left","right",target)
    }else{
      next
    }
    match <- gsub("banc","manc",match)
    line <- c(p, match, target)
    lines <- c(lines, paste(line, collapse = " "))
  }
  cat(lines, file = outfile)
}

# Write master regisration.sh
overall.outfile <- file.path(deform.folder,"registration.sh")
p <- "paramNeuronMesh.xml"
lines <- c("#!/bin/bash\n", "/Applications/deformetrica-2.1/deformetrica/bin/sparseMatching3 paramDiffeos.xml")
for (i in 1:length(banc.filelist.all)){
  target <- banc.filelist.all[i]
  p <- params[i]
  if(grepl("right",target)){
    match <- gsub("right","left",target)
  }else if(grepl("left",target)){
    match <- gsub("left","right",target)
  }else{
    next
  }
  match <- gsub("banc","manc",match)
  line <- c(p, match, target)
  lines <- c(lines, paste(line, collapse = " "))
}
cat(lines, file = overall.outfile)
