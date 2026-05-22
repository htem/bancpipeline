########################
### BANC-FAFB NBLAST ###
########################
source("banc/banc-startup.R")
redo <- FALSE
version <- banc.nblast.version

# Direct us to the BANC dataset
banc.nblast.fafb.swc.save.path.version <- file.path(banc.nblast.fafb.swc.save.path,version)
banc.nblast.fafb.obj.save.path.version <- file.path(banc.nblast.fafb.obj.save.path,version)

# Deformetrica folder
deform.folder <- file.path(banc.deform.save.path,'fafb')
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

#######################
### Matched neurons ###
#######################

# Add some matched nuclei
bc.meta <- banctable_query()
bc.meta.matches <- bc.meta %>%
  dplyr::filter(!is.na(cell_type),
                cell_type!="NA",
                !is.na(fafb_match),
                ! fafb_match %in% c("NA",""," "),
                side %in% c("right","left")) %>%
  dplyr::group_by(cell_type) %>%
  dplyr::filter(all(c("right","left")%in%side), proofread==TRUE) %>%
  dplyr::mutate(nct = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(nct==2) %>%
  dplyr::arrange(cell_type,side) %>%
  dplyr::distinct(root_id, .keep_all = TRUE) %>%
  dplyr::distinct(root_id, fafb_match, side, cell_type)

################################
### Matched nuclei positions ###
################################

# get potential soma matches
bc.meta.matches.nulcei <- bc.meta %>%
  dplyr::filter(!is.na(cell_type),
                cell_type!="NA",
                !is.na(fafb_match),
                ! fafb_match %in% c("NA",""," "),
                !is.na(nucleus_position_nm),
                !nucleus_position_nm%in%c(" ","","NA"),
                side %in% c("right","left")) %>%
  dplyr::group_by(cell_type) %>%
  dplyr::filter(all(c("right","left")%in%side), proofread==TRUE) %>%
  dplyr::mutate(nct = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(nct==2) %>%
  dplyr::arrange(cell_type,side) %>%
  dplyr::distinct(root_id, .keep_all = TRUE) %>%
  dplyr::distinct(root_id, fafb_match, side, cell_type, nucleus_position_nm)

# Get fafb meta data
fw.meta <- suppressWarnings(readr::read_csv(file = file.path(banc.meta.save.path,"flywire_meta.csv"), 
                                            col_types = hemibrainr:::sql_col_types))
banc.fw.meta <- fw.meta %>%
  dplyr::filter(!is.na(soma_x), !is.na(soma_y), !is.na(soma_z), !is.na(nucleus_id), nucleus_id!="0") %>%
  dplyr::select(root_783, soma_x, soma_y, soma_z) %>%
  dplyr::inner_join(bc.meta.matches.nulcei, by = c("root_783"="fafb_match"))

# put soma in BANC space
fw.soma <- nat::xyzmatrix(banc.fw.meta[,c("soma_x","soma_y","soma_z")])
fw.positions.xyz <- flywire_raw2nm(nat::xyzmatrix(fw.soma))
fw.positions.xyz.jrc2018f <- nat.templatebrains::xform_brain(fw.positions.xyz,
                                                             sample = "FAFB14",
                                                             reference = "JRC2018F")
fw.positions.xyz.jrc2018f.m <- nat.templatebrains::mirror_brain(fw.positions.xyz.jrc2018f,
                                                                brain = nat.flybrains::JRC2018F,
                                                                transform = "flip")
fw.positions.xyz.banc <- bancr::banc_to_JRC2018F(fw.positions.xyz.jrc2018f.m, region = "brain", method = "tpsreg", inverse = TRUE)
nat::xyzmatrix(fw.soma) <- fw.positions.xyz.banc
banc.fw.meta[,c("soma_x","soma_y","soma_z")] <- fw.positions.xyz.banc
banc.soma <- nat::xyzmatrix(banc.fw.meta$nucleus_position_nm)

# Remove things that are too far from each other
distances <- euclidean_distances(banc.soma, fw.soma)
keep <- distances<=100000
banc.soma <- banc.soma[keep,]
fw.soma <- fw.soma[keep,]
banc.fw.meta <- banc.fw.meta[keep,]
bc.meta.matches <- subset(bc.meta.matches, bc.meta.matches$root_id %in% banc.fw.meta$root_id) %>%
  dplyr::filter(side %in% c("right","left"))

# Save as .vtk for deformetrica
deformetricar::write.vtk(banc.soma/1000, filename = file.path(deform.folder.points,"matched_nuclei_banc.vtk"))
deformetricar::write.vtk(fw.soma/1000, filename = file.path(deform.folder.points,"matched_nuclei_fafb.vtk"))

#################################
### Matched neuropil surfaces ###
#################################

fafb.jrc2018f <- xform_brain(elmr::FAFB14.surf, reference = "JRC2018F", sample="FAFB14")
fafb.jrc2018f <- nat.templatebrains::mirror_brain(fafb.jrc2018f,
                                                  brain = nat.flybrains::JRC2018F,
                                                  transform = "flip")
fafb.b <- banc_to_JRC2018F(fafb.jrc2018f, region="brain", method="tpsreg", banc.units = "nm", inverse = TRUE)
bancr:::write_mesh3d_to_vtk(mesh = as.mesh3d(banc_brain_neuropil.surf)/1000, filename = file.path(deform.folder.surface,"brain_neuropil_banc.vtk"),
                            simplify = TRUE, percent = 0.25)
bancr:::write_mesh3d_to_vtk(mesh = as.mesh3d(fafb.b)/1000, filename = file.path(deform.folder.surface,"brain_neuropil_fafb.vtk"),
                            simplify = TRUE, percent = 0.25)

#######################
### Matched neurons ###
#######################

# Read data
flywire.nuclei <- banc.fw.meta %>%
  dplyr::select(-root_id) %>%
  dplyr::rename(root_id=root_783) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(nucleus_position_nm=paste(c(soma_x,soma_y,soma_z),collapse=","),
                nucleus_id='1') %>%
  dplyr::ungroup() %>%
  dplyr::distinct(root_id,nucleus_position_nm,nucleus_id)
write.mesh <- TRUE
for(i in 1:nrow(bc.meta.matches)){
  
  try({
    
    # Get IDs
    banc.id <- bc.meta.matches[i,]$root_id
    fafb.id <- bc.meta.matches[i,]$fafb_match
    ct <- bc.meta.matches[i,]$cell_type
    side <- bc.meta.matches[i,]$side
    if(is.na(ct)|grepl("ISSUE|NA",ct)|ct%in%c(""," ","NA")){
      next
    }
    ct <- gsub(" ","_",ct)
    ct.fafb <- paste0(ct,"_",side,"_fafb")
    ct.banc <- paste0(ct,"_",side,"_banc")
    ct.fafb.skel <- paste0("skeleton_",ct,"_", side,"_fafb")
    ct.banc.skel <- paste0("skeleton_",ct,"_",side,"_banc")
    message("Working on: ", ct,"_",side)
    make.plot <- 0
    
    ### BANC ###
    
    # BANC meshes
    banc.mesh.file <- file.path(deform.folder.meshes, paste0(ct.banc,".vtk"))
    if(!file.exists(banc.mesh.file)){
      if(write.mesh){
        #banc.mesh <- banc_read_neuron_meshes(banc.id)
        banc.mesh <- readobj::read.obj(file.path(banc.obj.save.path,paste0(banc.id,".obj")), convert.rgl=TRUE)
        banc.mesh <- banc_decapitate(banc.mesh[[1]], invert = TRUE, y.cut = 250000)
        banc.mesh.simp <- Rvcg::vcgQEdecim(banc.mesh, percent = 0.1)
        banc.mesh.simp.p <- bancr:::banc_in_neuropil.mesh3d(banc.mesh.simp, surf = banc_neuropil.surf, invert = TRUE)
        nat::xyzmatrix(banc.mesh.simp.p) <- round(nat::xyzmatrix(banc.mesh.simp.p))
        bancr:::write_mesh3d_to_vtk(mesh = banc.mesh.simp.p/1000, filename = banc.mesh.file, simplify = TRUE, percent = 1)
        make.plot <- make.plot + 1
      }
    }
    
    # Do we have the soma?
    banc.soma <- subset(bc.meta, root_id==banc.id & nucleus_id !='0')
    soma.exists <- nrow(banc.soma)
    fw.soma <- subset(flywire.nuclei, flywire.nuclei$root_id==fafb.id)
    if(soma.exists){
      soma.exists <- nrow(fw.soma)
    }
    
    # BANC skeletons
    banc.skel.file <- file.path(deform.folder.skeletons, paste0(ct.banc.skel,".vtk"))
    if(!file.exists(banc.skel.file)){
      # banc.skel <- banc_read_l2skel(banc.id)
      # banc.skel <- banc_reroot(banc.skel, id = banc.id, banc_nuclei = bc.meta)
      banc.skel <- nat::read.neuron(file.path(banc.l2swc.save.path,paste0(banc.id,".swc")))
      banc.skel <- banc_reroot(banc.skel,id = banc.id, banc_nuclei = banc.soma)
      banc.skel <- banc_decapitate(banc.skel, invert = TRUE, y.cut = 250000)
      banc.skel <- nat::resample(banc.skel, 1000)
      # banc.skel <- prune_strahler(banc.skel, orderstoprune = 1)
      banc.skel.p <- simplify_neuron(banc.skel, n = 1) #hemibrainr:::primary_neurite(banc.skel)
      if(soma.exists){
        nat::xyzmatrix(banc.skel.p) <- round(nat::xyzmatrix(banc.skel.p))
        write_neuron_to_vtk_paired(banc.skel.p/1000, 
                                   file = banc.skel.file)
      }
      make.plot <- make.plot + 1
    }
    
    
    ### FAFB ###
    
    #  FAFB meshes
    fafb.mesh.file <- file.path(deform.folder.meshes, paste0(ct.fafb,".vtk"))
    if(!file.exists(fafb.mesh.file)){
      if(write.mesh){
        # choose_segmentation("flywire31")
        # fafb.mesh <- read_cloudvolume_meshes(fafb.id,
        #                                      cloudvolume.url = "graphene://https://prod.flywire-daf.com/segmentation/table/fly_v31")
        # 
        # Transform into JRCVNC2918F
        # fafb.mesh.jrc2018f <- xform_brain(fafb.mesh, reference = "JRC2018F", sample="FAFB14")
        # fafb.mesh.jrc2018f.m <- nat.templatebrains::mirror_brain(fafb.mesh.jrc2018f,
        #                                                           brain = nat.flybrains::JRC2018F,
        #                                                           transform = "flip")
        # fafb.mesh.b <- banc_to_JRC2018F(fafb.mesh.jrc2018f.m, region="brain", method="tpsreg", banc.units = "nm", inverse = TRUE)
        fafb.mesh.b <- readobj::read.obj(file.path(banc.nblast.fafb.obj.save.path.version,paste0(fafb.id,".obj")), convert.rgl=TRUE)
        fafb.mesh.b.simp <- banc_decapitate(fafb.mesh.b[[1]], invert = TRUE, y.cut = 250000)
        fafb.mesh.b.simp <- Rvcg::vcgQEdecim(fafb.mesh.b.simp, percent = 0.1)
        fafb.mesh.b.simp.p <- bancr:::banc_in_neuropil.mesh3d(fafb.mesh.b.simp, surf = banc_neuropil.surf, invert = TRUE)
        nat::xyzmatrix(fafb.mesh.b.simp.p) <- round(nat::xyzmatrix(fafb.mesh.b.simp.p))
        bancr:::write_mesh3d_to_vtk(mesh = fafb.mesh.b.simp.p/1000, filename = fafb.mesh.file, simplify = TRUE, percent = 1)
        make.plot <- make.plot + 1
      }
    }
    
    # # Transform into the BANC
    fafb.skel.file <- file.path(deform.folder.skeletons, paste0(ct.fafb.skel,".vtk"))
    if(!file.exists(fafb.skel.file)){
      # fafb.neuron <- NULL
      # fafb.neuron <- try(read_l2skel(fafb.id,
      #                                cloudvolume.url = "graphene://https://prod.flywire-daf.com/segmentation/table/fly_v31"))
      # if(!length(fafb.neuron)){
      #   next
      # }
      # fafb.neuron.jrc2018f <- xform_brain(fafb.neuron, reference = "JRC2018F", sample="FAFB14")
      # fafb.neuron.jrc2018f.m <- nat.templatebrains::mirror_brain(fafb.neuron.jrc2018f,
      #                                                          brain = nat.flybrains::JRC2018F,
      #                                                          transform = "flip")
      # fafb.neuron.b <- banc_to_JRC2018F(fafb.neuron.jrc2018f.m, region="brain", method="tpsreg", banc.units = "nm", inverse = TRUE)
      # fafb.neuron.b <- nat::resample(fafb.neuron.b, 1000)
      # fafb.neuron.b <- banc_reroot(fafb.neuron.b, id = fafb.id, banc_nuclei = flywire.nuclei)
      fafb.neuron.b <- nat::read.neuron(file.path(banc.nblast.fafb.swc.save.path.version,paste0(fafb.id,".swc")))
      fafb.neuron.b <- banc_decapitate(fafb.neuron.b, invert = TRUE, y.cut = 250000)
      fafb.neuron.b <- nat::resample(fafb.neuron.b, 1000)
      # fafb.neuron.b <- prune_strahler(fafb.neuron.b, orderstoprune = 1)
      fafb.neuron.b.p <- simplify_neuron(fafb.neuron.b, n = 1) #hemibrainr:::primary_neurite(fafb.neuron.b)
      if(soma.exists){
        nat::xyzmatrix(fafb.neuron.b.p) <- round(nat::xyzmatrix(fafb.neuron.b.p))
        write_neuron_to_vtk_paired(fafb.neuron.b.p/1000, 
                                   file = fafb.skel.file)
      } 
      make.plot <- make.plot + 1
    }
    
    # Save 2D comparison plot
    if(make.plot==4){
      rotation_matrix <- bancr:::banc_rotation_matrices[["main"]]
      g <- ggplot2::ggplot() +
        geom_neuron(x=banc.mesh.simp.p, rotation_matrix = rotation_matrix, cols = "blue", alpha = 0.5, linewidth = 0.3) +
        #geom_neuron(x=banc.mesh.simp, rotation_matrix = rotation_matrix, cols = "blue", alpha = 0.1, linewidth = 0.3) +
        geom_neuron(x=banc.skel.p, rotation_matrix = rotation_matrix, cols = "cyan", alpha = 1, linewidth = 0.3) +
        geom_neuron(x=fafb.mesh.b.simp.p, rotation_matrix = rotation_matrix, cols = "red", alpha = 0.5, linewidth = 0.3) +
        geom_neuron(x=fafb.mesh.b.simp, rotation_matrix = rotation_matrix, cols = "red", alpha = 0.1, linewidth = 0.3) +
        geom_neuron(x=fafb.neuron.b.p, rotation_matrix = rotation_matrix, cols = "darkorange", alpha = 1, linewidth = 0.3) +
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
    }
  })
}

########################
### registration.xml ###
########################

# Create registration.sh
sample_files <- c()
subject_ids <- list()
deformable_object_types <- c()
noise_stds <- c()
kernel_widths <- c()
for(dd in deform.dirs){
  
  # object type
  if(basename(dd)=="Points"){
    deformable_object_type <- "LandMark"
    kernel_width <- 25
    noise_std <- 0.001
  }else if(basename(dd)=="NeuronMesh"){
    deformable_object_type <- "SurfaceMesh"
    kernel_width <- 25
    noise_std <- 0.1
  }else if(basename(dd)=="Neuron"){
    deformable_object_type <- "PolyLine"
    kernel_width <- 25
    noise_std <- 0.01
  }else if(basename(dd)=="Surface"){
    deformable_object_type <- "SurfaceMesh"
    kernel_width <- 25
    noise_std <- 1
  }

  # Get paired files
  banc.filelist <- list.files(dd, pattern = "banc.vtk", full.names = TRUE)
  flipped.filelist <- file.path(dd,gsub("banc","fafb",basename(banc.filelist)))
  keep <- file.exists(flipped.filelist)
  banc.filelist <- banc.filelist[keep]
  flipped.filelist <- flipped.filelist[keep]
  sample_files <- c(sample_files,flipped.filelist)
  subject_ids[[basename(dd)]] <- banc.filelist
  
  # Types
  deformable_object_types <- c(deformable_object_types,rep(deformable_object_type,length(banc.filelist)))
  noise_stds <- c(noise_stds,rep(noise_std,length(banc.filelist)))
  kernel_widths <- c(kernel_widths,rep(kernel_width,length(banc.filelist)))
  
  # The dataset .xml
  create_dataset_xml(data_folder="",
                     subject_files = banc.filelist,
                     output_file = file.path(dd,"data_set.xml"),
                     dataset="_banc")
  
  # The dataset .xml
  create_model_xml(flipped.filelist,
                   output_file = file.path(dd,"model.xml"),
                   deformable_object_type = deformable_object_type,
                   dataset="_fafb",
                   kernel_type = "torch",
                   attachment_type = "Varifold",
                   noise_std = noise_std,
                   kernel_width = kernel_width,
                   kernel_device = "cpu",
                   number_of_timepoints = 10)
  
  # The optimisation .xml
  create_optimization_xml(output_file = file.path(dd, "optimization_parameters.xml"),
                          optimization_method_type = "GradientAscent",
                          max_iterations = 100,
                          save_every_n_iters = 10,
                          print_every_n_iters = 1,
                          convergence_tolerance = 1e-4,
                          initial_step_size = 1e-3,
                          use_cuda = "On",
                          number_of_processes = 1,
                          freeze_template = "On")
  
  # source activate deformetrica
  # deformetrica estimate model.xml data_set.xml --p optimization_parameters.xml
  
  # file.edit('test.R')
  
}

# Advanced usage with custom subject and visit IDs
create_dataset_xml(data_folder="", 
                   output_file = file.path(deform.folder,"data_set.xml"),
                   subject_files = unlist(subject_ids),
                   #subject_ids = NULL, 
                   #visit_ids = subject_ids,
                   dataset="_banc")

# Advanced usage with custom parameters
create_model_xml(sample_files, 
                 output_file = file.path(deform.folder,"model.xml"),
                 deformable_object_type = deformable_object_types,
                 kernel_type = "torch",
                 attachment_type = "varifold",
                 noise_std = noise_stds,
                 kernel_width = kernel_widths,
                 kernel_device = "cpu",
                 number_of_timepoints = 10,
                 dataset="_fafb")

# create optimization .xml
create_optimization_xml(output_file = file.path(deform.folder, "optimization_parameters.xml"),
                        optimization_method_type = "GradientAscent",
                        max_iterations = 100,
                        save_every_n_iters = 10,
                        print_every_n_iters = 1,
                        convergence_tolerance = 1e-4,
                        initial_step_size = 1e-3,
                        use_cuda = "Off",
                        number_of_processes = 1,
                        freeze_template = "On")

