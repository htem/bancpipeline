### -----------------------------------------------------------------
### Legacy exploratory analysis by a member of the Wilson Lab.
### Historical record only — not part of the bancpipeline release.
### -----------------------------------------------------------------

######################
### load libraries ###
######################

# Load packages
library(bancr)
library(tidyverse)
library(nat.nblast)
library(dendextend)

# Make sure all functions query BANC and not FAFB
# choose_banc()

################
### get data ###
################

# Load neuronlist back if already gotten
# l2 <- readRDS("l2_an.rds")
# l2dps <- readRDS("l2dps_an.rds")
# nb.all.full<- readRDS("nb.all.full_an.rds")


# Get meta data
# make sure you run: banctable_set_token(user = "", pwd = "") to generate banc token
bc <- banctable_query()
meta <- bc %>%
  dplyr::mutate(cell_type = dplyr::case_when(
    is.na(cell_type) ~ gsub("auto\\:","",fafb_cell_type),
    is.na(cell_type) ~ gsub("auto\\:","",manc_cell_type),
    TRUE ~ cell_type
  ),
  top_nt = gsub("auto\\:","",top_nt),
  cell_sub_class = gsub("auto\\:","",cell_sub_class),
  cell_class = gsub("auto\\:","",cell_class),
  super_class = gsub("auto\\:","",super_class)) %>%
  dplyr::filter(region%in%c("neck_connective"),
                !is.na(side), 
                super_class=="ascending")
# Read L2 skeletons for NBLASTing (inaccurate but fast to acquire)

# Get all ids
choose_banc()
ids <- unique(meta$root_id)
ids <- na.omit(ids)
ids <- ids[ids!="0"]
ids <- banc_latestid(ids)
ids.left <- unique(subset(meta, side=="left")$root_id)
ids.left <- banc_latestid(ids.left)
ids.left <- na.omit(ids.left)
ids.left <- ids.left[ids.left!="0"]
l2dps <- banc_read_l2dp(ids)
# l2dps[ids.left] <- banc_mirror(l2dps[ids.left])

# Get skeletons for visualisation
# Get skeletons for visualization (units in nanometers)
l2 <- banc_read_l2skel(ids)
# l2[ids.left] <- banc_mirror(l2[ids.left])

# mirror neurons to right hemisphere (mirror require units in nanometers)
l2[ids.left] <- banc_mirror(l2[ids.left])

# remove probelm
# Assuming your neuronlist is named 'l2'
l2 <- remove(l2, "720575941628870438")

# Make dotprops for NBLASTing (inaccurate but fast to acquire) (units in microns)
l2dps <- dotprops(l2/1000)
 
# Re-root to soma where this is known
# banc.roots <- banc_roots()
banc.roots <- bancr:::banc_roots()
l2 <- banc_reroot(l2, roots = banc.roots, estimate = TRUE)

# Save the state of this R session, to avoid doing the step above next time you load this file
save.image()

########################
### visualizing data ###
########################

# plot all neurons (mirrored to right side)
open3d(windowRect = c(20, 30, 1000, 1000))
banc_view()
plot3d(banc_neuropil.surf, alpha = 0.1, col = "lightgrey")
plot3d(l2, lwd=0.5, soma = 5000)

figure.dir <- str_c(here::here(), "/figures/")
snapshot3d(str_c(figure.dir, "allANs_skt.png"))
close3d()

# plot example neuron (DNg34)
AN_GNG_175 <- subset(meta, cell_type == "auto:AN_GNG_175")
AN_GNG_175.ids <- AN_GNG_175$root_id

open3d(windowRect = c(20, 30, 1000, 1000))
wire3d(banc_neuropil.surf/1e3, alpha = 0.1, col='lightgrey')
plot3d(l2dps[AN_GNG_175.ids[1]], lwd=3, soma = 5000, col = 'red')
plot3d(l2dps[AN_GNG_175.ids[2]], lwd=100, soma = 5000, col = 'blue')

# Check NBLAST distance for example neuron
# DNg34 <- l2dps[DNg34.ids]
# nb.all.sub <- nat.nblast::nblast_allbyall(DNg34, normalisation = "mean")

###############
### NBLAST  ###
###############

# Run NBLAST
nb.all.full <- nat.nblast::nblast_allbyall(l2dps, normalisation = "mean")

# Save neuronlist using saveRDS()
# saveRDS(l2, file = "l2_an.rds")
# saveRDS(l2dps, file = "l2dps_an.rds")
# saveRDS(nb.all.full, file = "nb.all.full_an.rds")

# Cluster
hc <- nhclust(scoremat = nb.all.full, method = "ward.D2")

# Flatten
hc$height=sqrt(hc$height)
#hc$height = sqrt(hc$height)

# Color clusters
cut.height <- 1.9 # needs optimisation
# Split data into clusters and color clusters (each with a different color)
cut.height <- 1.9 # needs optimization
hc.col <- dendroextras::color_clusters(hc, h = cut.height)
hc.col <- dendextend::set_labels(hc.col,  paste0(meta[match(labels(hc.col),meta$root_id),"cell_type"],"_",meta[match(labels(hc.col),meta$root_id),"side"]))

# Plot
plot(hc.col, labels = T)
abline(h=cut.height)
# Set custom labels
hc.col <- dendextend::set_labels(hc.col, 
                                 paste0(meta[match(labels(hc.col), meta$root_id), "cell_type"], 
                                        "_", meta[match(labels(hc.col), meta$root_id), "side"]))

# See neurons
clear3d()
banc_view()
plot3d(hc, db = l2, h = cut.height, soma = 5000)
##############################
### visualizing clustering ###
##############################

# figure settings
width_in_inches <- 20
height_in_inches <- 20
dpi <- 100
width <- width_in_inches * dpi
height <- height_in_inches * dpi

# 1) plot and save dendrogram without labels
figure.dir <- str_c(here::here(), "/figures/")
png(filename = str_c(figure.dir, "AllDNs_clus_fulldend.png"), width = width, height = height, res = dpi)
plot(hc.col, ylab = "Height", leaflab = "none")
# to add labels run:
plot(hc.col, , labels = T)
abline(h = cut.height)
dev.off()

# 2) custom plot clusters (select one in red and the rest in black)
hc.col2 <- as.dendrogram(hc)

# Get groups/clusters
groups <- cutree(hc.col2, h = cut.height)

other_color <- "#000000"
target_color <- "#FF0000"
target_group <- 1
temp_col <- ifelse(groups == target_group, target_color, other_color)
temp_col <- temp_col[order.dendrogram(hc.col2)]
temp_col <- factor(temp_col, unique(temp_col))

hc.col2 <- hc.col2 %>% 
  color_branches(hc.col2, clusters = as.numeric(temp_col), col = levels(temp_col)) %>% 
  set("labels_colors", as.character(temp_col))

# Get groups
# Set custom labels
hc.col2 <- dendextend::set_labels(hc.col2, 
                                  paste0(meta[match(labels(hc.col2), meta$root_id), "cell_type"], 
                                         "_", meta[match(labels(hc.col2), meta$root_id), "side"]))

# plot and save dendrogram without labels
figure.dir <- str_c(here::here(), "/figures/")
png(filename = str_c(figure.dir, "AllANs_clus_", as.character(target_group), "_fulldend.png"), 
    width = width, height = height, res = dpi)
plot(hc.col, ylab = "Height", leaflab = "none")
# to add labels run:
plot(hc.col, , labels = T)
abline(h = cut.height)
dev.off()

##############################################################################
### Split into clusters and plot each cluster (ordered from left to right) ###
##############################################################################

# Get groups/clusters
groups <- cutree(hc, h = cut.height)
table(groups)

###########################
### Manual Cell typing  ###
###########################
# Reorder clusters to have the leftmost as cluster 1
# Get the order of the labels in the dendrogram
order_of_labels <- order.dendrogram(hc.col2)

# Get a vector of group labels ordered by the dendrogram
ordered_group_labels <- groups[order_of_labels]
table(ordered_group_labels)

# Find unique groups in their left-to-right order
unique_groups <- unique(ordered_group_labels)

choices <- list()
other_color <- "#000000"
target_color <- "#FF0000"
figure.dir <- str_c(here::here(), "/figures/AN_nblast/")

if (!dir.exists(figure.dir)) {
  dir.create(figure.dir, recursive = TRUE)
}



# clear an_type column
# bc$an_type <- ""
an_types <- subset(bc, super_class=="ascending")

choices <- list()
#(1:3)
next_clust <- length(unique(bc$an_type))
for(cluster in unique(groups)){
  # get ordered cluster number (from left to right)
  cluster.ordered <- unique_groups[cluster]
  
  # Get members 
  # members <- names(groups[groups==cluster])
  members <- names(groups[groups == cluster.ordered])
  neurons <- l2[members]
  
  
  # add to an_type
  #bc$auto_an_type[bc$root_id %in% members] <- as.character(cluster.ordered)
  
  # plot selected cluster in full dendrogram
  temp_hc <- as.dendrogram(hc)
  temp_col <- ifelse(groups == cluster.ordered, target_color, other_color)
  temp_col <- temp_col[order.dendrogram(temp_hc)]
  temp_col <- factor(temp_col, unique(temp_col))
  temp_hc <- temp_hc %>% 
    color_branches(temp_hc, clusters = as.numeric(temp_col), col = levels(temp_col)) %>% 
    set("labels_colors", as.character(temp_col))
  temp_hc <- dendextend::set_labels(temp_hc, 
                                    paste0(meta[match(labels(temp_hc), meta$root_id), "cell_type"], 
                                           "_", meta[match(labels(temp_hc), meta$root_id), "side"]))
  
  width_in_inches <- 20
  height_in_inches <- 10
  dpi <- 100
  width <- width_in_inches * dpi
  height <- height_in_inches * dpi
  output_file <- str_c(figure.dir, 'DN_nblast_clus_', cluster, "_fulldend.png")
  png(filename = output_file, width = width, height = height, res = dpi)
  
  plot(temp_hc, ylab = "Height", leaflab = "none")
  abline(h=cut.height)
  dev.off()
  
  # Plot new sub-dendrogram
  nb.all.sub <- nb.all.full[members,members]
  nb.all.sub <- nb.all.full[members, members]
  hc.sub <- nhclust(scoremat = nb.all.sub, method = "ward.D2")
  hc.sub$height <- sqrt(hc.sub$height)
  cut.height.sub <- 1.4 
  hc.sub.col <- dendroextras::color_clusters(hc.sub, h = cut.height.sub)
  hc.sub.col <- dendextend::set_labels(hc.sub.col, 
                                       paste0(meta[match(labels(hc.sub.col), meta$root_id), "cell_type"], 
                                              "_", meta[match(labels(hc.sub.col), meta$root_id), "side"]))
  plot(hc.sub.col, labels = T)
  abline(h=cut.height.sub)
    
  dpi <- 300
  width <- width_in_inches * dpi
  height <- height_in_inches * dpi
  output_file <- str_c(figure.dir, 'AN_nblast_clus_', cluster, "_subdend.png")
  png(filename = output_file, width = width, height = height, res = dpi)
  
  par(mar = c(3, 1, 1, 10))
  plot(hc.sub.col, labels = T, horiz = TRUE)
  par(mar = c(5, 4, 4, 2) + 0.1)
  dev.off()
  
  # Plot neurons
  clear3d()
  open3d(windowRect = c(20, 30, 1000, 1000))
  banc_view()
  plot3d(banc_neuropil.surf, alpha = 0.1, col = "lightgrey")
  plot3d(hc.sub, db=neurons, lwd=0.5, soma = 2500, h=cut.height.sub)
  
  output_file <- str_c(figure.dir, 'AN_nblast_clus_', cluster, "_skt.png")
  snapshot3d(output_file)
  close3d()
  
  # Get choices
  db <- nlscan(neurons, col = "black", lwd = 3, soma = 5000)
  choices[[cluster]] <- db
  if (length(choices[[cluster]])> 0 ){
     next_clust <- next_clust + 1
     bc$an_type[bc$root_id %in% choices[[cluster]]] <- as.character(next_clust)
   }

  # # Pause for user
  message(paste(choices,collapse=", "))
  readline("prompt: press any key to continue to the next cluster ")

}

# add to unique groups
addons <- c(42:next_clust)
unique_groups_add <- c(unique_groups, addons)

# Record manually entries here
# skipped cluster with 90 cells - 15
# skipped cluster last one - 41









# all cluster 1 good
# weird group 2 720575941454113773 - maybe L-R flip error
# all cluster 3 good
# split cluster 4 into 2 (this is one) c("720575941475559274", "720575941687777499", "720575941504240407", "720575941595780800", "720575941592508453", "720575941605024365", "720575941453611315")
# cluster 5 into at least 2 c("720575941600940473", "720575941547007624", "720575941545246396", "720575941446655252", "720575941550685063", "720575941533863194", "720575941492153214", "720575941450406825", "720575941505656866", "720575941501922498", "720575941516666584", "720575941720197690", "720575941558884078", "720575941602726134"),
# cluster 6 into at least 2 c("720575941477208426", "720575941623204298", "720575941474958401", "720575941700378074", "720575941548344637", "720575941475546673", "720575941529838821", "720575941576989173", "720575941574020534", "720575941612928019", "720575941532995521", "720575941731066923", "720575941438814043", "720575941590564203")
# cluster 7 into at least 2 c("720575941439043279", "720575941438594395", "720575941668513201", "720575941487747125", "720575941445166100", "720575941429400535", "720575941529998821", "720575941413115281", "720575941566804218", "720575941629121004", "720575941492467070", "720575941474315294", "720575941513720979", "720575941551121179", "720575941505059397", "720575941593486806", "720575941655872536", "720575941573970831", "720575941513587779", "720575941551181595")
# cluster 8 into at least 2 c("720575941553575681", "720575941521385291", "720575941550649599", "720575941501806786", "720575941652038417", "720575941446713108", "720575941438012479", "720575941533541720", "720575941689321112", "720575941584446598", "720575941546678152", "720575941514790137", "720575941539124861", "720575941562991511", "720575941525280166", "720575941502550359", "720575941630674848", "720575941437864511", "720575941424975956", "720575941430392078", "720575941544537012", "720575941523892324", "720575941612537241", 
# "720575941512515907")
# 9 split c("720575941617513855", "720575941546765601", "720575941608989869") ???
# 10 split c("720575941528671179", "720575941614995157", "720575941601256889", "720575941549208063", "720575941649660145", "720575941512656593")
# 11 split c("720575941526290086", "720575941514355802", "720575941452381747", "720575941521998411", "720575941614138265", "720575941533807338", "720575941689278939")
# 12 split c("720575941505559621", "720575941574492086", "720575941455137005", "720575941571900232", "720575941450363770", "720575941602242750", "720575941613446809", "720575941491803056", "720575941534977562", "720575941428996856", "720575941460382352", "720575941521613559", "720575941665374451", "720575941614894037", "720575941438936155", "720575941552480257", "720575941576989009", "720575941521521227", "720575941637686389")
# 13 none
# 14 split c("720575941593751077", "720575941547360829", "720575941515109267", "720575941497484257", "720575941534809098", "720575941562512018", "720575941515899820", "720575941451837865", "720575941610664621", "720575941582936927", "720575941668474545", "720575941552836762")
# 15 need to deal with later
# 16 no split
# 17 split c("720575941529735883", "720575941533250738", "720575941553307364", "720575941438781135", "720575941602404598", "720575941491974320", "720575941543533756", "720575941543838140", "720575941566550522", "720575941400300835")
# 18 no split
# 19 no split
# 20 split c("720575941501031511", "720575941443454634", "720575941403642992", "720575941622910716")
# 21 split need to deal with later
# 22 split c("720575941516028844", "720575941567097588", "720575941575462504", "720575941495978569", "720575941534751681", "720575941492454270", "720575941668617649", "720575941406964142", "720575941638173471", "720575941668877735", "720575941448837725", "720575941667034355", "720575941540978268", "720575941424935508", "720575941554156883", "720575941468939072", "720575941638028661", "720575941551318043", "720575941687317004", "720575941602232766", "720575941472526176")
# 23 split c("720575941474428980", "720575941587840324", "720575941437865023", "720575941540971356", "720575941639214367", "720575941474562334", "720575941547338557", "720575941571681693", "720575941472322912", "720575941536941461", "720575941528624075", "720575941532374360", "720575941573319945", "720575941539581532", "720575941684989935", "720575941679005885", "720575941572120136")
# 24 no split
# 25 c("720575941405855918", "720575941578442741", "720575941537692365", "720575941557128996", "720575941571120030", "720575941652049881", "720575941661382780", "720575941537872253", "720575941558552591", "720575941428247800", "720575941565317422", "720575941612174044", "720575941588672652", "720575941536268437", "720575941429366369", "720575941438669653", "720575941549063167", "720575941436503359", "720575941515403667", "720575941606245229", "720575941477603296", "720575941429071753", "720575941650621461")
# 26 no split
# 27 need to deal with later
# 28 c("720575941484283654", "720575941593607126", "720575941479644067", "720575941638071669", "720575941599321385", "720575941469844855", "720575941594043173", "720575941546884641", "720575941649949937")
# 29 c("720575941426287700", "720575941548674932", "720575941478958722", "720575941526154628", "720575941533857818", "720575941500465547", "720575941572713467", "720575941476052276", "720575941603458750", "720575941455227885", "720575941484963394", "720575941448941902", "720575941573099677", "720575941413008276", "720575941602408950", "720575941534431754", "720575941490333927", "720575941655935256")
# 30 keep for now
# 31 c("720575941642581960", "720575941652414737", "720575941623094730", "720575941529729227", "720575941468212189", "720575941502082027", "720575941653794009", "720575941610678445", "720575941520133367", "720575941478248015", "720575941474197278", "720575941390437270", "720575941536219107", "720575941466022422", "720575941498208162", "720575941496290889", "720575941495195008", "720575941597314496", "720575941729786411", "720575941667205553", "720575941485128770", "720575941651061489", "720575941437905727", 
# "720575941439389887", "720575941657709136", "720575941458480083", "720575941390400662", "720575941581527903", "720575941621749244", "720575941525340838", "720575941482336698")
# 32 keep fo now
# 33 c("720575941429620617", "720575941514107385", "720575941517626028", "720575941429604449", "720575941525368834", "720575941564118423", "720575941546861448", "720575941720487482")
# 34 c("720575941589176830", "720575941643643256", "720575941552806810", "720575941480387523", "720575941507028770", "720575941479272783", "720575941668465063", "720575941433883790", "720575941721899322", "720575941413640732", "720575941642933586")
# 35 c("720575941534989082", "720575941517277612", "720575941391860118", "720575941411659412", "720575941610588380", "720575941564190871", "720575941549086463", "720575941428464047")
# 36 c("720575941489638094", "720575941495189888", "720575941552815514", "720575941731063595", "720575941637843743", "720575941442194346", "720575941505698839", "720575941594043429")
# 37 keep for now - messy
# 38 no split
# 39 c("720575941429273559", "720575941493524167", "720575941622780412", "720575941537676925", "720575941668452007", "720575941657553744")
# 40 keep
# 41 keep

# round 2

an_types <- subset(an, super_class=="ascending")
i=1
for(cluster in unique_groups_add[39:68]){
  
  cluster_plot <- subset(bc, an_type==as.character(cluster))
  members <- cluster_plot$root_id
  neurons <- l2[unique(members)]
  
  clear3d()
  open3d(windowRect = c(20, 30, 1000, 1000))
  banc_view()
  plot3d(banc_neuropil.surf, alpha = 0.1, col = "lightgrey")
  plot3d(neurons)
  
  output_file <- str_c(figure.dir, 'AN_nblast_clus_', i, "_skt.png")
  snapshot3d(output_file)
  close3d()
  i=i+1
  
  # # Get choices
  # db <- nlscan(neurons, col = "black", lwd = 3, soma = 5000)
  # choices[[cluster]] <- db
  # if (length(choices[[cluster]])> 0 ){
  #   next_clust <- next_clust + 1
  #   #bc$an_type[bc$root_id %in% choices[[cluster]]] <- as.character(next_clust)
  # }
  
  # # Pause for user
  # message(paste(choices,collapse=", "))
  # readline("prompt: press any key to continue to the next cluster ")
}




# Upload to Seatable! -----------------------------------------------------
# 
#subset bc to just ANs before uploading to seatable
an_types <- subset(bc, super_class=="ascending")
an_types <- as.data.frame(an_types)
an_types[is.na(an_types)] = ""
an_type_cleaned <- subset(an_types, `_id` != "deru3mgEShCuanWx8Xcx8w MRRFPVRITvCT3D6mgjeCig")
banctable_update_rows(df = an_type_cleaned[,c("_id", "root_id", "an_type")], table = "banc_meta", base = "banc_meta")


# investigate clusters ----------------------------------------------------
types <- as.numeric(unique(bc$an_type))
types <- na.omit(types)

an1.left <- banc_latestid(an1.left)
an1.right <- banc_latestid(an1.right)



choices <- list()
pick = c(19)
#merge 7 and 62
#bc$an_type[bc$an_type == 49] <- 101

# 68, 45, 25 are similar class different legs

for(cluster in types){
  # get ordered cluster number (from left to right)
  members <- subset(bc, an_type == cluster)
  member_names <- members$root_id
  member_names <- banc_latestid(member_names)
  neurons <- l2[member_names]
  neurons <- neurons[sapply(neurons, nat::is.neuron)]
  neurons <- neurons[!duplicated(names(neurons))]
  # 
  # members2 <- subset(bc, an_type == 15)
  # member_names2 <- members2$root_id
  # member_names2 <- banc_latestid(member_names2)
  # neurons2 <- l2[member_names2]
  # neurons2 <- neurons2[sapply(neurons2, nat::is.neuron)]
  # neurons2 <- neurons2[!duplicated(names(neurons2))]
  # 
  # members3 <- subset(bc, an_type == 16)
  # member_names3 <- members3$root_id
  # member_names3 <- banc_latestid(member_names3)
  # neurons3 <- l2[member_names3]
  # neurons3 <- neurons3[sapply(neurons3, nat::is.neuron)]
  # neurons3 <- neurons3[!duplicated(names(neurons3))]
  
  # Plot neurons
  clear3d()
  open3d(windowRect = c(20, 30, 1000, 1000))
  banc_view()
  plot3d(banc_neuropil.surf, alpha = 0.1, col = "lightgrey")
  plot3d(neurons, lwd=0.5, soma = 2500, col="red")
  #plot3d(neurons2, lwd=0.5, soma = 2500, col="blue")
  #plot3d(neurons3, lwd=0.5, soma = 2500, col="green")

  output_file <- str_c(figure.dir, 'AN_nblast_clus_', cluster, "_skt.png")
  snapshot3d(output_file)
  close3d()
  
  # # Get choices
  # db <- nlscan(neurons, col = "black", lwd = 3, soma = 5000)
  # choices[[cluster]] <- dbx
  # 
  # # # Pause for user
  # message(paste(choices,collapse=", "))
  # readline("prompt: press any key to continue to the next cluster ")
  
}






# Plot new dendrogram -----------------------------------------------------
# recent nb.all.full names
row_names <- rownames(nb.all.full)
# Apply the 'banc_id_isrecent()' function to each row name
row_names <- banc_latestid(row_names)

# Rename the rows of the dataframe
rownames(nb.all.full) <- row_names
colnames(nb.all.full) <- row_names

# Function to remove duplicate rows or columns
remove_duplicates <- function(mat) {
  # Get unique row/column names
  unique_names <- unique(rownames(mat))
  
  # Keep only the first occurrence of each unique row/column
  mat <- mat[match(unique_names, rownames(mat)), ]
  mat <- mat[, match(unique_names, colnames(mat))]
  
  return(mat)
}

# Remove duplicate rows and columns
similarity_matrix <- remove_duplicates(nb.all.full)

# Print the dimensions of the new matrix
print(dim(similarity_matrix))

an_working <- subset(bc, root_id %in% row_names)
an_working <- an_working[!duplicated(an_working$root_id), ]
#an_working$an_type
an_working <- an_working[, c("root_id", "an_type")]

nb_df <- as.data.frame(similarity_matrix)
nb_df$root_id <- rownames(similarity_matrix)

merged_df <- nb_df %>% left_join(an_working, by = "root_id")
collapsed_df <- merged_df %>%
  group_by(an_type) %>%
  summarise(across(where(is.numeric), max, na.rm = TRUE),
            root_id = first(root_id)) %>%
  ungroup()


# Assuming 'an_working' contains the mapping of root_id to an_type
neuron_to_type <- an_working %>% select(root_id, an_type)


# Remove root_id if it exists
collapsed_df <- collapsed_df %>% select(-any_of("root_id"))
library(tidyr)

type_by_type <- collapsed_df %>%
  # Convert to long format
  pivot_longer(cols = -an_type, names_to = "neuron_id", values_to = "value") %>%
  # Join with the neuron_to_type mapping
  left_join(neuron_to_type, by = c("neuron_id" = "root_id")) %>%
  # Group by row type and column type
  group_by(an_type.x, an_type.y) %>%
  # Sum the values
  summarise(total_value = max(value, na.rm = TRUE), .groups = "drop") %>%
  # Spread back to wide format
  pivot_wider(names_from = an_type.y, values_from = total_value, values_fill = 0) %>%
  # Rename the row type column
  rename(an_type = an_type.x)


type_by_type <- subset(type_by_type, !is.na(an_type))

# Assuming type_by_type is your current dataframe
type_by_type <- type_by_type %>%
  column_to_rownames(var = "an_type")
type_by_type <- type_by_type[, !grepl("NA", colnames(type_by_type), ignore.case = TRUE)]


# now should make dendrogram

# just dendrogram
# Perform clustering
dist_matrix <- dist(type_by_type)
hc <- hclust(dist_matrix, method = "complete")  # You can change the method as needed

# Plot the dendrogram
plot(hc, main = "Dendrogram of Type-by-Type Matrix, Mean", 
     xlab = "", sub = "", 
     hang = -1)  # hang = -1 makes the labels appear at the same level
library(pheatmap)

# Create a heatmap with dendrograms
pheatmap(type_by_type,
         main = "Clustered Heatmap of Type-by-Type Matrix",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         fontsize = 8,
         cellwidth = 10,
         cellheight = 10,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete")