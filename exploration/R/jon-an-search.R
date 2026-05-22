### -----------------------------------------------------------------
### Legacy exploratory analysis by a member of the Wilson Lab.
### Historical record only — not part of the bancpipeline release.
### -----------------------------------------------------------------

# AVLP JO and AN convergence


library(dplyr)
library(bit64)
# Load packages
library(bancr)
library(tidyverse)
library(nat.nblast)
library(dendextend)

# Make sure all functions query BANC and not FAFB
choose_banc()

# Initialize parameters---------------------------------------------------------

# banc data
bc <- banctable_query()

# ans
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

# fetch all DNs
# fetch all ANs
an_meta = flytable_meta('AN.*')
job_meta = flytable_meta('type:JO-B')
joa_meta = flytable_meta('type:JO-A')
jo_meta = rbind(joa_meta, job_meta)
nJO = dim(jo_meta)[1]
# remove problem JONs

# Initialize an empty list to store root_ids that cause errors
remove_jo <- c()

for (d in 1:nJO) {
  tryCatch(
    {
      # Fetch JON
      thisJON <- jo_meta[d,]
      
      print(paste('Working on', thisJON$root_id, '...'))
      # Find all outputs of JON
      outputs_JON <- flywire_partner_summary2(thisJON$root_id, partners = 'out', threshold = min_syn)
    },
    error = function(e) {
      print("An error occurred with root_id:")
      print(thisJON$root_id)
      
      # Add the root_id to the remove_jo list
      remove_jo <<- c(remove_jo, thisJON$root_id)
    }
  )
}

# After the loop, remove all root_ids in remove_jo from jo_meta
jo_meta <- jo_meta[!jo_meta$root_id %in% remove_jo, ]




# ensure only descending neurons are included
an_meta = an_meta[grepl("ascending",an_meta$super_class),]
an_omit = c('720575940604705481', "720575940628415951", "720575940644598324", "720575940635837786", "720575940615900189", "720575940627430906", '720575940639264590', "720575940618092965", "720575940612294311", "720575940629292842", '720575940645264260', '720575940621949675', "720575940650261494", '720575940620262621', "720575940621978109")
for (o in 1:length(an_omit)){
  an_meta = an_meta[!grepl(an_omit[o],an_meta$root_id),]
}




# Find cell types that input to BOTH copies of a given DN-----------------------
print('Finding cell types downstream of both JONs and ANs...')


min_syn = 3;

outputs_ans <- data.frame()
for (l in 1:nAN){
  print(paste('Working on', an_meta[l,]$root_id, '...'))
  output_an = flywire_partner_summary2(an_meta[l,]$root_id, partners = 'out', threshold = min_syn)
  outputs_ans = rbind(outputs_ans, output_an)
}

# initialize
shared_tracker = c()
downstream_info <- data.frame()
dn_name = c()

nJO <- length(unique(jo_meta$root_id))

# for each JO
shared_tracker <- c()
d<-1
for (d in 1:nJO){ #
  # fetch JON
  thisJON = jo_meta[d,]
  # number of ANs
  nAN = dim(an_meta)[1]
  
  print(paste('Working on',thisJON$root_id,'...'))
  # find all outputs of jon
  outputs_JON = flywire_partner_summary2(thisJON$root_id, partners = 'out', threshold = min_syn)
  # for each left copy (if there is one)
  ids <- unique(outputs_ans$pre_pt_root_id)
  for (l in 1:nAN){
    # find all outputs of ANs
    # outputs_an = flywire_partner_summary2(an_meta[l,]$root_id, partners = 'out', threshold = min_syn)
    this_AN <- subset(outputs_ans, pre_pt_root_id == an_meta[l,]$root_id)
    
    # find shared inputs
    shared_outputs = nat::intersect(as.character(outputs_JON$post_pt_root_id),as.character(this_AN$post_pt_root_id))
    shared_output_int <- as.integer64(shared_outputs)
    nShared = length(shared_outputs)
    # if there were shared inputs, store
    if (nShared>0){
      shared_inputs_all_info <- subset(outputs_JON, post_pt_root_id %in% shared_output_int)
      shared_inputs_an <- subset(this_AN, post_pt_root_id %in% shared_output_int)
      shared_inputs_all_info$jon_weight <- shared_inputs_all_info$weight
      shared_inputs_all_info$an_weight <- shared_inputs_an$weight
      shared_inputs_all_info$an_name <- shared_inputs_an$pre_pt_root_id
      #shared_tracker = c(shared_tracker,t(shared_inputs))
      downstream_info <- rbind(downstream_info, shared_inputs_all_info)
    }
  }
}

# analyze big fry data ----------------------------------------------------

# count instances of single IDs
downstream_freq <- count(downstream_info, post_pt_root_id)
downstream_freq <- downstream_freq[order(-downstream_freq$n),]

# sort upstream info by synapse counts
downstream_info <- downstream_info[order(-downstream_info$weight),]
# remove 3-5 synapse ANs
downstream_info_crop <- subset(downstream_info, an_weight >10)


#collapse by cell ids - may want to better titrate how data is summarized
collapsed_downstream <- downstream_info_crop %>%
  group_by(pre_pt_root_id) %>%
  summarize(
    cell_type=first(cell_type),
    across(!cell_type , ~ toString(.))
  )




#collapse by cell type
collapsed_downstream2 <- collapsed_downstream %>%
  group_by(cell_type) %>%
  summarize(
    across(, ~ toString(.))
  )


# pull ANs specifically
total_ans <- unique(downstream_info_crop$an_name)

# now look at who those ANs are in the stream of things
jo_ans <- subset(an_meta, an_meta$root_id %in% total_ans)

# collapse over cell type
jo_ans <- jo_ans %>%
  group_by(cell_type) %>%
  summarize(
    across(, ~ toString(.))
  )