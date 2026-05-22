### -----------------------------------------------------------------
### Legacy exploratory analysis by a member of the Wilson Lab.
### Historical record only — not part of the bancpipeline release.
### -----------------------------------------------------------------

# Code to explore morphological classes of DNs

# 1) Load banctable and save DN metadata, 
#   DN skeletons & dotprops as "meta_dn.rds", "l2_dn.rds", and "l2dps_dn.rds"
#   all-to-all nblast scores as "nblast_score_all_dn.rds"
dn_gen_save_data.R

#   synapses per neuropil, abs number, norm number, number of neurons
dn_synapse_per_neuropil_save_data.R

# 2) Cluster DNs using morphological similarity NBLAST scores

# use various distances ([3, 2.4, 1.9, 1.5]) and explore the pattern of splitting
dn_nblast_clusters.R

# assign cluster IDs to DNs
dn_nblast_cluster_assign.R

# test that assigment in seatable is correct (match plots)
dn_nblast_cluster_assiggn_test.R

# plot full dendrogram and neurons colorcoded by cluster (using a arbitrary initial distance of 1.9)
dn_nblast_full_dendrogram.R

# plot full dendrogram and neurons colorcoded by cluster (using a arbitrary initial distance of 1.9)
#   plus neurotransmitter info
dn_nblast_cluster_1-9_with_nt.R

# plot full dendrogram and neurons colorcoded by cluster (using a arbitrary initial distance of 1.9)
#   and add plots of subclusters
dn_nblast_cluster_1-9_and_subclusters.R

# start comparing NBLAST vs connectivity information
#   load cosine distance and match with morphological dn types
dn_connectivity_save_data.R

#   plot cosine distance vs morphological dn types
dn_connectivity_plot.R