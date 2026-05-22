#!/bin/bash
###########################################################
### Publish per-version syn_id -> neuropil/region/side
### lookups to GCS.
###
### Builds slim parquets for v1, v2, v3 from on-disk sources
### and uploads to gs://lee-lab.../synapses/v{1.1,2.0,3.0}/.
###
### Usage: sbatch o2_banc_publish_lookups.sh
###########################################################
#SBATCH -c 2
#SBATCH -t 0-01:00
#SBATCH -p short
#SBATCH --mem=128G
#SBATCH -o /home/ab714/bancpipeline/jobs/publish_lookups_%j.out
#SBATCH -e /home/ab714/bancpipeline/jobs/publish_lookups_%j.err

cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0

Rscript banc/share/banc-publish-synapse-lookups.R
