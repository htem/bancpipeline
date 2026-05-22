#!/bin/bash
#SBATCH -c 10                              # Request cores
#SBATCH -t 4-00:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                          # Partition to run in
#SBATCH --mem-per-cpu=12G                  # Memory per core (bumped from 8G after 38901939 OOM 2026-05-09)
#SBATCH -o /home/ab714/bancpipeline/jobs/banc_update_%j.out         # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e /home/ab714/bancpipeline/jobs/banc_update_%j.err         # File to which STDERR will be written, including job ID (%j)

echo "RUNNING BANC UPDATE"

start=`date +%s`

cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0  # suppress core dumps from R shutdown crashes

# Data acquisition
echo "sorting results"
Rscript banc/meta/banc-sort-folders.R

echo "fetching banc root IDs"
Rscript banc/update/banc-ids.R

echo "deleting outdated BANC data"
Rscript banc/update/banc-delete.R

echo "saving banc l2 skeletons"
Rscript banc/metrics/banc-l2.R

# Calculate metrics
echo "calculating root positions"
Rscript banc/metrics/banc-calculate-root-positions.R

echo "calculating regions"
Rscript banc/metrics/banc-calculate-regions.R

echo "calculating synapses"
Rscript banc/metrics/banc-calculate-synapses.R

echo "reading l2 metrics"
Rscript banc/metrics/banc-calculate-l2-metrics.R

echo "calculating volumes"
Rscript banc/metrics/banc-calculate-volumes.R

# Skeletonisation and splitting
echo "getting BANC meshes"
Rscript banc/metrics/banc-obj.R

echo "skeletonising BANC neurons"
Rscript banc/metrics/banc-calculate-skeletons.R

echo "splitting BANC neurons"
Rscript banc/metrics/banc-calculate-split.R

# NBLAST data compilation
echo "handling wrong matches"
Rscript banc/nblast/banc-nblast-wrong-matches.R

echo "compiling NBLAST results"
Rscript banc/nblast/banc-nblast-compile.R

echo "sharing NBLAST to GCS"
Rscript banc/share/banc-nblast-share-gcs.R

echo "building banc meta data"
Rscript banc/meta/banc-meta.R

# Push to seatable
echo "pushing metrics to seatable"
Rscript banc/update/banc-update-metrics.R

echo "updating match columns"
Rscript banc/update/banc-update-matches.R

echo "updating cell types"
Rscript banc/update/banc-update-celltypes.R

echo "updating status flags"
Rscript banc/update/banc-update-status.R

echo "updating seatable (identity + new rows)"
Rscript banc/update/banc-update-seatable.R

echo "Re-plot our progress"
Rscript banc/utilities/banc-plot.R
Rscript banc/nblast/banc-nblast-plot.R

# Git operations removed — run manually after reviewing job output

# Report
end=`date +%s`
runtime=$((end-start))
echo "script completed in: "
echo $runtime

