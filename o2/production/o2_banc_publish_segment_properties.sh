#!/bin/bash
#SBATCH -c 2
#SBATCH -t 0-01:00
#SBATCH -p short
#SBATCH --mem=8G
#SBATCH -o /home/ab714/bancpipeline/jobs/banc_publish_segprop_%j.out
#SBATCH -e /home/ab714/bancpipeline/jobs/banc_publish_segprop_%j.err
###############################################################################
# Standalone publisher for neuroglancer segment_properties.
#
# Posts tag-rich segment_properties/info to each per-dataset GCS layer using
# the helpers in banc/transforms/banc-ngl-upload.R. Skips mesh uploads.
###############################################################################
set -euo pipefail
cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0
echo "=== publish segment_properties to GCS ==="
Rscript banc/transforms/banc-publish-segment-properties.R
echo "=== Done ==="
