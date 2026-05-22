#!/bin/bash
###########################################################
### Full v850 data rebuild for BANC pipeline
###
### Produces all versioned data files and pushes to GCS.
### Prerequisite: v850 synapse CSV must be on GCS at:
###   gs://lee-lab_brain-and-nerve-cord-fly-connectome/v850/synapses_v2_human_readable.csv.gz
###
### Usage:
###   sbatch -c 10 -t 0-24:00 -p medium --mem-per-cpu=8G \
###          -o jobs/v850_rebuild_%j.out -e jobs/v850_rebuild_%j.err \
###          o2/o2_banc_v850_rebuild.sh
###
### Steps:
###   1. Update SeaTable IDs (root_id, root_626, root_850, supervoxel_id)
###   2. Download L2 skeletons for any new/updated neurons
###   3. Calculate v850 connectivity (synapse parquet, edgelist, basic meta)
###  3b. Extract version-agnostic synapse neuropil lookup (reuse prior assignments)
###   4. Calculate neuropil inclusion for synapses (lookup-accelerated)
###   5. Calculate NT predictions from v850 synapse data
###   6. Calculate L2 metrics (cable length, node count)
###   7. Calculate root positions and regions
###   8. Calculate volumes
###   9. Calculate synapse split tables
###  10. Compile versioned data (meta, enriched synapses, edgelists)
###  11. Compile NBLAST feathers (with root_850 column)
###  12. Push NT predictions to SeaTable
###  13. Push NBLAST feathers to GCS
###  14. Clean up stale files
###########################################################

#SBATCH -c 1
#SBATCH -t 0-600:00
#SBATCH -p priority
#SBATCH --mem-per-cpu=250G
#SBATCH -o jobs/v850_rebuild_%j.out
#SBATCH -e jobs/v850_rebuild_%j.err

set -e
start=$(date +%s)
cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0  # suppress core dumps from R shutdown crashes
mkdir -p jobs

echo "=========================================="
echo "  BANC v850 rebuild started: $(date)"
echo "=========================================="

# Step 1: Update SeaTable IDs (root_id, root_626, root_850, supervoxel_id)
echo "=== Step 1: Updating SeaTable IDs ==="
Rscript banc/update/banc-ids.R
echo "  Done: $(date)"

# Step 2: Download L2 skeletons
echo "=== Step 2: Downloading L2 skeletons ==="
Rscript banc/metrics/banc-l2.R
echo "  Done: $(date)"

# Step 3: Calculate v850 connectivity
echo "=== Step 3: Calculating v850 connectivity ==="
Rscript banc/metrics/banc-calculate-connectivity.R
echo "  Done: $(date)"

# Step 3b: Extract synapse neuropil lookup from prior version data
echo "=== Step 3b: Extracting synapse neuropil lookup ==="
Rscript banc/metrics/banc-extract-synapse-lookups.R
echo "  Done: $(date)"

# Step 4: Calculate neuropil inclusion (uses lookup for known synapses)
echo "=== Step 4: Calculating neuropil inclusion ==="
Rscript banc/metrics/banc-calculate-neuropil-inclusion.R
echo "  Done: $(date)"

# Step 5: Calculate NT predictions
echo "=== Step 5: Calculating NT predictions ==="
Rscript banc/metrics/banc-calculate-ntpred.R
echo "  Done: $(date)"

# Step 6: Calculate L2 metrics
echo "=== Step 6: Calculating L2 metrics ==="
Rscript banc/metrics/banc-calculate-l2-metrics.R
echo "  Done: $(date)"

# Step 7: Calculate root positions and regions
echo "=== Step 7: Calculating root positions ==="
Rscript banc/metrics/banc-calculate-root-positions.R
echo "  Done: $(date)"

echo "=== Step 7b: Calculating regions ==="
Rscript banc/metrics/banc-calculate-regions.R
echo "  Done: $(date)"

# Step 8: Calculate volumes
echo "=== Step 8: Calculating volumes ==="
Rscript banc/metrics/banc-calculate-volumes.R
echo "  Done: $(date)"

# Step 9: Calculate synapse split tables
echo "=== Step 9: Calculating synapse splits ==="
Rscript banc/metrics/banc-calculate-split.R
echo "  Done: $(date)"

# Step 10: Compile versioned data
echo "=== Step 10: Compiling versioned data (banc_850_*.feather/parquet) ==="
Rscript banc/share/banc-data.R
echo "  Done: $(date)"

# Step 11: Compile NBLAST feathers (with root_850)
echo "=== Step 11: Compiling NBLAST feathers ==="
Rscript banc/nblast/banc-nblast-compile.R
echo "  Done: $(date)"

# Step 12: Push NT predictions to SeaTable
echo "=== Step 12: Pushing NT predictions to SeaTable ==="
Rscript banc/update/banc-update-ntpred.R
echo "  Done: $(date)"

# Step 13: Push NBLAST feathers to GCS
echo "=== Step 13: Pushing NBLAST feathers to GCS ==="
Rscript banc/share/banc-nblast-share-gcs.R
echo "  Done: $(date)"

# Step 14: Push v850 versioned data to GCS
echo "=== Step 14: Pushing v850 data to GCS ==="
BANC_CONNECTIVITY="/n/data1/hms/neurobio/wilson/banc/connectivity"
GCS_BASE="gs://lee-lab_brain-and-nerve-cord-fly-connectome"

# Meta feather to compiled_data/<dataset>_<ver>/
gsutil cp "${BANC_CONNECTIVITY}/banc_850_meta.feather" \
  "${GCS_BASE}/compiled_data/banc_850/banc_850_meta.feather"

# Versioned connectivity files (legacy `brain_and_nerve_cord/connectivity/`
# subdir was dropped at the new bucket; goes alongside the meta).
for f in banc_850_meta.feather banc_850_edgelist_simple.feather \
         banc_850_edgelist_split.feather banc_850_synapses_enriched.parquet; do
  src="${BANC_CONNECTIVITY}/${f}"
  if [ -f "$src" ]; then
    gsutil cp "$src" "${GCS_BASE}/compiled_data/banc_850/${f}"
    echo "  Pushed: ${f}"
  else
    echo "  WARNING: ${src} not found, skipping"
  fi
done
echo "  Done: $(date)"

# Step 15: Clean up stale files
echo "=== Step 15: Cleaning up stale files ==="
Rscript banc/update/banc-delete.R
echo "  Done: $(date)"

end=$(date +%s)
echo "=========================================="
echo "  BANC v850 rebuild complete: $(date)"
echo "  Total time: $(( (end - start) / 60 )) minutes"
echo "=========================================="
