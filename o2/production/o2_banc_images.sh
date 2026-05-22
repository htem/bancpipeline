#!/bin/bash
#SBATCH -c 1                              # Request cores
#SBATCH -t 1-12:00                       # Runtime in D-HH:MM format (actual ~27h)
#SBATCH -p medium                         # Partition to run in
#SBATCH --mem=48G                         # Total memory (bumped from 12G after recurring OOMs 2026-05-06 → 05-10)
#SBATCH -o /home/ab714/bancpipeline/jobs/banc_images_%j.out         # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e /home/ab714/bancpipeline/jobs/banc_images_%j.err         # File to which STDERR will be written, including job ID (%j)

# Why 48G: 4+ OOM events 2026-05-06 → 05-10. Image-build steps load NBLAST
# tables + neuron meshes; peak working set exceeded 12G. sort-folders.R is
# cheap but inherits whatever R sessions leaked before it. Trace which step
# OOMs via per-stage timestamps below.

echo "RUNNING BANC UPDATE"

start=`date +%s`

cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0  # suppress core dumps from R shutdown crashes

run_stage() {
  local label=$1; shift
  local stage_start=$(date +%s)
  echo "==================================================="
  echo "[$(date)] STAGE: $label"
  "$@" || echo "  STAGE $label exited non-zero (continuing)"
  local stage_end=$(date +%s)
  echo "[$(date)] STAGE $label DONE in $((stage_end-stage_start))s"
}

echo "building banc match images"
run_stage manc-images       Rscript banc/nblast/banc-manc-nblast-images.R
run_stage fafb-images       Rscript banc/nblast/banc-fafb-nblast-images.R
run_stage hemibrain-images  Rscript banc/nblast/banc-hemibrain-nblast-images.R
run_stage fanc-images       Rscript banc/nblast/banc-fanc-nblast-images.R
run_stage lr-images         Rscript banc/nblast/banc-lr-nblast-images.R

echo "sorting results"
run_stage sort-folders      Rscript banc/nblast/banc-sort-folders.R

# Report
end=`date +%s`
runtime=$((end-start))
echo "script completed in: "
echo $runtime

