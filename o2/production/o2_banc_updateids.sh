#!/bin/bash
#SBATCH -c 1                             # Request cores
#SBATCH -t 0-01:00                       # Runtime in D-HH:MM format (was 10m, hit timeout 2026-04-28 — banc-ids.R needs ~13m per orchestrator history)
#SBATCH -p short                         # Partition to run in
#SBATCH --mem-per-cpu=8G                 # Memory per core (was 3G — bumping for safety)
#SBATCH -o /home/ab714/bancpipeline/jobs/banc_updateids_%j.out    # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e /home/ab714/bancpipeline/jobs/banc_updateids_%j.err    # File to which STDERR will be written, including job ID (%j)

start=`date +%s`

cd /home/ab714/bancpipeline/
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0  # suppress core dumps from R shutdown crashes

echo "update banc ids in seatable"
Rscript banc/update/banc-updateids.R

end=`date +%s`
runtime=$((end-start))
echo "script completed in: "
echo $runtime

