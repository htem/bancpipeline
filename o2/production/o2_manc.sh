#!/bin/bash
#SBATCH -c 1                          # Request cores
#SBATCH -t 00-120:00                  # Runtime in D-HH:MM format
#SBATCH -p medium                     # Partition to run in
#SBATCH --mem-per-cpu=250G            # Memory per core
#SBATCH -o jobs/manc_%j.out           # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e jobs/manc_%j.err           # File to which STDERR will be written, including job ID (%j)

start=`date +%s`

cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh

echo "get MANC metadata"
Rscript manc/manc-meta.R

echo "do MANC split"
Rscript manc/manc-split.R

echo "make MANC data"
Rscript manc/manc-data.R

end=`date +%s`
runtime=$((end-start))
echo "script completed in: "
echo $runtime

