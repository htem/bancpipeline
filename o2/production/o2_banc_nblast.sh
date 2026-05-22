#!/bin/bash
#SBATCH -c 10                           # Request cores
#SBATCH -t 2-00:00                      # Runtime in D-HH:MM format. Bumped to 2d
                                        # (was 18h, still TIMEOUT on 39983439 +
                                        # 40106803, 2026-05-14/15). Cron fires
                                        # daily; with 2d walltime overlapping
                                        # runs are acceptable while we clear
                                        # the post-v888 NBLAST backlog. Once
                                        # caught up, daily runs should finish
                                        # within 24h again.
#SBATCH -p medium                       # Partition to run in (short caps at 12h; this workload regularly exceeds it)
#SBATCH --mem-per-cpu=5G                # Memory per core (50G total)
#SBATCH -o jobs/banc_nblast_%j.out      # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e jobs/banc_nblast_%j.err      # File to which STDERR will be written, including job ID (%j)

start=`date +%s`

cd /home/ab714/bancpipeline/
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0  # suppress core dumps from R shutdown crashes

echo "NBLASTing BANC skeletons against FAFB skeletons"
Rscript banc/nblast/banc-fafb-nblast.R

echo "NBLASTing BANC skeletons against maleCNS skeletons"
Rscript banc/nblast/banc-malecns-nblast.R

echo "NBLASTing BANC skeletons against MANC skeletons"
Rscript banc/nblast/banc-manc-nblast.R

echo "NBLASTing BANC skeletons against FANC skeletons"
Rscript banc/nblast/banc-fanc-nblast.R

echo "NBLASTing BANC skeletons against hemibrain skeletons"
Rscript banc/nblast/banc-hemibrain-nblast.R

echo "NBLASTing BANC skeletons left-right"
Rscript banc/nblast/banc-nblast-lr.R

echo "Removing invalidated NBLAST matches"
Rscript banc/nblast/banc-nblast-wrong-matches.R

echo "Compiling NBLAST results into feather files"
Rscript banc/nblast/banc-nblast-compile.R

echo "Syncing NBLAST matches to CAVE tables"
Rscript banc/nblast/banc-nblast-cave.R

echo "Update saved results"
Rscript banc/share/banc-nblast-share-gcs.R

echo "saving banc l2 skeletons"
Rscript banc/metrics/banc-l2.R

echo "deleting outdated BANC data"
Rscript banc/update/banc-delete.R

end=`date +%s`
runtime=$((end-start))
echo "script completed in: "
echo $runtime
