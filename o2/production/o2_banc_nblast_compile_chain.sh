#!/bin/bash
#SBATCH -J banc_nblast_compile_chain
#SBATCH -c 4
#SBATCH -t 0-23:00
#SBATCH -p medium
#SBATCH --mem=128G
#SBATCH -o jobs/banc_nblast_compile_chain_%j.out
#SBATCH -e jobs/banc_nblast_compile_chain_%j.err

# Post-fleet finalisation: dedupe wrong matches → compile feathers → push to
# CAVE → push to GCS. Submit with `--dependency=afterany:<jobid>:<jobid>:...`
# so it queues until the upstream NBLAST jobs finish (afterany = on any exit
# status, so a single TIMEOUT doesn't block the chain).
#
# Prior standalone compile-only (40811652) TIMEOUT'd at 4h inside the root_id
# refresh — bumped to 23h on medium with 128G headroom.

set -uo pipefail
start=$(date +%s)
cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0

echo "=== Stage 1: wrong matches @ $(date) ==="
Rscript banc/nblast/banc-nblast-wrong-matches.R || echo "  WARN: wrong-matches exited non-zero"

echo "=== Stage 2: compile feathers @ $(date) ==="
Rscript banc/nblast/banc-nblast-compile.R || echo "  WARN: compile exited non-zero"

echo "=== Stage 3: sync to CAVE @ $(date) ==="
Rscript banc/nblast/banc-nblast-cave.R || echo "  WARN: cave sync exited non-zero"

echo "=== Stage 4: push to GCS @ $(date) ==="
Rscript banc/share/banc-nblast-share-gcs.R || echo "  WARN: gcs share exited non-zero"

end=$(date +%s)
echo "compile chain completed in $((end-start))s"
