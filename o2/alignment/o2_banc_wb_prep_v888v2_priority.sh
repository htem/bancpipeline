#!/bin/bash
#SBATCH -c 8
#SBATCH -t 0-04:00
#SBATCH -p priority
#SBATCH --mem-per-cpu=12G
#SBATCH -o /home/ab714/bancpipeline/jobs/banc_wb_prep_v888v2_priority_%j.out
#SBATCH -e /home/ab714/bancpipeline/jobs/banc_wb_prep_v888v2_priority_%j.err
###############################################################################
# Priority-partition copy of o2_banc_wb_prep_v888v2.sh. Rebuilds
# data/whole_brain_alignment_v888v2/banc_brain_both_* using freshly-pushed
# banc_888_meta.feather + banc_888_edgelist_simple_v2.feather + latest NBLAST.
#
# Decoupled from the production alignment so that align + NTAC can each
# depend on this and run in parallel on the next priority slot.
###############################################################################
set -uo pipefail
cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh

WB_DIR=data/whole_brain_alignment_v888v2

start=$(date +%s)
echo "=== prep v888 + v2 (priority) ==="
BANC_WB_OUTPUT_DIR="$WB_DIR" \
BANC_SYN_SOURCE=v2 \
  Rscript alignment/banc-alignment-prep.R --region whole-brain both --source local || \
  echo "  prep R exited non-zero — verifying outputs..."

missing=0
for f in banc_brain_both_meta.csv fafb_brain_both_meta.csv \
         banc_brain_both_edgelist.feather fafb_brain_both_edgelist.feather \
         banc_brain_both_seeds.csv banc_fafb_brain_both_nblast.csv \
         fafb_type_capacity.csv; do
  path="$WB_DIR/$f"
  if [[ ! -s "$path" ]]; then
    echo "  MISSING: $path" >&2
    missing=$((missing+1))
  else
    mtime=$(stat -c %Y "$path")
    if (( mtime < start )); then
      echo "  STALE (not refreshed): $path" >&2
      missing=$((missing+1))
    fi
  fi
done
if (( missing > 0 )); then
  echo "FATAL: $missing prep output(s) missing or stale" >&2
  exit 2
fi
echo "  all prep outputs verified fresh"

# Regenerate seeds_production.csv (all is_holdout=FALSE) for the align job.
scripts/seeds_to_production.sh "$WB_DIR"
echo "  seeds_production.csv regenerated"

echo "=== Done ($(( ($(date +%s) - start)/60 )) min) ==="
