#!/bin/bash
#SBATCH -c 8
#SBATCH -t 0-04:00
#SBATCH -p short
#SBATCH --mem-per-cpu=12G
#SBATCH -o /home/ab714/bancpipeline/jobs/banc_wb_prep_v888v2_%j.out
#SBATCH -e /home/ab714/bancpipeline/jobs/banc_wb_prep_v888v2_%j.err
###############################################################################
# Re-prep WB alignment inputs for v888 + v2 — STANDALONE (no align phase).
#
# Pulls fresh meta from SeaTable (banctable_query() — uncached) and rebuilds
# the alignment cache in data/whole_brain_alignment_v888v2/ using:
#   - latest banc_888 + v2 edgelist on disk (refreshed by 38847335 data push)
#   - 2026-05-04 NBLAST compile feathers under banc.meta.save.path
#
# Outputs (overwrites in place):
#   banc_brain_both_meta.csv, fafb_brain_both_meta.csv
#   banc_brain_both_edgelist.feather, fafb_brain_both_edgelist.feather
#   banc_brain_both_seeds.csv
#   banc_fafb_brain_both_nblast.csv
#   fafb_type_capacity.csv
#
# After this finishes:
#   1. Re-run scripts/seeds_to_production.sh to regenerate seeds_production.csv
#      from the freshly-written seeds.csv.
#   2. Wait for ensemble 38762177 to finish, pick winner, run production.
###############################################################################
set -uo pipefail
cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh

WB_DIR=data/whole_brain_alignment_v888v2

start=$(date +%s)
echo "=== prep v888 + v2 (fresh SeaTable meta + 2026-05-04 NBLAST feathers) ==="
BANC_WB_OUTPUT_DIR="$WB_DIR" \
BANC_SYN_SOURCE=v2 \
  Rscript alignment/banc-alignment-prep.R --region whole-brain both --source local || \
  echo "  prep R exited non-zero — verifying outputs..."

# Job 38870493 hit a `free(): invalid pointer` on R cleanup AFTER writing
# every output — same reticulate/arrow destructor race that bites
# nblast-compile, banc-data, banc-ids. Verify outputs by mtime instead
# of trusting the exit code.
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

echo "=== Done ==="
