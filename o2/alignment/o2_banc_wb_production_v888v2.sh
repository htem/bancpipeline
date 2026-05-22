#!/bin/bash
#SBATCH -c 16
#SBATCH -t 2-12:00
#SBATCH -p medium
#SBATCH --mem-per-cpu=15G
#SBATCH -o /home/ab714/bancpipeline/jobs/banc_wb_production_v888v2_%j.out
#SBATCH -e /home/ab714/bancpipeline/jobs/banc_wb_production_v888v2_%j.err
###############################################################################
# v888 whole-brain production alignment — all tier-1 anchored, no holdout.
#
# Best holdout params from sweep + ensemble runs (64.7% @ iter 57 in 38762177):
#   ensemble metric, ensemble_blend=0.3, hop2_weight=1.0, nblast_threshold=0.15,
#   tau_start=4.0, alpha 0.05→0.95, nt_weight=0, ind_weight=0.3.
# Production switch (vs eval): `--manual-labels all` anchors every neuron with
# a known cell_type, OR equivalently we swap seeds.csv for seeds_production.csv
# (all is_holdout flipped to FALSE).  This script does both — belt and braces.
#
# Chain (intended dependencies):
#   38607783 (proofread-redo NBLAST) → o2_banc_update.sh (fresh meta + nblast
#   compile) → this script. Submit with --dependency=afterok:<update_id>.
#
# Wall: ensemble v888v2 took 1d 21h with 16c×15G. Production has more anchors
# → fewer unknowns to iterate, so similar or faster. Bumped to 2d 12h slack.
###############################################################################
set -euo pipefail
ulimit -c 0
cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"

ALIGN_PY=/home/ab714/.local/share/r-miniconda/envs/r-reticulate/bin/python
WB_DIR=data/whole_brain_alignment_v888v2

#-----------------------------------------------------------------------------
# Stage 1: fresh prep (reads latest NBLAST feathers + SeaTable meta)
#-----------------------------------------------------------------------------
echo "=== Stage 1: re-prep v888 + v2 with latest NBLAST + meta ==="
prep_start=$(date +%s)
BANC_WB_OUTPUT_DIR="$WB_DIR" \
BANC_SYN_SOURCE=v2 \
  Rscript alignment/banc-alignment-prep.R --region whole-brain both --source local || \
  echo "  prep R exited non-zero — verifying outputs..."

# Verify outputs refreshed (reticulate/arrow destructor race may exit non-zero
# AFTER writes complete — same pattern as 38870493).
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
    if (( mtime < prep_start )); then
      echo "  STALE (not refreshed): $path" >&2
      missing=$((missing+1))
    fi
  fi
done
if (( missing > 0 )); then
  echo "FATAL: $missing prep output(s) missing or stale" >&2
  exit 2
fi
echo "  all prep outputs refreshed"

#-----------------------------------------------------------------------------
# Stage 2: convert seeds.csv → seeds_production.csv (all tier-1 anchored)
#-----------------------------------------------------------------------------
echo ""
echo "=== Stage 2: convert seeds.csv → seeds_production.csv ==="
scripts/seeds_to_production.sh "$WB_DIR"

# Stash original seeds.csv, swap in production seeds so the align script picks
# them up (banc-alignment-run.py reads banc_brain_both_seeds.csv by name).
cp -p "$WB_DIR/banc_brain_both_seeds.csv" \
      "$WB_DIR/banc_brain_both_seeds.precapture.csv"
cp -p "$WB_DIR/banc_brain_both_seeds_production.csv" \
      "$WB_DIR/banc_brain_both_seeds.csv"
echo "  seeds.csv now points at the all-tier-1-anchored variant"
trap 'cp -p "$WB_DIR/banc_brain_both_seeds.precapture.csv" "$WB_DIR/banc_brain_both_seeds.csv" || true' EXIT

#-----------------------------------------------------------------------------
# Stage 3: ensemble alignment, all tier-1 anchored
#-----------------------------------------------------------------------------
echo ""
echo "=== Stage 3: ensemble + ind_weight=0.3 + --manual-labels all (80 iters) ==="
"$ALIGN_PY" alignment/banc-alignment-run.py \
  --side both --bilateral \
  --data-dir "$WB_DIR" --file-prefix brain \
  --metric ensemble \
  --ensemble-blend 0.3 \
  --hop2-weight 1.0 \
  --nblast-threshold 0.15 \
  --tau-start 4.0 --alpha-start 0.05 --alpha-end 0.95 \
  --nt-weight 0 \
  --ind-weight 0.3 \
  --manual-labels all \
  --max-iter 80 \
  --checkpoint-every 5 \
  --forbidden-matches "$WB_DIR/forbidden-matches.csv" \
  --output-suffix production_v888v2_all_tier1

echo ""
echo "=== Done ==="
