#!/bin/bash
#SBATCH -c 16                                # 16 cores (OMP-parallel sparse ops in align.py)
#SBATCH -t 2-00:00                           # 2 days
#SBATCH -p medium                            # medium partition (5d max, fits comfortably)
#SBATCH --mem=128G                           # whole-brain cosine+bilateral peaks ~25-30G; headroom
#SBATCH -o /home/ab714/bancpipeline/jobs/banc_alignment_%j.out
#SBATCH -e /home/ab714/bancpipeline/jobs/banc_alignment_%j.err
###############################################################################
# BANC cross-dataset alignment — O2 driver
#
# Mirrors alignment/run-scheduled-alignments.sh but uses local O2 v850
# data (edgelists + NBLAST) instead of GCS. Runs prep for optic lobe and whole
# brain, refreshes the forbidden-matches CSV from SeaTable, then runs three
# alignment passes:
#
#   1. Optic both, ensemble metric, ALL GT held out       (validation)
#   2. Optic both, ensemble metric, ALL GT as anchors     (production)
#   3. Whole brain both, cosine metric, ALL GT as anchors (production, focus)
#
# Whole brain is the primary target going forward; the two optic runs remain
# for regression tracking. To skip optic runs entirely, set SKIP_OPTIC=1.
#
# Usage:
#   sbatch o2_banc_alignment.sh
#   SKIP_OPTIC=1 sbatch o2_banc_alignment.sh
###############################################################################
set -euo pipefail

cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0  # suppress core dumps from R shutdown crashes

# Keep OMP parallelism in sync with the core allocation
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"

SKIP_OPTIC="${SKIP_OPTIC:-0}"

echo "=== BANC alignment (O2) ==="
echo "  Job ID:        ${SLURM_JOB_ID:-n/a}"
echo "  Cores:         ${SLURM_CPUS_PER_TASK:-?}"
echo "  Skip optic:    ${SKIP_OPTIC}"
start=$(date +%s)

###############################################################################
# Phase 1 — data prep (all from local O2 v850 paths)
###############################################################################
echo
echo "--- Phase 1: data prep ---"

# Tolerate Rscript shutdown crashes (free(): invalid pointer in libc teardown
# after all output files are written) — gate on file existence instead.
require_files() {
  for f in "$@"; do
    if [[ ! -f "$f" ]]; then
      echo "FATAL: missing prep output $f" >&2
      exit 4
    fi
  done
}

if [[ "$SKIP_OPTIC" != "1" ]]; then
  echo "[prep_optic]"
  Rscript alignment/banc-alignment-prep.R --region optic-lobe both --source local || true
  require_files \
    data/optic_lobe/banc_optic_both_meta.csv \
    data/optic_lobe/fafb_optic_both_meta.csv \
    data/optic_lobe/banc_optic_both_edgelist.feather \
    data/optic_lobe/fafb_optic_both_edgelist.feather \
    data/optic_lobe/banc_optic_both_seeds.csv \
    data/optic_lobe/fafb_type_capacity.csv
fi

echo "[prep_wb]"
Rscript alignment/banc-alignment-prep.R --region whole-brain both --source local || true
require_files \
  data/whole_brain_alignment/banc_brain_both_meta.csv \
  data/whole_brain_alignment/fafb_brain_both_meta.csv \
  data/whole_brain_alignment/banc_brain_both_edgelist.feather \
  data/whole_brain_alignment/fafb_brain_both_edgelist.feather \
  data/whole_brain_alignment/banc_brain_both_seeds.csv \
  data/whole_brain_alignment/fafb_type_capacity.csv

echo "[refresh_forbidden]"
Rscript alignment/banc-alignment-false-positives.R || true
require_files \
  alignment/presets/optic-lobe/forbidden-matches.csv \
  alignment/presets/whole-brain/forbidden-matches.csv

OPTIC_FORBIDDEN=alignment/presets/optic-lobe/forbidden-matches.csv
WB_FORBIDDEN=alignment/presets/whole-brain/forbidden-matches.csv

###############################################################################
# Phase 2 — alignment runs
###############################################################################
echo
echo "--- Phase 2: alignment runs ---"

if [[ "$SKIP_OPTIC" != "1" ]]; then
  echo "[align_optic_holdout]"
  python alignment/banc-alignment-run.py \
    --side both \
    --bilateral \
    --metric ensemble \
    --hop2-weight 1.0 \
    --nblast-threshold 0.15 \
    --tau-start 4.0 \
    --alpha-start 0.05 \
    --alpha-end 0.95 \
    --ensemble-blend 0.3 \
    --nt-weight 0 \
    --ind-weight 0.5 \
    --max-iter 60 \
    --forbidden-matches "$OPTIC_FORBIDDEN" \
    --output-suffix with_rules_holdout

  echo "[align_optic_anchored]"
  python alignment/banc-alignment-run.py \
    --side both \
    --bilateral \
    --metric ensemble \
    --hop2-weight 1.0 \
    --nblast-threshold 0.15 \
    --tau-start 4.0 \
    --alpha-start 0.05 \
    --alpha-end 0.95 \
    --ensemble-blend 0.3 \
    --nt-weight 0 \
    --ind-weight 0.5 \
    --max-iter 60 \
    --manual-labels \
    --forbidden-matches "$OPTIC_FORBIDDEN" \
    --output-suffix with_rules_anchored
fi

echo "[align_wb_anchored]"
# ind_weight=0.0 for WB: OOMs at 127k×35k scale otherwise. Ensemble metric on.
python alignment/banc-alignment-run.py \
  --side both \
  --bilateral \
  --data-dir data/whole_brain_alignment \
  --file-prefix brain \
  --metric ensemble \
  --hop2-weight 1.0 \
  --nblast-threshold 0.15 \
  --tau-start 4.0 \
  --alpha-start 0.05 \
  --alpha-end 0.95 \
  --ensemble-blend 0.3 \
  --nt-weight 0 \
  --ind-weight 0.0 \
  --max-iter 80 \
  --manual-labels \
  --forbidden-matches "$WB_FORBIDDEN" \
  --output-suffix with_rules_anchored

end=$(date +%s)
runtime=$((end-start))
echo
echo "=== Done (${runtime}s total) ==="
echo "Outputs:"
if [[ "$SKIP_OPTIC" != "1" ]]; then
  echo "  data/optic_lobe/banc_optic_both_alignment_with_rules_holdout.csv"
  echo "  data/optic_lobe/banc_optic_both_alignment_with_rules_anchored.csv"
fi
echo "  data/whole_brain_alignment/banc_brain_both_alignment_with_rules_anchored.csv"
