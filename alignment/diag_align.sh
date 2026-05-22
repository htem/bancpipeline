#!/bin/bash
# Quick diagnostic for the iteration-degradation issue in alignment.
# Runs banc-alignment-run.py with cosine metric, 10 iters, 50% holdout,
# varying only the alpha schedule. Outputs go to data/optic_lobe/ with the
# diag_<case> suffix for easy comparison.
#
# Usage: sbatch alignment/diag_align.sh {baseline|gentle|reverse|pure_nblast}
#SBATCH -c 8
#SBATCH -p short
#SBATCH -t 0-01:00
#SBATCH --mem=64G
set -euo pipefail
cd /home/ab714/bancpipeline
source o2/o2_env.sh
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
ALIGN_PY=/home/ab714/.local/share/r-miniconda/envs/r-reticulate/bin/python
OPTIC_DIR=data/optic_lobe
FORBIDDEN=$OPTIC_DIR/forbidden-matches.csv

CASE="${1:?test case required}"

COMMON=(
  --side both
  --data-dir "$OPTIC_DIR" --file-prefix optic
  --metric cosine
  --hop2-weight 1.0
  --nblast-threshold 0.15
  --tau-start 4.0 --tau-end 0.5
  --nt-weight 0 --ind-weight 0.5
  --max-iter 10
  --stratified-holdout 0.5
  --manual-labels
  --forbidden-matches "$FORBIDDEN"
)
# Common args WITHOUT --stratified-holdout / --manual-labels for "designed" cases
# that use the script's intended Tier 1 = full holdout setup.
COMMON_DESIGNED=(
  --side both
  --data-dir "$OPTIC_DIR" --file-prefix optic
  --metric cosine
  --hop2-weight 1.0
  --nblast-threshold 0.15
  --tau-start 4.0 --tau-end 0.5
  --nt-weight 0 --ind-weight 0.5
  --forbidden-matches "$FORBIDDEN"
)
# Keep --bilateral as a default for most cases; cases override by NOT including it.
BILAT=(--bilateral)

case "$CASE" in
  baseline)    ARGS=(--alpha-start 0.05  --alpha-end 0.95 --output-suffix diag_baseline) ;;
  gentle)      ARGS=(--alpha-start 0.0  --alpha-end 0.5  --output-suffix diag_gentle) ;;
  reverse)     ARGS=(--alpha-start 0.95 --alpha-end 0.5  --output-suffix diag_reverse) ;;
  pure_nblast) ARGS=(--alpha-start 0.0  --alpha-end 0.0  --output-suffix diag_pure_nblast) ;;
  # Round 2 — investigating other potential causes (cosine, fast):
  no_stage1)   ARGS=(--alpha-start 0.05  --alpha-end 0.95 --skip-stage1
                     --output-suffix diag_no_stage1) ;;
  no_confgate) ARGS=(--alpha-start 0.05  --alpha-end 0.95 --no-conf-gate
                     --output-suffix diag_no_confgate) ;;
  low_tau)     ARGS=(--alpha-start 0.05  --alpha-end 0.95 --tau-start 2.0
                     --output-suffix diag_low_tau) ;;
  flat_alpha)  ARGS=(--alpha-start 0.05  --alpha-end 0.5
                     --output-suffix diag_flat_alpha) ;;
  # Momentum runs — classic softmax momentum on top of baseline schedule:
  mom_03)      ARGS=(--alpha-start 0.05  --alpha-end 0.95 --momentum 0.3
                     --output-suffix diag_mom_03) ;;
  mom_05)      ARGS=(--alpha-start 0.05  --alpha-end 0.95 --momentum 0.5
                     --output-suffix diag_mom_05) ;;
  mom_07)      ARGS=(--alpha-start 0.05  --alpha-end 0.95 --momentum 0.7
                     --output-suffix diag_mom_07) ;;
  # Combined: morphology-first ramp + momentum + skip stage1 (most-different-from-current)
  combined)    ARGS=(--alpha-start 0.2  --alpha-end 0.7 --momentum 0.5 --skip-stage1
                     --output-suffix diag_combined) ;;
  # Round 3 — strip individual "new" features. Each removes ONE thing from baseline
  # (α 0.5→0.95, conf_gate, stage1, ind_weight, bilateral, soma_rule, capacity).
  no_ind)      ARGS=(--alpha-start 0.05  --alpha-end 0.95 --ind-weight 0
                     --output-suffix diag_no_ind) ;;
  no_bilat)    ARGS=(--alpha-start 0.05  --alpha-end 0.95
                     --output-suffix diag_no_bilat) ;;
  no_soma)     ARGS=(--alpha-start 0.05  --alpha-end 0.95 --no-soma-rule
                     --output-suffix diag_no_soma) ;;
  big_cap)     ARGS=(--alpha-start 0.05  --alpha-end 0.95 --capacity-scale 100
                     --output-suffix diag_big_cap) ;;
  # Stripped: remove everything that didn't exist in v1 except α blending
  v1_like)     ARGS=(--alpha-start 0.05  --alpha-end 0.95
                     --skip-stage1 --no-conf-gate --ind-weight 0 --no-soma-rule
                     --capacity-scale 100 --tau-start 2.0
                     --output-suffix diag_v1_like) ;;
  # Round 4: remove --manual-labels too (the only "always on" common flag).
  no_manual)   ARGS=(--alpha-start 0.05  --alpha-end 0.95
                     --output-suffix diag_no_manual) ;;
  v1_no_manual) ARGS=(--alpha-start 0.05 --alpha-end 0.95
                      --skip-stage1 --no-conf-gate --ind-weight 0 --no-soma-rule
                      --capacity-scale 100 --tau-start 2.0
                      --output-suffix diag_v1_no_manual) ;;
  no_nbthresh)  ARGS=(--alpha-start 0.05 --alpha-end 0.95 --nblast-threshold 0.0
                      --output-suffix diag_no_nbthresh) ;;
  # Round 5: isolate algorithm's true Mi1 quality with no anchoring + α=0
  pure_nb_no_manual) ARGS=(--alpha-start 0.0 --alpha-end 0.0
                           --output-suffix diag_pure_nb_no_manual) ;;
  # And the same with v1-like simplifications (also no manual_labels)
  pure_nb_v1_like)   ARGS=(--alpha-start 0.0 --alpha-end 0.0
                           --skip-stage1 --no-conf-gate --ind-weight 0 --no-soma-rule
                           --capacity-scale 100 --tau-start 2.0
                           --output-suffix diag_pure_nb_v1_like) ;;
  # Round 6: USE THE DESIGNED EVAL — no stratified-holdout, no manual-labels.
  # Tier 1 = full holdout (40k typed gold), Tier 2 = soft NBLAST anchors.
  # This is how v1 was actually tested.
  designed)         ARGS=(--alpha-start 0.05 --alpha-end 0.95 --max-iter 30
                          --output-suffix diag_designed) ;;
  designed_v1like)  ARGS=(--alpha-start 0.05 --alpha-end 0.95 --max-iter 30
                          --skip-stage1 --no-conf-gate --ind-weight 0 --no-soma-rule
                          --capacity-scale 100 --tau-start 2.0
                          --output-suffix diag_designed_v1like) ;;
  *) echo "unknown case $CASE" >&2; exit 2 ;;
esac

echo "=== diag $CASE ==="
# Cases that need --bilateral OFF: omit BILAT.
# Cases that need --manual-labels OFF: rebuild COMMON without it.
NO_MANUAL_CASES="no_manual v1_no_manual pure_nb_no_manual pure_nb_v1_like"
NO_BILAT_CASES="no_bilat v1_like v1_no_manual pure_nb_v1_like designed_v1like"

CMD_COMMON=("${COMMON[@]}")
if [ "$CASE" = "designed" ] || [ "$CASE" = "designed_v1like" ]; then
  CMD_COMMON=("${COMMON_DESIGNED[@]}")
elif echo " $NO_MANUAL_CASES " | grep -q " $CASE "; then
  # Strip --manual-labels from CMD_COMMON
  CMD_COMMON=()
  for tok in "${COMMON[@]}"; do
    [ "$tok" = "--manual-labels" ] && continue
    CMD_COMMON+=("$tok")
  done
fi
if echo " $NO_BILAT_CASES " | grep -q " $CASE "; then
  "$ALIGN_PY" alignment/banc-alignment-run.py "${CMD_COMMON[@]}" "${ARGS[@]}"
else
  "$ALIGN_PY" alignment/banc-alignment-run.py "${CMD_COMMON[@]}" "${BILAT[@]}" "${ARGS[@]}"
fi
