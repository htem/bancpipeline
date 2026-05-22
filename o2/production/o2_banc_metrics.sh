#!/bin/bash
#SBATCH -c 16                                      # cores
#SBATCH -t 6-23:00                                 # Runtime D-HH:MM (just under 1 week)
#SBATCH -p long                                    # 30-day cap; we fit comfortably
#SBATCH --mem-per-cpu=8G                           # 128G total
#SBATCH -o /home/ab714/bancpipeline/jobs/banc_metrics_%j.out
#SBATCH -e /home/ab714/bancpipeline/jobs/banc_metrics_%j.err

# Densify per-neuron metrics for the v888 proofread / roughly_proofread set.
# Each per-metric calc script is idempotent — it only processes neurons not
# already in its feather — so looping picks up CAVE transients, new IDs added
# during the run, and any neuron whose first-pass calc failed. Each pass ends
# with a SeaTable push and a "live" GCS upload of the combined feather; the
# final pass after the loop runs banc-data.R for the full versioned data push.

set -uo pipefail

cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0  # suppress core dumps from R shutdown crashes

mkdir -p jobs

start=$(date +%s)
# Stop kicking off a new pass after this point so the final banc-data.R has
# ~6h to finish before SLURM walltime kills us.
deadline=$((start + 6*24*3600 + 12*3600))

run_step () {
  local label="$1"; shift
  echo "[pass ${n_pass}] ${label} :: $(date)"
  if ! "$@"; then
    echo "[pass ${n_pass}] WARN: ${label} returned non-zero — continuing"
  fi
}

n_pass=0
while (( $(date +%s) < deadline )); do
  n_pass=$((n_pass+1))
  echo
  echo "==================== PASS ${n_pass} @ $(date) ===================="

  run_step "update root_ids"        Rscript banc/update/banc-ids.R
  run_step "L2 skeletons"           Rscript banc/metrics/banc-l2.R
  run_step "root positions"         Rscript banc/metrics/banc-calculate-root-positions.R
  run_step "regions"                Rscript banc/metrics/banc-calculate-regions.R
  run_step "synapses"               Rscript banc/metrics/banc-calculate-synapses.R
  run_step "L2 metrics"             Rscript banc/metrics/banc-calculate-l2-metrics.R
  run_step "volumes"                Rscript banc/metrics/banc-calculate-volumes.R
  run_step "push metrics SeaTable"  Rscript banc/update/banc-update-metrics.R

  # Live GCS upload of the combined feather so downstream consumers see
  # progress between full data pushes.
  COMBINED=/n/data1/hms/neurobio/wilson/banc/banc_metrics.feather
  if [ -s "$COMBINED" ]; then
    GCS_DST="gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/banc_888/banc_metrics_live.feather"
    echo "[pass ${n_pass}] live GCS push -> ${GCS_DST}"
    if ! gsutil -q cp "$COMBINED" "$GCS_DST"; then
      echo "[pass ${n_pass}] WARN: live GCS push failed"
    fi
  fi

  # CAVE/SeaTable cooldown before the next pass.
  echo "[pass ${n_pass}] cooldown 30 min before next pass"
  sleep 1800
done

echo
echo "==================== FINAL: banc-data.R @ $(date) ===================="
# Section 2 emits banc_888_metrics.feather from the now-dense combined feather;
# Section 9 rsyncs the whole banc_888/ tree to GCS.
Rscript banc/share/banc-data.R --source v3 || echo "WARN: banc-data.R returned non-zero"

end=$(date +%s)
echo
echo "completed ${n_pass} pass(es) in $((end-start))s"
