#!/bin/bash
###########################################################
### Full v890 data rebuild for BANC pipeline
###
### Cloned from o2_banc_v888_rebuild.sh with three structural changes:
###   - banc.version = "890" via banc-startup.R (single source of truth)
###   - v3 connectivity + completion read CSV.gz, not parquet (curated v3
###     parquet is no longer staged at GCS as of v890)
###   - Spectral clustering + betweenness moved BEFORE the GCS push so their
###     outputs ship in the same compiled_data/banc_890/ rsync batch
###
### Prerequisites at v890 GCS (gs://brain-and-nerve-cord_exports/brain_and_nerve_cord/v890/):
###   - synapses_v2_human_readable.csv.gz                (✓ verified 2026-05-12)
###   - synapses_v2_human_readable_id_size_prerootid_postrootid_prex_prey_prez.parquet (✓)
###   - synapses_v1_human_readable_id_size_prerootid_postrootid_prex_prey_prez_neuropil.parquet (✓)
###   - synapses_v3_human_readable.csv.gz                (✓ verified 2026-05-08)
###
### Prerequisite on SeaTable:
###   - banc_meta must have a `root_890` text column. Step 1b (banc-ids890.R)
###     aborts with a clear message if it's missing.
###
### Usage:
###   sbatch o2/o2_banc_v890_rebuild.sh
###
### Steps:
###   1. Update SeaTable IDs (banc-ids.R — runs full path including bc890 if root_890 exists)
###  1b. Backfill root_890 column in SeaTable (safety net; idempotent)
###   2. Download L2 skeletons for any new/updated neurons
###   3. Calculate v890 connectivity (v2 synapse parquet + _v2 edgelist, size>=5)
###  3b. Calculate v890 v3 edgelist (reads CSV.gz from GCS, CAVE-curated root_ids)
###  3c. Extract version-agnostic synapse neuropil lookup (reuse prior assignments)
###   4. Calculate neuropil inclusion for synapses (lookup-accelerated)
###  4b. Calculate completion metrics (v1/v2/v3 capture rates at thresholds 0/5/10)
###   5. Calculate NT predictions from v890 v2 synapse data
###   6. Calculate L2 metrics (cable length, node count)
###   7. Calculate root positions and regions
###   8. Calculate volumes
###   9. Calculate synapse split tables
###  10. Compile versioned data (meta, _v2 + _v3 enriched synapses, _v2 + _v3 edgelists)
###  11. Compile NBLAST feathers (with root_890 column)
###  12. Push NT predictions to SeaTable
###  13. Push NBLAST feathers to GCS
###  14a. Spectral clustering v3 + v2 (parallel)               ← MOVED UP from old Step 15
###  14b. Betweenness centrality v3 + v2 (parallel)            ← MOVED UP from old Step 16
###  14c. Push banc_890/ versioned dir to GCS (compiled_data)
###  14d. Push BANC-project/data/cns_network/ → compiled_data/banc_890/cns_network/
###  14e. Push BANC-project/data/betweeness/ → compiled_data/banc_890/betweeness/
###  15. Clean up stale files
###  16. Chain influence (sbatch array + aggregate, runs AFTER this job exits)
###
### NOT included (run separately):
###   - Alignment: SKIPPED for v890 — reusing v888v2 ensemble production results
###     (per 2026-05-13 plan; v888→v890 root drift is small for alignment purposes).
###########################################################

#SBATCH -c 10
#SBATCH -t 3-00:00
#SBATCH -p priority
#SBATCH --mem-per-cpu=25G
#SBATCH -o /home/ab714/bancpipeline/jobs/v890_rebuild_%j.out
#SBATCH -e /home/ab714/bancpipeline/jobs/v890_rebuild_%j.err

start=$(date +%s)
cd /home/ab714/bancpipeline
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0

run_step() {
  local step="$1"; shift
  echo "=== Step ${step}: $* ==="
  "$@" || echo "  WARNING: Step ${step} exited with code $?"
  echo "  Done: $(date)"
}

# Parallel variant: runs the step in the background, tagging output with the
# step label so interleaved lines remain readable. Used by parallel groups
# below; the caller is responsible for calling `wait` at the end of the group.
run_step_bg() {
  local step="$1"; shift
  local log="/tmp/v890_rebuild_step_${step}_$$_$(date +%s%N).log"
  (
    echo "=== Step ${step}: $* ==="
    "$@" || echo "  WARNING: Step ${step} exited with code $?"
    echo "  Done Step ${step}: $(date)"
  ) > "$log" 2>&1 &
  # record (pid, step, log) for tail after wait
  BG_PIDS+=("$!")
  BG_STEPS+=("$step")
  BG_LOGS+=("$log")
}

# After a parallel group: wait for all, then stream each step's log (in
# submission order) into the main job log and tidy up.
flush_parallel_group() {
  local group_name="$1"
  echo ""
  echo "--- waiting for parallel group: ${group_name} (${#BG_PIDS[@]} steps) ---"
  local had_error=0
  for pid in "${BG_PIDS[@]}"; do
    wait "$pid" || had_error=1
  done
  for i in "${!BG_LOGS[@]}"; do
    echo ""
    echo "----- [group ${group_name}] Step ${BG_STEPS[$i]} output -----"
    cat "${BG_LOGS[$i]}"
    rm -f "${BG_LOGS[$i]}"
  done
  echo ""
  echo "--- group ${group_name} done (exit-errors: ${had_error}) ---"
  BG_PIDS=(); BG_STEPS=(); BG_LOGS=()
}

# Parallel-group accumulator arrays.
BG_PIDS=(); BG_STEPS=(); BG_LOGS=()

echo "=========================================="
echo "  BANC v890 rebuild started: $(date)"
echo "  cores=${SLURM_CPUS_PER_TASK:-?}  mem_per_cpu=${SLURM_MEM_PER_CPU:-?}"
echo "=========================================="

# -----------------------------------------------------------------------------
# SeaTable mutations — strictly sequential (later steps read what these write)
# -----------------------------------------------------------------------------
run_step 1  Rscript banc/update/banc-ids.R
# Step 1b: safety-net bc890 backfill. Idempotent — runs even if Step 1
# partial-failed. Aborts cleanly if root_890 column is missing from SeaTable.
run_step 1b Rscript banc/update/banc-ids890.R

# -----------------------------------------------------------------------------
# GROUP A (parallel): L2 download + v2 connectivity + v3 connectivity + lookup.
# All 4 are IO-bound (GCS pulls) or moderate-CPU; read-only on each other's
# outputs. Peak combined memory ~120G; fits in the 250G allocation.
# -----------------------------------------------------------------------------
run_step_bg 2   Rscript banc/metrics/banc-l2.R
run_step_bg 3   Rscript banc/metrics/banc-calculate-connectivity.R --source v2
run_step_bg 3b  Rscript banc/metrics/banc-calculate-connectivity.R --source v3
run_step_bg 3c  Rscript banc/metrics/banc-extract-synapse-lookups.R
flush_parallel_group "A (producers)"

# -----------------------------------------------------------------------------
# Sequential — real dependency chain
# -----------------------------------------------------------------------------
run_step 4   Rscript banc/metrics/banc-calculate-neuropil-inclusion.R
run_step 4b  Rscript banc/metrics/banc-calculate-completion.R
run_step 5   Rscript banc/metrics/banc-calculate-ntpred.R --source v2
run_step 5b  Rscript banc/metrics/banc-calculate-ntpred.R --source v3
run_step 6   Rscript banc/metrics/banc-calculate-l2-metrics.R
run_step 7   Rscript banc/metrics/banc-calculate-root-positions.R
run_step 7b  Rscript banc/metrics/banc-calculate-regions.R
run_step 8   Rscript banc/metrics/banc-calculate-volumes.R
run_step 9   Rscript banc/metrics/banc-calculate-split.R
run_step 10  Rscript banc/share/banc-data.R --source v3
run_step 10b Rscript banc/share/banc-data.R --source v2
run_step 11  Rscript banc/nblast/banc-nblast-compile.R

# -----------------------------------------------------------------------------
# GROUP C (parallel): NT writeback to SeaTable + NBLAST share to GCS.
# Both lightweight, independent writers.
# -----------------------------------------------------------------------------
run_step_bg 12  Rscript banc/update/banc-update-ntpred.R --source v2
run_step_bg 13  Rscript banc/share/banc-nblast-share-gcs.R
flush_parallel_group "C (writebacks)"

# 12b runs AFTER 12 (not in parallel) — both push to SeaTable; serializing
# avoids a write race on the same banc_meta rows.
run_step 12b Rscript banc/update/banc-update-ntpred.R --source v3

# -----------------------------------------------------------------------------
# GROUP B1 (parallel): spectral clustering v3 + v2 concurrently.
# Moved BEFORE the GCS push (Step 14c) so the outputs ship in the same batch.
# -----------------------------------------------------------------------------
run_step_bg 14a   Rscript banc/clustering/banc-spectral-clustering.R --source v3
run_step_bg 14a-v2  Rscript banc/clustering/banc-spectral-clustering.R --source v2
flush_parallel_group "B1 (clustering)"

# -----------------------------------------------------------------------------
# GROUP B2 (parallel): betweenness centrality v3 + v2 concurrently.
# Moved BEFORE the GCS push (Step 14c).
# -----------------------------------------------------------------------------
run_step_bg 14b   Rscript banc/betweenness/banc-betweenness-run.R both --source v3
run_step_bg 14b-v2  Rscript banc/betweenness/banc-betweenness-run.R both --source v2
flush_parallel_group "B2 (betweenness)"

# Step 14c: Push v890 versioned data to GCS
echo "=== Step 14c: Pushing v890 data to GCS ==="
BANC_CONNECTIVITY="/n/data1/hms/neurobio/wilson/banc/connectivity"
BANC_VERSIONED="/n/data1/hms/neurobio/wilson/connectomes/banc/banc_890"
GCS_NEW="gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/banc_890"
BANC_PROJECT_DATA="${BANC_PROJECT_DATA_DIR:-/n/data1/hms/neurobio/wilson/banc/BANC-project/data}"

# Single source of truth: rsync the full versioned directory to the new bucket.
# banc-data.R writes here (Sections 1-5 etc.), so this picks up everything in
# one shot. Exclude *.sqlite — the influencer intermediate
# (banc_<ver>_proofread_influence.sqlite, ~690 MB) is local-only; nothing
# downstream reads it from GCS. v888 leaked one to GCS by accident.
gsutil -m rsync -r -x '.*\.sqlite$' "${BANC_VERSIONED}" "${GCS_NEW}"
echo "  Pushed banc_890/ → ${GCS_NEW} (excluding *.sqlite)"
echo "  Done: $(date)"

# Step 14d: spectral clustering CSVs (BANC-project producer writes here)
echo "=== Step 14d: Pushing spectral clustering outputs to GCS ==="
if [ -d "${BANC_PROJECT_DATA}/cns_network" ]; then
  gsutil -m rsync -r "${BANC_PROJECT_DATA}/cns_network" "${GCS_NEW}/cns_network"
  echo "  Pushed cns_network/ → ${GCS_NEW}/cns_network/"
else
  echo "  WARNING: ${BANC_PROJECT_DATA}/cns_network not found; skipping"
fi
echo "  Done: $(date)"

# Step 14e: betweenness outputs
echo "=== Step 14e: Pushing betweenness outputs to GCS ==="
if [ -d "${BANC_PROJECT_DATA}/betweeness" ]; then
  gsutil -m rsync -r "${BANC_PROJECT_DATA}/betweeness" "${GCS_NEW}/betweeness"
  echo "  Pushed betweeness/ → ${GCS_NEW}/betweeness/"
else
  echo "  WARNING: ${BANC_PROJECT_DATA}/betweeness not found; skipping"
fi
echo "  Done: $(date)"

# Step 15: Clean up stale files
run_step 15 Rscript banc/update/banc-delete.R

# -----------------------------------------------------------------------------
# Step 16: chain all-to-all influence after the rebuild body finishes.
# Submits the influence array immediately; aggregate+sync gated on afterok.
# Writes go to banc_890/influence/ (version-aware via banc-startup.R).
# -----------------------------------------------------------------------------
EDGELIST_V3="/n/data1/hms/neurobio/wilson/connectomes/banc/banc_${BANC_VERSION:-890}/banc_${BANC_VERSION:-890}_edgelist_simple_v3.feather"
if [ -f "$EDGELIST_V3" ]; then
  echo "=== Step 16: submitting chained all-to-all influence ==="
  INF_ARRAY_JID=$(sbatch --parsable o2_banc_influence_array.sh) && \
    echo "  Influence array submitted: ${INF_ARRAY_JID}"
  if [ -n "${INF_ARRAY_JID:-}" ]; then
    INF_AGG_JID=$(sbatch --parsable \
      --dependency=afterok:${INF_ARRAY_JID} \
      o2_banc_aggregate_influence.sh) && \
      echo "  Aggregate+sync submitted: ${INF_AGG_JID} (waits on ${INF_ARRAY_JID})"
  else
    echo "  WARNING: array submission failed; skipping aggregate chain" >&2
  fi
  echo "  Done: $(date)"
else
  echo "=== Step 16: SKIPPED — ${EDGELIST_V3} not found ==="
fi

end=$(date +%s)
echo "=========================================="
echo "  v890 rebuild complete in $((end-start))s"
echo "=========================================="
