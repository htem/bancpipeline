#!/bin/bash
###########################################################
### Full v888 data rebuild for BANC pipeline
###
### Produces all versioned data files and pushes to GCS.
### Prerequisites (must already exist on GCS at v888):
###   - synapses_v2_human_readable.csv.gz
###   - synapses_v3_human_readable_id_size_prerootid_postrootid_prex_prey_prez.parquet
###     (CAVE-ingested export — used by Step 3b for v3 edgelist root_ids,
###     and by Step 4b for v3 capture rates)
###
### As of 2026-04-20: locally-produced banc_<ver>_synapses_v3.parquet is no
### longer required as a prerequisite. The v3 root_ids in that file are
### corrupt (~30% spurious "0"s from svid_cache contamination); we now read
### root_ids straight from the GCS export. Local v3 batches still provide
### spatial neuropil annotations to the completion script.
###
### Usage:
###   sbatch o2/o2_banc_v888_rebuild.sh
###
### Steps:
###   1. Update SeaTable IDs (root_id, banc_match, banc_png_match, banc_nblast_match)
###      KNOWN ISSUE (2026-04-19): banc-ids.R crashes early at
###      banc_proofreading_notes(live=TRUE) with pt_position type mismatch.
###      Step 1b below patches the root_888 gap.
###  1b. Backfill root_888 column in SeaTable (runs even if Step 1 partial-failed)
###   2. Download L2 skeletons for any new/updated neurons
###   3. Calculate v888 connectivity (v2 synapse parquet + _v2 edgelist, size>=5)
###  3b. Calculate v888 v3 edgelist (reads GCS v3 export — CAVE-curated root_ids)
###  3c. Extract version-agnostic synapse neuropil lookup (reuse prior assignments)
###   4. Calculate neuropil inclusion for synapses (lookup-accelerated)
###  4b. Calculate completion metrics (v1/v2/v3 capture rates at thresholds 0/5/10)
###   5. Calculate NT predictions from v888 v2 synapse data
###   6. Calculate L2 metrics (cable length, node count)
###   7. Calculate root positions and regions
###   8. Calculate volumes
###   9. Calculate synapse split tables
###  10. Compile versioned data (meta, _v2 + _v3 enriched synapses, _v2 + _v3 edgelists)
###  11. Compile NBLAST feathers (with root_888 column)
###  12. Push NT predictions to SeaTable
###  13. Push NBLAST feathers to GCS (v888)
###  14. Push v888 versioned data to GCS (banc_888/)
###  15. Clean up stale files
###
### NOT included (run separately when ready):
###   - Influence: re-run banc-build-influence.R + banc-aggregate-influence.R
###   - Alignment: re-run banc-alignment-prep.R --region whole-brain + sweep with v888 root IDs
###   - Spectral clustering + betweenness: re-run with --source v2|v3
###########################################################

#SBATCH -c 10
#SBATCH -t 3-00:00
#SBATCH -p medium
#SBATCH --mem-per-cpu=25G
#SBATCH -o /home/ab714/bancpipeline/jobs/v888_rebuild_%j.out
#SBATCH -e /home/ab714/bancpipeline/jobs/v888_rebuild_%j.err

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
  local log="/tmp/v888_rebuild_step_${step}_$$_$(date +%s%N).log"
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
echo "  BANC v888 rebuild started: $(date)"
echo "  cores=${SLURM_CPUS_PER_TASK:-?}  mem_per_cpu=${SLURM_MEM_PER_CPU:-?}"
echo "=========================================="

# -----------------------------------------------------------------------------
# SeaTable mutations — strictly sequential (later steps read what these write)
# -----------------------------------------------------------------------------
run_step 1  Rscript banc/update/banc-ids.R
# Step 1b: safety-net bc888 backfill. Idempotent — runs even if Step 1
# partial-failed at the (now patched) banc_proofreading_notes call.
run_step 1b Rscript banc/update/banc-ids888.R

# -----------------------------------------------------------------------------
# GROUP A (parallel): L2 download + v2 connectivity + v3 connectivity + lookup.
# All 4 are IO-bound (GCS pulls) or moderate-CPU; read-only on each other's
# outputs. Peak combined memory ~120G; fits in the 250G allocation.
# Saves ~max(t3, t3b, t3c) + t2 vs serial; typically ~1h.
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
# Both lightweight, independent writers. Saves the length of one.
# -----------------------------------------------------------------------------
run_step_bg 12  Rscript banc/update/banc-update-ntpred.R --source v2
run_step_bg 13  Rscript banc/share/banc-nblast-share-gcs.R
flush_parallel_group "C (writebacks)"

# 12b runs AFTER 12 (not in parallel) — both push to SeaTable; serializing
# avoids a write race on the same banc_meta rows. Writes the v3-source
# per-neuron NT prediction to neurotransmitter_predicted_v3.
run_step 12b Rscript banc/update/banc-update-ntpred.R --source v3

# Step 14: Push v888 versioned data to GCS
echo "=== Step 14: Pushing v888 data to GCS ==="
BANC_CONNECTIVITY="/n/data1/hms/neurobio/wilson/banc/connectivity"
BANC_VERSIONED="/n/data1/hms/neurobio/wilson/connectomes/banc/banc_888"
GCS_NEW="gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/banc_888"

# Single source of truth: rsync the full versioned directory to the new bucket.
# banc-data.R writes here (Sections 1-5 etc.), so this picks up everything in one shot.
gsutil -m rsync -r "${BANC_VERSIONED}" "${GCS_NEW}"
echo "  Pushed banc_888/ → ${GCS_NEW}"
echo "  Done: $(date)"

# -----------------------------------------------------------------------------
# GROUP B1 (parallel): spectral clustering v3 + v2 concurrently. Each reads
# its own edgelist, independent memory. Saves one full clustering run.
# -----------------------------------------------------------------------------
run_step_bg 15   Rscript banc/clustering/banc-spectral-clustering.R --source v3
run_step_bg 15b  Rscript banc/clustering/banc-spectral-clustering.R --source v2
flush_parallel_group "B1 (clustering)"

# -----------------------------------------------------------------------------
# GROUP B2 (parallel): betweenness centrality v3 + v2 concurrently. Same
# independence as clustering. Saves one full betweenness run.
# -----------------------------------------------------------------------------
run_step_bg 16   Rscript banc/betweenness/banc-betweenness-run.R both --source v3
run_step_bg 16b  Rscript banc/betweenness/banc-betweenness-run.R both --source v2
flush_parallel_group "B2 (betweenness)"

# Step 17: Clean up stale files
run_step 17 Rscript banc/update/banc-delete.R

# -----------------------------------------------------------------------------
# Step 18: chain all-to-all influence after the rebuild body finishes.
# Submits the 12-shard array immediately; the aggregate+sync job is queued
# with --dependency=afterok:<array_jid> so it only fires if every shard
# returns exit 0. Writes go to banc_888/influence/ (version-aware via
# banc-startup.R). Skipped if the v3 edgelist isn't on disk.
# -----------------------------------------------------------------------------
EDGELIST_V3="/n/data1/hms/neurobio/wilson/connectomes/banc/banc_${BANC_VERSION:-888}/banc_${BANC_VERSION:-888}_edgelist_simple_v3.feather"
if [ -f "$EDGELIST_V3" ]; then
  echo "=== Step 18: submitting chained all-to-all influence ==="
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
  echo "=== Step 18: SKIPPED — ${EDGELIST_V3} not found ==="
fi

end=$(date +%s)
echo "=========================================="
echo "  BANC v888 rebuild complete: $(date)"
echo "  Total time: $(( (end - start) / 60 )) minutes"
echo "=========================================="
