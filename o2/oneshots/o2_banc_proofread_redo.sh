#!/bin/bash
###########################################################
### Proofread-restricted NBLAST redo across all datasets
###
### Goal: lift FAFB / MaleCNS / MANC / FANC / hemibrain / LR
### / native NBLAST coverage of BANC v888 by re-running every
### `proofread == TRUE` neuron with redo=TRUE, then compile
### and push to GCS.
###
### Mechanism:
### - banc/nblast/banc-make-proofread-ids.R writes the
###   proofread root_ids to a flat file.
### - BANC_TEST_IDS_FILE points each NBLAST script at that
###   file so banc.test.ids restricts nblast.todo.
### - BANC_NBLAST_REDO=TRUE forces re-NBLAST on existing CSVs
###   (otherwise the per-neuron file.exists() check skips).
### - Each NBLAST stage is bounded by `timeout` so a slow
###   stage can't starve the rest. Adjust budgets if needed.
###
### Wall: 7 days on the long partition.
###########################################################
#SBATCH -c 20
#SBATCH -t 7-00:00
#SBATCH -p long
#SBATCH --mem-per-cpu=8G  # 160G total — 5G/cpu accumulated 42 oom_kill events in 38607783
#SBATCH -o jobs/banc_proofread_redo_%j.out
#SBATCH -e jobs/banc_proofread_redo_%j.err
#SBATCH -J banc_proofread_redo

set -uo pipefail
start_total=$(date +%s)

cd /home/ab714/bancpipeline/
source /home/ab714/bancpipeline/o2/o2_env.sh
ulimit -c 0

PROOFREAD_IDS=/n/data1/hms/neurobio/wilson/banc/meta/banc_proofread_ids.txt
export BANC_TEST_IDS_FILE="$PROOFREAD_IDS"
export BANC_NBLAST_REDO=TRUE
export BANC_PROOFREAD_IDS_FILE="$PROOFREAD_IDS"

run_stage () {
  local label=$1; shift
  local budget_secs=$1; shift
  local cmd="$*"
  local s=$(date +%s)
  echo "==================================================="
  echo "[$(date)] STAGE: $label  (budget ${budget_secs}s)"
  echo "  cmd: $cmd"
  timeout "${budget_secs}s" bash -c "$cmd" || true
  local elapsed=$(( $(date +%s) - s ))
  echo "[$(date)] STAGE $label DONE in ${elapsed}s"
}

# --- Stage 0: refresh proofread id list ---
run_stage "make-proofread-ids" 600 "Rscript banc/nblast/banc-make-proofread-ids.R"

# --- Stage 1: refresh L2 skeletons (banc.l2.save.path) ---
# banc-l2.R fetches missing L2 SWCs; proofread filter doesn't apply here, but
# any newly-proofread neurons' SWCs must exist before NBLAST runs.
run_stage "banc-l2"             21600 "Rscript banc/metrics/banc-l2.R"

# --- Stages 2-8: per-target NBLAST, redo=TRUE, restricted to proofread ids ---
# Budgets allocate roughly:  FAFB 1.5d, MaleCNS 1d, MANC 8h, FANC 6h,
# hemibrain 1d, LR mirror 1d, native 1d. Total ~6.5d, leaves ~12h for
# compile + share + slack.
run_stage "fafb-nblast"        129600 "Rscript banc/nblast/banc-fafb-nblast.R"        # 36h
run_stage "malecns-nblast"      86400 "Rscript banc/nblast/banc-malecns-nblast.R"     # 24h
run_stage "manc-nblast"         28800 "Rscript banc/nblast/banc-manc-nblast.R"        # 8h
run_stage "fanc-nblast"         21600 "Rscript banc/nblast/banc-fanc-nblast.R"        # 6h
run_stage "hemibrain-nblast"    86400 "Rscript banc/nblast/banc-hemibrain-nblast.R"   # 24h
run_stage "lr-nblast"           86400 "Rscript banc/nblast/banc-nblast-lr.R"          # 24h
run_stage "native-nblast"       86400 "Rscript banc/nblast/banc-nblast-native.R"      # 24h

# --- Stage 9: drop wrong matches before compile ---
run_stage "wrong-matches"        7200 "Rscript banc/nblast/banc-nblast-wrong-matches.R"

# --- Stage 10: compile every per-target NBLAST into feathers ---
# (compile doesn't need the proofread filter; clear the env var so no scope
#  bleed.)
unset BANC_TEST_IDS_FILE
unset BANC_NBLAST_REDO
run_stage "compile"             43200 "Rscript banc/nblast/banc-nblast-compile.R"     # 12h

# --- Stage 11: push compiled feathers + reviewed-match CSVs to GCS ---
run_stage "share-gcs"            7200 "Rscript banc/share/banc-nblast-share-gcs.R"   # 2h

end_total=$(date +%s)
echo "==================================================="
echo "[$(date)] WORKFLOW COMPLETE in $(( end_total - start_total ))s"
