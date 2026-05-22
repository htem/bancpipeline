#!/usr/bin/env bash
###############################################################################
# Scheduled BANC alignment runs
#
# Re-pulls fresh data from GCS / SeaTable (meta, connectivity edgelists, NBLAST
# scores, forbidden-matches), then runs three alignments sequentially:
#
#   1. Optic both, ensemble metric, ALL GT held out (no manual labels)
#   2. Optic both, ensemble metric, ALL GT as anchors (--manual-labels)
#   3. Whole brain both, cosine metric, ALL GT as anchors (--manual-labels)
#
# All runs use:
#   --bilateral
#   --forbidden-matches <path>   (reviewed false positives, refreshed at start)
#   soma rule ON                 (default; vetoes nucleated -> sensory/L1-L5/R7/R8)
#
# Resource gate: bails if 5-min load is high or RAM is tight, retrying every
# 30 minutes up to MAX_RETRIES times before giving up.
#
# Usage:
#   bash alignment/run-scheduled-alignments.sh
#
# Logs land in data/scheduled_runs/.
###############################################################################
set -uo pipefail

REPO=/Users/abates/projects/flyconnectome/bancpipeline
PYTHON=/opt/miniconda3/bin/python3
RSCRIPT=/usr/local/bin/Rscript

cd "$REPO" || { echo "FATAL: cannot cd to $REPO" >&2; exit 1; }

LOG_DIR="$REPO/data/scheduled_runs"
mkdir -p "$LOG_DIR"
STAMP=$(date +%Y%m%d_%H%M%S)
MAIN_LOG="$LOG_DIR/run_${STAMP}.log"

log() {
  printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" | tee -a "$MAIN_LOG"
}

###############################################################################
# Resource gate
#
# Machine: 10 cores, 64 GB RAM (M-series).
#   - LOAD_THRESHOLD 4.0  → 5-min loadavg must be < 4.0 (>=6 cores idle)
#   - RAM_THRESHOLD_GB 25 → free+inactive RAM must exceed 25 GB
#                            (whole-brain cosine bilateral peaks ~25-30 GB)
# Retry up to MAX_RETRIES times at RETRY_INTERVAL seconds before bailing.
###############################################################################
LOAD_THRESHOLD=${LOAD_THRESHOLD:-8.0}
RAM_THRESHOLD_GB=${RAM_THRESHOLD_GB:-25}
MAX_RETRIES=${MAX_RETRIES:-4}
RETRY_INTERVAL=${RETRY_INTERVAL:-1800}   # 30 min

check_resources() {
  local loadavg5 free_pages inactive_pages free_gb pagesize
  loadavg5=$(uptime | sed -E 's/.*load averages?: //' | awk '{print $2}')
  pagesize=$(sysctl -n hw.pagesize)
  free_pages=$(vm_stat | awk '/Pages free/ {gsub(/\./,"",$3); print $3}')
  inactive_pages=$(vm_stat | awk '/Pages inactive/ {gsub(/\./,"",$3); print $3}')
  free_gb=$(awk -v fp="$free_pages" -v ip="$inactive_pages" -v ps="$pagesize" \
    'BEGIN { printf "%.1f", (fp+ip)*ps/1024/1024/1024 }')
  log "  Resources: loadavg5m=${loadavg5} (need <${LOAD_THRESHOLD}), free+inactive RAM=${free_gb}GB (need >${RAM_THRESHOLD_GB})"
  awk -v l="$loadavg5" -v lt="$LOAD_THRESHOLD" -v g="$free_gb" -v gt="$RAM_THRESHOLD_GB" \
    'BEGIN { exit !(l < lt && g > gt) }'
}

wait_for_resources() {
  local attempt=0
  while (( attempt < MAX_RETRIES )); do
    if check_resources; then
      log "Resources OK, proceeding."
      return 0
    fi
    attempt=$((attempt + 1))
    if (( attempt < MAX_RETRIES )); then
      log "Resources insufficient — retry ${attempt}/${MAX_RETRIES} in $((RETRY_INTERVAL/60)) min"
      sleep "$RETRY_INTERVAL"
    fi
  done
  log "FATAL: resources still insufficient after ${MAX_RETRIES} retries; aborting."
  return 1
}

###############################################################################
# Step runner: command → log file, continue past failures.
###############################################################################
run_step() {
  local label="$1"; shift
  local logfile="$LOG_DIR/${STAMP}_${label}.log"
  log "==> ${label}"
  log "    cmd: $*"
  log "    log: ${logfile}"
  local t0=$(date +%s)
  if "$@" > "$logfile" 2>&1; then
    local dt=$(( $(date +%s) - t0 ))
    log "    OK (${label}, ${dt}s)"
    return 0
  else
    local rc=$?
    local dt=$(( $(date +%s) - t0 ))
    log "    FAILED (${label}, rc=${rc}, ${dt}s) — see ${logfile}"
    return $rc
  fi
}

###############################################################################
# Force-fresh GCS pulls
#
# The R prep scripts cache GCS downloads in tempdir(), which is per-Rscript
# session — every Rscript invocation already starts fresh. NBLAST is only
# pulled from GCS if the O2 path is missing (which it is locally), so it gets
# re-downloaded too. Nothing to clear here, but log it explicitly.
###############################################################################
log "=== Scheduled BANC alignments ==="
log "  Repo: ${REPO}"
log "  Stamp: ${STAMP}"
log "  Main log: ${MAIN_LOG}"
log "  PATH: $PATH"

if ! wait_for_resources; then
  exit 1
fi

###############################################################################
# Phase 1 — fresh data prep
###############################################################################
log "--- Phase 1: re-pull data + refresh forbidden matches ---"

run_step prep_optic        "$RSCRIPT" alignment/banc-alignment-prep.R --region optic-lobe  both --source gcs || true
run_step prep_wb           "$RSCRIPT" alignment/banc-alignment-prep.R --region whole-brain both --source gcs || true
run_step refresh_forbidden "$RSCRIPT" alignment/banc-alignment-false-positives.R           || true

OPTIC_FORBIDDEN="alignment/presets/optic-lobe/forbidden-matches.csv"
WB_FORBIDDEN="alignment/presets/whole-brain/forbidden-matches.csv"

if [[ ! -f "$OPTIC_FORBIDDEN" ]]; then
  log "WARNING: ${OPTIC_FORBIDDEN} missing after refresh; optic runs will skip --forbidden-matches"
fi
if [[ ! -f "$WB_FORBIDDEN" ]]; then
  log "WARNING: ${WB_FORBIDDEN} missing after refresh; WB run will skip --forbidden-matches"
fi

###############################################################################
# Phase 2 — three alignment runs
###############################################################################
log "--- Phase 2: alignment runs ---"

# Run 1: optic both, ensemble, ALL GT held out
run_step align_optic_holdout \
  "$PYTHON" alignment/banc-alignment-run.py \
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

# Run 2: optic both, ensemble, ALL GT anchored
run_step align_optic_anchored \
  "$PYTHON" alignment/banc-alignment-run.py \
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

# Run 3: whole brain both, cosine, ALL GT anchored
run_step align_wb_anchored \
  "$PYTHON" alignment/banc-alignment-run.py \
    --side both \
    --bilateral \
    --data-dir data/whole_brain_alignment \
    --file-prefix brain \
    --metric cosine \
    --hop2-weight 1.0 \
    --nblast-threshold 0.15 \
    --tau-start 4.0 \
    --alpha-start 0.05 \
    --alpha-end 0.95 \
    --nt-weight 0 \
    --ind-weight 0.0 \
    --max-iter 80 \
    --manual-labels \
    --forbidden-matches "$WB_FORBIDDEN" \
    --output-suffix with_rules_anchored

log "=== Done ==="
log "  Outputs:"
log "    data/optic_lobe/banc_optic_both_alignment_with_rules_holdout.csv"
log "    data/optic_lobe/banc_optic_both_alignment_with_rules_anchored.csv"
log "    data/whole_brain_alignment/banc_brain_both_alignment_with_rules_anchored.csv"
log "  Per-step logs in ${LOG_DIR}/${STAMP}_*.log"
