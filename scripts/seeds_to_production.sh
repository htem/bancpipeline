#!/bin/bash
###############################################################################
# Convert WB alignment seeds.csv → seeds_production.csv (no holdout).
#
# Flips every is_holdout=TRUE → FALSE so all GT cells become anchors. Use
# this for the final 100% GT production alignment run after a re-prep.
#
# is_holdout is the last column; cell_type may contain commas inside quoted
# strings. Anchor on `,TRUE$` end-of-line to dodge that.
#
# Usage:
#   scripts/seeds_to_production.sh [WB_DIR]
# WB_DIR defaults to data/whole_brain_alignment_v888v2/.
###############################################################################
set -euo pipefail

WB_DIR="${1:-data/whole_brain_alignment_v888v2}"
SRC="$WB_DIR/banc_brain_both_seeds.csv"
DST="$WB_DIR/banc_brain_both_seeds_production.csv"

if [[ ! -f "$SRC" ]]; then
  echo "FATAL: $SRC not found" >&2
  exit 1
fi

n_true_in=$(grep -c ',TRUE$'  "$SRC")
n_false_in=$(grep -c ',FALSE$' "$SRC")

sed 's/,TRUE$/,FALSE/' "$SRC" > "$DST"

n_true_out=$(grep -c ',TRUE$'  "$DST" || true)
n_false_out=$(grep -c ',FALSE$' "$DST")

echo "Source : $SRC"
echo "  is_holdout=TRUE : $n_true_in"
echo "  is_holdout=FALSE: $n_false_in"
echo ""
echo "Production: $DST"
echo "  is_holdout=TRUE : $n_true_out  (should be 0)"
echo "  is_holdout=FALSE: $n_false_out"

if (( n_true_out != 0 )); then
  echo "FATAL: production still has TRUE rows" >&2
  exit 2
fi

n_anchors=$(awk -F, 'NR>1 && $NF=="FALSE"' "$DST" | wc -l)
echo ""
echo "Total anchors in production: $n_anchors"
