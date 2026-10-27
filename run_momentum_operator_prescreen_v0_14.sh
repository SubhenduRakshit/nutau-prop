#!/usr/bin/env bash
set -euo pipefail
PYTHON="${PYTHON:-python3}"
NPROC="${NUTAU_NPROC:-6}"
mkdir -p results/v0.14
OUT="results/v0.14/run_momentum_operator_prescreen_v0_14_summary.txt"
exec > >(tee "$OUT") 2>&1

echo "==== v0.14 momentum-operator prescreen started $(date '+%Y-%m-%d %H:%M:%S %Z') ===="
"$PYTHON" scripts/run_momentum_operator_prescreen_v0_14.py --workers "$NPROC"
echo "==== v0.14 momentum-operator prescreen completed $(date '+%Y-%m-%d %H:%M:%S %Z') ===="
echo "Summary: $OUT"
