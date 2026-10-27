#!/usr/bin/env bash
set -euo pipefail
PYTHON="${PYTHON:-python3}"
NPROC="${NUTAU_NPROC:-8}"
mkdir -p results data/icecube_2015
touch PROGRESS.log
exec > >(tee -a PROGRESS.log) 2>&1

echo
echo "---- v0.4 detector stage started $(date '+%Y-%m-%d %H:%M:%S %Z') ----"
echo "[D1/3] combined nu+nubar BSM transport/cross-section scan"
$PYTHON scripts/run_bsm_scan_combined.py --workers "$NPROC"
echo "[D2/3] official IceCube effective-area fold"
$PYTHON scripts/run_icecube_event_weighting.py
echo "[D3/3] detector-stage outputs"
echo "  results/event_weighted_counts.csv"
echo "  results/event_weighted_summary.csv"
echo "  results/event_weighted_summary.json"
echo "  results/event_weighted_best_model.png"
echo "---- v0.4 detector stage completed $(date '+%Y-%m-%d %H:%M:%S %Z') ----"
