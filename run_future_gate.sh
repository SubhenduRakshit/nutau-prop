#!/usr/bin/env bash
set -euo pipefail
PYTHON="${PYTHON:-python3}"
NPROC="${NUTAU_NPROC:-8}"
mkdir -p results/future_current results/future
touch PROGRESS.log
exec > >(tee -a PROGRESS.log) 2>&1

echo
echo "---- v0.5 future publishability gate started $(date '+%Y-%m-%d %H:%M:%S %Z') ----"

echo "[F1/3] expanded forward-blind / charge-radius scan"
$PYTHON scripts/run_bsm_scan_combined.py \
  --config config/future_scan.yaml \
  --workers "$NPROC" \
  --out results/bsm_scan_future.csv

echo "[F2/3] current-IceCube event-weighted fold for expanded scan"
$PYTHON scripts/run_icecube_event_weighting.py \
  --scan results/bsm_scan_future.csv \
  --config config/detector.yaml \
  --outdir results/future_current

echo "[F3/3] Gen2 + generic future-tau projections"
$PYTHON scripts/run_future_projection.py

echo "Outputs:"
echo "  results/future/future_sensitivity_summary.csv"
echo "  results/future/future_sensitivity_energy_bands.csv"
echo "  results/future/generic_tau_telescope_scaling.csv"
echo "  results/future/future_publishability_gate.json"
echo "---- v0.5 future publishability gate completed $(date '+%Y-%m-%d %H:%M:%S %Z') ----"
