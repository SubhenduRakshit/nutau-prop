#!/usr/bin/env bash
set -euo pipefail
PYTHON="${PYTHON:-python3}"
NPROC="${NUTAU_NPROC:-8}"
mkdir -p results/arca data/arca_zenodo
touch PROGRESS.log
exec > >(tee -a PROGRESS.log) 2>&1

echo
echo "---- v0.6 KM3NeT/ARCA gate started $(date '+%Y-%m-%d %H:%M:%S %Z') ----"

if [ ! -f results/bsm_scan_future.csv ]; then
  echo "[A1/3] future-scan transport table"
  $PYTHON scripts/run_bsm_scan_combined.py \
    --config config/future_scan.yaml \
    --workers "$NPROC" \
    --out results/bsm_scan_future.csv
else
  echo "[A1/3] reusing results/bsm_scan_future.csv"
fi

echo "[A2/3] ARCA official-record acquisition + track/shower likelihood"
$PYTHON scripts/run_arca_projection.py

echo "[A3/3] outputs"
echo "  results/arca/arca_sensitivity_summary.csv"
echo "  results/arca/arca_event_bins.csv"
echo "  results/arca/arca_publishability_gate.json"
echo "  results/arca/zenodo_inventory.json"

echo "---- v0.6 KM3NeT/ARCA gate completed $(date '+%Y-%m-%d %H:%M:%S %Z') ----"
