#!/usr/bin/env bash
set -euo pipefail
PYTHON="${PYTHON:-python3}"; NPROC="${NUTAU_NPROC:-8}"
mkdir -p results/historical cache; touch PROGRESS.log
exec > >(tee -a PROGRESS.log) 2>&1

echo
echo "---- v0.3 historical validation started $(date '+%Y-%m-%d %H:%M:%S %Z') ----"
echo "workers=$NPROC"
echo "[H1/5] historical neutrino table"
$PYTHON scripts/precompute_xsec.py --config config/historical.yaml --workers "$NPROC" \
  --out cache/historical_nu.npz --kind sm
echo "[H2/5] historical antineutrino table"
$PYTHON scripts/precompute_xsec.py --config config/historical.yaml --workers "$NPROC" \
  --out cache/historical_antinu.npz --kind sm --antinu
echo "[H3/5] Rakshit-Reya Z2/Z3 + direct transport"
$PYTHON scripts/validate_rakshit_reya.py --workers "$NPROC"
echo "[H4/5] Dutta-Reno-Sarcevic secondary channel"
$PYTHON scripts/validate_dutta_reno_sarcevic.py --workers "$NPROC"
echo "[H5/5] combined historical gate"
$PYTHON scripts/run_historical_gate.py
echo "---- historical validation completed $(date '+%Y-%m-%d %H:%M:%S %Z') ----"
