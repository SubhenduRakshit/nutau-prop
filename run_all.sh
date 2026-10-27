#!/usr/bin/env bash
set -euo pipefail
PYTHON="${PYTHON:-python3}"; NPROC="${NUTAU_NPROC:-8}"
mkdir -p results/historical cache; touch PROGRESS.log
exec > >(tee -a PROGRESS.log) 2>&1

echo "================================================================"
echo "nutau-prop v0.3 started $(date '+%Y-%m-%d %H:%M:%S %Z') | workers=$NPROC"
echo "================================================================"
echo "[1/10] environment"; $PYTHON scripts/check_environment.py
echo "[2/10] modern SM neutrino table"; $PYTHON scripts/precompute_xsec.py --workers "$NPROC" --out cache/sm.npz --kind sm
echo "[3/10] modern SM validation"; $PYTHON scripts/validate_sm_xsec.py | tee results/sm_xsec_validation.txt
echo "[4/10] historical neutrino table"; $PYTHON scripts/precompute_xsec.py --config config/historical.yaml --workers "$NPROC" --out cache/historical_nu.npz --kind sm
echo "[5/10] historical antineutrino table"; $PYTHON scripts/precompute_xsec.py --config config/historical.yaml --workers "$NPROC" --out cache/historical_antinu.npz --kind sm --antinu
echo "[6/10] Rakshit-Reya"; $PYTHON scripts/validate_rakshit_reya.py --workers "$NPROC"
echo "[7/10] Dutta-Reno-Sarcevic"; $PYTHON scripts/validate_dutta_reno_sarcevic.py --workers "$NPROC"
echo "[8/10] MANDATORY historical gate"; $PYTHON scripts/run_historical_gate.py
echo "[9/10] BSM scan"; $PYTHON scripts/run_bsm_scan.py
echo "[10/10] BSM gate"; $PYTHON scripts/run_publishability_gate.py | tee results/gate.txt
echo "Completed $(date '+%Y-%m-%d %H:%M:%S %Z')"
