#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p results/v0.84_eq4p3_M50_threshold_refine
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export NUTAU_NPROC="${NUTAU_NPROC:-14}"
echo "v0.84: Eq4.3 M50 threshold-edge refinement; no MG5/Pythia/Delphes/CheckMATE"
python scripts/pkr_eq4p3_M50_threshold_refine_v0_84.py | tee results/v0.84_eq4p3_M50_threshold_refine/run.txt
git add results/v0.84_eq4p3_M50_threshold_refine
git commit -m "results: v0.84 Eq4.3 M50 threshold-edge refinement" || true
git push || true
