#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export NUTAU_NPROC="${NUTAU_NPROC:-14}"
mkdir -p results/v0.76_eq4p4_mass_prescreen
python scripts/pkr_eq4p4_mass_prescreen_v0_76.py 2>&1 | tee results/v0.76_eq4p4_mass_prescreen/run.txt
git add results/v0.76_eq4p4_mass_prescreen
git commit -m "results: v0.76 Eq4.4 mass-dependent UHE/dijet prescreen" || true
git push || true
