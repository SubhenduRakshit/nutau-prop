#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
export NUTAU_NPROC="${NUTAU_NPROC:-14}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p results/v0.74_eq4p4_survival_envelope
python scripts/pkr_eq4p4_survival_envelope_v0_74.py | tee results/v0.74_eq4p4_survival_envelope/run.txt
git add results/v0.74_eq4p4_survival_envelope
git commit -m "results: v0.74 Eq4.4 M50 survival envelope" || true
git push || true
