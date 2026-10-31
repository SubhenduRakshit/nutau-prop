#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BROWSER=true
export PKR_CM_ENV="${PKR_CM_ENV:-nutau-checkmate-v051}"
export PKR_NORM_NEVENTS="${PKR_NORM_NEVENTS:-20000}"
export PKR_MG5_NCORE="${PKR_MG5_NCORE:-14}"
mkdir -p results/v0.78_eq4p4_full_envelope results/v0.79_eq4p4_mass_acceptance_pilot cache/v0.79_eq4p4_mass_acceptance_pilot

if [ ! -f results/v0.78_eq4p4_full_envelope/mass_summary.csv ]; then
  python scripts/pkr_eq4p4_full_envelope_v0_78.py | tee results/v0.78_eq4p4_full_envelope/run.txt
fi

eval "$(conda shell.bash hook)"
if ! conda env list | awk '{print $1}' | grep -qx "$PKR_CM_ENV" || [ ! -f cache/v0.51/env.sh ]; then
  bash setup_pkr_checkmate_v0_51.sh
fi
set +u; conda activate "$PKR_CM_ENV"; set -u
source cache/v0.51/env.sh
export BROWSER=true

echo "v0.79: 50 GeV reused; 100/220 GeV pilots, events=$PKR_NORM_NEVENTS, MG5 workers=$PKR_MG5_NCORE"
python scripts/pkr_eq4p4_mass_acceptance_pilot_v0_79.py 2>&1 | tee results/v0.79_eq4p4_mass_acceptance_pilot/run.txt

git add results/v0.78_eq4p4_full_envelope results/v0.79_eq4p4_mass_acceptance_pilot
git commit -m "results: v0.78 full envelope and v0.79 mass acceptance pilot" || true
git push || true
