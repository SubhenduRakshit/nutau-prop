#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BROWSER=true
export PKR_CM_ENV="${PKR_CM_ENV:-nutau-checkmate-v051}"
export PKR_NORM_NEVENTS="${PKR_NORM_NEVENTS:-100000}"
export PKR_MG5_NCORE="${PKR_MG5_NCORE:-14}"
mkdir -p cache/v0.75_eq4p4_bestpoint_checkmate results/v0.75_eq4p4_bestpoint_checkmate

eval "$(conda shell.bash hook)"
if ! conda env list | awk '{print $1}' | grep -qx "$PKR_CM_ENV" || [ ! -f cache/v0.51/env.sh ]; then
  bash setup_pkr_checkmate_v0_51.sh
fi
set +u; conda activate "$PKR_CM_ENV"; set -u
source cache/v0.51/env.sh
export BROWSER=true

echo "v0.75: exact single best-point benchmark; MG5 workers=$PKR_MG5_NCORE, events=$PKR_NORM_NEVENTS"
set +e
python scripts/pkr_eq4p4_bestpoint_checkmate_v0_75.py 2>&1 | tee results/v0.75_eq4p4_bestpoint_checkmate/run.log
RC=${PIPESTATUS[0]}
set -e

git add results/v0.75_eq4p4_bestpoint_checkmate
git commit -m "results: v0.75 exact Eq4.4 marginal-window CheckMATE benchmark" || true
git push || true
exit "$RC"
