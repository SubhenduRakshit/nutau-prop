#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PKR_CM_NPROC="${PKR_CM_NPROC:-6}"
export PKR_CM_NEVENTS="${PKR_CM_NEVENTS:-100000}"
export PKR_CM_ENV="${PKR_CM_ENV:-nutau-checkmate-v051}"
mkdir -p cache/v0.51 results/v0.51

eval "$(conda shell.bash hook)"
if ! conda env list | awk '{print $1}' | grep -qx "$PKR_CM_ENV" || [ ! -f cache/v0.51/env.sh ]; then
  bash setup_pkr_checkmate_v0_51.sh
fi
# Conda's macOS Fortran activation hook is not nounset-safe.
set +u
conda activate "$PKR_CM_ENV"
set -u
source cache/v0.51/env.sh

python scripts/pkr_run_checkmate_v0_51.py | tee results/v0.51/run.log

git add results/v0.51
git commit -m "results: v0.51 PKR MG5 Pythia CheckMATE monojet gate" || true
git push

echo ""
echo "v0.51 complete; shell remains open."
