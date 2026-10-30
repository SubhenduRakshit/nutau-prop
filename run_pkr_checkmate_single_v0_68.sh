#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

# Keep BLAS/OpenMP single-threaded; MG5 itself uses process-level multicore.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BROWSER=true
export PKR_CM_ENV="${PKR_CM_ENV:-nutau-checkmate-v051}"
export PKR_NORM_NEVENTS="${PKR_NORM_NEVENTS:-100000}"
export PKR_MG5_NCORE="${PKR_MG5_NCORE:-14}"
mkdir -p cache/v0.68_checkmate_single results/v0.68_checkmate_single

eval "$(conda shell.bash hook)"
if ! conda env list | awk '{print $1}' | grep -qx "$PKR_CM_ENV" || [ ! -f cache/v0.51/env.sh ]; then
  bash setup_pkr_checkmate_v0_51.sh
fi
set +u
conda activate "$PKR_CM_ENV"
set -u
source cache/v0.51/env.sh
export BROWSER=true

echo "v0.68: MG5 event generation will use $PKR_MG5_NCORE cores"
set +e
python scripts/pkr_checkmate_single_v0_68_multicore.py 2>&1 | tee results/v0.68_checkmate_single/run.log
RC=${PIPESTATUS[0]}
set -e

git add results/v0.68_checkmate_single
git commit -m "results: v0.68 single conservative CT18 CheckMATE benchmark" || true
git push || true

if [ "$RC" -eq 0 ]; then
  echo ""
  echo "v0.68 SINGLE CHECKMATE BENCHMARK: COMPLETE"
  echo "Inspect multibin_r before making any collider conclusion."
else
  echo ""
  echo "v0.68 SINGLE CHECKMATE BENCHMARK: FAIL (rc=$RC)"
  echo "Stop here; do not run a scan or infer a collider verdict."
fi
exit "$RC"
