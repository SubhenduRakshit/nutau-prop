#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PKR_CM_ENV="${PKR_CM_ENV:-nutau-checkmate-v051}"
export PKR_NORM_NEVENTS="${PKR_NORM_NEVENTS:-5000}"
mkdir -p cache/v0.61_mg5_norm results/v0.61_mg5_norm

eval "$(conda shell.bash hook)"
if ! conda env list | awk '{print $1}' | grep -qx "$PKR_CM_ENV" || [ ! -f cache/v0.51/env.sh ]; then
  bash setup_pkr_checkmate_v0_51.sh
fi
set +u
conda activate "$PKR_CM_ENV"
set -u
source cache/v0.51/env.sh

set +e
python scripts/pkr_certify_mg5_normalization_v0_61.py 2>&1 | tee results/v0.61_mg5_norm/run.log
RC=${PIPESTATUS[0]}
set -e

git add results/v0.61_mg5_norm
git commit -m "results: v0.61 MG5 normalization certification" || true
git push || true

if [ "$RC" -eq 0 ]; then
  echo ""
  echo "v0.61 MG5 NORMALIZATION: PASS"
  echo "Next step: exactly one CheckMATE benchmark; do not start a full campaign."
else
  echo ""
  echo "v0.61 MG5 NORMALIZATION: FAIL (rc=$RC)"
  echo "Stop here; no collider exclusion/survival interpretation is allowed."
fi
exit "$RC"
