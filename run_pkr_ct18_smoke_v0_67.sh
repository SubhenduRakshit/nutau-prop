#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BROWSER=true
export PKR_CM_ENV="${PKR_CM_ENV:-nutau-checkmate-v051}"
export PKR_NORM_NEVENTS="${PKR_NORM_NEVENTS:-1000}"
mkdir -p cache/v0.67_ct18_smoke results/v0.67_ct18_smoke

eval "$(conda shell.bash hook)"
set +u
conda activate "$PKR_CM_ENV"
set -u
source cache/v0.51/env.sh
export BROWSER=true

set +e
python scripts/pkr_ct18_smoke_v0_67.py 2>&1 | tee results/v0.67_ct18_smoke/run.log
RC=${PIPESTATUS[0]}
set -e

git add results/v0.67_ct18_smoke
git commit -m "results: v0.67 CT18 compatibility smoke test" || true
git push || true

if [ "$RC" -eq 0 ]; then
  echo "v0.67 CT18NLO SMOKE: PASS"
  echo "Next: exactly one CT18NLO CheckMATE benchmark; no scan."
else
  echo "v0.67 CT18NLO SMOKE: FAIL (rc=$RC)"
  echo "Stop here; do not run CheckMATE."
fi
exit "$RC"
