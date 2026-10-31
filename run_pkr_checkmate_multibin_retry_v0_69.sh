#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BROWSER=true
export PKR_CM_ENV="${PKR_CM_ENV:-nutau-checkmate-v051}"
mkdir -p cache/v0.69_multibin_retry results/v0.69_multibin_retry

eval "$(conda shell.bash hook)"
set +u
conda activate "$PKR_CM_ENV"
set -u
source cache/v0.51/env.sh
export BROWSER=true

python - <<'PY'
import spey, spey_pyhf, pyhf
print("stats imports: spey/spey_pyhf/pyhf OK")
PY

set +e
python scripts/pkr_checkmate_multibin_retry_v0_69.py 2>&1 | tee results/v0.69_multibin_retry/run.log
RC=${PIPESTATUS[0]}
set -e

git add results/v0.69_multibin_retry
git commit -m "results: v0.69 CheckMATE multibin retry" || true
git push || true
exit "$RC"
