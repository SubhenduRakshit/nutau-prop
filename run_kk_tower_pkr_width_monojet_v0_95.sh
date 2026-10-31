#!/usr/bin/env bash
set -euo pipefail

ROOT="${HOME}/Research/nutau-prop"
REPO="https://github.com/SubhenduRakshit/nutau-prop.git"

if [[ ! -d "${ROOT}/.git" ]]; then
  echo "ERROR: canonical repository is missing: ${ROOT}" >&2
  echo "Expected ${REPO}" >&2
  exit 2
fi

cd "${ROOT}"
origin="$(git remote get-url origin 2>/dev/null || true)"
case "${origin}" in
  https://github.com/SubhenduRakshit/nutau-prop|https://github.com/SubhenduRakshit/nutau-prop.git|git@github.com:SubhenduRakshit/nutau-prop.git) ;;
  *) echo "ERROR: unexpected git origin: ${origin}" >&2; exit 3;;
esac

export PYTHONPATH="${ROOT}/src:${PYTHONPATH:-}"
export NUTAU_NPROC="${NUTAU_NPROC:-14}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

mkdir -p results/v0.95_kk_pkr_width_monojet
python scripts/kk_tower_pkr_width_monojet_v0_95.py \
  --nproc "${NUTAU_NPROC}" \
  "$@" \
  2>&1 | tee results/v0.95_kk_pkr_width_monojet/run.log

cat results/v0.95_kk_pkr_width_monojet/PKR_KK_TOWER_WIDTH_MONOJET_v0.95.md
