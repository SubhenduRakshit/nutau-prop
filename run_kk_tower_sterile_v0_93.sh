#!/usr/bin/env bash
set -euo pipefail

# Canonical compute location for this project.  Do not run from historical
# nutau-prop-cert or nutau-prop-v0.* working copies.
ROOT="${HOME}/Research/nutau-prop"
REPO="https://github.com/SubhenduRakshit/nutau-prop.git"

if [[ ! -d "${ROOT}/.git" ]]; then
  echo "ERROR: canonical repository is missing: ${ROOT}" >&2
  echo "Clone ${REPO} into ${ROOT} first." >&2
  exit 2
fi

cd "${ROOT}"

origin="$(git remote get-url origin 2>/dev/null || true)"
case "${origin}" in
  https://github.com/SubhenduRakshit/nutau-prop|https://github.com/SubhenduRakshit/nutau-prop.git|git@github.com:SubhenduRakshit/nutau-prop.git)
    ;;
  *)
    echo "ERROR: ${ROOT} is not the expected SubhenduRakshit/nutau-prop clone." >&2
    echo "origin=${origin}" >&2
    exit 3
    ;;
esac

export NUTAU_NPROC="${NUTAU_NPROC:-14}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"

mkdir -p results/v0.93_kk_tower
python scripts/kk_tower_sterile_screen_v0_93.py \
  --nproc "${NUTAU_NPROC}" \
  "$@" \
  2>&1 | tee results/v0.93_kk_tower/run.log

cat results/v0.93_kk_tower/KK_TOWER_STERILE_SCREEN_v0.93.md
