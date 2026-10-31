#!/usr/bin/env bash
set -eo pipefail
cd "$(dirname "$0")"
mkdir -p results/v0.88b_cgc_allflavor cache/v0.88b_cgc_allflavor

eval "$(conda shell.bash hook)"
set +u
conda activate "${NUTAU_ENV:-nutau-checkmate-v051}"
[ -f cache/v0.51/env.sh ] && source cache/v0.51/env.sh
set -u

export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
export NUTAU_NPROC="${NUTAU_NPROC:-14}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "v0.88b: CGC all-flavour PREM propagation with ${NUTAU_NPROC} workers"
python scripts/cgc_earth_allflavor_v0_88b.py | tee results/v0.88b_cgc_allflavor/run.txt

git add results/v0.88b_cgc_allflavor
git commit -m "results: v0.88b CGC all-flavour Earth propagation" || true
git push || true
