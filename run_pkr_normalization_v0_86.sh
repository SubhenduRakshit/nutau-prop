#!/usr/bin/env bash
set -euo pipefail
PKR_ROOT="$(cd "$(dirname "$0")" && pwd)"
PKR_PYTHON="${PKR_PYTHON:-/Users/subhendurakshit/miniforge3/bin/python}"
export PYTHONPATH="$PKR_ROOT/src"
export NUTAU_NPROC=14
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd "$PKR_ROOT"
mkdir -p results/v0.86_normalization
"$PKR_PYTHON" scripts/check_pkr_traces_v0_86.py > results/v0.86_normalization/exact_trace_run.txt
"$PKR_PYTHON" scripts/pkr_normalization_fold_v0_86.py --stage initial > results/v0.86_normalization/initial_run.txt
"$PKR_PYTHON" scripts/pkr_normalization_fold_v0_86.py --stage verify > results/v0.86_normalization/verify_run.txt
"$PKR_PYTHON" scripts/pkr_refine_screens_v0_86.py > results/v0.86_normalization/refinement_run.txt
# Avoid the pre-existing unrelated root py.py shadowing pytest's dependency.
(cd /private/tmp && "$PKR_PYTHON" -m pytest "$PKR_ROOT/tests/test_pkr_resolved.py" -q) > results/v0.86_normalization/tests.txt
"$PKR_PYTHON" scripts/pkr_provenance_v0_86.py
cd manuscript/v0.86
TEXINPUTS="../v0.9:${TEXINPUTS:-}" latexmk -pdf -interaction=nonstopmode -halt-on-error normalization.tex
printf '%s\n' 'v0.86 complete: analytic checks, dynamic-Q folds, conditional screens, and PDF.'
