#!/usr/bin/env bash
set -euo pipefail
REF=${1:-results/v0.29/gluon_d8_C8pi.yoda}
OUT=${2:-results/v0.29/contur_scan}
CREF=25.132741228718345
mkdir -p "$OUT"
for C in 25.132741 24 22 21 20.8 20 18 16 12; do
  S=$(python - <<PY
c=$C; r=$CREF
print((c/r)**2)
PY
)
  TAG=$(printf '%05.1f' "$C" | tr '.' 'p')
  Y="$OUT/C${TAG}.yoda"
  yodascale -c ".* x $S" -o "$Y" "$REF"
  mkdir -p "$OUT/C${TAG}"
  contur "$Y" -o "$OUT/C${TAG}" --nomultip 2>&1 | tee "$OUT/C${TAG}.log"
done
printf 'Scan complete.  5%% UHE survival requires |C_G| >= 20.8 at Lambda=4 TeV.\n'
