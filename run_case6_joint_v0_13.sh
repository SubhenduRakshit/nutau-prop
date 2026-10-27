#!/usr/bin/env bash
set -euo pipefail
PYTHON="${PYTHON:-python3}"
mkdir -p results/v0.13
OUT="results/v0.13/run_case6_joint_v0_13_summary.txt"
exec > >(tee "$OUT") 2>&1

echo "==== Case-6 joint propagation+detector fold started $(date '+%Y-%m-%d %H:%M:%S %Z') ===="

for f in \
  results/v0.11/tables/case6_L1p5_nu.npz \
  results/v0.12/case6_earth_transport_v0_12.csv \
  data/uhe_official/icecube_ehe_effective_area.csv \
  data/uhe_official/icecube_ehe_differential_limit.csv \
  data/uhe_official/km3.json \
  data/uhe_official/flux.json; do
  [[ -f "$f" ]] || { echo "MISSING prerequisite: $f"; exit 2; }
done

echo "[J1/2] end-to-end ARCA + IceCube + KM3NeT public-response fold"
"$PYTHON" scripts/run_case6_joint_event_fold_v0_13_cached.py

echo "[J2/2] numerical bookkeeping gate"
"$PYTHON" scripts/validate_case6_joint_v0_13.py

echo "==== Case-6 joint fold completed $(date '+%Y-%m-%d %H:%M:%S %Z') ===="
echo "Summary: $OUT"
