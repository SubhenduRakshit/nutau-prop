#!/usr/bin/env bash
set -euo pipefail
PYTHON="${PYTHON:-python3}"
NPROC="${NUTAU_NPROC:-14}"
mkdir -p results/v0.11 results/v0.12
OUT="results/v0.12/run_case6_v0_11_v0_12_summary.txt"
exec > >(tee "$OUT") 2>&1

echo "==== Case-6 ordered gate started $(date '+%Y-%m-%d %H:%M:%S %Z') | workers=$NPROC ===="

# Preserve the already-established historical validation. Re-run only when the
# machine-local gate file is absent; this prevents silently proceeding after a
# clean checkout that lost generated validation outputs.
if [[ ! -f results/historical_gate.json ]]; then
  echo "[H] historical RR/DRS gate outputs absent locally -> rerunning validation"
  NUTAU_NPROC="$NPROC" ./run_historical_validation.sh
else
  "$PYTHON" - <<'PY'
import json
q=json.load(open('results/historical_gate.json'))
if not q.get('PASS'):
    raise SystemExit('historical RR/DRS gate is not PASS')
print('[H] historical RR/DRS combined gate PASS')
PY
fi

echo "[3a] build exact dSigma/dy Case-6 tables"
NUTAU_NPROC="$NPROC" "$PYTHON" scripts/build_case6_tables_v0_11.py --workers "$NPROC"

echo "[3b] differential-kernel regression gate"
"$PYTHON" scripts/validate_case6_differential_v0_11.py

echo "[3c] channel-resolved detector NSI fold (no modified Earth propagation yet)"
"$PYTHON" scripts/run_case6_detector_v0_11.py

echo "[1] full flavor-aware Earth transport with the same dSigma/dy kernel"
NUTAU_NPROC="$NPROC" "$PYTHON" scripts/run_case6_transport_v0_12.py --workers "$NPROC"

echo "==== Case-6 ordered gate completed $(date '+%Y-%m-%d %H:%M:%S %Z') ===="
echo "Summary: $OUT"
