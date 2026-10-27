#!/usr/bin/env bash
set -euo pipefail
PYTHON="${PYTHON:-python3}"
NPROC="${NUTAU_NPROC:-14}"
mkdir -p results/v0.15
OUT="results/v0.15/run_momentum_transport_v0_15_summary.txt"
exec > >(tee "$OUT") 2>&1

echo "==== v0.15 momentum-operator transport started $(date '+%Y-%m-%d %H:%M:%S %Z') | workers=$NPROC ===="
echo "[A] exact transport-ready PDF tables for Cases 5/7"
NUTAU_NPROC="$NPROC" "$PYTHON" scripts/build_momentum_tables_v0_15.py --workers "$NPROC"
echo "[B] Case-5/7 chirality-aware transport brackets"
"$PYTHON" scripts/run_momentum_transport_v0_15.py
echo "[C] Case-6 external-state robustness audit (reuses v0.11 tables)"
"$PYTHON" scripts/run_case6_chirality_audit_v0_15.py
echo "==== v0.15 completed $(date '+%Y-%m-%d %H:%M:%S %Z') ===="
echo "Summary: $OUT"
