#!/usr/bin/env bash
set -euo pipefail
PYTHON="${PYTHON:-python3}"
NPROC="${NUTAU_NPROC:-12}"

echo "==== v0.18 corrected transport rerun started $(date) ===="
echo "workers=$NPROC"

echo "[0/2] v0.17.1 refinement reassessment (no transport rerun)"
"$PYTHON" scripts/reassess_tau_transport_v0_17_1.py

echo "[1/2] Cases 5/6/7 corrected transport + absolute transmitted spectral weights"
# The Python production script internally completes contacts first.
# It then reuses the existing v0.16 resolved-mediator cross-section tables.
"$PYTHON" scripts/run_corrected_transport_v0_18.py --workers "$NPROC"

echo "[2/2] completed; detector weighting intentionally deferred"
echo "Outputs:"
echo "  results/v0.18/contact_transport_corrected_v0_18.csv"
echo "  results/v0.18/contact_transport_corrected_summary_v0_18.csv"
echo "  results/v0.18/resolved_mediator_transport_corrected_v0_18.csv"
echo "  results/v0.18/resolved_mediator_transport_corrected_summary_v0_18.csv"
echo "==== v0.18 corrected transport rerun completed $(date) ===="
