#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p results/v0.16
printf '==== v0.16 resolved mediator started %s ====\n' "$(date)"
python scripts/build_resolved_mediator_tables_v0_16.py --workers "${NUTAU_NPROC:-14}"
python scripts/run_resolved_mediator_transport_v0_16.py
printf '==== v0.16 resolved mediator completed %s ====\n' "$(date)"
printf 'Cross-section gate: results/v0.16/resolved_mediator_shape_gate_v0_16.csv\n'
printf 'Transport gate: results/v0.16/resolved_mediator_transport_summary_v0_16.csv\n'
