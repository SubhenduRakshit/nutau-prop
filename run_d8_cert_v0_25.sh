#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p results/v0.25

echo '=== v0.25 analytic normalization/unitarity ==='
python scripts/check_d8_analytic_v0_25.py | tee results/v0.25/analytic_gate.txt

echo '=== v0.25 exact dynamic-Q CT18NLO UHE fold: strict unitary-to-Lambda ==='
python scripts/d8_uhe_cert_v0_25.py --pdf CT18NLO --lambda-tev 4 --unitarity-mode to_lambda \
  --xis 1,0.5,0.25 --energies-pev 1,10,100,300 | tee results/v0.25/uhe_gate.txt

echo '=== v0.25 HighPT Run-II di-tau recast ==='
if command -v wolframscript >/dev/null 2>&1; then
  wolframscript -file scripts/highpt_d8_tautau_v0_25.wls | tee results/v0.25/highpt_gate.txt
else
  echo 'SKIP HighPT: wolframscript not found' | tee results/v0.25/highpt_gate.txt
fi

echo '=== v0.25 certification runner finished ==='
