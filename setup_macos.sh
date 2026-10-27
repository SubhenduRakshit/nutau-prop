#!/usr/bin/env bash
set -euo pipefail

cat <<'EOF'
Recommended Apple-Silicon setup (Miniforge/Mambaforge):
  1. Install Miniforge if you do not already have conda/mamba.
  2. Then this script creates an isolated environment.

LHAPDF itself is NOT bundled in this archive.
EOF

if command -v mamba >/dev/null 2>&1; then
  SOLVER=mamba
elif command -v conda >/dev/null 2>&1; then
  SOLVER=conda
else
  echo "Neither mamba nor conda found."
  echo "Install Miniforge, then rerun this script."
  exit 2
fi

$SOLVER create -y -n nutau-prop -c conda-forge python=3.12 lhapdf numpy scipy pyyaml matplotlib pandas pytest
echo
echo "Activate with:"
echo "  conda activate nutau-prop"
echo "Then install this package:"
echo "  python -m pip install -e ."
echo
echo "Install the PDF set (after activation):"
echo "  lhapdf install CT18NLO"
echo
echo "If 'lhapdf install' is not on PATH, use lhapdf-config --datadir"
echo "and install CT18NLO with the LHAPDF data installer supplied by your conda package."
