#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

NPROC="${PKR_CM_NPROC:-6}"
ENV_NAME="${PKR_CM_ENV:-nutau-checkmate-v051}"
WORK="$PWD/cache/v0.51"
mkdir -p "$WORK" "$WORK/bin" "$WORK/models"

if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda/miniforge not found in PATH" >&2
  return 2 2>/dev/null || false
fi

eval "$(conda shell.bash hook)"
if ! conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  conda create -y -n "$ENV_NAME" -c conda-forge \
    python=3.11 root=6.38.04 hepmc2=2.06.11 pythia8=8.312 lhapdf=6.5.6 \
    numpy scipy pandas pyhf jsonpatch future autoconf automake libtool make cmake \
    pkg-config compilers fortran-compiler git curl unzip
fi
# Conda's macOS Fortran activation hook may read GFORTRAN before defining it.
# Temporarily disable nounset only for environment activation, then restore it.
set +u
conda activate "$ENV_NAME"
set -u
python -m pip install -q --upgrade spey spey-pyhf || true

if command -v lhapdf >/dev/null 2>&1; then
  lhapdf install CT18NLO >/dev/null 2>&1 || true
fi

MODEL="$WORK/models/DMsimp_s_spin1"
if [ ! -f "$MODEL/particles.py" ]; then
  rm -rf "$WORK/model_unpack" "$MODEL"
  mkdir -p "$WORK/model_unpack"
  ZIP="$WORK/DMsimp_s_spin1.zip"
  if curl -fsSL "https://feynrules.irmp.ucl.ac.be/raw-attachment/wiki/DMsimp/DMsimp_s_spin1.zip" -o "$ZIP" && unzip -q -o "$ZIP" -d "$WORK/model_unpack"; then
    SRC_MODEL="$(find "$WORK/model_unpack" -type f -name particles.py -print -quit | xargs dirname)"
    cp -R "$SRC_MODEL" "$MODEL"
  else
    rm -rf "$WORK/semi-visible-jets"
    git clone -q --depth 1 https://github.com/smsharma/semi-visible-jets.git "$WORK/semi-visible-jets"
    cp -R "$WORK/semi-visible-jets/DMsimp_s_spin1" "$MODEL"
  fi
fi

MG5="$WORK/MG5_aMC"
if [ ! -x "$MG5/bin/mg5_aMC" ]; then
  rm -rf "$MG5"
  git clone -q --depth 1 https://github.com/mg5amcnlo/mg5amcnlo.git "$MG5"
fi

DELPHES="$WORK/Delphes-3.5.1"
if [ ! -f "$DELPHES/libDelphes.so" ] && [ ! -f "$DELPHES/libDelphes.dylib" ]; then
  rm -rf "$DELPHES"
  git clone -q --depth 1 --branch 3.5.1 https://github.com/delphes/delphes.git "$DELPHES" || {
    git clone -q --depth 1 https://github.com/delphes/delphes.git "$DELPHES"
  }
  make -C "$DELPHES" -j"$NPROC"
fi

CMSRC="$WORK/checkmate2"
CMINST="$WORK/checkmate-install"
if [ ! -x "$CMINST/bin/CheckMATE" ]; then
  # Preserve a partially compiled tree so reruns continue rather than starting over.
  if [ ! -d "$CMSRC/.git" ]; then
    rm -rf "$CMSRC"
    git clone -q --depth 1 https://github.com/CheckMATE2/checkmate2.git "$CMSRC"
  fi
  (
    cd "$CMSRC"
    if [ ! -f Makefile ]; then
      autoreconf -i
      ./configure \
        --prefix="$CMINST" \
        --with-delphes="$DELPHES" \
        --with-hepmc="$CONDA_PREFIX" \
        --with-rootsys="$CONDA_PREFIX" \
        --with-python="$CONDA_PREFIX/bin/python"
    fi

    # Apple Clang rejects two legacy ATLAS helper returns that narrow double -> float
    # inside list initialization. This analysis is unrelated to our CMS gate, but
    # CheckMATE builds all analyses. Add explicit casts; numerical values are unchanged.
    ATLAS_SRC="tools/analysis/src/analyses/ATLAS_13TeV/atlas_2411_02040.cc"
    if [ -f "$ATLAS_SRC" ]; then
      python - "$ATLAS_SRC" <<'PY'
from pathlib import Path
import sys
p = Path(sys.argv[1])
s = p.read_text()
old = "return {1.5*lambdas[0], 1.5*(lambdas[0] + lambdas[1]), 2.*lambdas[1]/(lambdas[1] + lambdas[2])};"
new = "return {static_cast<float>(1.5*lambdas[0]), static_cast<float>(1.5*(lambdas[0] + lambdas[1])), static_cast<float>(2.*lambdas[1]/(lambdas[1] + lambdas[2]))};"
if old in s:
    n = s.count(old)
    p.write_text(s.replace(old, new))
    print(f"Applied Apple-Clang narrowing fix to {p} ({n} occurrence(s)).")
elif new in s:
    print(f"Apple-Clang narrowing fix already present in {p}.")
else:
    raise SystemExit(f"ERROR: expected narrowing pattern not found in {p}; refusing blind patch")
PY
    fi

    BUILDLOG="$WORK/checkmate_build_serial.log"
    echo "Continuing CheckMATE build serially; full log: $BUILDLOG"
    set +e
    make -j1 V=1 >"$BUILDLOG" 2>&1
    RC=$?
    set -e
    if [ "$RC" -ne 0 ]; then
      echo ""
      echo "=== CheckMATE decisive build diagnostics ==="
      grep -nE "error:|fatal error:|undefined symbols|Undefined symbols|ld: |linker command failed|No such file|cannot find" "$BUILDLOG" | tail -n 60 || true
      echo ""
      echo "=== Last 100 build-log lines ==="
      tail -n 100 "$BUILDLOG"
      return "$RC" 2>/dev/null || exit "$RC"
    fi

    make install >>"$BUILDLOG" 2>&1 || {
      echo "CheckMATE make install failed; tail follows:"
      tail -n 100 "$BUILDLOG"
      return 5 2>/dev/null || exit 5
    }
  )
fi

# CheckMATE internally names the CMS EXO-20-004 / arXiv:2107.13021 analysis
# cms_2107_13201 (digits transposed). Use the package's actual identifier.
CM_ANALYSIS="cms_2107_13201"
if ! grep -Rqs "$CM_ANALYSIS" "$CMSRC/data" "$CMINST/share" 2>/dev/null; then
  echo "ERROR: current CheckMATE build does not contain $CM_ANALYSIS" >&2
  return 3 2>/dev/null || false
fi

SHOWER="$WORK/bin/pythia_lhe_to_hepmc_v0_51"
SRC="$PWD/scripts/pythia_lhe_to_hepmc_v0_51.cc"
if [ ! -x "$SHOWER" ] || [ "$SRC" -nt "$SHOWER" ]; then
  CXX_BIN="${CXX:-c++}"
  PY_CXX="$(pythia8-config --cxxflags 2>/dev/null || true)"
  PY_LIB="$(pythia8-config --libs 2>/dev/null || true)"
  HEP_CXX="$(HepMC-config --cflags 2>/dev/null || true)"
  HEP_LIB="$(HepMC-config --libs 2>/dev/null || true)"
  "$CXX_BIN" -std=c++17 -O2 -I"$CONDA_PREFIX/include" $PY_CXX $HEP_CXX "$SRC" \
    -L"$CONDA_PREFIX/lib" $PY_LIB $HEP_LIB -lpythia8 -lHepMC -o "$SHOWER"
fi

cat > "$WORK/env.sh" <<EOF
export PKR_MG5="$MG5/bin/mg5_aMC"
export PKR_CHECKMATE="$CMINST/bin/CheckMATE"
export PKR_CM_MODEL="$MODEL"
export PKR_CM_ANALYSIS="$CM_ANALYSIS"
export DYLD_LIBRARY_PATH="$DELPHES:$CONDA_PREFIX/lib:\${DYLD_LIBRARY_PATH:-}"
export LD_LIBRARY_PATH="$DELPHES:$CONDA_PREFIX/lib:\${LD_LIBRARY_PATH:-}"
EOF

echo "v0.51 collider environment ready"
echo "  env:       $ENV_NAME"
echo "  MG5:       $MG5/bin/mg5_aMC"
echo "  CheckMATE: $CMINST/bin/CheckMATE"
echo "  analysis:  $CM_ANALYSIS (CMS EXO-20-004 / arXiv:2107.13021)"
