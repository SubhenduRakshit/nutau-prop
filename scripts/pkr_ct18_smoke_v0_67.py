#!/usr/bin/env python3
"""v0.67: deterministic CT18NLO/LHAPDF smoke test.

v0.66 still hit AlphaS_Lambda4 at runtime.  This version does not assume which
MG5 LHAPDF bridge variant is selected.  It patches BOTH pdf_lhapdf6.cc and
pdf_lhapdf62.cc in the generated process, replacing every legacy Lambda4/5
metadata lookup by the fallback value already encoded by MG5 (0.0 in common
blocks, -1.0 in getLam4/getLam5).  It then requires zero remaining
AlphaS_Lambda4/5 strings, explicitly rebuilds libpdf, and only then runs the
single BR=0.45 CT18NLO sample.  No CheckMATE run occurs.
"""
from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = ROOT / "scripts" / "pkr_ct18_smoke_v0_66.py"
spec = importlib.util.spec_from_file_location("pkr_ct18_v066", BASE_PATH)
if spec is None or spec.loader is None:
    raise RuntimeError(f"cannot import {BASE_PATH}")
v066 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v066)

v065 = v066.v065
v065.WORK = ROOT / "cache" / "v0.67_ct18_smoke"
v065.OUT = ROOT / "results" / "v0.67_ct18_smoke"


def patch_one_source(src: Path) -> int:
    if not src.exists():
        return 0
    lines = src.read_text(errors="replace").splitlines(keepends=True)
    out = []
    changed = 0
    for raw in lines:
        line = raw
        if "AlphaS_Lambda4" in line:
            indent = line[: len(line) - len(line.lstrip())]
            nl = "\n" if line.endswith("\n") else ""
            if "w50512_.qcdl4" in line:
                line = f"{indent}w50512_.qcdl4 = 0.0;{nl}"
            elif "lhapdfr_.qcdlha4" in line:
                line = f"{indent}lhapdfr_.qcdlha4 = 0.0;{nl}"
            elif "qcdl4 =" in line:
                line = f"{indent}qcdl4 = -1.0;{nl}"
            elif "return" in line:
                line = f"{indent}return -1.0;{nl}"
            else:
                raise RuntimeError(f"unrecognised AlphaS_Lambda4 lookup in {src}: {raw.strip()}")
            changed += 1
        if "AlphaS_Lambda5" in line:
            indent = line[: len(line) - len(line.lstrip())]
            nl = "\n" if line.endswith("\n") else ""
            if "w50512_.qcdl5" in line:
                line = f"{indent}w50512_.qcdl5 = 0.0;{nl}"
            elif "lhapdfr_.qcdlha5" in line:
                line = f"{indent}lhapdfr_.qcdlha5 = 0.0;{nl}"
            elif "qcdl5 =" in line:
                line = f"{indent}qcdl5 = -1.0;{nl}"
            elif "return" in line:
                line = f"{indent}return -1.0;{nl}"
            else:
                raise RuntimeError(f"unrecognised AlphaS_Lambda5 lookup in {src}: {raw.strip()}")
            changed += 1
        out.append(line)
    text = "".join(out)
    if "AlphaS_Lambda4" in text or "AlphaS_Lambda5" in text:
        raise RuntimeError(f"Lambda metadata lookup survived patch in {src}")
    src.write_text(text)
    return changed


def explicit_rebuild_pdf(procdir: Path) -> None:
    pdfdir = procdir / "Source" / "PDF"
    log = v065.OUT / "pdf_rebuild.log"
    v065.OUT.mkdir(parents=True, exist_ok=True)
    with log.open("w") as f:
        for cmd in (["make", "clean"], ["make", "-j1"]):
            p = subprocess.run(cmd, cwd=pdfdir, stdout=f, stderr=subprocess.STDOUT, text=True)
            if p.returncode != 0:
                tail = "\n".join(log.read_text(errors="replace").splitlines()[-80:])
                raise RuntimeError(f"explicit PDF rebuild failed rc={p.returncode}: {' '.join(cmd)}\n{tail}")
    libpdf = procdir / "lib" / "libpdf.a"
    if not libpdf.exists():
        raise RuntimeError(f"explicit PDF rebuild produced no {libpdf}")
    print(f"{procdir.name}: explicitly rebuilt patched libpdf.a")


def patch_lhapdf_lambda_fallback_v067(procdir: Path) -> None:
    pdfdir = procdir / "Source" / "PDF"
    sources = [pdfdir / "pdf_lhapdf6.cc", pdfdir / "pdf_lhapdf62.cc"]
    counts = {p.name: patch_one_source(p) for p in sources}
    total = sum(counts.values())
    if total == 0:
        raise RuntimeError(f"no AlphaS_Lambda4/5 lookups found in generated PDF bridge: {counts}")
    # Both source variants must be clean if present, so MG5 cannot select an
    # unpatched bridge based on its LHAPDF version bookkeeping.
    for p in sources:
        if p.exists():
            t = p.read_text(errors="replace")
            if "AlphaS_Lambda4" in t or "AlphaS_Lambda5" in t:
                raise RuntimeError(f"unpatched Lambda metadata remains in {p}")
    for stale in [pdfdir / "pdf_lhapdf6.o", pdfdir / "pdf_lhapdf62.o",
                  procdir / "lib" / "libpdf.a", procdir / "lib" / "libpdf.dylib"]:
        stale.unlink(missing_ok=True)
    print(f"{procdir.name}: removed all Lambda4/5 metadata reads from both bridge variants: {counts}")
    explicit_rebuild_pdf(procdir)


v065.patch_lhapdf_lambda_fallback = patch_lhapdf_lambda_fallback_v067

# Retain v0.66's process-local Safari/notification suppression.
v065.base.generate_events = v066.generate_events_v066

if __name__ == "__main__":
    rc = v065.main()
    print(f"v0.67 CT18NLO smoke result: {'PASS' if rc == 0 else 'FAIL'}")
    sys.exit(rc)
