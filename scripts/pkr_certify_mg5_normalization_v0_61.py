#!/usr/bin/env python3
"""v0.61 deterministic MG5 normalization certification.

v0.60 reached the actual MadEvent compile stage and exposed a macOS/arm64
linker mismatch: the generated MG5 Source/make_opts inherited

    STDLIB=-lstdc++

while the conda-forge LHAPDF library in this environment is built with
Apple/Clang libc++ (the undefined std::__1 symbols in the linker failure).

v0.61 keeps every physics/card/normalization check unchanged and patches only
the freshly generated process-local Source/make_opts on Darwin:

    STDLIB=-lstdc++  ->  STDLIB=-lc++

The patch is strict: if neither the expected old nor already-correct setting is
present, the run stops rather than editing blindly.  No CheckMATE run occurs.
"""
from __future__ import annotations

import importlib.util
import platform
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = ROOT / "scripts" / "pkr_certify_mg5_normalization_v0_60.py"
spec = importlib.util.spec_from_file_location("pkr_norm_v060", BASE_PATH)
if spec is None or spec.loader is None:
    raise RuntimeError(f"cannot import {BASE_PATH}")
v060 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v060)

base = v060.base
previous_generate_events = base.generate_events

base.WORK = ROOT / "cache" / "v0.61_mg5_norm"
base.OUT = ROOT / "results" / "v0.61_mg5_norm"
base.PASS_MARKER = base.WORK / "NORMALIZATION_CERTIFIED.json"


def patch_macos_stdlib(procdir: Path) -> None:
    if platform.system() != "Darwin":
        return
    make_opts = procdir / "Source" / "make_opts"
    if not make_opts.exists():
        raise RuntimeError(f"generated MG5 make_opts missing: {make_opts}")
    text = make_opts.read_text(errors="replace")
    old_re = re.compile(r"(?m)^\s*STDLIB\s*=\s*-lstdc\+\+\s*$")
    new_re = re.compile(r"(?m)^\s*STDLIB\s*=\s*-lc\+\+\s*$")
    n_old = len(old_re.findall(text))
    n_new = len(new_re.findall(text))
    if n_old == 1 and n_new == 0:
        text = old_re.sub("STDLIB=-lc++", text, count=1)
        make_opts.write_text(text)
        print(f"{procdir.name}: patched MG5 STDLIB -lstdc++ -> -lc++ for macOS/LHAPDF")
    elif n_old == 0 and n_new == 1:
        print(f"{procdir.name}: MG5 STDLIB already uses -lc++")
    else:
        nearby = [ln.strip() for ln in text.splitlines() if "STDLIB" in ln][:8]
        raise RuntimeError(
            f"refusing blind make_opts patch: old={n_old}, new={n_new}, nearby={nearby}"
        )


def generate_events_v061(procdir: Path, log: Path) -> None:
    patch_macos_stdlib(procdir)
    previous_generate_events(procdir, log)


base.generate_events = generate_events_v061

if __name__ == "__main__":
    sys.exit(base.main())
