#!/usr/bin/env python3
"""v0.62 deterministic MG5 normalization certification.

v0.61 confirmed that the generated process-local Source/make_opts already had
STDLIB=-lc++, but the MadEvent gensym link still omitted the C++ runtime.

The upstream MG5 LO makefile links gensym with $(LINKLIBS) $(LDFLAGS), while
Source/make_opts only initializes LDFLAGS from $(STDLIB) when LDFLAGS is
undefined.  In the conda environment LDFLAGS is already defined, so STDLIB
never reaches the linker.

v0.62 keeps every physics/card/normalization check unchanged and adds one
Darwin-only process-local make_opts rule:

    LDFLAGS += $(STDLIB)

This preserves any conda linker flags while guaranteeing that the libc++
runtime required by LHAPDF is present in the final Fortran-driven link.
No CheckMATE run occurs.
"""
from __future__ import annotations

import importlib.util
import platform
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = ROOT / "scripts" / "pkr_certify_mg5_normalization_v0_61.py"
spec = importlib.util.spec_from_file_location("pkr_norm_v061", BASE_PATH)
if spec is None or spec.loader is None:
    raise RuntimeError(f"cannot import {BASE_PATH}")
v061 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v061)

base = v061.base
previous_generate_events = base.generate_events

base.WORK = ROOT / "cache" / "v0.62_mg5_norm"
base.OUT = ROOT / "results" / "v0.62_mg5_norm"
base.PASS_MARKER = base.WORK / "NORMALIZATION_CERTIFIED.json"

MARKER = "# v0.62 Darwin/LHAPDF: force C++ runtime into final link"


def patch_macos_ldflags(procdir: Path) -> None:
    if platform.system() != "Darwin":
        return
    make_opts = procdir / "Source" / "make_opts"
    if not make_opts.exists():
        raise RuntimeError(f"generated MG5 make_opts missing: {make_opts}")
    text = make_opts.read_text(errors="replace")
    if MARKER in text:
        print(f"{procdir.name}: MG5 Darwin LDFLAGS runtime patch already present")
        return
    if "STDLIB=-lc++" not in text and "STDLIB = -lc++" not in text:
        nearby = [ln.strip() for ln in text.splitlines() if "STDLIB" in ln][:8]
        raise RuntimeError(
            f"expected libc++ STDLIB before LDFLAGS patch; nearby={nearby}"
        )
    text += (
        "\n" + MARKER + "\n"
        "ifeq ($(UNAME),Darwin)\n"
        "LDFLAGS += $(STDLIB)\n"
        "endif\n"
    )
    make_opts.write_text(text)
    print(f"{procdir.name}: appended $(STDLIB) to MG5 LDFLAGS on Darwin")


def generate_events_v062(procdir: Path, log: Path) -> None:
    # v0.61 first certifies/patches STDLIB itself; v0.62 ensures it actually
    # reaches the gensym/madevent link even when conda predefines LDFLAGS.
    v061.patch_macos_stdlib(procdir)
    patch_macos_ldflags(procdir)
    # Call v0.60 implementation directly to avoid applying v0.61 twice.
    v061.previous_generate_events(procdir, log)


base.generate_events = generate_events_v062

if __name__ == "__main__":
    sys.exit(base.main())
