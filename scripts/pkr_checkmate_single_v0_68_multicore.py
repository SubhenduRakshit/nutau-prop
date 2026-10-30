#!/usr/bin/env python3
"""Multicore launcher for the v0.68 single CheckMATE benchmark.

Only the MG5/MadEvent event-generation stage is parallelised.  The underlying
v0.68 physics, CT18NLO compatibility patch, Pythia shower and CheckMATE logic
are unchanged.
"""
from __future__ import annotations

import importlib.util
import os
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = ROOT / "scripts" / "pkr_checkmate_single_v0_68.py"
spec = importlib.util.spec_from_file_location("pkr_cm_v068", BASE_PATH)
if spec is None or spec.loader is None:
    raise RuntimeError(f"cannot import {BASE_PATH}")
v068 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v068)

NCORE = int(os.environ.get("PKR_MG5_NCORE", "14"))
if NCORE < 1:
    raise RuntimeError(f"PKR_MG5_NCORE must be >=1, got {NCORE}")

previous_generate_events = v068.base.generate_events


def set_config_option(text: str, key: str, value: str) -> str:
    pat = re.compile(rf"(?m)^\s*#?\s*{re.escape(key)}\s*=.*$")
    line = f"{key} = {value}"
    if pat.search(text):
        return pat.sub(line, text, count=1)
    return text.rstrip() + "\n" + line + "\n"


def configure_mg5_multicore(procdir: Path) -> None:
    cfg = procdir / "Cards" / "me5_configuration.txt"
    if not cfg.exists():
        raise RuntimeError(f"process-local MG5 configuration missing: {cfg}")
    text = cfg.read_text(errors="replace")
    text = set_config_option(text, "run_mode", "2" if NCORE > 1 else "0")
    text = set_config_option(text, "nb_core", str(NCORE))
    cfg.write_text(text)

    verify = cfg.read_text(errors="replace")
    if NCORE > 1:
        if not re.search(r"(?m)^\s*run_mode\s*=\s*2\s*$", verify):
            raise RuntimeError("failed to persist MG5 run_mode=2")
        if not re.search(rf"(?m)^\s*nb_core\s*=\s*{NCORE}\s*$", verify):
            raise RuntimeError(f"failed to persist MG5 nb_core={NCORE}")
    print(f"{procdir.name}: MG5 multicore enabled with {NCORE} workers")


def generate_events_multicore(procdir: Path, log: Path) -> None:
    configure_mg5_multicore(procdir)
    previous_generate_events(procdir, log)


v068.base.generate_events = generate_events_multicore

if __name__ == "__main__":
    sys.exit(v068.main())
