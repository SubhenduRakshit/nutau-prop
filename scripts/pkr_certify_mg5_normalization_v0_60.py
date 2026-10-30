#!/usr/bin/env python3
"""v0.60 deterministic MG5 normalization certification.

v0.59 reached process-local event generation, but invoked

    bin/generate_events -f

without an explicit run name.  The robust MadEvent batch form is

    bin/generate_events run_01 -f

v0.60 changes only that invocation and adds a compact Events-tree diagnostic
if the expected LHE is still absent.  All card edits and normalization gates
remain unchanged.  No detector simulation is run here.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = ROOT / "scripts" / "pkr_certify_mg5_normalization_v0_59.py"
spec = importlib.util.spec_from_file_location("pkr_norm_v059", BASE_PATH)
if spec is None or spec.loader is None:
    raise RuntimeError(f"cannot import {BASE_PATH}")
v059 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v059)

# v059.base is the v0.55 deterministic implementation; its main() is reused.
base = v059.base

base.WORK = ROOT / "cache" / "v0.60_mg5_norm"
base.OUT = ROOT / "results" / "v0.60_mg5_norm"
base.PASS_MARKER = base.WORK / "NORMALIZATION_CERTIFIED.json"


def generate_events_v060(procdir: Path, log: Path) -> None:
    run_name = "run_01"
    gen = procdir / "bin" / "generate_events"
    madevent = procdir / "bin" / "madevent"
    if gen.exists():
        base.run([gen, run_name, "-f"], log, cwd=procdir)
    elif madevent.exists():
        base.run([madevent, "generate_events", run_name, "-f"], log, cwd=procdir)
    else:
        raise RuntimeError(f"no process-local event generator found below {procdir}/bin")

    # Fail with useful local structure rather than another opaque '0 LHE' error.
    expected_gz = procdir / "Events" / run_name / "unweighted_events.lhe.gz"
    expected = procdir / "Events" / run_name / "unweighted_events.lhe"
    if not expected_gz.exists() and not expected.exists():
        evdir = procdir / "Events"
        listing = []
        if evdir.exists():
            for p in sorted(evdir.rglob("*")):
                if p.is_file():
                    listing.append(str(p.relative_to(procdir)))
                if len(listing) >= 40:
                    break
        tail = ""
        if log.exists():
            tail = "\n".join(log.read_text(errors="replace").splitlines()[-40:])
        raise RuntimeError(
            f"generate_events returned success but no LHE appeared at {run_name}; "
            f"Events files={listing}\n--- generator log tail ---\n{tail}"
        )


base.generate_events = generate_events_v060

if __name__ == "__main__":
    sys.exit(base.main())
