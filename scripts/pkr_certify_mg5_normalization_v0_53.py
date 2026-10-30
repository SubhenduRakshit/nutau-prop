#!/usr/bin/env python3
"""v0.53 MG5 normalization certification.

v0.52 correctly diagnosed that all persisted cards and event samples were the
UFO defaults.  The cause was execution semantics: `launch ... -f` forced the
run immediately, so the subsequent launch-time `set` commands were never
applied.  This wrapper reuses the strict v0.52 certification logic but changes
only the launch command to plain `launch`, as required for scripted card edits.

No CheckMATE run is performed here.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = ROOT / "scripts" / "pkr_certify_mg5_normalization_v0_52.py"

spec = importlib.util.spec_from_file_location("pkr_norm_v052", BASE_PATH)
if spec is None or spec.loader is None:
    raise RuntimeError(f"cannot import {BASE_PATH}")
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)

# Preserve v0.52 outputs; v0.53 gets a fresh cache/results namespace.
base.WORK = ROOT / "cache" / "v0.53_mg5_norm"
base.OUT = ROOT / "results" / "v0.53_mg5_norm"
base.PASS_MARKER = base.WORK / "NORMALIZATION_CERTIFIED.json"


def make_card_v053(procdir: Path, gamma_total: float, gchi: float) -> str:
    sets = [
        f"set MY1 {base.M:.16g}",
        f"set MXd {base.MCHI:.16g}",
        f"set WY1 {gamma_total:.16g}",
        f"set gVXd {gchi:.16g}",
        "set gAXd 0",
        "set gVXc 0",
    ]
    for p in ["gVd11", "gVu11", "gVd22", "gVu22", "gVd33", "gVu33"]:
        sets.append(f"set {p} {base.GQ:.16g}")
    for p in [
        "gAd11", "gAu11", "gAd22", "gAu22", "gAd33", "gAu33",
        "gVl11", "gVl22", "gVl33", "gAl11", "gAl22", "gAl33",
        "gnu11", "gnu22", "gnu33", "gVu31", "gAu31", "gVd31", "gAd31", "gVh",
    ]:
        sets.append(f"set {p} 0")
    runsets = [
        f"set nevents {base.NEVENTS}",
        "set ebeam1 6500",
        "set ebeam2 6500",
        "set pdlabel lhapdf",
        f"set lhaid {base.PDF_LHAID}",
        f"set ptj {base.PTJMIN:.16g}",
        "set etaj 5.0",
        "set maxjetflavor 5",
        "set ickkw 0",
    ]
    return "\n".join([
        "set automatic_html_opening False",
        f"import model {base.MODEL}",
        "generate p p > xd xd~ j",
        f"output {procdir} -f",
        # IMPORTANT: no -f here.  Plain launch opens the launch-time card-edit
        # stage, so the scripted set commands below are actually applied.
        f"launch {procdir}",
        *runsets,
        *sets,
        "done",
        "",
    ])


base.make_card = make_card_v053


def relabel_outputs() -> None:
    summary_path = base.OUT / "normalization_certification.json"
    if summary_path.exists():
        d = json.loads(summary_path.read_text())
        d["version"] = "v0.53_mg5_norm"
        d["v0.52_diagnosis"] = (
            "v0.52 used `launch ... -f`; the forced launch completed before the "
            "subsequent scripted set commands could edit run_card/param_card."
        )
        d["launch_mode"] = "plain launch followed by scripted set commands"
        summary_path.write_text(json.dumps(d, indent=2) + "\n")
        if base.PASS_MARKER.exists():
            base.PASS_MARKER.write_text(json.dumps(d, indent=2) + "\n")


if __name__ == "__main__":
    rc = base.main()
    relabel_outputs()
    sys.exit(rc)
