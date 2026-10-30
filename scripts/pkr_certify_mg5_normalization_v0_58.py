#!/usr/bin/env python3
"""v0.58 deterministic MG5 normalization certification.

v0.57 correctly added beam-specific PDF labels (pdlabel1/pdlabel2) but failed
before touching MG5 because the wrapper looked for inherited configuration
constants directly on the v0.55 module.  In v0.55 those constants live in its
imported v0.52 module as ``base.base``.

v0.58 changes only that namespace lookup.  All normalization checks and MG5
card compatibility logic are otherwise unchanged.  No detector simulation is
run here.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = ROOT / "scripts" / "pkr_certify_mg5_normalization_v0_56.py"
spec = importlib.util.spec_from_file_location("pkr_norm_v056", BASE_PATH)
if spec is None or spec.loader is None:
    raise RuntimeError(f"cannot import {BASE_PATH}")
v056 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v056)

# v056.base is the v0.55 module; v0.55.base is the original v0.52 config module.
base = v056.base
cfg = base.base

base.WORK = ROOT / "cache" / "v0.58_mg5_norm"
base.OUT = ROOT / "results" / "v0.58_mg5_norm"
base.PASS_MARKER = base.WORK / "NORMALIZATION_CERTIFIED.json"


def configure_run_card_v058(path: Path, seed: int) -> None:
    vals = {
        "nevents": str(cfg.NEVENTS),
        "iseed": str(seed),
        "ebeam1": "6500",
        "ebeam2": "6500",
        "pdlabel1": "lhapdf",
        "pdlabel2": "lhapdf",
        "lhaid": str(cfg.PDF_LHAID),
        "ptj": f"{cfg.PTJMIN:.16g}",
        "etaj": "5.0",
        "maxjetflavor": "5",
        "ickkw": "0",
    }
    for k, v in vals.items():
        base.set_run_card_value(path, k, v)


base.configure_run_card = configure_run_card_v058

if __name__ == "__main__":
    sys.exit(base.main())
