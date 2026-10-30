#!/usr/bin/env python3
"""v0.57 deterministic MG5 normalization certification.

v0.56 reached the generated run card and showed that this MG5 version uses
beam-specific PDF labels (`pdlabel1`, `pdlabel2`) rather than the historical
single `pdlabel` key.  v0.57 changes only that compatibility point.

No detector simulation is run here.
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
base = v056.base

base.WORK = ROOT / "cache" / "v0.57_mg5_norm"
base.OUT = ROOT / "results" / "v0.57_mg5_norm"
base.PASS_MARKER = base.WORK / "NORMALIZATION_CERTIFIED.json"


def configure_run_card_v057(path: Path, seed: int) -> None:
    vals = {
        "nevents": str(base.NEVENTS),
        "iseed": str(seed),
        "ebeam1": "6500",
        "ebeam2": "6500",
        "pdlabel1": "lhapdf",
        "pdlabel2": "lhapdf",
        "lhaid": str(base.PDF_LHAID),
        "ptj": f"{base.PTJMIN:.16g}",
        "etaj": "5.0",
        "maxjetflavor": "5",
        "ickkw": "0",
    }
    for k, v in vals.items():
        base.set_run_card_value(path, k, v)


base.configure_run_card = configure_run_card_v057

if __name__ == "__main__":
    sys.exit(base.main())
