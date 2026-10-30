#!/usr/bin/env python3
"""v0.59 deterministic MG5 normalization certification.

v0.58 reached the modern MG5 run card successfully but found that this card
format has no `ickkw` entry.  Since this certification uses a single hard-jet
matrix element with no MLM/CKKW merging, forcing a nonexistent matching key is
unnecessary.

v0.59 changes only that compatibility point: all required run-card settings
remain strict, but `ickkw` is not written.  All parameter-card, hash,
MG5-vs-LHE, monotonicity, and sigma/BR tests are unchanged.  No detector
simulation is run here.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = ROOT / "scripts" / "pkr_certify_mg5_normalization_v0_58.py"
spec = importlib.util.spec_from_file_location("pkr_norm_v058", BASE_PATH)
if spec is None or spec.loader is None:
    raise RuntimeError(f"cannot import {BASE_PATH}")
v058 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v058)

# v058.base is the v0.55 deterministic implementation; cfg is the v0.52
# configuration module holding NEVENTS/PDF_LHAID/PTJMIN.
base = v058.base
cfg = v058.cfg

base.WORK = ROOT / "cache" / "v0.59_mg5_norm"
base.OUT = ROOT / "results" / "v0.59_mg5_norm"
base.PASS_MARKER = base.WORK / "NORMALIZATION_CERTIFIED.json"


def configure_run_card_v059(path: Path, seed: int) -> None:
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
    }
    for k, v in vals.items():
        base.set_run_card_value(path, k, v)


base.configure_run_card = configure_run_card_v059

if __name__ == "__main__":
    sys.exit(base.main())
