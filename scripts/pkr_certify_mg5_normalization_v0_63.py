#!/usr/bin/env python3
"""v0.63 deterministic MG5 normalization certification.

v0.62 fixed the macOS linker path and reached genuine event generation, but
CT18NLO then failed inside MG5's legacy LHAGLUE bridge because that bridge asks
for obsolete AlphaS_Lambda4 metadata which CT18NLO does not provide.

That LHAPDF compatibility problem is orthogonal to the present certification:
we only need to test whether, at fixed production inputs, the generated
pp -> Xd Xd~ j cross section scales with BR_inv as required by the NWA.
Therefore v0.63 uses MG5's built-in nn23lo1 PDF for BOTH beams only for this
normalization gate.  The same PDF is used for every BR point, so the
sigma/BR consistency test is unaffected by this choice.

CT18NLO remains mandatory for the later physical collider benchmark and is NOT
claimed fixed here.  No CheckMATE run occurs.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = ROOT / "scripts" / "pkr_certify_mg5_normalization_v0_60.py"
spec = importlib.util.spec_from_file_location("pkr_norm_v060", BASE_PATH)
if spec is None or spec.loader is None:
    raise RuntimeError(f"cannot import {BASE_PATH}")
v060 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v060)

# v060.base is the deterministic v0.55 implementation.  Its nested base module
# is the original v0.52 configuration holder.
base = v060.base
cfg = base.base

base.WORK = ROOT / "cache" / "v0.63_mg5_norm"
base.OUT = ROOT / "results" / "v0.63_mg5_norm"
base.PASS_MARKER = base.WORK / "NORMALIZATION_CERTIFIED.json"


def configure_run_card_v063(path: Path, seed: int) -> None:
    vals = {
        "nevents": str(cfg.NEVENTS),
        "iseed": str(seed),
        "ebeam1": "6500",
        "ebeam2": "6500",
        # Deliberately use MG5's built-in PDF for this normalization-only gate.
        "pdlabel1": "nn23lo1",
        "pdlabel2": "nn23lo1",
        "ptj": f"{cfg.PTJMIN:.16g}",
        "etaj": "5.0",
        "maxjetflavor": "5",
    }
    for k, v in vals.items():
        base.set_run_card_value(path, k, v)


base.configure_run_card = configure_run_card_v063


def relabel_summary() -> None:
    p = base.OUT / "normalization_certification.json"
    if not p.exists():
        return
    d = json.loads(p.read_text())
    d["version"] = "v0.63_mg5_norm"
    d["pdf_for_certification"] = "MG5 built-in nn23lo1 on both beams"
    d["scope_note"] = (
        "This certifies generator normalization/NWA BR scaling only. "
        "CT18NLO/LHAPDF compatibility remains unresolved and is required "
        "before the physical collider benchmark."
    )
    p.write_text(json.dumps(d, indent=2) + "\n")
    if base.PASS_MARKER.exists():
        base.PASS_MARKER.write_text(json.dumps(d, indent=2) + "\n")


if __name__ == "__main__":
    rc = base.main()
    relabel_summary()
    sys.exit(rc)
