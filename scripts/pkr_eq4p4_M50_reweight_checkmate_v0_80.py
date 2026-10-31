#!/usr/bin/env python3
"""v0.80 — exact same-sample CheckMATE reweight for the least-visible M=50 Eq4.4 point.

No MG5, no Pythia, no new events. Reuse the v0.75 100k HepMC sample at M=50 GeV.
At fixed mediator mass in the narrow-width surrogate, the accepted shape is held
fixed and only the on-shell normalization changes as gq^2*BRinv.  Rescale the
v0.75 signal cross section to the v0.78 minimum-coupling-factor representative
point, rerun CheckMATE/Spey, and classify only with the predeclared 0.8--1.2 band.

This is a surrogate monojet certification of the currently screened M=50 slice,
not a statement about the exact PKR derivative-UFO acceptance or all mediator masses.
"""
from __future__ import annotations

import importlib.util
import json
import os
import re
import shutil
import subprocess
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "v0.80_eq4p4_M50_reweight"
WORK = ROOT / "cache" / "v0.80_eq4p4_M50_reweight"
SUM78 = ROOT / "results" / "v0.78_eq4p4_full_envelope" / "mass_summary.csv"
REF75 = ROOT / "results" / "v0.75_eq4p4_bestpoint_checkmate" / "bestpoint_checkmate.json"
HEPMC75 = ROOT / "cache" / "v0.75_eq4p4_bestpoint_checkmate" / "hepmc" / "eq4p4_derivative__v074_best.hepmc"
ANALYSIS = "cms_2107_13201"

P51 = ROOT / "scripts" / "pkr_run_checkmate_v0_51.py"
s51 = importlib.util.spec_from_file_location("cm51", P51)
if s51 is None or s51.loader is None:
    raise RuntimeError(f"cannot import {P51}")
cm51 = importlib.util.module_from_spec(s51); s51.loader.exec_module(cm51)

P69 = ROOT / "scripts" / "pkr_checkmate_multibin_retry_v0_69.py"
s69 = importlib.util.spec_from_file_location("v069", P69)
if s69 is None or s69.loader is None:
    raise RuntimeError(f"cannot import {P69}")
v069 = importlib.util.module_from_spec(s69); s69.loader.exec_module(v069)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)
    if not SUM78.exists() or not REF75.exists() or not HEPMC75.exists():
        raise RuntimeError("missing v0.78 summary, v0.75 JSON, or v0.75 HepMC")

    s = pd.read_csv(SUM78)
    q = s[s.M_GeV == 50.0]
    if len(q) != 1:
        raise RuntimeError("expected one M=50 row in v0.78 summary")
    row = q.iloc[0]
    gq_rep = float(row.representative_gq)
    br_rep = float(row.representative_BR_inv)
    a_rep = float(row.representative_a)
    P_rep = float(row.representative_monojet_coupling_factor)

    ref = json.loads(REF75.read_text())
    gq_ref = float(ref["benchmark"]["gq"])
    br_ref = float(ref["benchmark"]["BR_invisible"])
    P_ref = gq_ref*gq_ref*br_ref
    xsec_ref = float(ref["generator"]["cross_section_pb"])
    scale = P_rep/P_ref
    xsec_rep = xsec_ref*scale

    for mod in ("spey", "spey_pyhf", "pyhf"):
        __import__(mod)
    v069.patch_optional_jax()

    cmout = WORK / "checkmate" / "M50_rep"
    card = WORK / "checkmate_M50_rep.in"
    card.write_text(f"""[Parameters]\nBackend: numpy\nName: pkr_v080_M50_rep\nAnalyses: {ANALYSIS}\nOutputDirectory: {cmout}\nOutputExists: Overwrite\nSkipParamCheck: True\nEffTab: True\nMultibin: detailed\nExpected: True\nUplim: True\nModel: simple\nInvisiblePIDs: 5000521\nRandomSeed: 8001\nTotalResultFileColumns: analysis, sr, o, b, db, s, ds, s95obs, s95exp, robscons, rexpcons\nBestPerAnalysisResultFileColumns: analysis, sr, o, b, db, s, ds, s95obs, s95exp, robscons, rexpcons\n\n[signal]\nEvents: {HEPMC75}\nXSect: {xsec_rep:.12g} PB\nXSectErr: 10 %%\n""")
    if cmout.exists():
        shutil.rmtree(cmout)
    checkmate = Path(os.environ.get("PKR_CHECKMATE", ROOT/"cache"/"v0.51"/"checkmate-install"/"bin"/"CheckMATE"))
    p = subprocess.run([str(checkmate), str(card)], text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (OUT/"checkmate.log").write_text(p.stdout)

    mu95, rmb, cls = cm51.parse_multibin(cmout)
    rsingle = cm51.parse_singlebin(cmout)
    if mu95 is None:
        vals = re.findall(r"Upper limit:\s*([0-9.eE+-]+)", p.stdout)
        if vals:
            mu95 = float(vals[-1]); rmb = 1.0/mu95 if mu95 > 0 else None
    if rmb is None:
        raise RuntimeError(f"multibin result unavailable; CheckMATE rc={p.returncode}")

    if rmb < 0.8:
        status = "ROBUST_SURVIVAL_THRESHOLD"
    elif rmb > 1.2:
        status = "ROBUST_EXCLUSION_THRESHOLD"
    else:
        status = "MARGINAL_BOUNDARY_UNRESOLVED"

    out = {
        "version": "v0.80_eq4p4_M50_reweight",
        "source": "same 100k v0.75 HepMC; no MG5/Pythia/new events",
        "representative": {"M_GeV": 50.0, "gq": gq_rep, "a_GeV_inv": a_rep, "BR_inv": br_rep,
                           "coupling_factor_gq2BR": P_rep},
        "reference_v075": {"gq": gq_ref, "BR_inv": br_ref, "coupling_factor_gq2BR": P_ref,
                            "cross_section_pb": xsec_ref},
        "normalization_scale": scale,
        "rescaled_cross_section_pb": xsec_rep,
        "checkmate": {"returncode": p.returncode, "mu95": mu95, "multibin_r": rmb,
                      "singlebin_robscons": rsingle, "status": status},
        "interpretation": "Use 0.8<r<1.2 as boundary-unresolved. Do not chase r=1 with more event generation."
    }
    (OUT/"M50_reweight.json").write_text(json.dumps(out, indent=2)+"\n")

    print("=== v0.80 Eq4.4 M50 same-sample reweight ===")
    print(f"P_rep/P_v075 = {scale:.6f}")
    print(f"xsec: {xsec_ref:.6g} -> {xsec_rep:.6g} pb")
    print(f"mu95={mu95}  multibin_r={rmb}  singlebin={rsingle}")
    print(f"status: {status}")
    print("MG5 rerun: no; Pythia rerun: no; new events: no")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
