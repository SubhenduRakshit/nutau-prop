#!/usr/bin/env python3
"""Timed smoke test for the v0.19 exact-interaction candidate.

Runs the previously pathological resolved M=3 TeV, vertical chord with the
Ehlert spectrum.  No certification claim is made.  The purpose is to determine
whether removing max_p subcycling actually removes the multi-hour straggler
before any broad validation is attempted.
"""
from __future__ import annotations
import json, time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

from nutau_prop.chirality_transport_expm import TauHelicityFlipTransportExpm

ROOT=Path(__file__).resolve().parents[1]
T=ROOT/"results"/"v0.16"/"tables"
UHE=ROOT/"data"/"uhe_official"
OUT=ROOT/"results"/"v0.19"

def load(s):
    p=T/f"resolved_M3_S1_{s}.npz"
    if not p.exists(): raise SystemExit(f"Missing {p}")
    z=np.load(p); return {k:z[k] for k in z.files}

def flux(E):
    q=pd.read_json(UHE/"flux.json")
    ef=q["Energy [GeV]"].to_numpy(float)
    sed=q["log10(E^2 F(E) [GeV.cm^-2.s^-1.sr^-1])"].to_numpy(float)
    v=interp1d(ef,sed,bounds_error=False,fill_value=-np.inf)(E)
    return (10.0**v/E)/6.0

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    tn,ta=load("nu"),load("antinu"); E=np.asarray(tn["energies"],float); p=flux(E)
    tr=TauHelicityFlipTransportExpm(tn,ta,beta_tau=0.85e-6)
    rows=[]
    for mode in ("dirac_sink","majorana_flip"):
        print(f"START {mode}",flush=True); t0=time.perf_counter()
        o=tr.propagate(1.0,p,p,mode=mode,use_np=True,max_step_km=12.0,track_secondary_numu=True)
        dt=time.perf_counter()-t0
        arrays=[o["nu_tau"],o["antinu_tau"],o["nu_mu_secondary"],o["antinu_mu_secondary"],o["tau_minus"],o["tau_plus"]]
        finite=all(np.all(np.isfinite(x)) for x in arrays); mn=min(float(np.min(x)) for x in arrays)
        rows.append(dict(mode=mode,seconds=dt,finite=finite,min_flux=mn))
        print(f"DONE {mode} seconds={dt:.3f} finite={finite} min_flux={mn:.6g}",flush=True)
        pd.DataFrame(rows).to_csv(OUT/"expm_smoke_v0_19.csv",index=False)
    status={"version":"v0.19-candidate","case":"resolved M=3 TeV, cosz=1, 12 km PREM slices",
            "purpose":"runtime smoke test only; not validation","results":rows}
    (OUT/"expm_smoke_status_v0_19.json").write_text(json.dumps(status,indent=2))
    print("STATUS",json.dumps(status,indent=2))
if __name__=="__main__": main()
