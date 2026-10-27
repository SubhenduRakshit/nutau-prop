#!/usr/bin/env python3
"""Case-6 external-state robustness audit.

The original v0.12 propagation treated the NP daughter as the same active
nu_tau species.  Dimension-7 chirality-flipping EFT language instead motivates
bracketing that result by:
  * Dirac sink: active -> RH sterile-to-SM daughter;
  * Majorana helicity flip: nu_tau <-> antinu_tau with the same NP y kernel.
This audit reuses the existing v0.11 Case-6 tables; no PDF rebuild is required.
"""
from __future__ import annotations

import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

from nutau_prop.chirality_transport import TauHelicityFlipTransport

ROOT=Path(__file__).resolve().parents[1]
TABLES=ROOT/"results"/"v0.11"/"tables"
OUT=ROOT/"results"/"v0.15"
UHE=ROOT/"data"/"uhe_official"

def tag(lam): return str(lam).replace(".","p")
def load(lam,species):
    z=np.load(TABLES/f"case6_L{tag(lam)}_{species}.npz")
    return {k:z[k] for k in z.files}
def incident(E):
    q=pd.read_json(UHE/"flux.json")
    ef=q["Energy [GeV]"].to_numpy(float); sed=q["log10(E^2 F(E) [GeV.cm^-2.s^-1.sr^-1])"].to_numpy(float)
    ss=interp1d(ef,sed,bounds_error=False,fill_value=-np.inf)(E)
    p=E*(10.0**ss/E**2); m=np.nanmax(p)
    return p/m if m>0 else p
def li(q,e,col):
    x=q.E_PeV.to_numpy(float); y=q[col].to_numpy(float); g=(x>0)&(y>0)&np.isfinite(y)
    return float(10**np.interp(np.log10(e),np.log10(x[g]),np.log10(y[g])))

def main():
    rows=[]
    for lam in (1.5,2.0,2.44,3.0):
        tn=load(lam,"nu"); ta=load(lam,"antinu"); E=tn["energies"]; p0=incident(E)
        tr=TauHelicityFlipTransport(tn,ta)
        for mu in (0.2,0.5,0.8,1.0):
            sm=tr.propagate(mu,p0,p0,mode="majorana_flip",use_np=False)
            smtau=sm["nu_tau"]+sm["antinu_tau"]
            for mode in ("dirac_sink","majorana_flip"):
                bs=tr.propagate(mu,p0,p0,mode=mode,use_np=True)
                bt=bs["nu_tau"]+bs["antinu_tau"]
                for i,e in enumerate(E):
                    rows.append(dict(Lambda_TeV=lam,cosz_abs=mu,mode=mode,E_GeV=e,E_PeV=e/1e6,
                                     tau_ratio=bt[i]/max(smtau[i],1e-300)))
    D=pd.DataFrame(rows); D.to_csv(OUT/"case6_chirality_audit_v0_15.csv",index=False)
    S=[]
    for keys,q in D.groupby(["Lambda_TeV","cosz_abs","mode"]):
        lam,mu,mode=keys; q=q.sort_values("E_PeV"); r=dict(Lambda_TeV=lam,cosz_abs=mu,mode=mode)
        for e in (1.,10.,100.,220.,300.): r[f"tau_ratio_{e:g}PeV"]=li(q,e,"tau_ratio")
        S.append(r)
    S=pd.DataFrame(S).sort_values(["Lambda_TeV","mode","cosz_abs"])
    S.to_csv(OUT/"case6_chirality_summary_v0_15.csv",index=False)
    print(S[S.cosz_abs.isin([0.5,1.0])].to_string(index=False,float_format=lambda x:f"{x:.4g}"))
    (OUT/"case6_chirality_status_v0_15.json").write_text(json.dumps({
      "version":"v0.15","purpose":"external-state robustness bracket for Case 6",
      "original_v0.12":"same-active-daughter NC-like propagation",
      "new_brackets":["Dirac active-to-RH sink","Majorana nu_tau<->antinu_tau helicity flip"],
      "interpretation":"If both brackets are at least as suppressive as v0.12 at UHE, the Case-6 negative conclusion is strengthened."
    },indent=2))
    print("WROTE",OUT)
if __name__=="__main__": main()
