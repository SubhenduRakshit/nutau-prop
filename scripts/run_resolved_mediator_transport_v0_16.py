#!/usr/bin/env python3
"""Earth-transport screen for v0.16 resolved derivative-Z' tables.

Uses the same chirality brackets as the Case-6 robustness audit:
  * dirac_sink: NP daughter is RH/sterile to SM;
  * majorana_flip: nu_tau <-> antinu_tau with NP energy degradation.

This stage asks only whether propagator softening improves the surviving
energy/zenith structure.  Detector folding is deferred until a viable shape
appears.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

from nutau_prop.chirality_transport import TauHelicityFlipTransport

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results"/"v0.16"
TABLES=OUT/"tables"
UHE=ROOT/"data"/"uhe_official"


def tag(x):
    return (f"{float(x):.6g}").replace(".","p")


def load(mass,scale,species):
    z=np.load(TABLES/f"resolved_M{tag(mass)}_S{tag(scale)}_{species}.npz")
    return {k:z[k] for k in z.files}


def incident(E):
    q=pd.read_json(UHE/"flux.json")
    ef=q["Energy [GeV]"].to_numpy(float)
    sed=q["log10(E^2 F(E) [GeV.cm^-2.s^-1.sr^-1])"].to_numpy(float)
    ss=interp1d(ef,sed,bounds_error=False,fill_value=-np.inf)(E)
    p=E*(10.0**ss/E**2)
    m=np.nanmax(p)
    return p/m if m>0 else p


def loginterp(q,target,col):
    x=q.E_PeV.to_numpy(float); y=q[col].to_numpy(float)
    good=np.isfinite(x)&np.isfinite(y)&(x>0)&(y>0)
    if not np.any(good) or target<x[good].min() or target>x[good].max():
        return np.nan
    return float(10**np.interp(np.log10(target),np.log10(x[good]),np.log10(y[good])))


def run_one(mass,scale,mu,mode):
    tn=load(mass,scale,"nu"); ta=load(mass,scale,"antinu")
    E=tn["energies"]; p0=incident(E)
    tr=TauHelicityFlipTransport(tn,ta)
    sm=tr.propagate(mu,p0,p0,mode=mode,use_np=False)
    bs=tr.propagate(mu,p0,p0,mode=mode,use_np=True)
    smtau=sm["nu_tau"]+sm["antinu_tau"]
    bstau=bs["nu_tau"]+bs["antinu_tau"]
    smsec=sm["nu_mu_secondary"]+sm["antinu_mu_secondary"]
    bssec=bs["nu_mu_secondary"]+bs["antinu_mu_secondary"]
    out=[]
    for i,e in enumerate(E):
        out.append(dict(M_Zprime_TeV=mass,lowq_scale=scale,cosz_abs=mu,mode=mode,
                        E_GeV=e,E_PeV=e/1e6,
                        tau_ratio=bstau[i]/max(smtau[i],1e-300),
                        secondary_mu_ratio=(bssec[i]/smsec[i]) if smsec[i]>0 else np.nan))
    return out


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--masses-tev",nargs="+",type=float,default=[0.1,0.2,0.5,1.0,2.0,3.0])
    ap.add_argument("--lowq-scale",type=float,default=1.0)
    ap.add_argument("--cosz",nargs="+",type=float,default=[0.2,0.5,0.8,1.0])
    a=ap.parse_args()
    rows=[]
    for mass in a.masses_tev:
        for mu in a.cosz:
            for mode in ("dirac_sink","majorana_flip"):
                rows.extend(run_one(mass,a.lowq_scale,mu,mode))
    D=pd.DataFrame(rows)
    D.to_csv(OUT/"resolved_mediator_transport_v0_16.csv",index=False)
    summaries=[]
    for keys,q in D.groupby(["M_Zprime_TeV","lowq_scale","cosz_abs","mode"]):
        mass,scale,mu,mode=keys; q=q.sort_values("E_PeV")
        s=dict(M_Zprime_TeV=mass,lowq_scale=scale,cosz_abs=mu,mode=mode)
        for e in (1.,10.,100.,220.,300.,1000.,10000.):
            s[f"tau_ratio_{e:g}PeV"]=loginterp(q,e,"tau_ratio")
            s[f"secondary_mu_ratio_{e:g}PeV"]=loginterp(q,e,"secondary_mu_ratio")
        summaries.append(s)
    S=pd.DataFrame(summaries).sort_values(["M_Zprime_TeV","mode","cosz_abs"])
    S.to_csv(OUT/"resolved_mediator_transport_summary_v0_16.csv",index=False)
    show=S[S.cosz_abs.isin([0.5,1.0])]
    print("\nResolved-mediator transport gate:")
    print(show.to_string(index=False,float_format=lambda x:f"{x:.4g}"))
    status={
      "version":"v0.16",
      "modes":["dirac_sink","majorana_flip"],
      "decision":"Look for mediator masses that retain a measurable 10-300 PeV distortion while making the 1-10 EeV suppression materially less catastrophic than contact Case 6.",
      "next_if_positive":"joint detector fold and current mediator-coupling constraint audit",
      "next_if_negative":"vary low-Q normalization; if no window, derivative-Zprime parent does not solve the structural opacity problem."
    }
    (OUT/"resolved_mediator_transport_status_v0_16.json").write_text(json.dumps(status,indent=2))
    print("WROTE",OUT)

if __name__=="__main__": main()
