#!/usr/bin/env python3
"""v0.15 transport screen for momentum-growing Cases 5 and 7.

Case 5 is bracketed by:
  * Dirac sink (chirality-flipped RH daughter is sterile to SM),
  * Majorana tau<->mu flavor/helicity transition, with reverse transition.
Case 7 is bracketed by:
  * Dirac sink,
  * Majorana nu_tau<->antinu_tau helicity flip.

The matched SM baseline uses the same tables with NP disabled.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

from nutau_prop.chirality_transport import TauHelicityFlipTransport, TauMuonTransitionTransport
from nutau_prop.momentum_operators import C5_BENCHMARK,C7_BENCHMARK,lambda_4pi_tev

ROOT=Path(__file__).resolve().parents[1]
TABLES=ROOT/"results"/"v0.15"/"tables"
OUT=ROOT/"results"/"v0.15"
UHE=ROOT/"data"/"uhe_official"


def tag(x): return (f"{float(x):.6g}").replace(".","p")

def load(op,lam,species):
    z=np.load(TABLES/f"case{op}_L{tag(lam)}_{species}.npz")
    return {k:z[k] for k in z.files}

def incident(E):
    q=pd.read_json(UHE/"flux.json")
    ef=q["Energy [GeV]"].to_numpy(float)
    sed=q["log10(E^2 F(E) [GeV.cm^-2.s^-1.sr^-1])"].to_numpy(float)
    ss=interp1d(ef,sed,bounds_error=False,fill_value=-np.inf)(E)
    p=E*(10.0**ss/E**2)
    m=np.nanmax(p)
    return p/m if m>0 else p

def li(E,q,target,col):
    x=q.E_PeV.to_numpy(float); y=q[col].to_numpy(float)
    good=np.isfinite(x)&np.isfinite(y)&(x>0)&(y>0)
    if not np.any(good) or target<x[good].min() or target>x[good].max(): return np.nan
    return float(10**np.interp(np.log10(target),np.log10(x[good]),np.log10(y[good])))

def run_case5(lam,mu,mode):
    tn=load(5,lam,"nu"); ta=load(5,lam,"antinu"); E=tn["energies"]; p0=incident(E)
    tr=TauMuonTransitionTransport(tn,ta)
    sm=tr.propagate(mu,p0,p0,p0,p0,mode=mode,use_np=False)
    bs=tr.propagate(mu,p0,p0,p0,p0,mode=mode,use_np=True)
    smtau=sm["nu_tau"]+sm["antinu_tau"]; bstau=bs["nu_tau"]+bs["antinu_tau"]
    smmu=sm["nu_mu"]+sm["antinu_mu"]; bsmu=bs["nu_mu"]+bs["antinu_mu"]
    out=[]
    for i,e in enumerate(E):
        out.append(dict(operator=5,Lambda_UV_TeV=lam,cosz_abs=mu,mode=mode,E_GeV=e,E_PeV=e/1e6,
                        tau_ratio=bstau[i]/max(smtau[i],1e-300),
                        mu_ratio=bsmu[i]/max(smmu[i],1e-300),
                        active_pair_ratio=(bstau[i]+bsmu[i])/max(smtau[i]+smmu[i],1e-300)))
    return out

def run_case7(lam,mu,mode):
    tn=load(7,lam,"nu"); ta=load(7,lam,"antinu"); E=tn["energies"]; p0=incident(E)
    tr=TauHelicityFlipTransport(tn,ta)
    sm=tr.propagate(mu,p0,p0,mode=mode,use_np=False)
    bs=tr.propagate(mu,p0,p0,mode=mode,use_np=True)
    smtau=sm["nu_tau"]+sm["antinu_tau"]; bstau=bs["nu_tau"]+bs["antinu_tau"]
    smsec=sm["nu_mu_secondary"]+sm["antinu_mu_secondary"]
    bssec=bs["nu_mu_secondary"]+bs["antinu_mu_secondary"]
    out=[]
    for i,e in enumerate(E):
        out.append(dict(operator=7,Lambda_UV_TeV=lam,cosz_abs=mu,mode=mode,E_GeV=e,E_PeV=e/1e6,
                        tau_ratio=bstau[i]/max(smtau[i],1e-300),
                        secondary_mu_ratio=bssec[i]/max(smsec[i],1e-300) if smsec[i]>0 else np.nan))
    return out

def main():
    l5=[3.0,lambda_4pi_tev(C5_BENCHMARK)]
    l7=[3.0,lambda_4pi_tev(C7_BENCHMARK)]
    mus=[0.2,0.5,0.8,1.0]
    rows=[]
    for lam in l5:
        for mu in mus:
            for mode in ("dirac_sink","majorana_taumu"):
                rows.extend(run_case5(lam,mu,mode))
    for lam in l7:
        for mu in mus:
            for mode in ("dirac_sink","majorana_flip"):
                rows.extend(run_case7(lam,mu,mode))
    D=pd.DataFrame(rows)
    D.to_csv(OUT/"momentum_transport_v0_15.csv",index=False)
    summaries=[]
    for keys,q in D.groupby(["operator","Lambda_UV_TeV","cosz_abs","mode"]):
        op,lam,mu,mode=keys; q=q.sort_values("E_PeV")
        s=dict(operator=op,Lambda_UV_TeV=lam,cosz_abs=mu,mode=mode)
        for e in (1.,10.,100.,220.,300.):
            s[f"tau_ratio_{e:g}PeV"]=li(None,q,e,"tau_ratio")
            if op==5:
                s[f"mu_ratio_{e:g}PeV"]=li(None,q,e,"mu_ratio")
                s[f"active_pair_ratio_{e:g}PeV"]=li(None,q,e,"active_pair_ratio")
            else:
                s[f"secondary_mu_ratio_{e:g}PeV"]=li(None,q,e,"secondary_mu_ratio")
        summaries.append(s)
    S=pd.DataFrame(summaries).sort_values(["operator","Lambda_UV_TeV","mode","cosz_abs"])
    S.to_csv(OUT/"momentum_transport_summary_v0_15.csv",index=False)
    show=S[S.cosz_abs.isin([0.5,1.0])]
    print(show.to_string(index=False,float_format=lambda x:f"{x:.4g}"))
    status={
      "version":"v0.15",
      "case5":"Dirac sink and Majorana tau<->mu active-transition brackets",
      "case7":"Dirac sink and Majorana nu_tau<->antinu_tau helicity-flip brackets",
      "important":"Majorana Case-5 transition includes reverse mu<->tau scattering; it is not treated as one-way flavor escape.",
      "decision_rule":"If neither contact operator gives a distinctive surviving propagation/flavor signal at its perturbative endpoint, advance to resolved mediator turnover rather than broadening contact scans."
    }
    (OUT/"momentum_transport_status_v0_15.json").write_text(json.dumps(status,indent=2))
    print("WROTE",OUT)

if __name__=="__main__": main()
