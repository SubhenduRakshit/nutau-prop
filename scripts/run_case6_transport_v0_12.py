#!/usr/bin/env python3
"""v0.12 full Case-6 Earth transport using the v0.11 dσ/dy tables.

For each cutoff, zenith chord and incident spectrum this script propagates
nu_tau and antinu_tau separately through PREM, including:
  * SM CC tau production;
  * SM+Case6 NC energy degradation with the actual combined y kernel;
  * tau energy loss and tau -> nu_tau regeneration;
  * polarized secondary nu_mu/antinu_mu production;
  * SM-only propagation of the secondary muon-flavor channel.

The matched SM baseline uses the same CT10 grid and the stored p_nc_sm kernel.
"""
from __future__ import annotations

import argparse
import json
import math
import multiprocessing as mp
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

from nutau_prop.transport_case6 import TauTransportFlavorAware
from nutau_prop.util import powerlaw_per_logE

ROOT=Path(__file__).resolve().parents[1]
TABLES=ROOT/"results"/"v0.11"/"tables"
OUT=ROOT/"results"/"v0.12"
UHE=ROOT/"data"/"uhe_official"
_CACHE={}


def tag(lam): return str(lam).replace(".","p")


def load_table(lam,species):
    key=(float(lam),species,"raw")
    if key not in _CACHE:
        z=np.load(TABLES/f"case6_L{tag(lam)}_{species}.npz")
        _CACHE[key]={k:z[k] for k in z.files}
    return _CACHE[key]


def sm_view(t):
    return dict(energies=t["energies"],y=t["y"],dy=t["dy"],
                sigma_cc=t["sigma_cc"],sigma_nc=t["sigma_nc_sm"],
                p_cc=t["p_cc"],p_nc=t["p_nc_sm"])


def bsm_view(t):
    return dict(energies=t["energies"],y=t["y"],dy=t["dy"],
                sigma_cc=t["sigma_cc"],sigma_nc=t["sigma_nc"],
                p_cc=t["p_cc"],p_nc=t["p_nc"])


def transport(lam,species,bsm):
    key=(float(lam),species,"bsm" if bsm else "sm")
    if key in _CACHE: return _CACHE[key]
    t=load_table(lam,species); base=sm_view(t)
    primary=bsm_view(t) if bsm else base
    P=-1.0 if species=="nu" else +1.0
    tr=TauTransportFlavorAware(primary,base,beta_tau=0.85e-6,
                               tau_polarization=P,secondary_branching=0.18)
    _CACHE[key]=tr
    return tr


def incident(E,name):
    if name=="E^-2": return powerlaw_per_logE(E,2.0)
    if name=="E^-2.46": return powerlaw_per_logE(E,2.46)
    if name=="Ehlert":
        q=pd.read_json(UHE/"flux.json")
        ef=q["Energy [GeV]"].to_numpy(float)
        sed=q["log10(E^2 F(E) [GeV.cm^-2.s^-1.sr^-1])"].to_numpy(float)
        ss=interp1d(ef,sed,bounds_error=False,fill_value=-np.inf)(E)
        p=E*(10.0**ss/E**2)
        m=np.nanmax(p)
        return p/m if m>0 else p
    raise ValueError(name)


def run_one(task):
    lam,mu,spec=task
    tn=load_table(lam,"nu"); E=tn["energies"]; p0=incident(E,spec)
    results={}
    for species in ("nu","antinu"):
        sm=transport(lam,species,False).propagate(mu,p0,max_step_km=12.0,max_p=0.02,track_secondary_numu=True)
        bs=transport(lam,species,True).propagate(mu,p0,max_step_km=12.0,max_p=0.02,track_secondary_numu=True)
        results[species]=(sm,bs)
    smtau=results["nu"][0]["nu"]+results["antinu"][0]["nu"]
    bstau=results["nu"][1]["nu"]+results["antinu"][1]["nu"]
    smmu=results["nu"][0]["numu"]+results["antinu"][0]["numu"]
    bsmu=results["nu"][1]["numu"]+results["antinu"][1]["numu"]
    rows=[]
    for i,e in enumerate(E):
        rows.append(dict(Lambda_TeV=lam,cosz_abs=mu,spectrum=spec,E_GeV=e,E_PeV=e/1e6,
                         incident_per_logE=2.0*p0[i],
                         SM_nutau=smtau[i],BSM_nutau=bstau[i],
                         nutau_ratio=bstau[i]/max(smtau[i],1e-300),
                         SM_secondary_numu=smmu[i],BSM_secondary_numu=bsmu[i],
                         secondary_numu_ratio=bsmu[i]/max(smmu[i],1e-300) if smmu[i]>0 else np.nan))
    return rows


def loginterp_target(q,target,col):
    x=q.E_PeV.to_numpy(float); y=q[col].to_numpy(float)
    good=np.isfinite(x)&np.isfinite(y)&(x>0)&(y>0)
    if target<x[good].min() or target>x[good].max(): return np.nan
    return float(10**np.interp(np.log10(target),np.log10(x[good]),np.log10(y[good])))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--workers",type=int,default=int(os.environ.get("NUTAU_NPROC","14")))
    ap.add_argument("--lambdas",nargs="+",type=float,default=[1.5,2.0,2.44,3.0])
    ap.add_argument("--cosz",nargs="+",type=float,default=[0.2,0.5,0.8,1.0])
    ap.add_argument("--spectra",nargs="+",default=["E^-2","E^-2.46","Ehlert"])
    a=ap.parse_args(); OUT.mkdir(parents=True,exist_ok=True)
    for lam in a.lambdas:
        for species in ("nu","antinu"):
            p=TABLES/f"case6_L{tag(lam)}_{species}.npz"
            if not p.exists(): raise SystemExit(f"Missing {p}; run v0.11 table builder first")
    tasks=[(lam,mu,spec) for lam in a.lambdas for mu in a.cosz for spec in a.spectra]
    ctx=mp.get_context("spawn")
    with ctx.Pool(processes=max(1,a.workers)) as pool:
        chunks=list(pool.imap_unordered(run_one,tasks,chunksize=1))
    D=pd.DataFrame([r for chunk in chunks for r in chunk])
    D.to_csv(OUT/"case6_earth_transport_v0_12.csv",index=False)

    summaries=[]
    for (lam,mu,spec),q in D.groupby(["Lambda_TeV","cosz_abs","spectrum"]):
        q=q.sort_values("E_PeV")
        physical=(q.E_PeV>=1.0)&(q.E_PeV<=300.0)&(q.SM_nutau>1e-12*np.nanmax(q.SM_nutau))
        dev=np.abs(q.loc[physical,"nutau_ratio"]-1.0)
        s=dict(Lambda_TeV=lam,cosz_abs=mu,spectrum=spec,
               max_abs_nutau_fractional_change_1_300PeV=float(dev.max()) if len(dev) else np.nan)
        for e in (1.0,10.0,100.0,220.0,300.0):
            s[f"nutau_ratio_{e:g}PeV"]=loginterp_target(q,e,"nutau_ratio")
            s[f"secondary_numu_ratio_{e:g}PeV"]=loginterp_target(q,e,"secondary_numu_ratio")
        summaries.append(s)
    S=pd.DataFrame(summaries).sort_values(["spectrum","Lambda_TeV","cosz_abs"])
    S.to_csv(OUT/"case6_earth_transport_summary_v0_12.csv",index=False)
    gate={
      "version":"v0.12",
      "physics":"full differential Case-6 NC energy degradation + tau regeneration + SM-only secondary muon flavor",
      "historical_transport_core":"PREM/tau-decay machinery inherited from the v0.3 RR/DRS-validated implementation",
      "submission_caution":"constant-beta tau energy loss and public detector responses remain approximations; inspect numerical output before any submission claim",
    }
    (OUT/"transport_status_v0_12.json").write_text(json.dumps(gate,indent=2))
    show=S[(S.spectrum=="Ehlert") & (S.cosz_abs.isin([0.5,1.0]))]
    print(show.to_string(index=False,float_format=lambda x:f"{x:.4g}"))
    print("WROTE",OUT)


if __name__=="__main__": main()
