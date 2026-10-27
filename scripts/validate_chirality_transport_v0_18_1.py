#!/usr/bin/env python3
"""v0.18.1 validation of the actual chirality-aware production solvers.

This closes the gap identified in the v0.18 audit: v0.17 validated the finite-beta
primitive and TauTransportFlavorAware, but production uses TauHelicityFlipTransport
and TauMuonTransitionTransport.  Here those exact classes are tested with the
Ehlert spectrum, including pointwise 1--300 PeV convergence.

No PDF integration is performed; existing v0.11/v0.15/v0.16 tables are reused.
"""
from __future__ import annotations

import argparse, json, multiprocessing as mp, os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

from nutau_prop.chirality_transport import TauHelicityFlipTransport, TauMuonTransitionTransport
from nutau_prop.momentum_operators import C5_BENCHMARK, C7_BENCHMARK, lambda_4pi_tev

ROOT = Path(__file__).resolve().parents[1]
T11 = ROOT / "results" / "v0.11" / "tables"
T15 = ROOT / "results" / "v0.15" / "tables"
T16 = ROOT / "results" / "v0.16" / "tables"
UHE = ROOT / "data" / "uhe_official"
OUT = ROOT / "results" / "v0.18.1"
TARGETS_PEV = (1.0, 10.0, 100.0, 220.0, 300.0)


def tag(x): return (f"{float(x):.6g}").replace(".", "p")

def _npz(path):
    if not path.exists(): raise FileNotFoundError(path)
    z=np.load(path); return {k:z[k] for k in z.files}

def load5or7(op,lam,s): return _npz(T15/f"case{op}_L{tag(lam)}_{s}.npz")
def load6(lam,s):
    for p in (T11/f"case6_L{str(float(lam)).replace('.', 'p')}_{s}.npz",
              T11/f"case6_L{tag(lam)}_{s}.npz"):
        if p.exists(): return _npz(p)
    raise FileNotFoundError(f"Case6 {lam} {s}")
def loadr(m,s): return _npz(T16/f"resolved_M{tag(m)}_S1_{s}.npz")


def ehlert_one_species(E):
    q=pd.read_json(UHE/"flux.json")
    ef=q["Energy [GeV]"].to_numpy(float)
    sed=q["log10(E^2 F(E) [GeV.cm^-2.s^-1.sr^-1])"].to_numpy(float)
    f=interp1d(ef,sed,bounds_error=False,fill_value=-np.inf)(E)
    return (10.0**f/E)/6.0


def point_interp(E,y,target_gev):
    E=np.asarray(E,float); y=np.asarray(y,float)
    j=np.searchsorted(E,target_gev)
    if j==0 or j>=len(E): return np.nan
    i=j-1
    if not (np.isfinite(y[i]) and np.isfinite(y[j])): return np.nan
    if y[i]>0 and y[j]>0:
        t=(np.log(target_gev)-np.log(E[i]))/(np.log(E[j])-np.log(E[i]))
        return float(np.exp(np.log(y[i])+t*(np.log(y[j])-np.log(y[i]))))
    t=(target_gev-E[i])/(E[j]-E[i])
    return float(y[i]+t*(y[j]-y[i]))


def err(a,b,scale):
    # Relative where reference is meaningful; otherwise normalized absolute error.
    floor=max(abs(scale)*1e-12,1e-300)
    return abs(a-b)/max(abs(b),floor)


def run_config(kind,param,mu,mode,max_p,max_step):
    if kind=="case5":
        tn=load5or7(5,param,"nu"); ta=load5or7(5,param,"antinu"); E=tn["energies"]
        p=ehlert_one_species(E); tr=TauMuonTransitionTransport(tn,ta)
        o=tr.propagate(mu,p,p,p,p,mode=mode,use_np=True,max_p=max_p,max_step_km=max_step)
        return E,p,dict(nut=o["nu_tau"],anut=o["antinu_tau"],num=o["nu_mu"],anum=o["antinu_mu"])
    if kind in ("case6","case7","resolved"):
        if kind=="case6": tn,ta=load6(param,"nu"),load6(param,"antinu")
        elif kind=="case7": tn,ta=load5or7(7,param,"nu"),load5or7(7,param,"antinu")
        else: tn,ta=loadr(param,"nu"),loadr(param,"antinu")
        E=tn["energies"]; p=ehlert_one_species(E); tr=TauHelicityFlipTransport(tn,ta)
        o=tr.propagate(mu,p,p,mode=mode,use_np=True,max_p=max_p,max_step_km=max_step,track_secondary_numu=True)
        return E,p,dict(nut=o["nu_tau"],anut=o["antinu_tau"],num=o["nu_mu_secondary"],anum=o["antinu_mu_secondary"])
    raise ValueError(kind)


def _one(task):
    kind,param,mu,mode=task
    # Production resolution vs a jointly refined reference.
    E,p,a=run_config(kind,param,mu,mode,0.02,12.0)
    E2,p2,b=run_config(kind,param,mu,mode,0.01,6.0)
    if not np.allclose(E,E2): raise RuntimeError("energy grids differ")
    rows=[]; max_point=0.0
    for ep in TARGETS_PEV:
        tg=ep*1e6
        inc=point_interp(E,p,tg)
        for sp in ("nut","anut","num","anum"):
            x=point_interp(E,a[sp],tg); y=point_interp(E,b[sp],tg)
            if not (np.isfinite(x) and np.isfinite(y) and np.isfinite(inc)): continue
            er=err(x,y,inc); max_point=max(max_point,er)
            rows.append(dict(kind=kind,param=float(param),cosz_abs=float(mu),mode=mode,
                             E_PeV=ep,species=sp,production=x,reference=y,
                             pointwise_error=er))
    arrays=list(a.values())+list(b.values())
    finite=all(np.all(np.isfinite(x)) for x in arrays)
    minflux=min(float(np.min(x)) for x in arrays)
    # Integrated L1 on bins with meaningful incident support and reference flux.
    integ=[]
    for sp in ("nut","anut","num","anum"):
        ref=b[sp]; prod=a[sp]
        mask=(p>1e-12*np.max(p)) & (ref>1e-12*max(np.max(ref),1e-300))
        val=0.0 if not np.any(mask) else float(np.sum(np.abs(prod[mask]-ref[mask]))/max(np.sum(np.abs(ref[mask])),1e-300))
        integ.append(val)
    return rows,dict(kind=kind,param=float(param),cosz_abs=float(mu),mode=mode,
                     max_pointwise_error=max_point,max_integrated_L1=max(integ),
                     min_flux=minflux,finite=finite)


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--workers",type=int,default=int(os.environ.get("NUTAU_NPROC","12")))
    a=ap.parse_args(); OUT.mkdir(parents=True,exist_ok=True)
    l5=lambda_4pi_tev(C5_BENCHMARK); l7=lambda_4pi_tev(C7_BENCHMARK)
    tasks=[]
    # Strong/end-point contact benchmarks, shallow and long chords, both physical brackets.
    for mu in (0.1,0.5,1.0):
        for mode in ("dirac_sink","majorana_taumu"): tasks.append(("case5",l5,mu,mode))
        for mode in ("dirac_sink","majorana_flip"): tasks.append(("case7",l7,mu,mode))
        for mode in ("dirac_sink","majorana_flip"): tasks.append(("case6",2.44,mu,mode))
    # Resolved-mediator extremes plus an intermediate mass.
    for m in (0.1,1.0,3.0):
        for mu in (0.1,1.0):
            for mode in ("dirac_sink","majorana_flip"): tasks.append(("resolved",m,mu,mode))
    ctx=mp.get_context("spawn")
    with ctx.Pool(processes=min(max(1,a.workers),len(tasks))) as pool:
        out=list(pool.imap_unordered(_one,tasks,chunksize=1))
    P=pd.DataFrame([r for rows,_ in out for r in rows])
    S=pd.DataFrame([s for _,s in out]).sort_values(["kind","param","cosz_abs","mode"])
    # 0.5% pointwise and integrated tolerance; positivity/finite exact gate.
    S["PASS_pointwise"]=S.max_pointwise_error<5e-3
    S["PASS_integrated"]=S.max_integrated_L1<5e-3
    S["PASS_positive_finite"]=S.finite & (S.min_flux>=-1e-14)
    S["PASS"]=S.PASS_pointwise & S.PASS_integrated & S.PASS_positive_finite
    P.to_csv(OUT/"chirality_pointwise_convergence_v0_18_1.csv",index=False)
    S.to_csv(OUT/"chirality_validation_summary_v0_18_1.csv",index=False)
    status={"version":"v0.18.1","solver_classes":["TauHelicityFlipTransport","TauMuonTransitionTransport"],
            "spectrum":"Ehlert all-flavor split equally among six active nu/antinu species",
            "production_resolution":{"max_p":0.02,"max_step_km":12.0},
            "reference_resolution":{"max_p":0.01,"max_step_km":6.0},
            "pointwise_targets_PeV":list(TARGETS_PEV),"tolerance":0.005,
            "all_gates_PASS":bool(S.PASS.all()),"n_failed":int((~S.PASS).sum())}
    (OUT/"chirality_validation_status_v0_18_1.json").write_text(json.dumps(status,indent=2))
    print(S.to_string(index=False,float_format=lambda x:f"{x:.5g}")); print("\nSTATUS",json.dumps(status,indent=2))
    if not status["all_gates_PASS"]: raise SystemExit("v0.18.1 chirality validation FAILED")
if __name__=="__main__": main()
