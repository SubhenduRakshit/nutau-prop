#!/usr/bin/env python3
"""Validation gate for the v0.19 sparse-exponential chirality transport.

Goals
-----
1. Compare the expm interaction evolution against the established small-max_p
   solver on deliberately tractable shallow-chord cases.
2. Demonstrate matter-slice convergence of the expm solver on strong/vertical
   cases, including the former multi-hour resolved-M=3 TeV straggler.
3. Check species-resolved pointwise errors at 1, 10, 100, 220, 300 PeV,
   integrated spectral agreement, positivity and finiteness.
4. Checkpoint every completed configuration so interruption never discards
   previous work.

No detector folding and no physics interpretation are performed here.
"""
from __future__ import annotations

import argparse, json, os, time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

from nutau_prop.chirality_transport import TauHelicityFlipTransport, TauMuonTransitionTransport
from nutau_prop.chirality_transport_expm import TauHelicityFlipTransportExpm, TauMuonTransitionTransportExpm
from nutau_prop.momentum_operators import C5_BENCHMARK, C7_BENCHMARK, lambda_4pi_tev

ROOT=Path(__file__).resolve().parents[1]
T11=ROOT/"results"/"v0.11"/"tables"
T15=ROOT/"results"/"v0.15"/"tables"
T16=ROOT/"results"/"v0.16"/"tables"
UHE=ROOT/"data"/"uhe_official"
OUT=ROOT/"results"/"v0.19"
CK=OUT/"checkpoints"
TARGETS=(1.,10.,100.,220.,300.)
POINT_TOL=5e-3
L1_TOL=5e-3


def tag(x): return (f"{float(x):.6g}").replace(".","p")
def npz(p):
    if not p.exists(): raise FileNotFoundError(p)
    z=np.load(p); return {k:z[k] for k in z.files}
def load57(op,lam,s): return npz(T15/f"case{op}_L{tag(lam)}_{s}.npz")
def load6(lam,s):
    for p in (T11/f"case6_L{str(float(lam)).replace('.', 'p')}_{s}.npz",
              T11/f"case6_L{tag(lam)}_{s}.npz"):
        if p.exists(): return npz(p)
    raise FileNotFoundError(f"case6 {lam} {s}")
def loadr(m,s): return npz(T16/f"resolved_M{tag(m)}_S1_{s}.npz")


def ehlert_one(E):
    q=pd.read_json(UHE/"flux.json")
    ef=q["Energy [GeV]"].to_numpy(float)
    sed=q["log10(E^2 F(E) [GeV.cm^-2.s^-1.sr^-1])"].to_numpy(float)
    f=interp1d(ef,sed,bounds_error=False,fill_value=-np.inf)(E)
    return (10.0**f/E)/6.0


def interp_point(E,y,tgev):
    E=np.asarray(E,float); y=np.asarray(y,float); j=int(np.searchsorted(E,tgev))
    if j==0 or j>=len(E): return np.nan
    i=j-1
    if not (np.isfinite(y[i]) and np.isfinite(y[j])): return np.nan
    if y[i]>0 and y[j]>0:
        a=(np.log(tgev)-np.log(E[i]))/(np.log(E[j])-np.log(E[i]))
        return float(np.exp(np.log(y[i])+a*(np.log(y[j])-np.log(y[i])))
    a=(tgev-E[i])/(E[j]-E[i]); return float(y[i]+a*(y[j]-y[i]))


def relerr(a,b,scale):
    return abs(a-b)/max(abs(b),abs(scale)*1e-12,1e-300)


def unpack(kind,o):
    if kind=="case5":
        return dict(nut=o["nu_tau"],anut=o["antinu_tau"],num=o["nu_mu"],anum=o["antinu_mu"])
    return dict(nut=o["nu_tau"],anut=o["antinu_tau"],num=o["nu_mu_secondary"],anum=o["antinu_mu_secondary"])


def load_model(kind,param):
    if kind=="case5": return load57(5,param,"nu"),load57(5,param,"antinu")
    if kind=="case6": return load6(param,"nu"),load6(param,"antinu")
    if kind=="case7": return load57(7,param,"nu"),load57(7,param,"antinu")
    if kind=="resolved": return loadr(param,"nu"),loadr(param,"antinu")
    raise ValueError(kind)


def run(kind,param,mu,mode,solver,max_step,max_p=None):
    tn,ta=load_model(kind,param); E=np.asarray(tn["energies"],float); p=ehlert_one(E)
    if kind=="case5":
        C=TauMuonTransitionTransportExpm if solver=="expm" else TauMuonTransitionTransport
        tr=C(tn,ta)
        kw=dict(mode=mode,use_np=True,max_step_km=max_step)
        if solver=="old": kw["max_p"]=max_p
        o=tr.propagate(mu,p,p,p,p,**kw)
    else:
        C=TauHelicityFlipTransportExpm if solver=="expm" else TauHelicityFlipTransport
        tr=C(tn,ta)
        kw=dict(mode=mode,use_np=True,max_step_km=max_step,track_secondary_numu=True)
        if solver=="old": kw["max_p"]=max_p
        o=tr.propagate(mu,p,p,**kw)
    return E,p,unpack(kind,o)


def compare(E,p,a,b):
    point=[]; mx=0.0
    for ep in TARGETS:
        inc=interp_point(E,p,ep*1e6)
        for sp in ("nut","anut","num","anum"):
            x=interp_point(E,a[sp],ep*1e6); y=interp_point(E,b[sp],ep*1e6)
            if not (np.isfinite(x) and np.isfinite(y) and np.isfinite(inc)): continue
            er=relerr(x,y,inc); mx=max(mx,er)
            point.append(dict(E_PeV=ep,species=sp,value_a=x,value_b=y,relative_error=er))
    l1=[]
    for sp in ("nut","anut","num","anum"):
        ref=b[sp]; val=a[sp]
        mask=(p>1e-12*np.max(p)) & (ref>1e-12*max(np.max(ref),1e-300))
        e=0.0 if not np.any(mask) else float(np.sum(np.abs(val[mask]-ref[mask]))/max(np.sum(np.abs(ref[mask])),1e-300))
        l1.append(e)
    arrays=list(a.values())+list(b.values())
    return point,mx,max(l1),min(float(np.min(x)) for x in arrays),all(np.all(np.isfinite(x)) for x in arrays)


def key(task):
    return "__".join(str(x).replace(".","p") for x in task)


def atomic_json(path,obj):
    tmp=path.with_suffix(path.suffix+".tmp"); tmp.write_text(json.dumps(obj,indent=2)); tmp.replace(path)


def do_old_compare(task):
    _,kind,param,mu,mode=task; t=time.perf_counter()
    E,p,new=run(kind,param,mu,mode,"expm",3.0)
    E2,p2,old=run(kind,param,mu,mode,"old",3.0,0.005)
    if not np.allclose(E,E2): raise RuntimeError("energy grids differ")
    pts,mx,l1,mn,fin=compare(E,p,new,old)
    return dict(test="old_vs_expm",kind=kind,param=param,cosz_abs=mu,mode=mode,
                pointwise=pts,max_pointwise_error=mx,max_integrated_L1=l1,
                min_flux=mn,finite=fin,seconds=time.perf_counter()-t,
                PASS=bool(mx<POINT_TOL and l1<L1_TOL and fin and mn>=-1e-14))


def do_slice(task):
    _,kind,param,mu,mode=task; t=time.perf_counter()
    E,p,a12=run(kind,param,mu,mode,"expm",12.0)
    _,_,a6=run(kind,param,mu,mode,"expm",6.0)
    _,_,a3=run(kind,param,mu,mode,"expm",3.0)
    p126,m126,l126,_,_=compare(E,p,a12,a6)
    p63,m63,l63,mn,fin=compare(E,p,a6,a3)
    # Require the production 12-km result to agree with 6 km within 0.5%,
    # and the refinement not to reveal a larger unresolved error.
    trend=(m63<=1.25*m126+1e-5) or (max(m126,m63)<1e-4)
    return dict(test="slice_convergence",kind=kind,param=param,cosz_abs=mu,mode=mode,
                pointwise_12_vs_6=p126,pointwise_6_vs_3=p63,
                max_pointwise_12_vs_6=m126,max_pointwise_6_vs_3=m63,
                max_L1_12_vs_6=l126,max_L1_6_vs_3=l63,min_flux=mn,finite=fin,
                refinement_trend=trend,seconds=time.perf_counter()-t,
                PASS=bool(m126<POINT_TOL and l126<L1_TOL and m63<POINT_TOL and l63<L1_TOL and trend and fin and mn>=-1e-14))


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--resume",action="store_true",default=True); a=ap.parse_args()
    OUT.mkdir(parents=True,exist_ok=True); CK.mkdir(parents=True,exist_ok=True)
    l5=float(lambda_4pi_tev(C5_BENCHMARK)); l7=float(lambda_4pi_tev(C7_BENCHMARK))

    # Old-solver comparison is intentionally shallow/tractable.
    old_tasks=[
        ("old","case5",l5,0.1,"dirac_sink"),("old","case5",l5,0.1,"majorana_taumu"),
        ("old","case6",1.5,0.1,"dirac_sink"),("old","case6",1.5,0.1,"majorana_flip"),
        ("old","case7",3.0,0.1,"majorana_flip"),("old","resolved",0.1,0.1,"majorana_flip"),
    ]
    # New-solver slice convergence includes the formerly pathological strong vertical case.
    slice_tasks=[
        ("slice","case5",l5,1.0,"majorana_taumu"),
        ("slice","case6",2.44,1.0,"majorana_flip"),
        ("slice","case7",l7,1.0,"majorana_flip"),
        ("slice","resolved",0.1,1.0,"majorana_flip"),
        ("slice","resolved",1.0,1.0,"majorana_flip"),
        ("slice","resolved",3.0,1.0,"dirac_sink"),
        ("slice","resolved",3.0,1.0,"majorana_flip"),
    ]
    # Expensive-first ordering for predictable progress.
    tasks=slice_tasks+old_tasks
    results=[]
    for i,t in enumerate(tasks,1):
        cp=CK/(key(t)+".json")
        if cp.exists():
            r=json.loads(cp.read_text()); print(f"RESUME {i}/{len(tasks)} {t} PASS={r.get('PASS')}",flush=True)
        else:
            print(f"START {i}/{len(tasks)} {t}",flush=True)
            r=do_slice(t) if t[0]=="slice" else do_old_compare(t)
            atomic_json(cp,r)
            print(f"DONE  {i}/{len(tasks)} seconds={r['seconds']:.3f} PASS={r['PASS']} maxpoint="
                  f"{r.get('max_pointwise_12_vs_6',r.get('max_pointwise_error')):.4g}",flush=True)
        results.append(r)
        if not r.get("PASS",False):
            print("FAIL-FAST: checkpoint preserved; stopping before further tests",flush=True)
            break

    R=pd.DataFrame([{k:v for k,v in r.items() if not isinstance(v,list)} for r in results])
    R.to_csv(OUT/"expm_validation_summary_v0_19.csv",index=False)
    all_done=len(results)==len(tasks)
    allpass=all_done and all(bool(r.get("PASS",False)) for r in results)
    status=dict(version="v0.19-candidate",all_tasks_completed=all_done,all_gates_PASS=allpass,
                n_completed=len(results),n_total=len(tasks),n_failed=sum(not bool(r.get("PASS",False)) for r in results),
                pointwise_tolerance=POINT_TOL,integrated_L1_tolerance=L1_TOL,
                note="Candidate only. Production remains forbidden until all gates pass.")
    atomic_json(OUT/"expm_validation_status_v0_19.json",status)
    print("STATUS",json.dumps(status,indent=2),flush=True)
    if not allpass: raise SystemExit("v0.19 candidate validation NOT certified")

if __name__=="__main__": main()
