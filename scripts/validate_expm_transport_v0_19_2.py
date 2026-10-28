#!/usr/bin/env python3
"""v0.19.2 validation gate for sparse-exponential chirality transport.

Changes from v0.19.1:
- keeps the strict 1.5-km reference; does NOT loosen the 0.5% tolerance;
- chooses the coarsest production PREM step among 12, 6, 3 km that agrees
  directly with the 1.5-km reference within both pointwise and integrated L1 gates;
- requires 3 km vs 1.5 km to pass as an independent fine-grid convergence gate;
- expands Phase A to 12 worst-case vertical configurations: both external-state
  modes for Cases 5, 6, 7 and resolved M=0.1,1,3 TeV;
- checkpoints every completed configuration atomically and uses dynamic scheduling.

No detector folding and no physics interpretation are performed here.
"""
from __future__ import annotations

import argparse, json, multiprocessing as mp, os, time
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
OUT=ROOT/"results"/"v0.19.2"
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
        return float(np.exp(np.log(y[i])+a*(np.log(y[j])-np.log(y[i]))))
    a=(tgev-E[i])/(E[j]-E[i])
    return float(y[i]+a*(y[j]-y[i]))


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
    tmp=path.with_suffix(path.suffix+".tmp")
    tmp.write_text(json.dumps(obj,indent=2))
    tmp.replace(path)


def do_old_compare(task):
    _,kind,param,mu,mode=task; t=time.perf_counter()
    E,p,new=run(kind,param,mu,mode,"expm",3.0)
    E2,_,old=run(kind,param,mu,mode,"old",3.0,0.005)
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
    _,_,a15=run(kind,param,mu,mode,"expm",1.5)

    p1215,m1215,l1215,_,_=compare(E,p,a12,a15)
    p615,m615,l615,_,_=compare(E,p,a6,a15)
    p315,m315,l315,mn,fin=compare(E,p,a3,a15)

    direct={
        "12":dict(max_pointwise=m1215,max_L1=l1215,pass_gate=bool(m1215<POINT_TOL and l1215<L1_TOL)),
        "6":dict(max_pointwise=m615,max_L1=l615,pass_gate=bool(m615<POINT_TOL and l615<L1_TOL)),
        "3":dict(max_pointwise=m315,max_L1=l315,pass_gate=bool(m315<POINT_TOL and l315<L1_TOL)),
    }
    chosen=None
    for step in (12.0,6.0,3.0):
        if direct[str(int(step))]["pass_gate"]:
            chosen=step; break

    fine_ok=direct["3"]["pass_gate"]
    passed=bool(chosen is not None and fine_ok and fin and mn>=-1e-14)
    return dict(test="adaptive_slice_convergence_v0.19.2",kind=kind,param=param,cosz_abs=mu,mode=mode,
                pointwise_12_vs_1p5=p1215,pointwise_6_vs_1p5=p615,pointwise_3_vs_1p5=p315,
                max_pointwise_12_vs_1p5=m1215,max_pointwise_6_vs_1p5=m615,max_pointwise_3_vs_1p5=m315,
                max_L1_12_vs_1p5=l1215,max_L1_6_vs_1p5=l615,max_L1_3_vs_1p5=l315,
                chosen_production_step_km=chosen,direct_gate_summary=direct,
                min_flux=mn,finite=fin,seconds=time.perf_counter()-t,
                gate_definition="Choose coarsest of 12/6/3 km agreeing directly with 1.5 km at <0.5% pointwise and L1; additionally 3 km vs 1.5 km must pass",
                PASS=passed)


def worker(task):
    return task, (do_slice(task) if task[0]=="slice" else do_old_compare(task))


def run_phase(name,tasks,nproc):
    done=[]; pending=[]
    for t in tasks:
        cp=CK/(key(t)+".json")
        if cp.exists():
            r=json.loads(cp.read_text())
            print(f"RESUME {name} {t} PASS={r.get('PASS')}",flush=True)
            done.append((t,r))
            if not r.get("PASS",False):
                return done,False
        else:
            pending.append(t)
    if not pending:
        return done,True

    nworkers=min(nproc,len(pending))
    print(f"PHASE {name}: {len(pending)} pending tasks, workers={nworkers}",flush=True)
    ctx=mp.get_context("spawn")
    pool=ctx.Pool(processes=nworkers)
    ok=True
    try:
        for t,r in pool.imap_unordered(worker,pending,chunksize=1):
            cp=CK/(key(t)+".json"); atomic_json(cp,r)
            metric=r.get("max_pointwise_3_vs_1p5",r.get("max_pointwise_error"))
            chosen=r.get("chosen_production_step_km","-")
            print(f"DONE {name} {t} seconds={r['seconds']:.3f} PASS={r['PASS']} fine={metric:.4g} chosen_step={chosen}",flush=True)
            done.append((t,r))
            if not r.get("PASS",False):
                ok=False; pool.terminate(); break
        if ok: pool.close()
    finally:
        pool.join()
    return done,ok


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--nproc",type=int,default=int(os.environ.get("NUTAU_NPROC","12"))); a=ap.parse_args()
    OUT.mkdir(parents=True,exist_ok=True); CK.mkdir(parents=True,exist_ok=True)
    l5=float(lambda_4pi_tev(C5_BENCHMARK)); l7=float(lambda_4pi_tev(C7_BENCHMARK))

    slice_tasks=[
        ("slice","case5",l5,1.0,"dirac_sink"),("slice","case5",l5,1.0,"majorana_taumu"),
        ("slice","case6",2.44,1.0,"dirac_sink"),("slice","case6",2.44,1.0,"majorana_flip"),
        ("slice","case7",l7,1.0,"dirac_sink"),("slice","case7",l7,1.0,"majorana_flip"),
        ("slice","resolved",0.1,1.0,"dirac_sink"),("slice","resolved",0.1,1.0,"majorana_flip"),
        ("slice","resolved",1.0,1.0,"dirac_sink"),("slice","resolved",1.0,1.0,"majorana_flip"),
        ("slice","resolved",3.0,1.0,"dirac_sink"),("slice","resolved",3.0,1.0,"majorana_flip"),
    ]
    old_tasks=[
        ("old","case5",l5,0.1,"dirac_sink"),("old","case5",l5,0.1,"majorana_taumu"),
        ("old","case6",1.5,0.1,"dirac_sink"),("old","case6",1.5,0.1,"majorana_flip"),
        ("old","case7",3.0,0.1,"dirac_sink"),("old","case7",3.0,0.1,"majorana_flip"),
        ("old","resolved",0.1,0.1,"dirac_sink"),("old","resolved",0.1,0.1,"majorana_flip"),
    ]

    print("v0.19.2: phase A adaptive direct fine-grid convergence",flush=True)
    ra,oka=run_phase("slice",slice_tasks,a.nproc); results=[r for _,r in ra]
    if oka:
        print("v0.19.2: phase B old-vs-expm tractable comparisons",flush=True)
        rb,okb=run_phase("old",old_tasks,a.nproc); results.extend(r for _,r in rb)
    else:
        okb=False

    R=pd.DataFrame([{k:v for k,v in r.items() if not isinstance(v,(list,dict))} for r in results])
    R.to_csv(OUT/"expm_validation_summary_v0_19_2.csv",index=False)
    steps={f"{r.get('kind')}:{r.get('param')}:{r.get('mode')}":r.get("chosen_production_step_km")
           for r in results if r.get("test","").startswith("adaptive_slice")}
    n_total=len(slice_tasks)+len(old_tasks); all_done=len(results)==n_total
    allpass=all_done and oka and okb and all(bool(r.get("PASS",False)) for r in results)
    status=dict(version="v0.19.2-candidate",all_tasks_completed=all_done,all_gates_PASS=allpass,
                n_completed=len(results),n_total=n_total,n_failed=sum(not bool(r.get("PASS",False)) for r in results),
                pointwise_tolerance=POINT_TOL,integrated_L1_tolerance=L1_TOL,
                production_step_selection="coarsest of 12/6/3 km passing direct comparison to 1.5 km; 3 km must also pass",
                chosen_steps_km=steps,checkpointed=True,nproc=a.nproc,
                note="Candidate only. Production remains forbidden until all gates pass.")
    atomic_json(OUT/"expm_validation_status_v0_19_2.json",status)
    print("STATUS",json.dumps(status,indent=2),flush=True)
    if not allpass: raise SystemExit("v0.19.2 candidate validation NOT certified")

if __name__=="__main__": main()
