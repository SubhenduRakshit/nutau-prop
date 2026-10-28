#!/usr/bin/env python3
"""v0.20.8 adaptive multicore coupled-transport gate.

Purpose: close the transport prerequisite without multi-hour brute-force runs.
The combined generator is exact within each constant-density slab, so the
remaining spatial error is PREM density discretisation.  Start coarse and
refine only when the 0.5% gate requires it.

Initial six independent full-chord jobs (run concurrently):
  shallow |cos z|=0.1: 48 km tau16, 24 km tau16, 24 km tau32
  vertical |cos z|=1.0: 192 km tau16, 96 km tau16, 192 km tau32

Thus both chords get a spatial check and a full-chord tau16->32 check.  If a
spatial comparison fails 0.5%, only that chord is escalated by another factor
of two (24->12 km shallow or 96->48 km vertical).  No tolerance is loosened.

All jobs checkpoint every 20 PREM segments.  Independent jobs use separate
processes; BLAS/OpenMP threads are kept at one by the launcher.
"""
from __future__ import annotations

import json, math, os, time, multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

from nutau_prop.chirality_transport_mixed_expm import TauHelicityFlipTransportMixedExpm
from nutau_prop.constants import R_EARTH_KM
from nutau_prop.prem import density_g_cm3

ROOT=Path(__file__).resolve().parents[1]
T15=ROOT/'results'/'v0.15'/'tables'
UHE=ROOT/'data'/'uhe_official'
OUT=ROOT/'results'/'v0.20.8'
CK=OUT/'checkpoints'
TARGETS=(1.,10.,100.,220.,300.)
INTERFACES=(1221.5,3480.,5701.,5771.,5971.,6151.,6346.6,6356.,6368.)
MODE='majorana_flip'
GATE=5e-3
CASES={
    'shallow':dict(mu=0.1,base_km=48.0),
    'vertical':dict(mu=1.0,base_km=192.0),
}


def load(spec):
    z=np.load(T15/f'case7_L3_{spec}.npz')
    return {k:z[k] for k in z.files}


def flux(E):
    q=pd.read_json(UHE/'flux.json')
    ef=q['Energy [GeV]'].to_numpy(float)
    sed=q['log10(E^2 F(E) [GeV.cm^-2.s^-1.sr^-1])'].to_numpy(float)
    f=interp1d(ef,sed,bounds_error=False,fill_value=-np.inf)(E)
    return (10.0**f/E)/6.0


def nested_profile(mu,base_km,refine):
    """PREM-interface-aligned profile, exactly nested for refine=1,2,4."""
    mu=float(mu); refine=int(refine); base_km=float(base_km)
    if refine not in (1,2,4): raise ValueError('refine must be 1,2,4')
    L=2*R_EARTH_KM*mu
    impact=R_EARTH_KM*np.sqrt(max(0.,1.-mu*mu)); half=L/2
    cuts=[-half,half]
    for r in INTERFACES:
        if r>impact:
            x=np.sqrt(max(0.,r*r-impact*impact))
            if 1e-12<x<half-1e-12: cuts.extend((-x,x))
    cuts=np.array(sorted(set(round(float(x),12) for x in cuts)),float)
    edges=[cuts[0]]
    for a,b in zip(cuts[:-1],cuts[1:]):
        n=max(1,int(np.ceil((b-a)/base_km)))
        macro=np.linspace(a,b,n+1)
        for u,v in zip(macro[:-1],macro[1:]):
            sub=np.linspace(u,v,refine+1)
            edges.extend(sub[1:].tolist())
    edges=np.asarray(edges,float)
    ds=np.diff(edges); mid=(edges[:-1]+edges[1:])/2
    r=np.sqrt(impact*impact+mid*mid)
    return dict(ds_km=ds,rho=np.asarray(density_g_cm3(r),float),r_km=r)


def split_state(x,Nn,Nt):
    i=0
    n=x[i:i+Nn]; i+=Nn
    a=x[i:i+Nn]; i+=Nn
    mn=x[i:i+Nn]; i+=Nn
    ma=x[i:i+Nn]; i+=Nn
    tm=x[i:i+Nt]; i+=Nt
    tp=x[i:i+Nt]
    return dict(nu_tau=n,antinu_tau=a,nu_mu_secondary=mn,
                antinu_mu_secondary=ma,tau_minus=tm,tau_plus=tp)


def tag(task):
    return f"{task['name']}_r{task['refine']}_tau{task['tau_refine']}"


def run_task(task):
    tn,ta=load('nu'),load('antinu')
    E=np.asarray(tn['energies'],float); inc=flux(E)
    tg=tag(task); npz=CK/f'{tg}.npz'; meta=CK/f'{tg}.json'; partial=CK/f'{tg}_partial.npz'
    if npz.exists() and meta.exists():
        return json.loads(meta.read_text())

    tr=TauHelicityFlipTransportMixedExpm(tn,ta,tau_refine=int(task['tau_refine']))
    prof=nested_profile(task['mu'],task['base_km'],task['refine'])
    Nn,Nt=tr.Nn,tr.Nt
    start=0
    if partial.exists():
        z=np.load(partial); x=z['x']; start=int(z['next_index'])
    else:
        x=np.concatenate([inc,inc,np.zeros(2*Nn),np.zeros(2*Nt)])
    t0=time.perf_counter()
    for k in range(start,len(prof['ds_km'])):
        G=tr.generator(float(prof['rho'][k]),MODE,True,True)
        x=tr._advance(G,x,float(prof['ds_km'][k])*1e5)
        if (k+1)%20==0:
            np.savez_compressed(partial,x=x,next_index=k+1)
    sec=time.perf_counter()-t0
    out=split_state(x,Nn,Nt)
    if not all(np.all(np.isfinite(v)) for v in out.values()):
        raise RuntimeError(f'nonfinite result in {tg}')
    neg=sum(float(np.sum(np.abs(v[v<0]))) for v in out.values())
    pos=sum(float(np.sum(np.abs(v[v>0]))) for v in out.values())
    rec=dict(tag=tg,name=task['name'],mu=task['mu'],base_km=task['base_km'],
             refine=task['refine'],nominal_max_step_km=task['base_km']/task['refine'],
             tau_refine=task['tau_refine'],segments=len(prof['ds_km']),seconds=sec,
             negative_mass_fraction=neg/max(pos,1e-300))
    np.savez_compressed(npz,**out)
    meta.write_text(json.dumps(rec,indent=2))
    if partial.exists(): partial.unlink()
    return rec


def load_result(name,refine,tau_refine):
    tg=f'{name}_r{refine}_tau{tau_refine}'
    npz=CK/f'{tg}.npz'; meta=CK/f'{tg}.json'
    if not (npz.exists() and meta.exists()): raise RuntimeError(f'missing {tg}')
    z=np.load(npz); return {k:z[k] for k in z.files},json.loads(meta.read_text())


def interp(E,y,e):
    y=np.asarray(y,float); j=int(np.searchsorted(E,e))
    if j<=0 or j>=len(E): return np.nan
    i=j-1
    if y[i]>0 and y[j]>0:
        f=(np.log(e)-np.log(E[i]))/(np.log(E[j])-np.log(E[i]))
        return float(np.exp(np.log(y[i])+f*(np.log(y[j])-np.log(y[i]))))
    f=(e-E[i])/(E[j]-E[i]); return float(y[i]+f*(y[j]-y[i]))


def metrics(E,inc,a,b):
    species=('nu_tau','antinu_tau','nu_mu_secondary','antinu_mu_secondary')
    maxima={s:0. for s in species}; rows=[]; l1={}
    for ep in TARGETS:
        ii=interp(E,inc,ep*1e6)
        for s in species:
            av=interp(E,a[s],ep*1e6); bv=interp(E,b[s],ep*1e6)
            if not (np.isfinite(av) and np.isfinite(bv)): raise RuntimeError('nonfinite metric')
            den=max(abs(bv),abs(ii)*1e-12,1e-300)
            d=(av-bv)/den; maxima[s]=max(maxima[s],abs(d))
            rows.append(dict(E_PeV=ep,species=s,signed_relative=d))
    for s in species:
        scale=max(float(np.max(np.abs(b[s]))),1e-300)
        mask=(inc>1e-12*np.max(inc))&(np.abs(b[s])>1e-12*scale)
        l1[s]=0. if not np.any(mask) else float(np.sum(np.abs(a[s][mask]-b[s][mask]))/max(np.sum(np.abs(b[s][mask])),1e-300))
    return dict(maxima=maxima,L1=l1,pointwise=rows,
                worst_primary=max(maxima['nu_tau'],maxima['antinu_tau']),
                worst_secondary=max(maxima['nu_mu_secondary'],maxima['antinu_mu_secondary']))


def run_parallel(tasks,workers):
    if not tasks: return
    ctx=mp.get_context('spawn')
    n=min(workers,len(tasks),os.cpu_count() or 1)
    print(f'PARALLEL workers={n} tasks={len(tasks)} internal_threads=1',flush=True)
    with ProcessPoolExecutor(max_workers=n,mp_context=ctx) as ex:
        futs={ex.submit(run_task,t):t for t in tasks}
        for fut in as_completed(futs):
            rec=fut.result()
            print(f"DONE {rec['tag']} segments={rec['segments']} sec={rec['seconds']:.1f}",flush=True)


def initial_tasks():
    out=[]
    for name,c in CASES.items():
        out.append(dict(name=name,mu=c['mu'],base_km=c['base_km'],refine=1,tau_refine=16))
        out.append(dict(name=name,mu=c['mu'],base_km=c['base_km'],refine=2,tau_refine=16))
        # full-chord accumulated tau-grid check; use fine shallow, coarse vertical
        rtau=2 if name=='shallow' else 1
        out.append(dict(name=name,mu=c['mu'],base_km=c['base_km'],refine=rtau,tau_refine=32))
    return out


def main():
    OUT.mkdir(parents=True,exist_ok=True); CK.mkdir(parents=True,exist_ok=True)
    tn=load('nu'); E=np.asarray(tn['energies'],float); inc=flux(E)
    workers=max(1,int(os.environ.get('NUTAU_NPROC','6')))

    tasks=initial_tasks()
    # Print segment counts before starting, so accidental over-refinement is visible.
    for t in tasks:
        n=len(nested_profile(t['mu'],t['base_km'],t['refine'])['ds_km'])
        print(f"PLAN {tag(t)} segments={n}",flush=True)
    run_parallel(tasks,workers)

    report={'version':'v0.20.8-adaptive-coupled-gate','benchmark':'case7_L3_majorana_flip',
            'gate':GATE,'workers_requested':workers,'cases':[],'adaptive_runs':[]}
    overall=True

    for name,c in CASES.items():
        r1,m1=load_result(name,1,16); r2,m2=load_result(name,2,16)
        spatial12=metrics(E,inc,r1,r2)
        spatial=spatial12; spatial_pair='r1_vs_r2'
        need_refine=(spatial12['worst_primary']>=GATE or spatial12['worst_secondary']>=GATE)
        if need_refine:
            t=dict(name=name,mu=c['mu'],base_km=c['base_km'],refine=4,tau_refine=16)
            print(f"ADAPT {name}: spatial gate failed; running {tag(t)} only",flush=True)
            run_parallel([t],workers)
            r4,m4=load_result(name,4,16)
            spatial=metrics(E,inc,r2,r4); spatial_pair='r2_vs_r4'
            report['adaptive_runs'].append(tag(t))
        else:
            m4=None

        rtau=2 if name=='shallow' else 1
        ra,ma=load_result(name,rtau,16); rb,mb=load_result(name,rtau,32)
        taugrid=metrics(E,inc,ra,rb)

        neg_ok=all(m['negative_mass_fraction']<1e-12 for m in (m1,m2,ma,mb) if m is not None)
        spatial_ok=(spatial['worst_primary']<GATE and spatial['worst_secondary']<GATE)
        tau_ok=(taugrid['worst_primary']<GATE and taugrid['worst_secondary']<GATE)
        passed=bool(neg_ok and spatial_ok and tau_ok); overall &= passed
        print(f"CASE {name} spatial({spatial_pair}) primary={spatial['worst_primary']:.6g} secondary={spatial['worst_secondary']:.6g}",flush=True)
        print(f"CASE {name} tau16->32 primary={taugrid['worst_primary']:.6g} secondary={taugrid['worst_secondary']:.6g} PASS={passed}",flush=True)
        report['cases'].append(dict(name=name,spatial_initial=spatial12,spatial_final=spatial,
                                    spatial_pair=spatial_pair,tau_grid=taugrid,PASS=passed))

    report['ALL_PASS']=bool(overall)
    (OUT/'coupled_gate_status_v0_20_8.json').write_text(json.dumps(report,indent=2))
    print('ALL_PASS='+str(bool(overall)),flush=True)
    print('STATUS wrote',OUT/'coupled_gate_status_v0_20_8.json',flush=True)
    raise SystemExit(0 if overall else 2)

if __name__=='__main__':
    main()
