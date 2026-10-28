#!/usr/bin/env python3
"""v0.20.9 reuse-first coupled-transport gate.

This validator never discards completed v0.20.7 work.  The v0.20.7 runs used
the same analytic mixed-grid solver and are finer than the v0.20.8 redesign,
so completed v0.20.7 checkpoints are reused directly as convergence evidence.

For each chord:
  spatial gate: prefer old v0.20.7 tau16 pair (3->1.5 km shallow,
                12->6 km vertical).  If unavailable, run only the cheap
                v0.20.8-style pair (48->24 km or 192->96 km).
  tau-grid gate: prefer old v0.20.7 matched-profile tau16->tau32 pair.
                 If unavailable, run only a cheap matched-profile pair
                 (24 km shallow or 192 km vertical).

No automatic escalation to finer grids.  If a reused or cheap pair fails the
0.5% gate, the script stops and reports NEEDS_FURTHER_REFINEMENT rather than
launching a surprise long job.

Completed/partial v0.20.8 coarse checkpoints are also reused.  New independent
jobs run concurrently; internal BLAS/OpenMP threading remains one per worker.
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
OUT=ROOT/'results'/'v0.20.9'
OLDCK=ROOT/'results'/'v0.20.7'/'checkpoints'
COARSECK=ROOT/'results'/'v0.20.8'/'checkpoints'  # reuse v0.20.8 partial/final work
TARGETS=(1.,10.,100.,220.,300.)
INTERFACES=(1221.5,3480.,5701.,5771.,5971.,6151.,6346.6,6356.,6368.)
MODE='majorana_flip'
GATE=5e-3
CASES={
    'shallow':dict(mu=0.1,base_km=48.0,old_coarse_km=3.0,old_fine_km=1.5),
    'vertical':dict(mu=1.0,base_km=192.0,old_coarse_km=12.0,old_fine_km=6.0),
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
    mu=float(mu); refine=int(refine); base_km=float(base_km)
    if refine not in (1,2): raise ValueError('refine must be 1 or 2')
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


def tg(name,refine,tau_refine):
    return f'{name}_r{int(refine)}_tau{int(tau_refine)}'


def final_exists(directory,name,refine,tau_refine):
    tag=tg(name,refine,tau_refine)
    return (directory/f'{tag}.npz').exists() and (directory/f'{tag}.json').exists()


def load_result(directory,name,refine,tau_refine):
    tag=tg(name,refine,tau_refine)
    npz=directory/f'{tag}.npz'; meta=directory/f'{tag}.json'
    if not (npz.exists() and meta.exists()): raise RuntimeError(f'missing {npz}')
    z=np.load(npz)
    return {k:z[k] for k in z.files},json.loads(meta.read_text())


def task_tag(task):
    return tg(task['name'],task['refine'],task['tau_refine'])


def run_task(task):
    tn,ta=load('nu'),load('antinu')
    E=np.asarray(tn['energies'],float); inc=flux(E)
    tag=task_tag(task)
    npz=COARSECK/f'{tag}.npz'; meta=COARSECK/f'{tag}.json'; partial=COARSECK/f'{tag}_partial.npz'
    if npz.exists() and meta.exists():
        rec=json.loads(meta.read_text()); rec['reused_existing']=True
        return rec
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
        raise RuntimeError(f'nonfinite result in {tag}')
    neg=sum(float(np.sum(np.abs(v[v<0]))) for v in out.values())
    pos=sum(float(np.sum(np.abs(v[v>0]))) for v in out.values())
    rec=dict(tag=tag,name=task['name'],mu=task['mu'],base_km=task['base_km'],
             refine=task['refine'],nominal_max_step_km=task['base_km']/task['refine'],
             tau_refine=task['tau_refine'],segments=len(prof['ds_km']),seconds=sec,
             negative_mass_fraction=neg/max(pos,1e-300),reused_existing=False)
    np.savez_compressed(npz,**out)
    meta.write_text(json.dumps(rec,indent=2))
    if partial.exists(): partial.unlink()
    return rec


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
    # de-duplicate exact tasks
    unique=[]; seen=set()
    for t in tasks:
        key=(t['name'],int(t['refine']),int(t['tau_refine']))
        if key not in seen: seen.add(key); unique.append(t)
    tasks=unique
    n=min(max(1,int(workers)),len(tasks),os.cpu_count() or 1)
    print(f'PARALLEL workers={n} tasks={len(tasks)} internal_threads=1',flush=True)
    ctx=mp.get_context('spawn')
    with ProcessPoolExecutor(max_workers=n,mp_context=ctx) as ex:
        futs={ex.submit(run_task,t):t for t in tasks}
        for fut in as_completed(futs):
            rec=fut.result()
            print(f"DONE {rec['tag']} segments={rec['segments']} sec={rec['seconds']:.1f} reused_existing={rec.get('reused_existing',False)}",flush=True)


def main():
    OUT.mkdir(parents=True,exist_ok=True); COARSECK.mkdir(parents=True,exist_ok=True)
    tn=load('nu'); E=np.asarray(tn['energies'],float); inc=flux(E)
    workers=max(1,int(os.environ.get('NUTAU_NPROC','6')))

    report={'version':'v0.20.9-reuse-first','gate':GATE,'benchmark':'case7_L3_majorana_flip',
            'old_checkpoint_dir':str(OLDCK),'coarse_checkpoint_dir':str(COARSECK),
            'reuse':{},'new_tasks':[],'cases':[]}
    tasks=[]; source={}

    for name,c in CASES.items():
        # Spatial: prefer completed old fine pair.
        if final_exists(OLDCK,name,1,16) and final_exists(OLDCK,name,2,16):
            source[(name,'spatial')]=('old',1,2)
            print(f'REUSE {name} spatial: v0.20.7 r1_tau16 + r2_tau16',flush=True)
        else:
            source[(name,'spatial')]=('coarse',1,2)
            for r in (1,2):
                tasks.append(dict(name=name,mu=c['mu'],base_km=c['base_km'],refine=r,tau_refine=16))
            print(f'NEW {name} spatial: cheap r1_tau16 + r2_tau16',flush=True)

        # Tau-grid: prefer completed old matched fine pair r2 tau16/tau32.
        if final_exists(OLDCK,name,2,16) and final_exists(OLDCK,name,2,32):
            source[(name,'tau')]=('old',2)
            print(f'REUSE {name} tau-grid: v0.20.7 r2 tau16 + tau32',flush=True)
        else:
            rtau=2 if name=='shallow' else 1
            source[(name,'tau')]=('coarse',rtau)
            for trf in (16,32):
                tasks.append(dict(name=name,mu=c['mu'],base_km=c['base_km'],refine=rtau,tau_refine=trf))
            print(f'NEW {name} tau-grid: cheap r{rtau} tau16 + tau32',flush=True)

    # Show exactly what would run, before any compute.
    unique=[]; seen=set()
    for t in tasks:
        key=(t['name'],int(t['refine']),int(t['tau_refine']))
        if key in seen: continue
        seen.add(key); unique.append(t)
    for t in unique:
        nseg=len(nested_profile(t['mu'],t['base_km'],t['refine'])['ds_km'])
        if nseg>200 and os.environ.get('NUTAU_ALLOW_LONG','0')!='1':
            raise RuntimeError(f'refusing unexpected long task {task_tag(t)} with {nseg} segments')
        print(f"PLAN_NEW {task_tag(t)} segments={nseg}",flush=True)
        report['new_tasks'].append(dict(tag=task_tag(t),segments=nseg))

    run_parallel(unique,workers)

    overall=True; needs_refine=False
    for name,c in CASES.items():
        ssrc=source[(name,'spatial')]
        if ssrc[0]=='old':
            a,ma=load_result(OLDCK,name,1,16); b,mb=load_result(OLDCK,name,2,16)
            sdesc=f"v0.20.7 {c['old_coarse_km']}km->{c['old_fine_km']}km"
        else:
            a,ma=load_result(COARSECK,name,1,16); b,mb=load_result(COARSECK,name,2,16)
            sdesc=f"coarse {c['base_km']}km->{c['base_km']/2}km"
        spatial=metrics(E,inc,a,b)

        tsrc=source[(name,'tau')]
        if tsrc[0]=='old':
            r=tsrc[1]; a16,m16=load_result(OLDCK,name,r,16); a32,m32=load_result(OLDCK,name,r,32)
            tdesc=f"v0.20.7 matched profile r{r}"
        else:
            r=tsrc[1]; a16,m16=load_result(COARSECK,name,r,16); a32,m32=load_result(COARSECK,name,r,32)
            tdesc=f"coarse matched profile r{r}"
        taugrid=metrics(E,inc,a16,a32)

        spatial_ok=spatial['worst_primary']<GATE and spatial['worst_secondary']<GATE
        tau_ok=taugrid['worst_primary']<GATE and taugrid['worst_secondary']<GATE
        neg_ok=all(m.get('negative_mass_fraction',0.0)<1e-12 for m in (ma,mb,m16,m32))
        passed=bool(spatial_ok and tau_ok and neg_ok)
        if not passed: needs_refine=True
        overall &= passed
        print(f"CASE {name} spatial[{sdesc}] primary={spatial['worst_primary']:.6g} secondary={spatial['worst_secondary']:.6g}",flush=True)
        print(f"CASE {name} tau16->32[{tdesc}] primary={taugrid['worst_primary']:.6g} secondary={taugrid['worst_secondary']:.6g} PASS={passed}",flush=True)
        report['cases'].append(dict(name=name,spatial_source=sdesc,tau_source=tdesc,
                                    spatial=spatial,tau_grid=taugrid,PASS=passed))

    report['ALL_PASS']=bool(overall)
    report['NEEDS_FURTHER_REFINEMENT']=bool(needs_refine)
    (OUT/'coupled_gate_status_v0_20_9.json').write_text(json.dumps(report,indent=2))
    print('ALL_PASS='+str(bool(overall)),flush=True)
    print('NEEDS_FURTHER_REFINEMENT='+str(bool(needs_refine)),flush=True)
    print('STATUS wrote',OUT/'coupled_gate_status_v0_20_9.json',flush=True)
    raise SystemExit(0 if overall else 2)

if __name__=='__main__':
    main()
