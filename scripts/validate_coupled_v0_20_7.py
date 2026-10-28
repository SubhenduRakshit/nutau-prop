#!/usr/bin/env python3
"""v0.20.7 bounded coupled-transport gate, multicore.

Purpose: close the transport prerequisite and stop numerical-method expansion.
Uses the analytic tau-decay projection already wired into the mixed-grid solver.

One known difficult benchmark only: Case 7, Lambda=3 TeV, Majorana flip.
Checks:
  * shallow |cos z|=0.1: interface-aligned 3 -> 1.5 km
  * vertical |cos z|=1.0: interface-aligned 12 -> 6 km
  * tau grid factor 16 -> 32 on each fine spatial profile

The six independent configurations are run in separate worker processes.
Each worker remains single-threaded internally to avoid BLAS/OpenMP
oversubscription.  Existing final and partial checkpoints are reused.
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
OUT=ROOT/'results'/'v0.20.7'
CK=OUT/'checkpoints'
TARGETS=(1.,10.,100.,220.,300.)
INTERFACES=(1221.5,3480.,5701.,5771.,5971.,6151.,6346.6,6356.,6368.)
MODE='majorana_flip'
GATE=5e-3
CASES=(
    dict(name='shallow',mu=0.1,coarse_km=3.0),
    dict(name='vertical',mu=1.0,coarse_km=12.0),
)


def load(spec):
    z=np.load(T15/f'case7_L3_{spec}.npz')
    return {k:z[k] for k in z.files}


def flux(E):
    q=pd.read_json(UHE/'flux.json')
    ef=q['Energy [GeV]'].to_numpy(float)
    sed=q['log10(E^2 F(E) [GeV.cm^-2.s^-1.sr^-1])'].to_numpy(float)
    f=interp1d(ef,sed,bounds_error=False,fill_value=-np.inf)(E)
    return (10.0**f/E)/6.0


def profile(mu,coarse_km,refine):
    """PREM-interface-aligned edge nesting; refine=1 or 2."""
    mu=float(mu); refine=int(refine)
    L=2*R_EARTH_KM*mu
    impact=R_EARTH_KM*np.sqrt(max(0.,1.-mu*mu)); half=L/2
    cuts=[-half,half]
    for r in INTERFACES:
        if r>impact:
            x=np.sqrt(max(0.,r*r-impact*impact))
            if 1e-12<x<half-1e-12: cuts.extend((-x,x))
    cuts=np.array(sorted(set(round(float(x),12) for x in cuts)))
    edges=[cuts[0]]
    for a,b in zip(cuts[:-1],cuts[1:]):
        n=max(1,int(np.ceil((b-a)/float(coarse_km))))
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


def run_config(tn,ta,inc,name,mu,coarse_km,refine,tau_refine):
    tag=f'{name}_r{refine}_tau{tau_refine}'
    npz=CK/f'{tag}.npz'; meta=CK/f'{tag}.json'; partial=CK/f'{tag}_partial.npz'
    if npz.exists() and meta.exists():
        z=np.load(npz)
        return {k:z[k] for k in z.files},json.loads(meta.read_text())
    tr=TauHelicityFlipTransportMixedExpm(tn,ta,tau_refine=tau_refine)
    prof=profile(mu,coarse_km,refine)
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
        if (k+1)%100==0:
            np.savez_compressed(partial,x=x,next_index=k+1)
    sec=time.perf_counter()-t0
    out=split_state(x,Nn,Nt)
    if not all(np.all(np.isfinite(v)) for v in out.values()):
        raise RuntimeError(f'nonfinite result in {tag}')
    neg=sum(float(np.sum(np.abs(v[v<0]))) for v in out.values())
    pos=sum(float(np.sum(np.abs(v[v>0]))) for v in out.values())
    rec=dict(tag=tag,name=name,mu=mu,coarse_km=coarse_km,refine=refine,
             nominal_max_step_km=coarse_km/refine,tau_refine=tau_refine,
             segments=len(prof['ds_km']),seconds=sec,
             negative_mass_fraction=neg/max(pos,1e-300))
    np.savez_compressed(npz,**out)
    meta.write_text(json.dumps(rec,indent=2))
    if partial.exists(): partial.unlink()
    return out,rec


def worker(task):
    # Load inside each spawned process: avoids pickling solver/input state.
    tn,ta=load('nu'),load('antinu')
    E=np.asarray(tn['energies'],float); inc=flux(E)
    out,rec=run_config(tn,ta,inc,**task)
    return task['name'],int(task['refine']),int(task['tau_refine']),rec


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
    primary=max(maxima['nu_tau'],maxima['antinu_tau'])
    secondary=max(maxima['nu_mu_secondary'],maxima['antinu_mu_secondary'])
    return dict(maxima=maxima,L1=l1,pointwise=rows,worst_primary=primary,worst_secondary=secondary)


def load_completed(name,refine,tau_refine):
    tag=f'{name}_r{refine}_tau{tau_refine}'
    npz=CK/f'{tag}.npz'; meta=CK/f'{tag}.json'
    if not (npz.exists() and meta.exists()):
        raise RuntimeError(f'missing completed checkpoint {tag}')
    z=np.load(npz)
    return {k:z[k] for k in z.files},json.loads(meta.read_text())


def main():
    OUT.mkdir(parents=True,exist_ok=True); CK.mkdir(parents=True,exist_ok=True)
    tn=load('nu'); E=np.asarray(tn['energies'],float); inc=flux(E)

    tasks=[]
    for c in CASES:
        base=dict(name=c['name'],mu=c['mu'],coarse_km=c['coarse_km'])
        tasks.extend([
            dict(**base,refine=1,tau_refine=16),
            dict(**base,refine=2,tau_refine=16),
            dict(**base,refine=2,tau_refine=32),
        ])

    requested=max(1,int(os.environ.get('NUTAU_NPROC','6')))
    workers=min(requested,len(tasks),os.cpu_count() or 1)
    print(f'PARALLEL workers={workers} tasks={len(tasks)} internal_threads=1',flush=True)

    ctx=mp.get_context('spawn')
    with ProcessPoolExecutor(max_workers=workers,mp_context=ctx) as ex:
        futures={ex.submit(worker,t):t for t in tasks}
        for fut in as_completed(futures):
            t=futures[fut]
            name,refine,tau_refine,rec=fut.result()
            print(f"DONE {name} r{refine} tau{tau_refine} segments={rec['segments']} sec={rec['seconds']:.1f}",flush=True)

    report={'version':'v0.20.7-coupled-gate-multicore',
            'benchmark':'case7_L3_majorana_flip','workers':workers,'cases':[]}
    overall=True
    for c in CASES:
        coarse,mc=load_completed(c['name'],1,16)
        fine,mf=load_completed(c['name'],2,16)
        fine32,m32=load_completed(c['name'],2,32)
        spatial=metrics(E,inc,coarse,fine)
        taugrid=metrics(E,inc,fine,fine32)
        passed=(spatial['worst_primary']<GATE and spatial['worst_secondary']<GATE and
                taugrid['worst_primary']<GATE and taugrid['worst_secondary']<GATE and
                mc['negative_mass_fraction']<1e-12 and mf['negative_mass_fraction']<1e-12 and m32['negative_mass_fraction']<1e-12)
        overall &= passed
        print(f"CASE {c['name']} spatial primary={spatial['worst_primary']:.6g} secondary={spatial['worst_secondary']:.6g}",flush=True)
        print(f"CASE {c['name']} tau16->32 primary={taugrid['worst_primary']:.6g} secondary={taugrid['worst_secondary']:.6g} PASS={passed}",flush=True)
        report['cases'].append(dict(config=c,coarse=mc,fine=mf,fine_tau32=m32,
                                    spatial=spatial,tau_grid=taugrid,PASS=passed))
    report['gate']=GATE; report['ALL_PASS']=bool(overall)
    (OUT/'coupled_gate_status_v0_20_7.json').write_text(json.dumps(report,indent=2))
    print('ALL_PASS='+str(bool(overall)),flush=True)
    print('STATUS wrote',OUT/'coupled_gate_status_v0_20_7.json',flush=True)
    raise SystemExit(0 if overall else 2)


if __name__=='__main__':
    main()
