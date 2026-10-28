#!/usr/bin/env python3
"""Optional v0.19.4 energy-grid sensitivity diagnostic.

Compares the original 91-point Case-7 Lambda=3 TeV tables with freshly built
181-point tables using the same full expm A-B-A solver and the same nested,
PREM-interface-aligned 0.375-km spatial grid.
"""
from __future__ import annotations
import json, time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from nutau_prop.chirality_transport_expm import TauHelicityFlipTransportExpm
from nutau_prop.constants import R_EARTH_KM
from nutau_prop.prem import density_g_cm3

ROOT=Path(__file__).resolve().parents[1]
COARSE=ROOT/'results'/'v0.15'/'tables'
FINE=ROOT/'results'/'v0.19.4'/'tables'
UHE=ROOT/'data'/'uhe_official'
OUT=ROOT/'results'/'v0.19.4'
MU=0.1; MODE='majorana_flip'; TARGETS=(1.,10.,100.,220.,300.)
PREM=(1221.5,3480.0,5701.0,5771.0,5971.0,6151.0,6346.6,6356.0,6368.0)


def load(p):
    z=np.load(p); return {k:z[k] for k in z.files}


def flux(E):
    q=pd.read_json(UHE/'flux.json'); ef=q['Energy [GeV]'].to_numpy(float)
    sed=q['log10(E^2 F(E) [GeV.cm^-2.s^-1.sr^-1])'].to_numpy(float)
    f=interp1d(ef,sed,bounds_error=False,fill_value=-np.inf)(E)
    return (10.0**f/E)/6.0


def profile(step=0.375):
    mu=MU; L=2*R_EARTH_KM*mu; b=R_EARTH_KM*np.sqrt(1-mu*mu); half=L/2
    cuts=[-half,half]
    for r in PREM:
        if r>b:
            x=np.sqrt(max(0.0,r*r-b*b))
            if 1e-12<x<half-1e-12: cuts.extend((-x,x))
    cuts=np.array(sorted(set(round(float(x),12) for x in cuts)),float)
    factor=int(round(3.0/step)); edges=[cuts[0]]
    for a,c in zip(cuts[:-1],cuts[1:]):
        nmacro=max(1,int(np.ceil((c-a)/3.0))); macro=np.linspace(a,c,nmacro+1)
        for u,v in zip(macro[:-1],macro[1:]): edges.extend(np.linspace(u,v,factor+1)[1:].tolist())
    edges=np.asarray(edges,float); ds=np.diff(edges); mids=.5*(edges[:-1]+edges[1:])
    r=np.sqrt(b*b+mids*mids); return ds,np.asarray(density_g_cm3(r),float)


def run(tn,ta):
    E=np.asarray(tn['energies'],float); p=flux(E); tr=TauHelicityFlipTransportExpm(tn,ta)
    H=tr._operator(MODE,True); N=len(E)
    n=p.copy(); a=p.copy(); mn=np.zeros_like(n); ma=np.zeros_like(n); tm=np.zeros_like(n); tp=np.zeros_like(n)
    ds,rho=profile(0.375)
    for d,r in zip(ds,rho):
        h=d*1e5; hh=.5*h
        tm,tp,n,a,mn,ma=tr._tau_half(tm,tp,n,a,mn,ma,r,hh,True)
        x=np.concatenate([n,a,mn,ma,tm,tp]); x=tr._advance_linear(H,x,r,h)
        n,a,mn,ma,tm,tp=[x[i*N:(i+1)*N] for i in range(6)]
        tm,tp,n,a,mn,ma=tr._tau_half(tm,tp,n,a,mn,ma,r,hh,True)
    return E,p,dict(nut=n,anut=a,num=mn,anum=ma)


def interp_log(E,y,x):
    E=np.asarray(E,float); y=np.asarray(y,float); x=np.asarray(x,float)
    out=np.zeros_like(x)
    pos=y>0
    if np.sum(pos)>=2:
        out=np.exp(np.interp(np.log(x),np.log(E[pos]),np.log(y[pos]),left=-np.inf,right=-np.inf))
        out[~np.isfinite(out)]=0.0
    else:
        out=np.interp(x,E,y,left=0.0,right=0.0)
    return out


def main():
    cn=COARSE/'case7_L3_nu.npz'; ca=COARSE/'case7_L3_antinu.npz'
    fn=FINE/'case7_L3_refined181_nu.npz'; fa=FINE/'case7_L3_refined181_antinu.npz'
    for p in (cn,ca,fn,fa):
        if not p.exists(): raise FileNotFoundError(p)
    t=time.perf_counter(); Ec,pc,oc=run(load(cn),load(ca)); tc=time.perf_counter()-t
    t=time.perf_counter(); Ef,pf,of=run(load(fn),load(fa)); tf=time.perf_counter()-t
    rows=[]; maxima={'primary':0.0,'secondary':0.0}
    for ep in TARGETS:
        x=ep*1e6; inc=float(interp_log(Ef,pf,np.array([x]))[0])
        for sp in ('nut','anut','num','anum'):
            c=float(interp_log(Ec,oc[sp],np.array([x]))[0]); f=float(interp_log(Ef,of[sp],np.array([x]))[0])
            den=max(abs(f),abs(inc)*1e-12,1e-300); d=(c-f)/den
            cls='primary' if sp in ('nut','anut') else 'secondary'; maxima[cls]=max(maxima[cls],abs(d))
            rows.append(dict(E_PeV=ep,species=sp,coarse91=c,fine181=f,signed_relative_difference=d,
                             absolute_relative_difference=abs(d)))
    common=np.geomspace(max(Ec[0],Ef[0]),min(Ec[-1],Ef[-1]),241); l1={}
    for sp in ('nut','anut','num','anum'):
        c=interp_log(Ec,oc[sp],common); f=interp_log(Ef,of[sp],common)
        mask=f>1e-12*max(float(np.max(f)),1e-300)
        l1[sp]=0.0 if not np.any(mask) else float(np.sum(np.abs(c[mask]-f[mask]))/max(np.sum(np.abs(f[mask])),1e-300))
    out=dict(version='v0.19.4-energy-grid-diagnostic',certification=False,coarse_nenergy=len(Ec),fine_nenergy=len(Ef),
             coarse_seconds=tc,fine_seconds=tf,max_errors=maxima,L1_common_log_grid=l1,pointwise=rows,
             note='91-vs-181 energy-grid sensitivity at fixed nested 0.375-km spatial grid; not a production gate.')
    (OUT/'energy_grid_diagnostic_v0_19_4.json').write_text(json.dumps(out,indent=2))
    pd.DataFrame(rows).to_csv(OUT/'energy_grid_pointwise_v0_19_4.csv',index=False)
    print(f"ENERGY GRID 91->181 primary={maxima['primary']:.6g} secondary={maxima['secondary']:.6g}",flush=True)
    print('WROTE',OUT/'energy_grid_diagnostic_v0_19_4.json',flush=True)

if __name__=='__main__': main()
