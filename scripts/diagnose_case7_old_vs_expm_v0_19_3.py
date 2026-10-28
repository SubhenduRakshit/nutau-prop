#!/usr/bin/env python3
"""Targeted diagnostic for the v0.19.2 Phase-B Case-7 discrepancy.

Compares the expm solver with the old one-interaction solver as max_p -> 0
for Case 7, Lambda=3 TeV, |cos z|=0.1, Majorana flip.  No production output.
Each old-solver max_p run is checkpointed independently.
"""
from __future__ import annotations
import json, multiprocessing as mp, os, time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from nutau_prop.chirality_transport import TauHelicityFlipTransport
from nutau_prop.chirality_transport_expm import TauHelicityFlipTransportExpm

ROOT=Path(__file__).resolve().parents[1]
T15=ROOT/'results'/'v0.15'/'tables'
UHE=ROOT/'data'/'uhe_official'
OUT=ROOT/'results'/'v0.19.3'
CK=OUT/'checkpoints'
TARGETS=(1.,10.,100.,220.,300.)
LAM=3.0; MU=0.1; MODE='majorana_flip'
MAXPS=(0.02,0.01,0.005,0.0025,0.00125)

def tag(x): return (f'{float(x):.6g}').replace('.','p')
def npz(p):
    z=np.load(p); return {k:z[k] for k in z.files}
def load(s): return npz(T15/f'case7_L{tag(LAM)}_{s}.npz')
def flux(E):
    q=pd.read_json(UHE/'flux.json')
    ef=q['Energy [GeV]'].to_numpy(float)
    sed=q['log10(E^2 F(E) [GeV.cm^-2.s^-1.sr^-1])'].to_numpy(float)
    f=interp1d(ef,sed,bounds_error=False,fill_value=-np.inf)(E)
    return (10.0**f/E)/6.0

def ip(E,y,t):
    j=int(np.searchsorted(E,t))
    if j==0 or j>=len(E): return np.nan
    i=j-1
    if y[i]>0 and y[j]>0:
        a=(np.log(t)-np.log(E[i]))/(np.log(E[j])-np.log(E[i]))
        return float(np.exp(np.log(y[i])+a*(np.log(y[j])-np.log(y[i]))))
    a=(t-E[i])/(E[j]-E[i]); return float(y[i]+a*(y[j]-y[i]))

def unpack(o):
    return dict(nut=o['nu_tau'], anut=o['antinu_tau'],
                num=o['nu_mu_secondary'], anum=o['antinu_mu_secondary'])

def run_expm():
    tn,ta=load('nu'),load('antinu'); E=tn['energies'].astype(float); p=flux(E)
    tr=TauHelicityFlipTransportExpm(tn,ta)
    return E,p,unpack(tr.propagate(MU,p,p,mode=MODE,use_np=True,max_step_km=1.5,track_secondary_numu=True))

def run_old(maxp):
    tn,ta=load('nu'),load('antinu'); E=tn['energies'].astype(float); p=flux(E)
    tr=TauHelicityFlipTransport(tn,ta)
    t=time.perf_counter()
    o=tr.propagate(MU,p,p,mode=MODE,use_np=True,max_step_km=1.5,max_p=maxp,track_secondary_numu=True)
    return time.perf_counter()-t,E,p,unpack(o)

def metrics(E,p,val,ref):
    rows=[]; maxima={'primary':0.0,'secondary':0.0,'all':0.0}
    for ep in TARGETS:
        inc=ip(E,p,ep*1e6)
        for sp in ('nut','anut','num','anum'):
            a=ip(E,val[sp],ep*1e6); b=ip(E,ref[sp],ep*1e6)
            den=max(abs(b),abs(inc)*1e-12,1e-300)
            er=abs(a-b)/den
            cls='primary' if sp in ('nut','anut') else 'secondary'
            maxima[cls]=max(maxima[cls],er); maxima['all']=max(maxima['all'],er)
            rows.append(dict(E_PeV=ep,species=sp,old=a,expm=b,relative_error=er,
                             ratio_old_over_expm=(a/b if b!=0 else np.nan)))
    l1={}
    for sp in ('nut','anut','num','anum'):
        mask=(p>1e-12*np.max(p)) & (np.abs(ref[sp])>1e-12*max(np.max(np.abs(ref[sp])),1e-300))
        l1[sp]=0.0 if not np.any(mask) else float(np.sum(np.abs(val[sp][mask]-ref[sp][mask]))/max(np.sum(np.abs(ref[sp][mask])),1e-300))
    return rows,maxima,l1

def worker(maxp):
    cp=CK/f'old_maxp_{tag(maxp)}.json'
    if cp.exists(): return json.loads(cp.read_text())
    E,p,ref=run_expm()
    sec,E2,p2,val=run_old(maxp)
    rows,mx,l1=metrics(E,p,val,ref)
    out=dict(max_p=maxp,seconds=sec,max_errors=mx,L1=l1,pointwise=rows)
    tmp=cp.with_suffix('.tmp'); tmp.write_text(json.dumps(out,indent=2)); tmp.replace(cp)
    return out

def main():
    OUT.mkdir(parents=True,exist_ok=True); CK.mkdir(parents=True,exist_ok=True)
    n=min(int(os.environ.get('NUTAU_NPROC','5')),len(MAXPS))
    print(f'Case7 diagnostic: Lambda={LAM} TeV cosz={MU} mode={MODE} workers={n}',flush=True)
    ctx=mp.get_context('spawn')
    with ctx.Pool(n) as pool:
        rr=list(pool.imap_unordered(worker,MAXPS,chunksize=1))
    rr=sorted(rr,key=lambda r:r['max_p'],reverse=True)
    for r in rr:
        m=r['max_errors']
        print(f"max_p={r['max_p']:.5g} sec={r['seconds']:.1f} primary={m['primary']:.6g} secondary={m['secondary']:.6g} all={m['all']:.6g}",flush=True)
    # empirical convergence ratios versus max_p
    summary=[]
    for r in rr:
        summary.append(dict(max_p=r['max_p'],seconds=r['seconds'],
                            max_primary=r['max_errors']['primary'],
                            max_secondary=r['max_errors']['secondary'],
                            max_all=r['max_errors']['all'],
                            **{f'L1_{k}':v for k,v in r['L1'].items()}))
    pd.DataFrame(summary).to_csv(OUT/'case7_old_to_expm_convergence.csv',index=False)
    (OUT/'case7_old_to_expm_convergence.json').write_text(json.dumps(rr,indent=2))
    print('DETAIL worst rows at smallest max_p:',flush=True)
    best=rr[-1]
    for row in sorted(best['pointwise'],key=lambda x:x['relative_error'],reverse=True)[:8]:
        print(json.dumps(row),flush=True)

if __name__=='__main__': main()
