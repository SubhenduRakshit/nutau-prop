#!/usr/bin/env python3
"""v0.20.3 mixed-grid full-generator smoke test.

No full chord.  Benchmarks one representative 3-km slab, checks the exact
matrix semigroup property, checks rectangular CC production against the
historical represented-support row sums, and rechecks the factor-16 tau-only
one-slab daughter accuracy with an independent continuous-loss reference.
"""
from __future__ import annotations

import json, math, time
from pathlib import Path
import numpy as np
import pandas as pd
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import interp1d
from scipy.sparse import bmat, csr_matrix, diags

from nutau_prop.chirality_transport_mixed_expm import TauHelicityFlipTransportMixedExpm, _deposit_row
from nutau_prop.kernels import tau_nutau_pdf, tau_secondary_numu_pdf
from nutau_prop.constants import M_TAU, C_TAU_CM
from nutau_prop.prem import chord_profile

ROOT=Path(__file__).resolve().parents[1]
T15=ROOT/'results'/'v0.15'/'tables'
UHE=ROOT/'data'/'uhe_official'
OUT=ROOT/'results'/'v0.20.3'
RHO=2.9; H=0.75e5; BR=0.18


def load(spec):
    z=np.load(T15/f'case7_L3_{spec}.npz')
    return {k:z[k] for k in z.files}


def flux(E):
    q=pd.read_json(UHE/'flux.json')
    ef=q['Energy [GeV]'].to_numpy(float)
    sed=q['log10(E^2 F(E) [GeV.cm^-2.s^-1.sr^-1])'].to_numpy(float)
    f=interp1d(ef,sed,bounds_error=False,fill_value=-np.inf)(E)
    return (10.0**f/E)/6.0


def rel_l1(a,b):
    return float(np.sum(np.abs(np.asarray(a)-np.asarray(b)))/max(np.sum(np.abs(b)),1e-300))


def exact_tau_reference(E0, Enu, Etau, sign):
    pol=-1.0 if sign=='minus' else +1.0
    beta=0.85e-6; b=beta*RHO; a=M_TAU/(C_TAU_CM*E0)
    U=(a/b)*np.expm1(b*H); surv=math.exp(-min(float(U),745.0)); Eend=E0*math.exp(-b*H)
    tau=np.zeros(len(Etau)); _deposit_row(tau,np.log(Etau),Eend,surv)
    p=np.zeros(len(Enu)); s=np.zeros(len(Enu)); ln=np.log(Enu)
    pdec=-np.expm1(-min(float(U),745.0))
    qx,qw=leggauss(48); tt=0.5*(qx+1.0)*pdec; ww=0.5*pdec*qw
    uu=-np.log1p(-tt); Edec=E0/(1.0+(b/a)*uu)
    z=np.linspace(1e-6,1-1e-6,800)
    fp=tau_nutau_pdf(z,pol); fp/=np.trapezoid(fp,z)
    fs=tau_secondary_numu_pdf(z,pol,BR); ns=np.trapezoid(fs,z)
    if ns>0: fs*=BR/ns
    dz=z[1]-z[0]
    for ed,wd in zip(Edec,ww):
        for zz,pp,ss in zip(z,fp,fs):
            ee=float(ed*zz)
            if pp: _deposit_row(p,ln,ee,float(wd*pp*dz))
            if ss: _deposit_row(s,ln,ee,float(wd*ss*dz))
    return tau,p,s,surv,Eend


def tau_only_generator(tr,sign):
    Nn,Nt=tr.Nn,tr.Nt
    Znn=csr_matrix((Nn,Nn)); Znt=csr_matrix((Nn,Nt)); Ztn=csr_matrix((Nt,Nn))
    L=RHO*tr._tau_loss_per_rho(); Q,Dm,Dp,Sm,Sp=tr._tau_decay_blocks(True)
    if sign=='minus': D,S=Dm,Sm
    else: D,S=Dp,Sp
    return bmat([[Znn,Znn,Znt],[Znn,Znn,Znt],[D,S,L+Q]],format='csr')


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    tn,ta=load('nu'),load('antinu')
    t0=time.perf_counter(); tr=TauHelicityFlipTransportMixedExpm(tn,ta,tau_refine=16); build=time.perf_counter()-t0
    E=tr.E; Nt=tr.Nt; Nn=tr.Nn; p=flux(E)
    print(f'BUILD mixed factor=16 Nnu={Nn} Ntau={Nt} seconds={build:.3f}',flush=True)

    # CC represented-support consistency against the old same-grid operator.
    oldn=np.asarray(tr.nu.cc_to_tau().sum(axis=1)).ravel()
    olda=np.asarray(tr.anu.cc_to_tau().sum(axis=1)).ravel()
    newn=np.asarray(tr.Tn.sum(axis=1)).ravel(); newa=np.asarray(tr.Ta.sum(axis=1)).ravel()
    cc_n=float(np.max(np.abs(newn-oldn)/np.maximum(np.abs(oldn),1e-300)))
    cc_a=float(np.max(np.abs(newa-olda)/np.maximum(np.abs(olda),1e-300)))
    print(f'CC-TRANSFER rowsum maxrel nu={cc_n:.6g} antinu={cc_a:.6g}',flush=True)

    # One representative full mixed-grid slab only.
    x=np.concatenate([p,p,np.zeros(2*Nn),np.zeros(2*Nt)])
    G=tr.generator(RHO,'majorana_flip',True,True)
    t0=time.perf_counter(); y=tr._advance(G,x,3e5); slab=time.perf_counter()-t0
    prof=chord_profile(0.1,max_step_km=3.0); nseg=len(prof['ds_km'])
    print(f'SLAB 3km seconds={slab:.3f} projected_serial_sec={slab*nseg:.1f} segments={nseg} finite={np.all(np.isfinite(y))} min={np.min(y):.6g}',flush=True)

    # Exact discrete semigroup.
    h=0.75e5; one=tr._advance(G,x,h); two=tr._advance(G,tr._advance(G,x,h/2),h/2)
    semi=rel_l1(two,one)
    print(f'MIXED-GENERATOR semigroup_L1={semi:.6g}',flush=True)

    rows=[]
    for sign in ('minus','plus'):
        Gt=tau_only_generator(tr,sign)
        for target in (1.,10.,100.,300.):
            jn=int(np.argmin(np.abs(np.log(E)-np.log(target*1e6))))
            E0=float(E[jn]); jt=int(np.argmin(np.abs(np.log(tr.Etau)-np.log(E0))))
            xx=np.zeros(2*Nn+Nt); xx[2*Nn+jt]=1.0
            yy=np.asarray(tr._advance(Gt,xx,H),float)
            prim=yy[:Nn]; sec=yy[Nn:2*Nn]; tau=yy[2*Nn:]
            rt,rp,rs,surv,Eend=exact_tau_reference(E0,E,tr.Etau,sign)
            nt=float(np.sum(tau)); prob=tau/max(nt,1e-300)
            w1=float(np.sum(prob*np.abs(np.log(tr.Etau)-math.log(Eend)))) if nt>0 else None
            rec=dict(sign=sign,E_requested_PeV=target,E0_grid_PeV=E0/1e6,
                     tau_number_relative_error=abs(nt-surv)/max(surv,1e-300),
                     W1_logE=w1,tau_L1=rel_l1(tau,rt),
                     primary_daughter_L1=rel_l1(prim,rp),secondary_daughter_L1=rel_l1(sec,rs))
            rows.append(rec)
            print('TAU-MIXED '+json.dumps(rec),flush=True)

    worstp=max(r['primary_daughter_L1'] for r in rows)
    worsts=max(r['secondary_daughter_L1'] for r in rows)
    status=dict(version='v0.20.3-smoke',tau_refine=16,N_nu=Nn,N_tau=Nt,
                build_seconds=build,one_slab_seconds=slab,projected_serial_seconds=slab*nseg,
                semigroup_L1=semi,cc_rowsum_maxrel_nu=cc_n,cc_rowsum_maxrel_antinu=cc_a,
                tau_benchmarks=rows,worst_primary_daughter_L1=worstp,
                worst_secondary_daughter_L1=worsts,
                note='Smoke only. Full chord remains forbidden until runtime and these gates are reviewed.')
    (OUT/'mixed_expm_smoke_status_v0_20_3.json').write_text(json.dumps(status,indent=2))
    print(f'SUMMARY worst_primary={worstp:.6g} worst_secondary={worsts:.6g}',flush=True)
    print('STATUS wrote',OUT/'mixed_expm_smoke_status_v0_20_3.json',flush=True)

if __name__=='__main__': main()
