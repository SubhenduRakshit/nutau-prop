#!/usr/bin/env python3
"""v0.20.4 normalization-corrected mixed-grid smoke validation.

This is diagnostic/certification work, not a production run.

Changes relative to v0.20.3
--------------------------
* decay kernels use consistent trapezoid quadrature weights;
* the continuous-loss reference uses independent split Gauss-Legendre z
  quadrature and a locally implemented daughter deposition routine;
* explicit numerical gates produce a nonzero exit code on failure.

No full Earth chord is run here.
"""
from __future__ import annotations

import json, math, sys, time
from pathlib import Path

import numpy as np
import pandas as pd
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import interp1d
from scipy.sparse import bmat, csr_matrix

from nutau_prop.chirality_transport_mixed_expm import TauHelicityFlipTransportMixedExpm
from nutau_prop.kernels import (
    tau_nutau_pdf, tau_secondary_numu_pdf, normalized_trapezoid_weights,
    MTAU, MPI, MRHO, MA1,
)
from nutau_prop.constants import M_TAU, C_TAU_CM
from nutau_prop.prem import chord_profile

ROOT=Path(__file__).resolve().parents[1]
T15=ROOT/'results'/'v0.15'/'tables'
UHE=ROOT/'data'/'uhe_official'
OUT=ROOT/'results'/'v0.20.4'
RHO=2.9
H=0.75e5
BR=0.18
TARGETS=(1.,10.,100.,300.)

GATES={
    'semigroup_L1': 1e-10,
    'cc_rowsum_maxrel': 1e-10,
    'quadrature_multiplicity_abs': 1e-12,
    'negative_mass_fraction': 1e-12,
    'worst_primary_daughter_L1': 5e-3,
    'worst_secondary_daughter_L1': 5e-3,
    'max_tau_survival_abs_error': 1e-4,
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


def rel_l1(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    return float(np.sum(np.abs(a-b))/max(np.sum(np.abs(b)),1e-300))


def ref_deposit(row, log_grid, eout, weight):
    """Independent local log-E linear deposition used only by the reference."""
    if weight==0.0 or eout<=0:
        return
    u=math.log(float(eout))
    if u<log_grid[0] or u>log_grid[-1]:
        return
    j=int(np.searchsorted(log_grid,u))
    if j==0:
        row[0]+=weight
    elif j>=len(log_grid):
        row[-1]+=weight
    else:
        frac=(u-log_grid[j-1])/(log_grid[j]-log_grid[j-1])
        row[j-1]+=weight*(1.0-frac)
        row[j]+=weight*frac


def split_gauss_z_weights(pol, secondary=False, order=64):
    """Independent z quadrature split at all PDF discontinuities."""
    bounds=[0.0,0.3,
            1.0-(MA1/MTAU)**2,
            1.0-(MRHO/MTAU)**2,
            1.0-(MPI/MTAU)**2,1.0]
    bounds=sorted(set(max(0.0,min(1.0,float(x))) for x in bounds))
    xg,wg=leggauss(int(order))
    zs=[]; ws=[]
    for a,b in zip(bounds[:-1],bounds[1:]):
        if b<=a:
            continue
        z=0.5*(xg+1.0)*(b-a)+a
        jac=0.5*(b-a)
        if secondary:
            f=tau_secondary_numu_pdf(z,pol,BR)
        else:
            f=tau_nutau_pdf(z,pol)
        zs.append(z); ws.append(jac*wg*f)
    z=np.concatenate(zs); w=np.concatenate(ws)
    target=BR if secondary else 1.0
    total=float(np.sum(w))
    if total<=0:
        raise RuntimeError('independent z quadrature has zero normalization')
    w*=target/total
    return z,w


def exact_tau_reference(E0, Enu, Etau, sign):
    pol=-1.0 if sign=='minus' else +1.0
    beta=0.85e-6
    b=beta*RHO
    a=M_TAU/(C_TAU_CM*E0)
    U=(a/b)*np.expm1(b*H) if b!=0 else a*H
    surv=math.exp(-min(float(U),745.0))
    Eend=E0*math.exp(-b*H)

    tau=np.zeros(len(Etau)); lnt=np.log(Etau)
    ref_deposit(tau,lnt,Eend,surv)
    p=np.zeros(len(Enu)); s=np.zeros(len(Enu)); lnn=np.log(Enu)

    pdec=-np.expm1(-min(float(U),745.0))
    qx,qw=leggauss(64)
    tt=0.5*(qx+1.0)*pdec
    ww=0.5*pdec*qw
    uu=-np.log1p(-tt)
    Edec=E0/(1.0+(b/a)*uu) if b!=0 else np.full_like(uu,E0)

    zp,wp=split_gauss_z_weights(pol,False,64)
    zs,ws=split_gauss_z_weights(pol,True,64)
    for ed,wd in zip(Edec,ww):
        for zz,wz in zip(zp,wp):
            ref_deposit(p,lnn,float(ed*zz),float(wd*wz))
        for zz,wz in zip(zs,ws):
            ref_deposit(s,lnn,float(ed*zz),float(wd*wz))
    return tau,p,s,surv,Eend


def tau_only_generator(tr,sign):
    Nn,Nt=tr.Nn,tr.Nt
    Znn=csr_matrix((Nn,Nn)); Znt=csr_matrix((Nn,Nt))
    L=RHO*tr._tau_loss_per_rho(); Q,Dm,Dp,Sm,Sp=tr._tau_decay_blocks(True)
    D,S=(Dm,Sm) if sign=='minus' else (Dp,Sp)
    return bmat([[Znn,Znn,Znt],[Znn,Znn,Znt],[D,S,L+Q]],format='csr')


def multiplicity_regression():
    z=np.linspace(1e-6,1-1e-6,800)
    out={}
    for pol,label in ((-1.0,'minus'),(+1.0,'plus')):
        wp=normalized_trapezoid_weights(z,tau_nutau_pdf(z,pol),1.0)
        ws=normalized_trapezoid_weights(z,tau_secondary_numu_pdf(z,pol,BR),BR)
        out[f'primary_{label}']=abs(float(np.sum(wp))-1.0)
        out[f'secondary_{label}']=abs(float(np.sum(ws))-BR)
    return out


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    tn,ta=load('nu'),load('antinu')
    t0=time.perf_counter()
    tr=TauHelicityFlipTransportMixedExpm(tn,ta,tau_refine=16)
    build=time.perf_counter()-t0
    E=tr.E; Nt=tr.Nt; Nn=tr.Nn; p=flux(E)
    print(f'BUILD v0.20.4 factor=16 Nnu={Nn} Ntau={Nt} seconds={build:.3f}',flush=True)

    mult=multiplicity_regression()
    multmax=max(mult.values())
    print('QUADRATURE multiplicity_abs '+json.dumps(mult,sort_keys=True),flush=True)

    oldn=np.asarray(tr.nu.cc_to_tau().sum(axis=1)).ravel()
    olda=np.asarray(tr.anu.cc_to_tau().sum(axis=1)).ravel()
    newn=np.asarray(tr.Tn.sum(axis=1)).ravel(); newa=np.asarray(tr.Ta.sum(axis=1)).ravel()
    cc_n=float(np.max(np.abs(newn-oldn)/np.maximum(np.abs(oldn),1e-300)))
    cc_a=float(np.max(np.abs(newa-olda)/np.maximum(np.abs(olda),1e-300)))
    ccmax=max(cc_n,cc_a)
    print(f'CC-TRANSFER rowsum maxrel nu={cc_n:.6g} antinu={cc_a:.6g}',flush=True)

    x=np.concatenate([p,p,np.zeros(2*Nn),np.zeros(2*Nt)])
    G=tr.generator(RHO,'majorana_flip',True,True)
    t0=time.perf_counter(); y=tr._advance(G,x,3e5); slab=time.perf_counter()-t0
    prof=chord_profile(0.1,max_step_km=3.0); nseg=len(prof['ds_km'])
    neg=float(np.sum(np.abs(y[y<0])))
    pos=float(np.sum(np.abs(y[y>0])))
    negfrac=neg/max(pos,1e-300)
    print(f'SLAB 3km seconds={slab:.3f} projected_serial_sec={slab*nseg:.1f} segments={nseg} finite={np.all(np.isfinite(y))} negative_mass_fraction={negfrac:.6g}',flush=True)

    h=0.75e5
    one=tr._advance(G,x,h)
    two=tr._advance(G,tr._advance(G,x,h/2),h/2)
    semi=rel_l1(two,one)
    print(f'MIXED-GENERATOR semigroup_L1={semi:.6g}',flush=True)

    rows=[]
    for sign in ('minus','plus'):
        Gt=tau_only_generator(tr,sign)
        for target in TARGETS:
            jn=int(np.argmin(np.abs(np.log(E)-np.log(target*1e6))))
            E0=float(E[jn]); jt=int(np.argmin(np.abs(np.log(tr.Etau)-np.log(E0))))
            xx=np.zeros(2*Nn+Nt); xx[2*Nn+jt]=1.0
            yy=np.asarray(tr._advance(Gt,xx,H),float)
            prim=yy[:Nn]; sec=yy[Nn:2*Nn]; tau=yy[2*Nn:]
            rt,rp,rs,surv,Eend=exact_tau_reference(E0,E,tr.Etau,sign)
            nt=float(np.sum(tau)); prob=tau/max(nt,1e-300)
            w1=float(np.sum(prob*np.abs(np.log(tr.Etau)-math.log(Eend)))) if nt>0 else None
            rec=dict(sign=sign,E_requested_PeV=target,E0_grid_PeV=E0/1e6,
                     tau_survival_absolute_error=abs(nt-surv),
                     tau_number_relative_error=abs(nt-surv)/max(surv,1e-300),
                     W1_logE=w1,tau_L1=rel_l1(tau,rt),
                     primary_daughter_L1=rel_l1(prim,rp),
                     secondary_daughter_L1=rel_l1(sec,rs))
            rows.append(rec)
            print('TAU-MIXED '+json.dumps(rec),flush=True)

    worstp=max(r['primary_daughter_L1'] for r in rows)
    worsts=max(r['secondary_daughter_L1'] for r in rows)
    max_tau_abs=max(r['tau_survival_absolute_error'] for r in rows)

    checks={
        'finite': bool(np.all(np.isfinite(y))),
        'semigroup': semi <= GATES['semigroup_L1'],
        'cc_rowsum': ccmax <= GATES['cc_rowsum_maxrel'],
        'quadrature_multiplicity': multmax <= GATES['quadrature_multiplicity_abs'],
        'nonnegative': negfrac <= GATES['negative_mass_fraction'],
        'primary_daughters': worstp <= GATES['worst_primary_daughter_L1'],
        'secondary_daughters': worsts <= GATES['worst_secondary_daughter_L1'],
        'tau_survival_absolute': max_tau_abs <= GATES['max_tau_survival_abs_error'],
    }
    passed=all(checks.values())
    status=dict(version='v0.20.4-smoke',tau_refine=16,N_nu=Nn,N_tau=Nt,
                build_seconds=build,one_slab_seconds=slab,
                projected_serial_seconds=slab*nseg,
                semigroup_L1=semi,cc_rowsum_maxrel_nu=cc_n,
                cc_rowsum_maxrel_antinu=cc_a,
                quadrature_multiplicity_abs_errors=mult,
                negative_mass_fraction=negfrac,tau_benchmarks=rows,
                worst_primary_daughter_L1=worstp,
                worst_secondary_daughter_L1=worsts,
                max_tau_survival_absolute_error=max_tau_abs,
                gates=GATES,checks=checks,all_gates_PASS=passed,
                note='Smoke only. Passing does not certify full-chord density/neutrino-grid convergence or Case 5.')
    (OUT/'mixed_expm_smoke_status_v0_20_4.json').write_text(json.dumps(status,indent=2))
    print(f'SUMMARY worst_primary={worstp:.6g} worst_secondary={worsts:.6g} max_tau_abs={max_tau_abs:.6g}',flush=True)
    print('CHECKS '+json.dumps(checks,sort_keys=True),flush=True)
    print('ALL_GATES_PASS='+str(passed),flush=True)
    print('STATUS wrote',OUT/'mixed_expm_smoke_status_v0_20_4.json',flush=True)
    if not passed:
        sys.exit(2)

if __name__=='__main__':
    main()
