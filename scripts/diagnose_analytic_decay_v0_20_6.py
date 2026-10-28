#!/usr/bin/env python3
"""v0.20.6 exact-z decay diagnostic.

Goals
-----
1. Validate the closed-form projected decay row against an independent
   grid-boundary-aware Gauss-Legendre integration.
2. Remove inner z-quadrature uncertainty from the continuous-loss reference.
3. Check convergence of the remaining outer decay-probability integral.
4. With z projection exact, measure tau-grid advection/decay error separately
   for refinement factors 8, 16, 32.

No full Earth chord and no production certification.
"""
from __future__ import annotations

import json, math, time
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.sparse import bmat, csr_matrix, diags, lil_matrix
from scipy.sparse.linalg import expm_multiply

from nutau_prop.analytic_decay_projection import (
    projected_decay_row, projected_decay_kernel,
)
from nutau_prop.chirality_transport_mixed_expm import _refined_log_grid
from nutau_prop.kernels import tau_nutau_pdf, tau_secondary_numu_pdf, MTAU, MPI, MRHO, MA1
from nutau_prop.constants import M_TAU, C_TAU_CM

ROOT=Path(__file__).resolve().parents[1]
T15=ROOT/'results'/'v0.15'/'tables'
OUT=ROOT/'results'/'v0.20.6'
RHO=2.9
H=0.75e5
BETA=0.85e-6
BR=0.18
TARGETS=(1.,10.,100.,300.)
FACTORS=(8,16,32)
OUTER_ORDERS=(32,64,128,256)


def load_energy():
    z=np.load(T15/'case7_L3_nu.npz')
    return np.asarray(z['energies'],float)


def rel_l1(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    if not (np.all(np.isfinite(a)) and np.all(np.isfinite(b))):
        raise FloatingPointError('non-finite spectrum in L1 metric')
    den=float(np.sum(np.abs(b)))
    if not np.isfinite(den):
        raise FloatingPointError('non-finite L1 denominator')
    return float(np.sum(np.abs(a-b))/max(den,1e-300))


def physical_bounds():
    b=[0.0,0.3,1.0-(MA1/MTAU)**2,1.0-(MRHO/MTAU)**2,
       1.0-(MPI/MTAU)**2,1.0]
    return sorted(set(max(0.0,min(1.0,float(x))) for x in b))


def deposit(row,lnE,eout,w):
    if w==0.0 or eout<=0.0:
        return
    u=math.log(float(eout))
    if u<lnE[0] or u>lnE[-1]:
        return
    j=int(np.searchsorted(lnE,u))
    if j==0: row[0]+=w
    elif j>=len(lnE): row[-1]+=w
    else:
        f=(u-lnE[j-1])/(lnE[j]-lnE[j-1])
        row[j-1]+=w*(1.0-f); row[j]+=w*f


def grid_aware_gl_row(Ep,Enu,pol,secondary=False,order=32):
    """Independent numerical row split at PDF AND daughter-grid kinks."""
    cuts=list(physical_bounds())
    for e in Enu:
        z=float(e)/float(Ep)
        if 0.0<z<1.0:
            cuts.append(z)
    cuts=sorted(set(cuts))
    xg,wg=leggauss(int(order)); row=np.zeros(len(Enu)); lnE=np.log(Enu)
    for a,b in zip(cuts[:-1],cuts[1:]):
        if b<=a: continue
        z=0.5*(xg+1.0)*(b-a)+a
        jac=0.5*(b-a)
        f=(tau_secondary_numu_pdf(z,pol,BR) if secondary else tau_nutau_pdf(z,pol))
        for zz,ww,ff in zip(z,wg,f):
            deposit(row,lnE,Ep*float(zz),float(jac*ww*ff))
    # Physical formulas integrate to BR / 1; no post-normalization is applied.
    return row


def validate_analytic_rows(Enu):
    rows=[]; wp=0.0; ws=0.0
    for pol,label in ((-1.0,'minus'),(+1.0,'plus')):
        for target in TARGETS:
            j=int(np.argmin(np.abs(np.log(Enu)-np.log(target*1e6))))
            Ep=float(Enu[j])
            for secondary in (False,True):
                a=projected_decay_row(Ep,Enu,pol,secondary,BR)
                n=grid_aware_gl_row(Ep,Enu,pol,secondary,32)
                err=rel_l1(a,n)
                rows.append(dict(sign=label,E_PeV=target,secondary=secondary,L1=err,
                                 analytic_sum=float(np.sum(a)),numeric_sum=float(np.sum(n))))
                if secondary: ws=max(ws,err)
                else: wp=max(wp,err)
    return rows,wp,ws


def continuous_reference(E0,Enu,sign,outer_order):
    """Continuous tau loss/decay; exact analytic z projection at each decay energy."""
    pol=-1.0 if sign=='minus' else +1.0
    b=BETA*RHO; a=M_TAU/(C_TAU_CM*E0)
    U=(a/b)*np.expm1(b*H) if b!=0 else a*H
    pdec=-np.expm1(-min(float(U),745.0))
    p=np.zeros(len(Enu)); s=np.zeros(len(Enu))
    xg,wg=leggauss(int(outer_order))
    t=0.5*(xg+1.0)*pdec; w=0.5*pdec*wg
    u=-np.log1p(-t)
    Edec=E0/(1.0+(b/a)*u) if b!=0 else np.full_like(u,E0)
    for ed,ww in zip(Edec,w):
        p += float(ww)*projected_decay_row(float(ed),Enu,pol,False,BR)
        s += float(ww)*projected_decay_row(float(ed),Enu,pol,True,BR)
    surv=math.exp(-min(float(U),745.0))
    Eend=E0*math.exp(-b*H)
    return p,s,surv,Eend


def outer_convergence(Enu):
    refs={}; rows=[]
    for sign in ('minus','plus'):
        for target in TARGETS:
            j=int(np.argmin(np.abs(np.log(Enu)-np.log(target*1e6))))
            E0=float(Enu[j])
            by={}
            for q in OUTER_ORDERS:
                by[q]=continuous_reference(E0,Enu,sign,q)
            refs[(sign,target)]=(E0,*by[256])
            for q1,q2 in zip(OUTER_ORDERS[:-1],OUTER_ORDERS[1:]):
                ep=rel_l1(by[q1][0],by[q2][0]); es=rel_l1(by[q1][1],by[q2][1])
                rows.append(dict(sign=sign,E_PeV=target,order_lo=q1,order_hi=q2,
                                 primary_L1=ep,secondary_L1=es))
    return refs,rows


def tau_loss_generator(Etau):
    l=np.log(Etau); dl=np.diff(l)
    rate=np.empty(len(Etau)); rate[0]=BETA*RHO/dl[0]; rate[1:]=BETA*RHO/dl
    Q=lil_matrix((len(Etau),len(Etau)),dtype=float); Q.setdiag(-rate)
    for j in range(1,len(Etau)): Q[j,j-1]=rate[j]
    return Q.tocsr()


def tau_only_generator(Etau,Enu,pol):
    Nt,Nn=len(Etau),len(Enu)
    Kp=projected_decay_kernel(Etau,Enu,pol,False,BR)
    Ks=projected_decay_kernel(Etau,Enu,pol,True,BR)
    d=M_TAU/(C_TAU_CM*Etau)
    T=tau_loss_generator(Etau)-diags(d,0,format='csr')
    D=diags(d,0,format='csr')@Kp
    S=diags(d,0,format='csr')@Ks
    Znn=csr_matrix((Nn,Nn)); Znt=csr_matrix((Nn,Nt))
    return bmat([[Znn,Znn,Znt],[Znn,Znn,Znt],[D,S,T]],format='csr')


def factor_test(Enu,refs,factor):
    Etau=_refined_log_grid(Enu,factor)
    rows=[]; worstp=0.0; worsts=0.0; maxabs=0.0
    t0=time.perf_counter()
    gens={sign:tau_only_generator(Etau,Enu,-1.0 if sign=='minus' else +1.0)
          for sign in ('minus','plus')}
    build=time.perf_counter()-t0
    for sign in ('minus','plus'):
        G=gens[sign]
        for target in TARGETS:
            E0,rp,rs,surv,Eend=refs[(sign,target)]
            jt=int(np.argmin(np.abs(np.log(Etau)-np.log(E0))))
            x=np.zeros(2*len(Enu)+len(Etau)); x[2*len(Enu)+jt]=1.0
            y=np.asarray(expm_multiply(G.T*H,x),float)
            if not np.all(np.isfinite(y)):
                raise FloatingPointError(f'non-finite tau-grid evolution factor={factor}')
            p=y[:len(Enu)]; s=y[len(Enu):2*len(Enu)]; tau=y[2*len(Enu):]
            ep=rel_l1(p,rp); es=rel_l1(s,rs)
            nt=float(np.sum(tau)); abs_surv=abs(nt-surv)
            prob=tau/max(nt,1e-300)
            w1=float(np.sum(prob*np.abs(np.log(Etau)-math.log(Eend)))) if nt>0 else None
            rows.append(dict(factor=factor,sign=sign,E_PeV=target,primary_L1=ep,
                             secondary_L1=es,tau_survival_absolute_error=abs_surv,
                             W1_logE=w1))
            worstp=max(worstp,ep); worsts=max(worsts,es); maxabs=max(maxabs,abs_surv)
    return rows,worstp,worsts,maxabs,build,len(Etau)


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    Enu=load_energy()

    rowchecks,arp,ars=validate_analytic_rows(Enu)
    print(f'ANALYTIC-ROW vs grid-aware GL32 worst_primary={arp:.3e} worst_secondary={ars:.3e}',flush=True)

    refs,outer=outer_convergence(Enu)
    for q1,q2 in zip(OUTER_ORDERS[:-1],OUTER_ORDERS[1:]):
        subset=[r for r in outer if r['order_lo']==q1 and r['order_hi']==q2]
        wp=max(r['primary_L1'] for r in subset); ws=max(r['secondary_L1'] for r in subset)
        print(f'OUTER GL{q1}->GL{q2} worst_primary={wp:.3e} worst_secondary={ws:.3e}',flush=True)

    factor_rows=[]; summaries=[]
    for f in FACTORS:
        rr,wp,ws,ma,build,Nt=factor_test(Enu,refs,f)
        factor_rows.extend(rr)
        summaries.append(dict(factor=f,N_tau=Nt,worst_primary_L1=wp,
                              worst_secondary_L1=ws,max_tau_survival_abs=ma,
                              build_seconds=build))
        print(f'TAU-GRID factor={f:2d} Ntau={Nt:4d} build_sec={build:.2f} '
              f'worst_primary={wp:.6g} worst_secondary={ws:.6g} max_tau_abs={ma:.3e}',flush=True)

    status=dict(version='v0.20.6-diagnostic',N_nu=len(Enu),tau_factors=list(FACTORS),
                analytic_row_checks=rowchecks,analytic_row_worst_primary=arp,
                analytic_row_worst_secondary=ars,outer_convergence=outer,
                reference='continuous tau loss/decay; analytic z projection; outer GL256',
                tau_grid_rows=factor_rows,tau_grid_summaries=summaries,
                note='Diagnostic only. No full chord. NaN/nonfinite values raise immediately.')
    (OUT/'analytic_decay_status_v0_20_6.json').write_text(json.dumps(status,indent=2))
    print('STATUS wrote',OUT/'analytic_decay_status_v0_20_6.json',flush=True)


if __name__=='__main__':
    main()
