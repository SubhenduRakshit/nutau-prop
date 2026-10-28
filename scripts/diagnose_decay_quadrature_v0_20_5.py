#!/usr/bin/env python3
"""v0.20.5 decay-z quadrature diagnostic for mixed-grid transport.

Hold the neutrino grid (91) and tau grid (factor 16, 1441) fixed.  Vary only
how the tau-decay z integral is discretised:
  trapezoid nz=800,1600,3200
  split Gauss-Legendre order=16,32,64
and compare the resulting one-slab daughter spectra against an independent
higher-order continuous-loss reference using split GL(256) in z and GL(128)
in decay optical depth.

Diagnostic only: no Earth chord and no production certification.
"""
from __future__ import annotations

import json, math, time
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.sparse import bmat, csr_matrix, diags

from nutau_prop.chirality_transport_mixed_expm import (
    TauHelicityFlipTransportMixedExpm, _deposit_row
)
from nutau_prop.kernels import (
    tau_nutau_pdf, tau_secondary_numu_pdf, normalized_trapezoid_weights,
    MTAU, MPI, MRHO, MA1,
)
from nutau_prop.constants import M_TAU, C_TAU_CM

ROOT=Path(__file__).resolve().parents[1]
T15=ROOT/'results'/'v0.15'/'tables'
OUT=ROOT/'results'/'v0.20.5'
RHO=2.9
H=0.75e5
BR=0.18
TARGETS=(1.,10.,100.,300.)
VARIANTS=(('trap',800),('trap',1600),('trap',3200),
          ('gl',16),('gl',32),('gl',64))


def load(spec):
    z=np.load(T15/f'case7_L3_{spec}.npz')
    return {k:z[k] for k in z.files}


def ref_deposit(row, log_grid, eout, weight):
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
        f=(u-log_grid[j-1])/(log_grid[j]-log_grid[j-1])
        row[j-1]+=weight*(1-f); row[j]+=weight*f


def split_bounds():
    b=[0.0,0.3,
       1.0-(MA1/MTAU)**2,
       1.0-(MRHO/MTAU)**2,
       1.0-(MPI/MTAU)**2,1.0]
    return sorted(set(max(0.0,min(1.0,float(x))) for x in b))


def gl_nodes_weights(pol, secondary, order):
    xg,wg=leggauss(int(order)); zs=[]; ws=[]
    for a,b in zip(split_bounds()[:-1],split_bounds()[1:]):
        if b<=a: continue
        z=0.5*(xg+1.0)*(b-a)+a
        jac=0.5*(b-a)
        f=(tau_secondary_numu_pdf(z,pol,BR) if secondary
           else tau_nutau_pdf(z,pol))
        zs.append(z); ws.append(jac*wg*f)
    z=np.concatenate(zs); w=np.concatenate(ws)
    target=BR if secondary else 1.0
    s=float(np.sum(w))
    if s<=0: raise RuntimeError('zero GL normalization')
    w*=target/s
    return z,w


def trap_nodes_weights(pol, secondary, nz):
    z=np.linspace(1e-6,1-1e-6,int(nz))
    f=(tau_secondary_numu_pdf(z,pol,BR) if secondary
       else tau_nutau_pdf(z,pol))
    target=BR if secondary else 1.0
    w=normalized_trapezoid_weights(z,f,target)
    return z,w


def rect_kernel(Etau, Enu, z, w):
    ln=np.log(np.asarray(Enu,float))
    K=np.zeros((len(Etau),len(Enu)),float)
    for j,E in enumerate(Etau):
        row=K[j]
        for zz,ww in zip(z,w):
            _deposit_row(row,ln,E*float(zz),float(ww))
    return csr_matrix(K)


def make_kernels(tr, kind, n):
    out=[]
    for pol,secondary in ((-1.0,False),(+1.0,False),(-1.0,True),(+1.0,True)):
        if kind=='trap': z,w=trap_nodes_weights(pol,secondary,n)
        else: z,w=gl_nodes_weights(pol,secondary,n)
        out.append(rect_kernel(tr.Etau,tr.E,z,w))
    return out


def tau_generator(tr,sign):
    Nn,Nt=tr.Nn,tr.Nt
    Znn=csr_matrix((Nn,Nn)); Znt=csr_matrix((Nn,Nt))
    L=RHO*tr._tau_loss_per_rho()
    d=M_TAU/(C_TAU_CM*tr.Etau)
    Q=-diags(d,0,format='csr')
    if sign=='minus': Kp,Ks=tr.Kdec_minus,tr.Ksec_minus
    else: Kp,Ks=tr.Kdec_plus,tr.Ksec_plus
    D=diags(d,0,format='csr')@Kp
    S=diags(d,0,format='csr')@Ks
    return bmat([[Znn,Znn,Znt],[Znn,Znn,Znt],[D,S,L+Q]],format='csr')


def exact_ref(E0,Enu,Etau,sign):
    pol=-1.0 if sign=='minus' else +1.0
    b=0.85e-6*RHO
    a=M_TAU/(C_TAU_CM*E0)
    U=(a/b)*np.expm1(b*H)
    pdec=-np.expm1(-min(float(U),745.0))
    p=np.zeros(len(Enu)); s=np.zeros(len(Enu)); ln=np.log(Enu)
    qx,qw=leggauss(128)
    tt=0.5*(qx+1.0)*pdec; ww=0.5*pdec*qw
    uu=-np.log1p(-tt); Edec=E0/(1.0+(b/a)*uu)
    zp,wp=gl_nodes_weights(pol,False,256)
    zs,ws=gl_nodes_weights(pol,True,256)
    for ed,wd in zip(Edec,ww):
        for zz,wz in zip(zp,wp):
            ref_deposit(p,ln,float(ed*zz),float(wd*wz))
        for zz,wz in zip(zs,ws):
            ref_deposit(s,ln,float(ed*zz),float(wd*wz))
    return p,s


def rel_l1(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    return float(np.sum(np.abs(a-b))/max(np.sum(np.abs(b)),1e-300))


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    tn,ta=load('nu'),load('antinu')
    tr=TauHelicityFlipTransportMixedExpm(tn,ta,tau_refine=16)
    refs={}
    for sign in ('minus','plus'):
        for target in TARGETS:
            j=int(np.argmin(np.abs(np.log(tr.E)-np.log(target*1e6))))
            E0=float(tr.E[j])
            refs[(sign,target)]=(E0,*exact_ref(E0,tr.E,tr.Etau,sign))

    rows=[]
    for kind,n in VARIANTS:
        t0=time.perf_counter()
        km,kp,sm,sp=make_kernels(tr,kind,n)
        build=time.perf_counter()-t0
        tr.Kdec_minus, tr.Kdec_plus = km,kp
        tr.Ksec_minus, tr.Ksec_plus = sm,sp
        worstp=0.0; worsts=0.0
        for sign in ('minus','plus'):
            G=tau_generator(tr,sign)
            for target in TARGETS:
                E0,rp,rs=refs[(sign,target)]
                jt=int(np.argmin(np.abs(np.log(tr.Etau)-np.log(E0))))
                x=np.zeros(2*tr.Nn+tr.Nt); x[2*tr.Nn+jt]=1.0
                y=np.asarray(tr._advance(G,x,H),float)
                p=y[:tr.Nn]; s=y[tr.Nn:2*tr.Nn]
                ep=rel_l1(p,rp); es=rel_l1(s,rs)
                worstp=max(worstp,ep); worsts=max(worsts,es)
                rows.append(dict(kind=kind,n=int(n),sign=sign,E_PeV=target,
                                 primary_L1=ep,secondary_L1=es,
                                 build_seconds=build))
        print(f'QUAD {kind}{n} build_sec={build:.2f} worst_primary={worstp:.6g} worst_secondary={worsts:.6g}',flush=True)

    status=dict(version='v0.20.5-diagnostic',tau_refine=16,N_nu=tr.Nn,N_tau=tr.Nt,
                rho_g_cm3=RHO,slab_km=H/1e5,variants=[list(v) for v in VARIANTS],
                reference='split GL256 in z + GL128 in decay optical depth',rows=rows,
                note='Diagnostic only. Choose decay quadrature only after observing convergence; no full chord.')
    (OUT/'decay_quadrature_status_v0_20_5.json').write_text(json.dumps(status,indent=2))
    print('STATUS wrote',OUT/'decay_quadrature_status_v0_20_5.json',flush=True)

if __name__=='__main__':
    main()
