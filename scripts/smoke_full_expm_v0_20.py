#!/usr/bin/env python3
"""v0.20 combined-generator smoke + tau benchmark.

No production run.  This deliberately benchmarks ONE representative 3-km slab
before any full chord is attempted, in accordance with the project runtime rule.
It also checks exact discrete semigroup composition and compares the new tau
finite-volume generator against the analytic one-slab CoupledTauEvolution map.
"""
from __future__ import annotations

import json, time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.sparse import bmat, csr_matrix

from nutau_prop.chirality_transport_full_expm import TauHelicityFlipTransportFullExpm
from nutau_prop.tau_evolution import CoupledTauEvolution
from nutau_prop.prem import chord_profile

ROOT=Path(__file__).resolve().parents[1]
T15=ROOT/'results'/'v0.15'/'tables'
UHE=ROOT/'data'/'uhe_official'
OUT=ROOT/'results'/'v0.20'


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


def tau_only_generator(tr,rho,secondary=True):
    N=tr.N; Z=csr_matrix((N,N))
    L=float(rho)*tr._tau_loss_per_rho()
    Q,Dm,Dp,Sm,Sp=tr._tau_decay_blocks(secondary)
    return bmat([
        [Z,Z,Z,Z,Z,Z], [Z,Z,Z,Z,Z,Z],
        [Z,Z,Z,Z,Z,Z], [Z,Z,Z,Z,Z,Z],
        [Dm,Z,Z,Sm,L+Q,Z],
        [Z,Dp,Sp,Z,Z,L+Q],
    ],format='csr')


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    tn,ta=load('nu'),load('antinu')
    tr=TauHelicityFlipTransportFullExpm(tn,ta)
    E=tr.E; N=tr.N; p=flux(E)
    prof=chord_profile(0.1,max_step_km=3.0)
    rho=float(np.median(prof['rho']))
    nseg=len(prof['ds_km'])

    # Runtime benchmark: one representative slab only.
    x=np.concatenate([p,p,np.zeros(4*N)])
    G=tr.generator(rho,'majorana_flip',True,True)
    t0=time.perf_counter(); y=tr._advance(G,x,3.0e5); sec=time.perf_counter()-t0
    print(f'SLAB 3km rho={rho:.6g} sec={sec:.3f} nsegments_3km={nseg} projected_serial_sec={sec*nseg:.1f}',flush=True)
    print(f'SLAB finite={np.all(np.isfinite(y))} min={np.min(y):.6g}',flush=True)

    # A true matrix exponential must compose as a semigroup to roundoff.
    h=0.75e5
    one=tr._advance(G,x,h)
    two=tr._advance(G,tr._advance(G,x,0.5*h),0.5*h)
    semi=rel_l1(two,one)
    print(f'FULL-GENERATOR semigroup_L1={semi:.6g}',flush=True)

    # Tau-only one-slab comparison against the analytic continuous solution.
    Te=CoupledTauEvolution(E,beta_tau=tr.beta)
    Gt=tau_only_generator(tr,rho,True)
    tau_rows=[]
    for epev in (1.0,10.0,100.0,300.0):
        j=int(np.argmin(np.abs(np.log(E)-np.log(epev*1e6))))
        tm=np.zeros(N); tm[j]=1.0
        z=np.zeros(N)
        xx=np.concatenate([z,z,z,z,tm,z])
        yy=tr._advance(Gt,xx,h)
        n,a,mn,ma,tm_new,tp=[yy[i*N:(i+1)*N] for i in range(6)]
        tau_exact,n_exact,ma_exact,_=Te.advance_with_kernels(tm,rho,h,tr.Kdec_minus,tr.Ksec_minus)
        rec=dict(E_requested_PeV=epev,E_grid_PeV=float(E[j]/1e6),
                 tau_L1=rel_l1(tm_new,tau_exact),
                 primary_daughter_L1=rel_l1(n,n_exact),
                 secondary_daughter_L1=rel_l1(ma,ma_exact),
                 tau_number_new=float(np.sum(tm_new)),tau_number_exact=float(np.sum(tau_exact)))
        tau_rows.append(rec)
        print('TAU-BENCH '+json.dumps(rec),flush=True)

    status=dict(version='v0.20-smoke',case='case7_L3_cosz0.1_majorana_flip',
                one_slab_seconds=sec,nsegments_3km=nseg,projected_serial_seconds=sec*nseg,
                semigroup_L1=semi,tau_one_slab_benchmarks=tau_rows,
                note='Smoke/diagnostic only. Do not run full production until runtime and energy-grid accuracy are assessed.')
    (OUT/'full_expm_smoke_status_v0_20.json').write_text(json.dumps(status,indent=2))
    print('STATUS wrote',OUT/'full_expm_smoke_status_v0_20.json',flush=True)

if __name__=='__main__': main()
