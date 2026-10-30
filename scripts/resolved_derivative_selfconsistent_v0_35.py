#!/usr/bin/env python3
"""v0.35 self-consistent resolved derivative-Z' UHE screen.

For every trial kappa = c2*gq/Lambda [GeV^-1], enforce the derivative-vertex
validity domain sqrt(shat) < 1/kappa INSIDE the PDF fold.  No fixed-cut
quadratic rescaling is used.

Kernel convention follows v0.32:
 d2sigma/dxdy = (4/pi) kappa^2 x s^2 y(1-y) sum_q[x f_q]
                 /(M_Z'^2 + Q^2)^2.

The code solves for the first kappa for which sigma_NP/sigma_SM = target,
using a logarithmic scan followed by log-space bisection.  If the function
never reaches target within the requested kappa range, no root is reported.
"""
from __future__ import annotations

import argparse, csv, json, math
from pathlib import Path
import numpy as np

from nutau_prop.case6 import Case6Config, Case6DIS
from nutau_prop.constants import M_N, GEV2_TO_CM2

PEV_TO_GEV = 1.0e6


def integrate_y(arr, dy):
    return float(np.sum(np.asarray(arr, dtype=float) * np.asarray(dy, dtype=float)))


def sm_total(calc: Case6DIS, E):
    _, y, dy = calc.y_grid()
    vals=[]
    for yy in y:
        cc1,nc1=calc.dsdy_sm(E,yy,antinu=False)
        cc2,nc2=calc.dsdy_sm(E,yy,antinu=True)
        vals.append(0.5*(cc1+nc1+cc2+nc2))
    return integrate_y(vals,dy)


def resolved_dsdy(calc,E,yval,M_TeV,kappa):
    # Self-consistent c=1 derivative-vertex domain: sqrt(shat) < 1/kappa.
    lambda_vertex_tev = 1.0/(float(kappa)*1.0e3)
    bounds=calc._x_bounds(E,yval,lambda_tev=lambda_vertex_tev,cut="shat")
    if bounds is None:
        return 0.0
    xmin,xmax=bounds
    s=2.0*M_N*E
    m2=(1e3*float(M_TeV))**2
    kap=float(kappa)
    def f(x):
        q2=x*yval*s
        if q2<calc.cfg.qmin_gev**2:
            return 0.0
        p=calc._pdf_components(x,q2)
        sumxf=sum(p.values())
        d2=(4.0/math.pi)*kap*kap*x*s*s*yval*(1.0-yval)*sumxf/(m2+q2)**2
        return d2*GEV2_TO_CM2
    return calc._integrate_logx(f,xmin,xmax)


def resolved_sigma(calc,E,M_TeV,kappa):
    _,y,dy=calc.y_grid()
    vals=[resolved_dsdy(calc,E,yy,M_TeV,kappa) for yy in y]
    sig=integrate_y(vals,dy)
    if sig<=0:
        return 0.0,float('nan')
    my=float(np.sum(np.asarray(vals)*np.asarray(y)*np.asarray(dy))/sig)
    return sig,my


def ratio(calc,E,M_TeV,kappa,sm):
    sig,my=resolved_sigma(calc,E,M_TeV,kappa)
    return sig/sm,my,sig


def first_target_root(calc,E,M_TeV,sm,target,kmin,kmax,nscan=36,niter=18):
    ks=np.geomspace(kmin,kmax,nscan)
    vals=[]
    prev=None
    bracket=None
    peak=(-1.0,None)
    for k in ks:
        r,my,sig=ratio(calc,E,M_TeV,float(k),sm)
        vals.append((float(k),r,my,sig))
        if r>peak[0]: peak=(r,float(k))
        if prev is not None and prev[1] < target <= r:
            bracket=(prev[0],float(k))
            break
        prev=(float(k),r)
    if bracket is None:
        return None,vals,peak
    lo,hi=bracket
    # Find the first crossing, bisecting in log k while recomputing the moving cutoff.
    for _ in range(niter):
        mid=math.sqrt(lo*hi)
        r,_,_=ratio(calc,E,M_TeV,mid,sm)
        if r>=target: hi=mid
        else: lo=mid
    kroot=hi
    r,my,sig=ratio(calc,E,M_TeV,kroot,sm)
    return dict(kappa_GeV_m1=kroot,ratio=r,mean_y=my,sigma_NP_cm2=sig,
                lambda_vertex_TeV=1.0/(kroot*1e3)),vals,peak


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--pdfs',default='CT18NLO,NNPDF31_nnlo_as_0118')
    ap.add_argument('--energies-pev',default='100,300')
    ap.add_argument('--masses-gev',default='5,10,20,40,70,100,150,220,300,500,750,1000,1500,2000')
    ap.add_argument('--target',type=float,default=0.05)
    ap.add_argument('--kmin',type=float,default=1e-6)
    ap.add_argument('--kmax',type=float,default=5e-3)
    ap.add_argument('--nx',type=int,default=112)
    ap.add_argument('--outdir',default='results/v0.35')
    a=ap.parse_args()
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    pdfs=[x.strip() for x in a.pdfs.split(',') if x.strip()]
    energies=[float(x) for x in a.energies_pev.split(',')]
    masses=[float(x) for x in a.masses_gev.split(',')]
    rows=[]; scans=[]
    for pdf in pdfs:
        calc=Case6DIS(Case6Config(pdf_set=pdf,nx=a.nx))
        for ep in energies:
            E=ep*PEV_TO_GEV
            sm=sm_total(calc,E)
            for mg in masses:
                mt=mg/1e3
                root,scan,peak=first_target_root(calc,E,mt,sm,a.target,a.kmin,a.kmax)
                for k,r,my,sig in scan:
                    scans.append(dict(pdf=pdf,E_PeV=ep,M_GeV=mg,kappa_GeV_m1=k,
                                      lambda_vertex_TeV=1.0/(k*1e3),ratio=r,
                                      mean_y=my,sigma_NP_cm2=sig,sigma_SM_cm2=sm))
                row=dict(pdf=pdf,E_PeV=ep,M_GeV=mg,target=a.target,sigma_SM_cm2=sm,
                         peak_ratio_scan=peak[0],peak_kappa_scan=peak[1] or 0.0)
                if root:
                    row.update(root); row['status']='root_found'
                else:
                    row.update(kappa_GeV_m1=0.0,ratio=0.0,mean_y=float('nan'),sigma_NP_cm2=0.0,
                               lambda_vertex_TeV=0.0,status='no_first_crossing')
                rows.append(row)
    with open(out/'resolved_derivative_kappa5_selfconsistent_v0_35.csv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    with open(out/'resolved_derivative_scan_selfconsistent_v0_35.csv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(scans[0].keys())); w.writeheader(); w.writerows(scans)
    with open(out/'resolved_derivative_kappa5_selfconsistent_v0_35.json','w') as f:
        json.dump({'meta':vars(a),'rows':rows},f,indent=2)
    lines=['# v0.35 self-consistent resolved derivative-Zprime summary','',
           'Every trial kappa recomputes the domain sqrt(shat)<1/kappa inside the PDF fold.','']
    for pdf in pdfs:
        lines += [f'## {pdf}','', '| M [GeV] | kappa_5%(100 PeV) | kappa_5%(300 PeV) |',
                  '|---:|---:|---:|']
        for mg in masses:
            vals=[]
            for ep in (100.0,300.0):
                rr=[r for r in rows if r['pdf']==pdf and r['E_PeV']==ep and r['M_GeV']==mg]
                if rr and rr[0]['status']=='root_found': vals.append(f"{rr[0]['kappa_GeV_m1']:.4e}")
                else: vals.append('none')
            lines.append(f'| {mg:g} | {vals[0]} | {vals[1]} |')
        lines.append('')
    (out/'resolved_derivative_kappa5_summary_v0_35.md').write_text('\n'.join(lines)+'\n')
    print('wrote',out,'roots',sum(r['status']=='root_found' for r in rows),'/',len(rows))

if __name__=='__main__':
    main()
