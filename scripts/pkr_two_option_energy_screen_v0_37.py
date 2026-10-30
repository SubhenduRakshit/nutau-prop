#!/usr/bin/env python3
"""v0.37 low-energy IceCube screen for PKR Eq.(4.4) derivative and Eq.(4.3) dipole Z'.

Uses the validated v0.35 derivative kernel and the exact Eq.(4.3)=2*Eq.(4.4)
resolved-kernel relation at equal kappa, with sqrt(shat)<1/kappa imposed inside
every fold. For each mediator mass, evaluate fixed benchmark kappas chosen as the
100-PeV self-consistent 5% thresholds from v0.35/v0.36 across 0.1--300 PeV.
"""
from __future__ import annotations
import csv, math
from pathlib import Path
import numpy as np
from nutau_prop.case6 import Case6Config, Case6DIS
from nutau_prop.constants import M_N, GEV2_TO_CM2

PEV_TO_GEV=1e6
K_D={5:3.0943e-5,40:5.4254e-5,100:7.9706e-5,220:1.2008e-4,500:2.1643e-4}
K_M={5:2.1880e-5,15:2.7932e-5,40:3.8364e-5,100:5.6361e-5,220:8.4594e-5,500:1.4390e-4,1000:3.3001e-4}
ENERGIES=[0.1,0.3,1,3,10,30,100,300]

def integrate_y(a,dy): return float(np.sum(np.asarray(a)*np.asarray(dy)))
def sm_total(c,E):
    _,y,dy=c.y_grid(); vals=[]
    for yy in y:
        cc1,nc1=c.dsdy_sm(E,yy,False); cc2,nc2=c.dsdy_sm(E,yy,True)
        vals.append(.5*(cc1+nc1+cc2+nc2))
    return integrate_y(vals,dy)

def dsdy(c,E,y,MGeV,kappa,factor):
    lam=1/(kappa*1e3)
    b=c._x_bounds(E,y,lambda_tev=lam,cut='shat')
    if b is None:return 0.
    xmin,xmax=b; s=2*M_N*E; m2=MGeV*MGeV
    def f(x):
        q2=x*y*s
        if q2<c.cfg.qmin_gev**2:return 0.
        sumxf=sum(c._pdf_components(x,q2).values())
        return factor*(4/math.pi)*kappa*kappa*x*s*s*y*(1-y)*sumxf/(m2+q2)**2*GEV2_TO_CM2
    return c._integrate_logx(f,xmin,xmax)

def ratio(c,E,M,k,factor):
    _,y,dy=c.y_grid(); vals=[dsdy(c,E,yy,M,k,factor) for yy in y]
    return integrate_y(vals,dy)/sm_total(c,E)

def main():
    out=Path('results/v0.37');out.mkdir(parents=True,exist_ok=True)
    c=Case6DIS(Case6Config(pdf_set='CT18NLO',nx=48))
    rows=[]
    for name,table,factor in [('derivative_eq4p4',K_D,1.0),('dipole_eq4p3',K_M,2.0)]:
        for M,k in table.items():
            for ep in ENERGIES:
                r=ratio(c,ep*PEV_TO_GEV,M,k,factor)
                rows.append(dict(model=name,M_GeV=M,kappa=k,E_PeV=ep,ratio=r,percent=100*r))
    with open(out/'pkr_two_option_energy_screen_v0_37.csv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
    lines=['# v0.37 PKR two-option energy screen','','Fixed kappa = each model/mass 100-PeV self-consistent 5% threshold.','']
    for name in ['derivative_eq4p4','dipole_eq4p3']:
        lines += [f'## {name}','', '| M [GeV] | 0.1 PeV | 1 PeV | 3 PeV | 10 PeV | 100 PeV | 300 PeV |','|---:|---:|---:|---:|---:|---:|---:|']
        masses=sorted(set(r['M_GeV'] for r in rows if r['model']==name))
        for M in masses:
            d={r['E_PeV']:r['percent'] for r in rows if r['model']==name and r['M_GeV']==M}
            lines.append(f"| {M:g} | {d[0.1]:.3f}% | {d[1]:.3f}% | {d[3]:.3f}% | {d[10]:.3f}% | {d[100]:.3f}% | {d[300]:.3f}% |")
        lines.append('')
    (out/'pkr_two_option_energy_screen_v0_37.md').write_text('\n'.join(lines)+'\n')
    print('wrote',out)
if __name__=='__main__': main()
