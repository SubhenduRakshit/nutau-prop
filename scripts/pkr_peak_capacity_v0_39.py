#!/usr/bin/env python3
"""v0.39 peak-capacity scan for PKR resolved kernels.
Scan the gq=1 derivative kernel over k'=1/Lambda_vert. Because for fixed gq,
physical ratio = factor*gq^2*R_derivative(k'), this single scan yields the
absolute maximum attainable ratio for any physical kappa at each gq.
"""
import csv,math
from pathlib import Path
import numpy as np
from nutau_prop.case6 import Case6Config,Case6DIS
from nutau_prop.constants import M_N,GEV2_TO_CM2
PEV=1e6; MASSES=[5,40,100,220,500,1000,2000]; GQS=[0.04,0.07,0.1,0.3,1.0]
def integ(a,d):return float(np.sum(np.asarray(a)*np.asarray(d)))
def sm(c,E):
 _,y,dy=c.y_grid();v=[]
 for z in y:
  a,b=c.dsdy_sm(E,z,False);d,e=c.dsdy_sm(E,z,True);v.append(.5*(a+b+d+e))
 return integ(v,dy)
def ratio(c,E,M,k):
 _,ys,dy=c.y_grid();s=2*M_N*E;m2=M*M;lam=1/(k*1e3);v=[]
 for y in ys:
  b=c._x_bounds(E,y,lambda_tev=lam,cut='shat')
  if b is None:v.append(0.);continue
  def f(x):
   q2=x*y*s
   if q2<c.cfg.qmin_gev**2:return 0.
   sf=sum(c._pdf_components(x,q2).values())
   return (4/math.pi)*k*k*x*s*s*y*(1-y)*sf/(m2+q2)**2*GEV2_TO_CM2
  v.append(c._integrate_logx(f,*b))
 return integ(v,dy)/sm(c,E)
def main():
 out=Path('results/v0.39');out.mkdir(parents=True,exist_ok=True);c=Case6DIS(Case6Config(pdf_set='CT18NLO',nx=24));E=100*PEV
 ks=np.geomspace(1e-6,1e-2,26);rows=[]
 for M in MASSES:
  vals=[(float(k),ratio(c,E,M,float(k))) for k in ks];kp,rp=max(vals,key=lambda z:z[1])
  for g in GQS:
   for model,factor in [('derivative_eq4p4',1.),('dipole_eq4p3',2.)]:
    rows.append(dict(model=model,M_GeV=M,gq=g,peak_percent=100*factor*g*g*rp,kprime_peak=kp,physical_kappa_peak=g*kp,vertex_scale_GeV=1/kp))
 with open(out/'pkr_peak_capacity_v0_39.csv','w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
 lines=['# v0.39 PKR peak-capacity screen','','For fixed gq, max physical ratio = factor*gq^2*max R(kprime), kprime=kappa/gq.','']
 for model in ['derivative_eq4p4','dipole_eq4p3']:
  lines += [f'## {model}','', '| M | gq=.04 | gq=.07 | gq=.10 | gq=.30 | gq=1 |','|---:|---:|---:|---:|---:|---:|']
  for M in MASSES:
   d={r['gq']:r['peak_percent'] for r in rows if r['model']==model and r['M_GeV']==M}
   lines.append(f"| {M:g} | {d[.04]:.3f}% | {d[.07]:.3f}% | {d[.1]:.3f}% | {d[.3]:.3f}% | {d[1.]:.2f}% |")
  lines.append('')
 (out/'pkr_peak_capacity_v0_39.md').write_text('\n'.join(lines)+'\n')
 print('wrote',out)
if __name__=='__main__':main()
