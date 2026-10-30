#!/usr/bin/env python3
"""v0.38 split-coupling validity screen for PKR Eq.(4.3)/(4.4).

For fixed physical kappa=c*gq/Lambda, impose the actual dimension-5 vertex
validity bound sqrt(shat)<Lambda/c=gq/kappa rather than the gq=1 proxy 1/kappa.
Evaluate the previously determined 100-PeV 5% kappas for representative gq.
"""
import csv,math
from pathlib import Path
import numpy as np
from nutau_prop.case6 import Case6Config,Case6DIS
from nutau_prop.constants import M_N,GEV2_TO_CM2
PEV=1e6
KD={5:3.0943e-5,40:5.4254e-5,100:7.9706e-5,220:1.2008e-4,500:2.1643e-4}
KM={5:2.1880e-5,40:3.8364e-5,100:5.6361e-5,220:8.4594e-5,500:1.4390e-4,1000:3.3001e-4}
GQS=[1.0,0.3,0.1,0.07,0.04]

def integ(a,dy):return float(np.sum(np.asarray(a)*np.asarray(dy)))
def sm(c,E):
 _,y,dy=c.y_grid();v=[]
 for z in y:
  a,b=c.dsdy_sm(E,z,False);d,e=c.dsdy_sm(E,z,True);v.append(.5*(a+b+d+e))
 return integ(v,dy)
def np_sig(c,E,M,k,gq,factor):
 _,ys,dy=c.y_grid();s=2*M_N*E; m2=M*M; lamtev=(gq/k)/1e3; vals=[]
 for y in ys:
  b=c._x_bounds(E,y,lambda_tev=lamtev,cut='shat')
  if b is None: vals.append(0.); continue
  xmin,xmax=b
  def f(x):
   q2=x*y*s
   if q2<c.cfg.qmin_gev**2:return 0.
   sf=sum(c._pdf_components(x,q2).values())
   return factor*(4/math.pi)*k*k*x*s*s*y*(1-y)*sf/(m2+q2)**2*GEV2_TO_CM2
  vals.append(c._integrate_logx(f,xmin,xmax))
 return integ(vals,dy)
def main():
 out=Path('results/v0.38');out.mkdir(parents=True,exist_ok=True)
 c=Case6DIS(Case6Config(pdf_set='CT18NLO',nx=32));rows=[]
 for model,tab,factor in [('derivative_eq4p4',KD,1.),('dipole_eq4p3',KM,2.)]:
  for M,k in tab.items():
   for g in GQS:
    for ep in [100.,300.]:
     E=ep*PEV;r=np_sig(c,E,M,k,g,factor)/sm(c,E)
     rows.append(dict(model=model,M_GeV=M,kappa=k,gq=g,vertex_scale_GeV=g/k,E_PeV=ep,ratio=r,percent=100*r))
 with open(out/'split_coupling_validity_v0_38.csv','w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
 lines=['# v0.38 split-coupling validity screen','','Actual validity scale is Lambda/c = gq/kappa. Percentages below use each model/mass previous 100-PeV kappa5 benchmark.','']
 for model in ['derivative_eq4p4','dipole_eq4p3']:
  lines += [f'## {model}','', '| M | gq | vertex scale | 100 PeV | 300 PeV |','|---:|---:|---:|---:|---:|']
  for r in rows:
   if r['model']!=model or r['E_PeV']!=100.:continue
   q=[z for z in rows if z['model']==model and z['M_GeV']==r['M_GeV'] and z['gq']==r['gq'] and z['E_PeV']==300.][0]
   lines.append(f"| {r['M_GeV']:g} | {r['gq']:.2g} | {r['vertex_scale_GeV']:.1f} GeV | {r['percent']:.3f}% | {q['percent']:.3f}% |")
  lines.append('')
 (out/'split_coupling_validity_v0_38.md').write_text('\n'.join(lines)+'\n')
 print('wrote',out)
if __name__=='__main__':main()
