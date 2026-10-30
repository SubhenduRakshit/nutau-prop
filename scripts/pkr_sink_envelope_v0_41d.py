#!/usr/bin/env python3
from __future__ import annotations
import math,json
from pathlib import Path
import numpy as np,pandas as pd
from nutau_prop.case6 import Case6Config,Case6DIS
from nutau_prop.constants import M_N,GEV2_TO_CM2,NA
from nutau_prop.prem import chord_profile
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results'/'v0.41d';OUT.mkdir(parents=True,exist_ok=True)
PDF='CT18NLO';GQ=.03;M=50.;MODELS={'eq4p4_derivative':1.,'eq4p3_dipole':2.};E0=1e8;ENERGIES=np.array([1,3,10,30,100,300.])*1e6;MUS=(.1,.3,.5,.7,1.)
class R(Case6DIS):
 def dsnp(self,E,y,k,f):
  s=2*M_N*E;xmin=max(self.cfg.qmin_gev**2/(s*y),1e-9);xmax=min(1-1e-10,(GQ/k)**2/s)
  if xmin>=xmax:return 0.
  def z(x):
   q2=x*y*s
   if q2<self.cfg.qmin_gev**2:return 0.
   p=self._pdf_components(x,q2);c=k/M**2
   return f*(4/math.pi)*c*c*x*s*s*y*(1-y)*sum(p.values())*(M*M/(M*M+q2))**2*GEV2_TO_CM2
  return self._integrate_logx(z,xmin,xmax)
C=R(Case6Config(pdf_set=PDF,nx=32,y_log_bins=12,y_linear_bins=12))
def sigmas(E,k,f):
 _,y,dy=C.y_grid();snp=ssm=0.
 for yy,d in zip(y,dy):
  cc,nc=C.dsdy_sm(E,yy);snp+=C.dsnp(E,yy,k,f)*d;ssm+=(cc+nc)*d
 return snp,ssm
def solve(f):
 ks=np.geomspace(2e-6,2e-3,48);rs=[]
 for k in ks:
  a,b=sigmas(E0,k,f);rs.append(a/b)
 rs=np.array(rs);ii=np.where((rs[:-1]<.05)&(rs[1:]>=.05))[0]
 if not len(ii):j=int(np.argmax(rs));return None,float(ks[j]),float(rs[j])
 lo,hi=ks[ii[0]],ks[ii[0]+1]
 for _ in range(22):
  md=math.sqrt(lo*hi);a,b=sigmas(E0,md,f)
  if a/b>=.05:hi=md
  else:lo=md
 a,b=sigmas(E0,hi,f);return float(hi),float(hi),float(a/b)
def column(mu):
 p=chord_profile(mu,max_step_km=12.);return float(np.sum(np.asarray(p['rho'])*np.asarray(p['ds_km'])*1e5))
def main():
 gates=[];rows=[]
 cols={mu:column(mu) for mu in MUS}
 for model,f in MODELS.items():
  k,kp,rp=solve(f);gates.append(dict(model=model,M_GeV=M,gq=GQ,kappa5=k,vertex_scale_GeV=GQ/k if k else np.nan,peak_kappa=kp,peak_ratio=rp,status='OPEN_5PCT' if k else 'CLOSED'))
  if not k:continue
  for E in ENERGIES:
   snp,ssm=sigmas(E,k,f)
   for mu in MUS:
    tau=NA*cols[mu]*snp
    rows.append(dict(model=model,M_GeV=M,gq=GQ,kappa=k,E_PeV=E/1e6,mu=mu,column_g_cm2=cols[mu],sigma_np_cm2=snp,sigma_sm_cm2=ssm,np_sm_ratio=snp/ssm,tau_np=tau,pure_sink_survival=math.exp(-tau),pure_sink_frac=math.exp(-tau)-1))
 pd.DataFrame(gates).to_csv(OUT/'benchmark_gate.csv',index=False);pd.DataFrame(rows).to_csv(OUT/'sink_envelope.csv',index=False)
 (OUT/'metadata.json').write_text(json.dumps(dict(pdf=PDF,gq=GQ,M_GeV=M,validity='sqrt(shat)<gq/kappa',meaning='pure BSM sterile-sink optical depth; SM regeneration not included'),indent=2))
 print(pd.DataFrame(gates).to_string(index=False));print(pd.DataFrame(rows).query('mu==1.0').to_string(index=False))
if __name__=='__main__':main()
