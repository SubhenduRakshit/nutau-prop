#!/usr/bin/env python3
from __future__ import annotations
import json, math
from pathlib import Path
import numpy as np
import pandas as pd
from nutau_prop.case6 import Case6Config, Case6DIS, probability_from_dsdy
from nutau_prop.constants import M_N, GEV2_TO_CM2
from nutau_prop.chirality_transport_mixed_expm import TauHelicityFlipTransportMixedExpm
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'results'/'v0.41'; OUT.mkdir(parents=True,exist_ok=True)
PDF='CT18NLO'; GQ=0.03; MASSES=(50.0,); MODELS={'eq4p4_derivative':1.0,'eq4p3_dipole':2.0}
TARGET_E=1e8; TARGET_RATIO=.05; ENERGIES=np.geomspace(1e4,2e10,81); ZENITHS=(.1,.3,.5,.7,1.0); GAMMAS=(2.0,2.5); PROBE=(1.,3.,10.,30.,100.,300.)
class R(Case6DIS):
 def _d2(self,E,x,y,M,k,f):
  s=2*M_N*E; q2=x*y*s
  if q2<self.cfg.qmin_gev**2:return 0.
  p=self._pdf_components(x,q2); c=k/M**2
  return f*(4/math.pi)*c*c*x*s*s*y*(1-y)*sum(p.values())*(M*M/(M*M+q2))**2*GEV2_TO_CM2
 def ds(self,E,y,M,k,g,f):
  s=2*M_N*E; xmin=max(self.cfg.qmin_gev**2/(s*y),1e-9); xmax=min(1-1e-10,(g/k)**2/s)
  if xmin>=xmax:return 0.
  return self._integrate_logx(lambda x:self._d2(E,x,y,M,k,f),xmin,xmax)
C=R(Case6Config(pdf_set=PDF,nx=48,y_log_bins=20,y_linear_bins=20))
def ratio(E,M,k,g,f):
 _,ys,dy=C.y_grid(); snp=ssm=0.
 for y,d in zip(ys,dy):
  cc,nc=C.dsdy_sm(E,y); snp+=C.ds(E,y,M,k,g,f)*d; ssm+=(cc+nc)*d
 return snp/max(ssm,1e-300)
def solve(M,g,f):
 ks=np.geomspace(2e-6,2e-3,48); rs=np.array([ratio(TARGET_E,M,k,g,f) for k in ks]); ii=np.where((rs[:-1]<.05)&(rs[1:]>=.05))[0]
 if not len(ii):
  j=int(np.argmax(rs)); return None,float(ks[j]),float(rs[j])
 lo,hi=ks[ii[0]],ks[ii[0]+1]
 for _ in range(24):
  m=math.sqrt(lo*hi)
  if ratio(TARGET_E,M,m,g,f)>=.05:hi=m
  else:lo=m
 return float(hi),float(hi),float(ratio(TARGET_E,M,hi,g,f))
def table(M,k,g,f,anti):
 _,ys,dy=C.y_grid(); ss=[]
 for E in ENERGIES:
  cc=np.zeros_like(ys);nc=np.zeros_like(ys);npv=np.zeros_like(ys)
  for i,y in enumerate(ys): cc[i],nc[i]=C.dsdy_sm(E,y,antinu=anti); npv[i]=C.ds(E,y,M,k,g,f)
  a,b=probability_from_dsdy(cc,dy);c,d=probability_from_dsdy(nc,dy);e,h=probability_from_dsdy(npv,dy);ss.append((a,c,e,b,d,h))
 return dict(energies=ENERGIES,y=ys,dy=dy,sigma_cc=np.array([r[0] for r in ss]),sigma_nc_sm=np.array([r[1] for r in ss]),sigma_np=np.array([r[2] for r in ss]),p_cc=np.array([r[3] for r in ss]),p_nc_sm=np.array([r[4] for r in ss]),p_np=np.array([r[5] for r in ss]))
def interp(y,e):return float(np.interp(np.log(e),np.log(ENERGIES),y))
def run(model,M,k,f):
 tr=TauHelicityFlipTransportMixedExpm(table(M,k,GQ,f,False),table(M,k,GQ,f,True),tau_refine=8); rows=[]
 for gam in GAMMAS:
  inc=(ENERGIES/1e6)**(-gam)
  for mu in ZENITHS:
   sm=tr.propagate(mu,inc,inc,mode='dirac_sink',use_np=False,max_step_km=48,track_secondary_numu=True); bs=tr.propagate(mu,inc,inc,mode='dirac_sink',use_np=True,max_step_km=48,track_secondary_numu=True)
   ts=sm['nu_tau']+sm['antinu_tau'];tb=bs['nu_tau']+bs['antinu_tau'];ms=sm['nu_mu_secondary']+sm['antinu_mu_secondary'];mb=bs['nu_mu_secondary']+bs['antinu_mu_secondary']
   for ep in PROBE:
    e=ep*1e6;a=interp(ts,e);b=interp(tb,e);c=interp(ms,e);d=interp(mb,e)
    rows.append(dict(model=model,M_GeV=M,gq=GQ,kappa=k,vertex_scale_GeV=GQ/k,gamma=gam,mu=mu,E_PeV=ep,tau_frac=(b-a)/max(a,1e-300),secondary_numu_frac=(d-c)/max(c,1e-300) if c>0 else np.nan))
 return rows
def main():
 gate=[];rows=[]
 for model,f in MODELS.items():
  for M in MASSES:
   k,kp,rp=solve(M,GQ,f);gate.append(dict(model=model,M_GeV=M,gq=GQ,kappa5=k,peak_kappa=kp,peak_ratio=rp,status='OPEN_5PCT' if k else 'CLOSED'))
   print(gate[-1],flush=True)
   if k:rows+=run(model,M,k,f)
 pd.DataFrame(gate).to_csv(OUT/'benchmark_gate_v0_41.csv',index=False);df=pd.DataFrame(rows);df.to_csv(OUT/'earth_propagation_v0_41.csv',index=False)
 if len(df):
  out=[]
  for (m,M,g),q in df.groupby(['model','M_GeV','gamma']):
   i=q.tau_frac.abs().idxmax();j=q.secondary_numu_frac.abs().idxmax();a=q.loc[i];b=q.loc[j];out.append(dict(model=m,M_GeV=M,gamma=g,max_abs_tau_frac=abs(a.tau_frac),tau_at_mu=a.mu,tau_at_E_PeV=a.E_PeV,max_abs_secondary_numu_frac=abs(b.secondary_numu_frac),numu_at_mu=b.mu,numu_at_E_PeV=b.E_PeV))
  pd.DataFrame(out).to_csv(OUT/'earth_propagation_summary_v0_41.csv',index=False)
 (OUT/'metadata_v0_41.json').write_text(json.dumps(dict(pdf=PDF,gq=GQ,validity='sqrt(shat)<gq/kappa',mode='dirac_sink',note='100/220 GeV skipped: v0.39 peak-capacity already below 5% at gq=0.03'),indent=2))
if __name__=='__main__':main()
