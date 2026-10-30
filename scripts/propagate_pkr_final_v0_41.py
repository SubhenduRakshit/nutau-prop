#!/usr/bin/env python3
from __future__ import annotations
import json,math
from pathlib import Path
import numpy as np,pandas as pd
from nutau_prop.case6 import Case6Config,Case6DIS,probability_from_dsdy
from nutau_prop.constants import M_N,GEV2_TO_CM2
from nutau_prop.chirality_transport_mixed_expm import TauHelicityFlipTransportMixedExpm
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results'/'v0.41b';OUT.mkdir(parents=True,exist_ok=True)
PDF='CT18NLO';GQ=.03;M=50.;MODELS={'eq4p4_derivative':1.,'eq4p3_dipole':2.};E0=1e8;EGRID=np.geomspace(1e4,2e10,51);PROBE=(1.,3.,10.,30.,100.,300.);GAMMAS=(2.,2.5)
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
C=R(Case6Config(pdf_set=PDF,nx=36,y_log_bins=14,y_linear_bins=14))
def ratio(k,f):
 _,y,dy=C.y_grid();a=b=0.
 for yy,d in zip(y,dy):
  cc,nc=C.dsdy_sm(E0,yy);a+=C.dsnp(E0,yy,k,f)*d;b+=(cc+nc)*d
 return a/max(b,1e-300)
def solve(f):
 ks=np.geomspace(2e-6,2e-3,40);rs=np.array([ratio(k,f) for k in ks]);ii=np.where((rs[:-1]<.05)&(rs[1:]>=.05))[0]
 if not len(ii):j=int(np.argmax(rs));return None,ks[j],rs[j]
 lo,hi=ks[ii[0]],ks[ii[0]+1]
 for _ in range(20):
  m=math.sqrt(lo*hi)
  if ratio(m,f)>=.05:hi=m
  else:lo=m
 return float(hi),float(hi),float(ratio(hi,f))
def tab(k,f,anti):
 _,y,dy=C.y_grid();q=[]
 for E in EGRID:
  cc=np.zeros_like(y);nc=np.zeros_like(y);npv=np.zeros_like(y)
  for i,yy in enumerate(y):cc[i],nc[i]=C.dsdy_sm(E,yy,antinu=anti);npv[i]=C.dsnp(E,yy,k,f)
  a,A=probability_from_dsdy(cc,dy);b,B=probability_from_dsdy(nc,dy);c,D=probability_from_dsdy(npv,dy);q.append((a,b,c,A,B,D))
 return dict(energies=EGRID,y=y,dy=dy,sigma_cc=np.array([r[0] for r in q]),sigma_nc_sm=np.array([r[1] for r in q]),sigma_np=np.array([r[2] for r in q]),p_cc=np.array([r[3] for r in q]),p_nc_sm=np.array([r[4] for r in q]),p_np=np.array([r[5] for r in q]))
def it(y,e):return float(np.interp(np.log(e),np.log(EGRID),y))
def main():
 gates=[];rows=[]
 for model,f in MODELS.items():
  k,kp,rp=solve(f);gates.append(dict(model=model,M_GeV=M,gq=GQ,kappa5=k,peak_kappa=kp,peak_ratio=rp,status='OPEN_5PCT' if k else 'CLOSED'));print(gates[-1],flush=True)
  if not k:continue
  tr=TauHelicityFlipTransportMixedExpm(tab(k,f,False),tab(k,f,True),tau_refine=4)
  for gam in GAMMAS:
   inc=(EGRID/1e6)**(-gam);sm=tr.propagate(1.,inc,inc,mode='dirac_sink',use_np=False,max_step_km=192,track_secondary_numu=True);bs=tr.propagate(1.,inc,inc,mode='dirac_sink',use_np=True,max_step_km=192,track_secondary_numu=True)
   ts=sm['nu_tau']+sm['antinu_tau'];tb=bs['nu_tau']+bs['antinu_tau'];ms=sm['nu_mu_secondary']+sm['antinu_mu_secondary'];mb=bs['nu_mu_secondary']+bs['antinu_mu_secondary']
   for ep in PROBE:
    e=ep*1e6;a=it(ts,e);b=it(tb,e);c=it(ms,e);d=it(mb,e);rows.append(dict(model=model,M_GeV=M,gq=GQ,kappa=k,vertex_scale_GeV=GQ/k,gamma=gam,mu=1.,E_PeV=ep,tau_frac=(b-a)/max(a,1e-300),secondary_numu_frac=(d-c)/max(c,1e-300) if c>0 else np.nan))
 pd.DataFrame(gates).to_csv(OUT/'benchmark_gate.csv',index=False);df=pd.DataFrame(rows);df.to_csv(OUT/'vertical_propagation.csv',index=False)
 out=[]
 for (m,g),q in df.groupby(['model','gamma']):
  i=q.tau_frac.abs().idxmax();j=q.secondary_numu_frac.abs().idxmax();a=q.loc[i];b=q.loc[j];out.append(dict(model=m,gamma=g,max_abs_tau_frac=abs(a.tau_frac),tau_E_PeV=a.E_PeV,max_abs_secondary_numu_frac=abs(b.secondary_numu_frac),numu_E_PeV=b.E_PeV))
 pd.DataFrame(out).to_csv(OUT/'summary.csv',index=False);(OUT/'metadata.json').write_text(json.dumps(dict(screen='vertical-first',pdf=PDF,gq=GQ,M_GeV=M,mode='dirac_sink',validity='sqrt(shat)<gq/kappa',tau_refine=4,max_step_km=192),indent=2))
if __name__=='__main__':main()
