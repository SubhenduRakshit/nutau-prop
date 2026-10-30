#!/usr/bin/env python3
"""v0.42: coupled PKR absorption-horizon transport + public-response detector weighting.

Centers the analysis on 1 <= tau_SM <= 3 rather than the vertical chord.
The detector weighting is a public-response projection, not an official ARCA likelihood.
"""
from __future__ import annotations
import json, math, os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import numpy as np
import pandas as pd
from nutau_prop.case6 import Case6Config, Case6DIS, probability_from_dsdy
from nutau_prop.constants import M_N, GEV2_TO_CM2, NA
from nutau_prop.prem import chord_profile
from nutau_prop.chirality_transport_mixed_expm import TauHelicityFlipTransportMixedExpm
from nutau_prop.arca import fig7_response, calibrate_fig7_response
from nutau_prop.config import load_config

ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'results'/'v0.42'; OUT.mkdir(parents=True,exist_ok=True)
TABLEDIR=OUT/'tables'; TABLEDIR.mkdir(parents=True,exist_ok=True)
PDF='CT18NLO'; MZ=50.; GQ=.03; TARGET_E=1e8; TARGET_RATIO=.05
MODELS={'eq4p4_derivative':1.0,'eq4p3_dipole':2.0}
TOPOLOGIES={'sterile_sink':'dirac_sink','active_redistribution_proxy':'majorana_flip'}
MUS=np.array([.08,.12,.18,.25,.35,.45,.55,.70,1.00])
EGRID=np.geomspace(1e4,3e8,61); GAMMA=2.5; SIGMA_LOG10E=.30; TAU_REFINE=8; MAX_STEP_KM=48.

class ResolvedPKRDIS(Case6DIS):
 def _d2(self,E,x,y,k,f):
  s=2*M_N*E; q2=x*y*s
  if q2<self.cfg.qmin_gev**2:return 0.
  p=self._pdf_components(x,q2); c=k/MZ**2
  return f*(4/math.pi)*c*c*x*s*s*y*(1-y)*sum(p.values())*(MZ*MZ/(MZ*MZ+q2))**2*GEV2_TO_CM2
 def dsnp(self,E,y,k,f):
  s=2*M_N*E; xmin=max(self.cfg.qmin_gev**2/(s*y),1e-9); xmax=min(1-1e-10,(GQ/k)**2/s)
  if xmin>=xmax:return 0.
  return self._integrate_logx(lambda x:self._d2(E,x,y,k,f),xmin,xmax)
CALC=ResolvedPKRDIS(Case6Config(pdf_set=PDF,nx=52,y_log_bins=22,y_linear_bins=22))

def integrated_sigmas(E,k,f,antinu=False):
 _,ys,dy=CALC.y_grid(); a=b=0.
 for y,d in zip(ys,dy):
  cc,nc=CALC.dsdy_sm(E,y,antinu=antinu); a+=CALC.dsnp(E,y,k,f)*d; b+=(cc+nc)*d
 return a,b

def ratio100(k,f):
 r=[]
 for anti in (False,True):
  a,b=integrated_sigmas(TARGET_E,k,f,anti); r.append(a/max(b,1e-300))
 return .5*(r[0]+r[1])

def solve_kappa(f):
 ks=np.geomspace(2e-6,2e-3,90); rs=np.array([ratio100(k,f) for k in ks]); ii=np.where((rs[:-1]<.05)&(rs[1:]>=.05))[0]
 if not len(ii):
  j=int(np.argmax(rs)); return None,float(ks[j]),float(rs[j])
 lo,hi=float(ks[ii[0]]),float(ks[ii[0]+1])
 for _ in range(28):
  m=math.sqrt(lo*hi)
  if ratio100(m,f)>=.05:hi=m
  else:lo=m
 return hi,hi,ratio100(hi,f)

def build_table(k,f,antinu=False):
 _,ys,dy=CALC.y_grid(); rows=[]
 for E in EGRID:
  cc=np.zeros_like(ys); nc=np.zeros_like(ys); npv=np.zeros_like(ys)
  for i,y in enumerate(ys): cc[i],nc[i]=CALC.dsdy_sm(E,y,antinu=antinu); npv[i]=CALC.dsnp(E,y,k,f)
  scc,pcc=probability_from_dsdy(cc,dy); snc,pnc=probability_from_dsdy(nc,dy); snp,pnp=probability_from_dsdy(npv,dy)
  rows.append((scc,snc,snp,pcc,pnc,pnp))
 return dict(energies=EGRID,y=ys,dy=dy,sigma_cc=np.asarray([r[0] for r in rows]),sigma_nc_sm=np.asarray([r[1] for r in rows]),sigma_np=np.asarray([r[2] for r in rows]),p_cc=np.asarray([r[3] for r in rows]),p_nc_sm=np.asarray([r[4] for r in rows]),p_np=np.asarray([r[5] for r in rows]))

def save_table(path,t):np.savez_compressed(path,**t)
def load_table(path):
 z=np.load(path); return {k:z[k] for k in z.files}
def incident(E):return (E/1e6)**(-GAMMA)
def column_depth(mu):
 p=chord_profile(float(mu),max_step_km=24.); return float(np.sum(np.asarray(p['rho'])*np.asarray(p['ds_km'])*1e5))
def sm_tau(mu,tnu,tanu):
 s=.5*(tnu['sigma_cc']+tnu['sigma_nc_sm']+tanu['sigma_cc']+tanu['sigma_nc_sm']); return NA*column_depth(mu)*s

def transport_task(arg):
 model,mu=arg; tnu=load_table(TABLEDIR/f'{model}_nu.npz'); tanu=load_table(TABLEDIR/f'{model}_antinu.npz')
 tr=TauHelicityFlipTransportMixedExpm(tnu,tanu,tau_refine=TAU_REFINE); inc=incident(EGRID)
 sm=tr.propagate(mu,inc,inc,mode='dirac_sink',use_np=False,max_step_km=MAX_STEP_KM,track_secondary_numu=True)
 tguide=sm_tau(mu,tnu,tanu); out=[]
 for label,mode in TOPOLOGIES.items():
  bs=tr.propagate(mu,inc,inc,mode=mode,use_np=True,max_step_km=MAX_STEP_KM,track_secondary_numu=True)
  ts=sm['nu_tau']+sm['antinu_tau']; tb=bs['nu_tau']+bs['antinu_tau']; ms=sm['nu_mu_secondary']+sm['antinu_mu_secondary']; mb=bs['nu_mu_secondary']+bs['antinu_mu_secondary']
  for i,E in enumerate(EGRID):
   out.append(dict(model=model,topology=label,mu=mu,E_GeV=E,tau_SM=tguide[i],tau_flux_SM=ts[i],tau_flux_BSM=tb[i],tau_frac=(tb[i]-ts[i])/max(ts[i],1e-300),sec_numu_SM=ms[i],sec_numu_BSM=mb[i],sec_numu_frac=(mb[i]-ms[i])/max(ms[i],1e-300) if ms[i]>0 else np.nan))
 return out

def lognormal_kernel(Etrue,Erec,sigma=SIGMA_LOG10E):
 K=np.exp(-.5*((np.log10(Erec)[None,:]-np.log10(Etrue)[:,None])/sigma)**2); K/=np.maximum(K.sum(axis=1,keepdims=True),1e-300); return K

def mu_width(mu):
 mids=.5*(MUS[:-1]+MUS[1:]); edges=np.r_[0.,mids,1.]
 i=int(np.argmin(np.abs(MUS-mu))); return float(max(edges[i+1]-edges[i],0.))

def detector_weight(df):
 cfg=load_config(ROOT/'config'/'arca.yaml'); sf_sh,_=calibrate_fig7_response(cfg); A=np.asarray(fig7_response(EGRID)['shower']['up_tcc'])*sf_sh
 dE=np.gradient(EGRID); Erec=np.geomspace(1e5,1e8,31); K=lognormal_kernel(EGRID,Erec); rows=[]
 for (model,topology,mu),q in df.groupby(['model','topology','mu']):
  q=q.sort_values('E_GeV'); sm=q.tau_flux_SM.to_numpy(float); bs=q.tau_flux_BSM.to_numpy(float)
  # propagated spectra already contain the incident E^-gamma flux: do NOT multiply by it again.
  base=2*math.pi*mu_width(float(mu))*A*dE
  nsm_true=base*sm; nbs_true=base*bs; nsm=nsm_true@K; nbs=nbs_true@K
  tau=q.tau_SM.to_numpy(float); te=np.zeros_like(Erec)
  for j in range(len(Erec)):
   w=nsm_true*K[:,j]; te[j]=np.sum(w*tau)/max(np.sum(w),1e-300)
  for j,er in enumerate(Erec):
   rows.append(dict(model=model,topology=topology,mu=float(mu),delta_mu=mu_width(float(mu)),Erec_GeV=float(er),tau_SM_eff=float(te[j]),N_SM_proxy=float(nsm[j]),N_BSM_proxy=float(nbs[j]),frac_distortion=float((nbs[j]-nsm[j])/max(nsm[j],1e-300)),in_horizon_band=bool(1<=te[j]<=3)))
 return pd.DataFrame(rows)

def main():
 gates=[]
 for model,f in MODELS.items():
  k,kp,rp=solve_kappa(f)
  if k is None:gates.append(dict(model=model,status='NO_5PCT',peak_kappa=kp,peak_ratio=rp));continue
  gates.append(dict(model=model,status='OPEN_CONDITIONAL',M_GeV=MZ,gq=GQ,kappa=k,vertex_scale_GeV=GQ/k,ratio100=rp))
  save_table(TABLEDIR/f'{model}_nu.npz',build_table(k,f,False)); save_table(TABLEDIR/f'{model}_antinu.npz',build_table(k,f,True))
 pd.DataFrame(gates).to_csv(OUT/'benchmark_gate.csv',index=False)
 tasks=[(m,float(mu)) for m in MODELS for mu in MUS if any(g.get('model')==m and g.get('status')=='OPEN_CONDITIONAL' for g in gates)]
 workers=min(max(1,int(os.environ.get('NUTAU_NPROC','4'))),len(tasks)); allrows=[]
 with ProcessPoolExecutor(max_workers=workers) as ex:
  futs={ex.submit(transport_task,t):t for t in tasks}
  for fut in as_completed(futs): allrows.extend(fut.result()); print('DONE',futs[fut],flush=True)
 df=pd.DataFrame(allrows); df.to_csv(OUT/'coupled_true_energy_map.csv',index=False)
 det=detector_weight(df); det.to_csv(OUT/'arca_public_response_map.csv',index=False)
 summary=[]
 for (model,topology),q in det.groupby(['model','topology']):
  h=q[q.in_horizon_band].copy()
  if not len(h):continue
  h['Sproxy']=np.abs(h.N_BSM_proxy-h.N_SM_proxy)/np.sqrt(np.maximum(h.N_SM_proxy,1e-300)); j=h.Sproxy.idxmax(); b=h.loc[j]
  w=h.N_SM_proxy.to_numpy(float); weighted=np.sum(w*np.abs(h.frac_distortion))/max(np.sum(w),1e-300)
  summary.append(dict(model=model,topology=topology,max_abs_horizon_distortion=float(h.frac_distortion.abs().max()),SM_event_weighted_abs_horizon_distortion=float(weighted),best_mu=float(b.mu),best_Erec_GeV=float(b.Erec_GeV),best_tau_SM_eff=float(b.tau_SM_eff),Sproxy=float(b.Sproxy)))
 pd.DataFrame(summary).to_csv(OUT/'horizon_summary.csv',index=False)
 (OUT/'metadata.json').write_text(json.dumps(dict(version='v0.42',pdf=PDF,M_GeV=MZ,gq=GQ,target_ratio=.05,validity='sqrt(shat)<gq/kappa inside PDF fold',gamma=GAMMA,mus=MUS.tolist(),tau_refine=TAU_REFINE,max_step_km=MAX_STEP_KM,topologies=TOPOLOGIES,sigma_log10E=SIGMA_LOG10E,detector='calibrated public ARCA Fig7 tau-CC upgoing energy response; zenith from PREM coupled transport; not official likelihood'),indent=2))
 print(pd.DataFrame(summary).to_string(index=False),flush=True)
if __name__=='__main__':main()
