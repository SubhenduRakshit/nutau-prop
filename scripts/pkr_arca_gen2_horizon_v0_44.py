#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, math
from functools import lru_cache
from pathlib import Path
import numpy as np, pandas as pd
from scipy.optimize import minimize_scalar
from scipy.special import ndtr
from nutau_prop.arca import calibrate_fig7_response, fig7_response
from nutau_prop.config import load_config
from nutau_prop.transport import NuOnlyTransport
ROOT=Path(__file__).resolve().parents[1]
V42=ROOT/'results'/'v0.42d'; OUT=ROOT/'results'/'v0.44'
MODELS=['eq4p3_dipole','eq4p4_derivative']; TOPOLOGIES=['active_redistribution_proxy','sterile_sink']
GAMMA_MAP=2.5
MU_EDGES=np.arange(0.25,0.55+1e-12,0.05); MU_CENTERS=.5*(MU_EDGES[:-1]+MU_EDGES[1:])
EREC_EDGES=np.array([1e5,3e5,1e6,3e6,1e7,3e7,1e8],float)
SIGMA_LOG10E=.30; ETRUE_MAX=1e8; TAU_MIN=1.; TAU_MAX=3.

def load_npz(path):
 z=np.load(path); return {k:z[k] for k in z.files}

def nuonly_table(tab):
 return {'energies':tab['energies'],'sigma_cc':tab['sigma_cc'],'sigma_nc':tab['sigma_nc_sm'],'y':tab['y'],'p_nc':tab['p_nc_sm']}

def migration_prob(E,lo,hi,sigma=SIGMA_LOG10E):
 x=np.log10(np.asarray(E,float)); return ndtr((math.log10(hi)-x)/sigma)-ndtr((math.log10(lo)-x)/sigma)

def bilinear(mus,energies,grid,mu,E):
 per=np.array([np.interp(np.log10(E),np.log10(energies),row,left=row[0],right=row[-1]) for row in grid])
 return np.array([np.interp(mu,mus,per[:,i]) for i in range(len(E))])

def track_parent_key(mu):
 if .25<=mu<.50:return 'cm025_m050'
 if .50<=mu<=.55:return 'cm050_m075'
 raise ValueError(mu)

@lru_cache(maxsize=None)
def make_maps(model,topology):
 df=pd.read_csv(V42/'coupled_true_energy_map.csv'); q=df[(df.model==model)&(df.topology==topology)].copy()
 if q.empty: raise RuntimeError(f'Missing map {model}/{topology}')
 mus=np.sort(q.mu.unique().astype(float)); efull=np.sort(q.E_GeV.unique().astype(float)); use=efull<=ETRUE_MAX; E=efull[use]
 def grid(col):
  return np.asarray([q[np.isclose(q.mu,m)].sort_values('E_GeV')[col].to_numpy(float)[use] for m in mus])
 tau_sm=grid('tau_flux_SM'); tau_bs=grid('tau_flux_BSM'); sec_sm=grid('sec_numu_SM'); sec_bs=grid('sec_numu_BSM'); tg=grid('tau_SM')
 rtau=np.divide(tau_bs,tau_sm,out=np.ones_like(tau_sm),where=tau_sm>0)
 tnu=load_npz(V42/'tables'/f'{model}_nu.npz'); tanu=load_npz(V42/'tables'/f'{model}_antinu.npz')
 trn=NuOnlyTransport(nuonly_table(tnu)); tra=NuOnlyTransport(nuonly_table(tanu)); inc=(tnu['energies']/1e6)**(-GAMMA_MAP)
 rmu=np.ones_like(rtau)
 for im,m in enumerate(mus):
  pn=trn.propagate(float(m),inc,max_step_km=96.,max_p=.05)['nu']; pa=tra.propagate(float(m),inc,max_step_km=96.,max_p=.05)['nu']
  primary=(pn+pa)[use]; sm=primary+sec_sm[im]; bs=primary+sec_bs[im]
  rmu[im]=np.divide(bs,sm,out=np.ones_like(sm),where=sm>0)
 snp=.5*(tnu['sigma_np']+tanu['sigma_np']); snc=.5*(tnu['sigma_nc_sm']+tanu['sigma_nc_sm'])
 rho=np.divide(snp,snc,out=np.zeros_like(snp),where=snc>0)[use]
 return dict(mus=mus,energies=E,rtau=rtau,rmu=rmu,tau_guide=tg,rho_np_nc=rho)

def flux(E,cfg,gamma,norm=1.):
 a=cfg['analysis_flux']; return norm*float(a['phi0_per_gev_cm2_s_sr'])*(E/float(a['pivot_gev']))**(-gamma)

def integ(E,y): return float(np.trapezoid(y,E))

def build_counts(model,topology,cfg,gamma,norm,years,direct_np=True):
 M=make_maps(model,topology); E=M['energies']; R=fig7_response(E); sf_sh,sf_tr=calibrate_fig7_response(cfg); ph=flux(E,cfg,gamma,norm)
 sec=years*365.25*86400.; rows=[]
 for lo,hi,mu in zip(MU_EDGES[:-1],MU_EDGES[1:],MU_CENTERS):
  domega=2*math.pi*(hi-lo); rt=bilinear(M['mus'],E,M['rtau'],float(mu),E); rm=bilinear(M['mus'],E,M['rmu'],float(mu),E); tg=bilinear(M['mus'],E,M['tau_guide'],float(mu),E)
  rho=M['rho_np_nc'] if direct_np else np.zeros_like(E)
  sh=R['shower']; Ae=np.asarray(sh['up_ecc'])*sf_sh; At=np.asarray(sh['up_tcc'])*sf_sh; Anc=np.asarray(sh['up_nc'])*sf_sh
  rnc=(1.+rm+rt*(1.+rho))/3.; Ash_sm=Ae+At+Anc; Ash_bs=Ae+At*rt+Anc*rnc
  key=track_parent_key(float(mu)); Atr_sm=np.asarray(R['track'][key])*sf_tr; Atr_bs=Atr_sm*rm
  for elo,ehi in zip(EREC_EDGES[:-1],EREC_EDGES[1:]):
   P=migration_prob(E,float(elo),float(ehi))
   for ch,Asm,Abs in [('shower',Ash_sm,Ash_bs),('track',Atr_sm,Atr_bs)]:
    base=sec*domega*1e4*ph*P; nsm=integ(E,base*Asm); nbs=integ(E,base*Abs); den=integ(E,base*Asm); te=integ(E,base*Asm*tg)/max(den,1e-300)
    rows.append(dict(model=model,topology=topology,channel=ch,mu_lo=float(lo),mu_hi=float(hi),mu=float(mu),Erec_lo_GeV=float(elo),Erec_hi_GeV=float(ehi),tau_SM_eff=float(te),in_horizon_band=bool(TAU_MIN<=te<=TAU_MAX),N_SM=nsm,N_BSM=nbs,frac=(nbs-nsm)/max(nsm,1e-300),direct_np_in_showers=bool(direct_np)))
 return pd.DataFrame(rows)

def asimov_q(n,m):
 n=np.asarray(n,float); m=np.maximum(np.asarray(m,float),1e-300); term=np.where(n>0,n*np.log(np.maximum(n,1e-300)/m),0.); return float(2*np.sum(m-n+term))

def profiled_Z(model,topology,cfg,years,channels,horizon_only=True,direct_np=True):
 g0=float(cfg['analysis_flux']['gamma']); truth=build_counts(model,topology,cfg,g0,1.,years,direct_np); mask=truth.channel.isin(channels)
 if horizon_only: mask &= truth.in_horizon_band
 truth=truth[mask].reset_index(drop=True)
 if truth.empty:return dict(Z=0.,q=0.,alpha=1.,gamma=g0,n_bins=0,N_SM=0.,N_BSM=0.),truth
 n=truth.N_BSM.to_numpy(float); gmin=float(cfg['profile']['gamma_min']); gmax=float(cfg['profile']['gamma_max'])
 def obj(g):
  sm=build_counts(model,topology,cfg,float(g),1.,years,direct_np)[mask].reset_index(drop=True); s=sm.N_SM.to_numpy(float); a=n.sum()/max(s.sum(),1e-300); return asimov_q(n,a*s)
 opt=minimize_scalar(obj,bounds=(gmin,gmax),method='bounded',options={'xatol':2e-4}); sm=build_counts(model,topology,cfg,float(opt.x),1.,years,direct_np)[mask].reset_index(drop=True); s=sm.N_SM.to_numpy(float); a=n.sum()/max(s.sum(),1e-300); q=asimov_q(n,a*s)
 return dict(Z=math.sqrt(max(q,0.)),q=q,alpha=float(a),gamma=float(opt.x),n_bins=len(truth),N_SM=float(truth.N_SM.sum()),N_BSM=float(truth.N_BSM.sum())),truth

def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--config',default='config/arca.yaml'); ap.add_argument('--outdir',default='results/v0.44'); ap.add_argument('--no-direct-np',action='store_true'); a=ap.parse_args()
 cfg=load_config(ROOT/a.config); out=ROOT/a.outdir; out.mkdir(parents=True,exist_ok=True)
 need=[V42/'coupled_true_energy_map.csv']
 for m in MODELS: need += [V42/'tables'/f'{m}_nu.npz',V42/'tables'/f'{m}_antinu.npz']
 miss=[str(p) for p in need if not p.exists()]
 if miss: raise RuntimeError('Missing v0.42d inputs:\n'+'\n'.join(miss))
 years=float(cfg['arca']['exposure_years']); direct=not a.no_direct_np; sums=[]; bins=[]
 for m in MODELS:
  for t in TOPOLOGIES:
   for lab,chs in [('showers_only',('shower',)),('tracks_only',('track',)),('tracks_plus_showers',('track','shower'))]:
    r,b=profiled_Z(m,t,cfg,years,chs,True,direct); sums.append(dict(model=m,topology=t,selection=lab,exposure_years=years,direct_np_in_showers=direct,**r)); qq=b.copy(); qq['selection_for_summary']=lab; bins.append(qq)
   r,_=profiled_Z(m,t,cfg,years,('track','shower'),False,direct); sums.append(dict(model=m,topology=t,selection='available_mu_all_tau',exposure_years=years,direct_np_in_showers=direct,**r))
 S=pd.DataFrame(sums); S.to_csv(out/'arca_horizon_summary.csv',index=False); pd.concat(bins,ignore_index=True).to_csv(out/'arca_horizon_bins.csv',index=False)
 scale=[]
 for _,r in S[S.selection=='tracks_plus_showers'].iterrows():
  z=float(r.Z)
  for mult in [1.,5.,10.]: scale.append(dict(model=r.model,topology=r.topology,acceptance_multiplier=mult,Z_scaled=z*math.sqrt(mult),interpretation='ARCA_10y' if mult==1 else 'Gen2-like generic scaling; not detector response'))
  scale += [dict(model=r.model,topology=r.topology,acceptance_multiplier=(2/z)**2 if z>0 else np.inf,Z_scaled=2.,interpretation='required multiplier for 2sigma'),dict(model=r.model,topology=r.topology,acceptance_multiplier=(3/z)**2 if z>0 else np.inf,Z_scaled=3.,interpretation='required multiplier for 3sigma')]
 G=pd.DataFrame(scale); G.to_csv(out/'gen2_acceptance_scaling.csv',index=False)
 gate=[]
 for _,r in S[S.selection=='tracks_plus_showers'].iterrows(): gate.append(dict(model=r.model,topology=r.topology,Z_ARCA_10y_upper_bound=float(r.Z),verdict='ARCA_CLOSED_BY_SIGNAL_ONLY_UPPER_BOUND' if r.Z<2 else 'ARCA_PROMISING_BUT_BACKGROUND_LIKELIHOOD_REQUIRED'))
 pd.DataFrame(gate).to_csv(out/'detector_gate.csv',index=False)
 meta={'version':'v0.44','input_transport':'results/v0.42d converged coupled transport','arca_response':'published Fig.7 coarse digitisation calibrated to exact published 220/year track and 60/year shower cosmic rates','important_absorption_fix':'ARCA upgoing Aeff already contains SM Earth absorption; use incident flux times SM Aeff and reweight only by BSM/SM transport ratios','zenith':{'mu_edges':MU_EDGES.tolist(),'limitation':'transport uses fine mu bins; public shower Aeff is hemisphere-averaged and track Aeff held constant inside published coarse parent bin'},'reco_energy_edges_GeV':EREC_EDGES.tolist(),'sigma_log10E':SIGMA_LOG10E,'horizon_band':[TAU_MIN,TAU_MAX],'profiled':['astrophysical normalization','spectral index'],'secondary_numu':'included in total primary+tau-decay-secondary muon-flavor modifier','direct_np_showers':'included with NC response proxy' if direct else 'disabled','background_status':'no bin-by-bin atmospheric background; ARCA Z is optimistic signal-only upper bound','gen2_status':'no fabricated response; x5/x10 generic acceptance/exposure scaling','map_gamma':GAMMA_MAP,'limitation_gamma':'transport modifiers fixed at gamma=2.5 while detector gamma is profiled'}
 (out/'metadata.json').write_text(json.dumps(meta,indent=2))
 print('\n=== v0.44 ARCA absorption-horizon summary ==='); print(S[S.selection.isin(['showers_only','tracks_plus_showers'])][['model','topology','selection','Z','n_bins','N_SM','N_BSM','alpha','gamma']].to_string(index=False)); print('\n=== v0.44 detector gate ==='); print(pd.DataFrame(gate).to_string(index=False)); print('\n=== v0.44 Gen2-like scaling ==='); print(G[G.interpretation.str.contains('ARCA_10y|Gen2-like',regex=True)].to_string(index=False))
if __name__=='__main__': main()
