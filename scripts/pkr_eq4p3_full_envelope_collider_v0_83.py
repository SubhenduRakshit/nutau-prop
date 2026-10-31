#!/usr/bin/env python3
"""v0.83 — PKR Eq.4.3 full-a UHE + collider screen with no new events.

Uses:
  * v0.76 CT18 Eq4.4 base envelopes; Eq4.3 ratio is exactly 2x at equal M,gq,a.
  * v0.82 physical O+h.c. dipole width Gamma_inv=a^2 M^3/(6 pi).
  * v0.79 mass-dependent CMS production x acceptance bin vectors.
  * v0.81 canonical M=50 fixed-shape multibin normalization and the v0.79
    M=100/220 pilot multibin normalizations.

No MG5, Pythia, Delphes or CheckMATE is run.  At fixed mediator mass and in the
NWA, changing the invisible operator/couplings only rescales the fully invisible
recoil sample by P=gq^2 BR_inv.  Therefore fixed-shape multibin r scales linearly
with P.  M=100/220 remain pilot-shape diagnostics, not final exclusions.

Dijet treatment is deliberately one-sided: gq*sqrt(BRqq)<0.03 is retained only
as a conservative global-minimum screening pass based on the published CMS
statement that observed universal-gq limits span about 0.03--0.13.  Failing this
screen is NOT called exclusion; exact HEPData mass-specific values are still
required near a boundary.
"""
from __future__ import annotations

import json, math
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results'/'v0.76_eq4p4_mass_prescreen'/'base_mass_envelopes.csv'
KBFILE=ROOT/'results'/'v0.79_eq4p4_mass_acceptance_pilot'/'bin_yields_and_Kb.csv'
RANK79=ROOT/'results'/'v0.79_eq4p4_mass_acceptance_pilot'/'mass_ranking.csv'
V81=ROOT/'results'/'v0.81_eq4p4_M50_same_seed_reweight'/'M50_same_seed_reweight.json'
OUT=ROOT/'results'/'v0.83_eq4p3_full_envelope_collider'
OUT.mkdir(parents=True,exist_ok=True)

TARGET=0.05
MASSES=(50.0,100.0,220.0)
GQS=np.arange(0.005,0.16001,0.001)
QUARK_MASS={'u':0.00216,'d':0.00467,'s':0.0934,'c':1.27,'b':4.18,'t':172.69}
GLOBAL_DIJET_FLOOR=0.03


def gamma_qq_one(M,gq,mq):
    if M<=2*mq: return 0.0
    r=(mq/M)**2
    return 3*gq*gq*M/(12*math.pi)*(1+2*r)*math.sqrt(1-4*r)


def log_cross(a0,a1,f0,f1,target):
    x0,x1=math.log(a0),math.log(a1)
    y0,y1=math.log(max(f0,1e-300)),math.log(max(f1,1e-300))
    if y1==y0: return math.sqrt(a0*a1)
    t=(math.log(target)-y0)/(y1-y0)
    return math.exp(x0+t*(x1-x0))


def intervals_above(a,f,target):
    a=np.asarray(a,float); f=np.asarray(f,float); mask=f>=target
    out=[]; i=0
    while i<len(a):
        if not mask[i]: i+=1; continue
        j=i
        while j+1<len(a) and mask[j+1]: j+=1
        lo=float(a[0]) if i==0 else log_cross(a[i-1],a[i],f[i-1],f[i],target)
        hi=float(a[-1]) if j==len(a)-1 else log_cross(a[j],a[j+1],f[j],f[j+1],target)
        out.append((lo,hi))
        i=j+1
    return out


def metrics(M,gq,a):
    gqwidth=sum(gamma_qq_one(M,gq,mq) for mq in QUARK_MASS.values())
    ginv=a*a*M**3/(6*math.pi)  # v0.82 physical O+h.c. Eq4.3 width
    gtot=gqwidth+ginv
    brinv=ginv/gtot; brqq=gqwidth/gtot
    P=gq*gq*brinv
    return dict(Gamma_q_GeV=gqwidth,Gamma_inv_GeV=ginv,Gamma_total_GeV=gtot,
                Gamma_over_M=gtot/M,BR_inv=brinv,BR_qq=brqq,
                gq_equiv_dijet=gq*math.sqrt(brqq),coupling_factor=P,
                vertex_scale_GeV=1/a,narrow_width_pass=bool(gtot/M<0.10),
                vertex_scale_pass=bool(1/a>M),
                dijet_global_floor_screen_pass=bool(gq*math.sqrt(brqq)<GLOBAL_DIJET_FLOOR))


def load_refs():
    rank=pd.read_csv(RANK79)
    refs={}
    for M in (100.0,220.0):
        q=rank[np.isclose(rank.M_GeV,M)].iloc[0]
        refs[M]=dict(P=float(q.coupling_factor),r=float(q.pilot_multibin_r),
                     mc=float(q.pilot_mc_scale),source='v0.79 20k pilot shape')
    j=json.loads(V81.read_text())
    q=rank[np.isclose(rank.M_GeV,50.0)].iloc[0]
    refs[50.0]=dict(P=float(j['representative']['coupling_factor_gq2BR']),
                    r=float(j['checkmate']['multibin_r']),mc=float(q.pilot_mc_scale),
                    source='v0.81 100k deterministic M50 shape')
    return refs


def mass_score_prefactor(kb,M):
    q=kb[np.isclose(kb.M_GeV,M)]
    if q.empty: raise RuntimeError(f'no K_b rows for M={M}')
    denom=np.maximum(q.b.to_numpy(float)+q.db.to_numpy(float)**2,1e-12)
    k=q.K_b.to_numpy(float)
    return float(math.sqrt(np.sum(k*k/denom)))


def main():
    for p in (BASE,KBFILE,RANK79,V81):
        if not p.exists(): raise RuntimeError(f'missing required input {p}')
    base=pd.read_csv(BASE); kb=pd.read_csv(KBFILE); refs=load_refs()
    probe_rows=[]; interval_rows=[]; summary=[]

    for M in MASSES:
        q=base[np.isclose(base.M_GeV,M)].sort_values('a_GeV_inv')
        agrid=q.a_GeV_inv.to_numpy(float)
        # Eq4.3 is exactly twice Eq4.4 at equal M,gq,a.
        f43=2.0*q.ratio100_at_gq1.to_numpy(float)
        pref=mass_score_prefactor(kb,M)
        ref=refs[M]
        best=None; best_floor=None; nint=0; nvalid=0

        for gq in GQS:
            target_base=TARGET/(gq*gq)
            ints=intervals_above(agrid,f43,target_base)
            for ii,(lo,hi) in enumerate(ints):
                nint+=1
                pts=np.geomspace(lo,hi,180) if hi>lo else np.asarray([lo])
                local=[]
                for av in pts:
                    m=metrics(M,float(gq),float(av))
                    if not (m['narrow_width_pass'] and m['vertex_scale_pass']):
                        continue
                    nvalid+=1
                    P=m['coupling_factor']
                    rpred=ref['r']*P/ref['P']
                    mc=ref['mc']
                    rlo=rpred*(1-mc); rhi=rpred*(1+mc)
                    if rhi<0.8: mclass='ROBUST_SURVIVAL_SCREEN'
                    elif rlo>1.2: mclass='ROBUST_EXCLUSION_SCREEN'
                    else: mclass='MARGINAL_SCREEN'
                    row=dict(M_GeV=M,gq=float(gq),a_GeV_inv=float(av),interval_index=ii,
                             target_ratio=TARGET,monojet_r_pred=rpred,monojet_r_mc_lo=rlo,
                             monojet_r_mc_hi=rhi,monojet_class=mclass,
                             diag_background_weighted_score=P*pref,**m)
                    probe_rows.append(row); local.append(row)
                    if best is None or rpred<best['monojet_r_pred']: best=row
                    if m['dijet_global_floor_screen_pass'] and (best_floor is None or rpred<best_floor['monojet_r_pred']):
                        best_floor=row
                interval_rows.append(dict(M_GeV=M,gq=float(gq),interval_index=ii,
                                          a_lo_GeV_inv=lo,a_hi_GeV_inv=hi,
                                          n_width_vertex_valid=len(local),
                                          has_global_floor_dijet_point=any(x['dijet_global_floor_screen_pass'] for x in local)))

        if best is None:
            status='NO_WIDTH_VERTEX_VALID_5PCT_POINT'
            rep={}
        else:
            rep=best_floor if best_floor is not None else best
            if best_floor is not None and rep['monojet_r_mc_hi']<0.8:
                status='ROBUST_MONOJET_PLUS_GLOBAL_DIJET_FLOOR_SCREEN_SURVIVOR'
            elif best_floor is not None and rep['monojet_r_pred']<1.2:
                status='MARGINAL_MONOJET_GLOBAL_DIJET_FLOOR_SCREEN_POINT'
            elif best_floor is not None:
                status='GLOBAL_DIJET_FLOOR_POINT_BUT_MONOJET_STRONG'
            else:
                status='NO_GLOBAL_DIJET_FLOOR_SCREEN_POINT'
        summary.append(dict(M_GeV=M,status=status,n_ge5pct_intervals=nint,
                            n_width_vertex_valid_probes=nvalid,
                            representative_gq=rep.get('gq',np.nan),
                            representative_a=rep.get('a_GeV_inv',np.nan),
                            representative_BR_inv=rep.get('BR_inv',np.nan),
                            representative_Gamma_over_M=rep.get('Gamma_over_M',np.nan),
                            representative_gq_equiv_dijet=rep.get('gq_equiv_dijet',np.nan),
                            representative_coupling_factor=rep.get('coupling_factor',np.nan),
                            representative_monojet_r=rep.get('monojet_r_pred',np.nan),
                            representative_r_mc_hi=rep.get('monojet_r_mc_hi',np.nan),
                            representative_score=rep.get('diag_background_weighted_score',np.nan),
                            monojet_reference=ref['source']))

    P=pd.DataFrame(probe_rows); I=pd.DataFrame(interval_rows); S=pd.DataFrame(summary)
    P.to_csv(OUT/'valid_probe_points.csv',index=False)
    I.to_csv(OUT/'full_intervals.csv',index=False)
    S.to_csv(OUT/'mass_summary.csv',index=False)

    candidates=S[S.status=='ROBUST_MONOJET_PLUS_GLOBAL_DIJET_FLOOR_SCREEN_SURVIVOR']
    cand=None
    if len(candidates):
        cand=float(candidates.loc[candidates.representative_monojet_r.idxmin()].M_GeV)
    meta={
      'version':'v0.83_eq4p3_full_envelope_collider',
      'width':'physical Eq4.3 O+h.c.: Gamma_inv=a^2 M^3/(6 pi), certified v0.82',
      'uhe':'Eq4.3 ratio = 2 x corrected Eq4.4 base envelope at equal M,gq,a; full connected a intervals retained',
      'monojet':'no-new-events NWA fixed-shape scaling using existing mass-dependent CMS bin vectors',
      'dijet':'gq*sqrt(BRqq)<0.03 is sufficient screening only; failure is not exclusion; exact HEPData remains required near boundary',
      'candidate_M_GeV':cand,
      'warning':'M100/M220 monojet shapes come from 20k v0.79 pilots and are screening diagnostics only.'
    }
    (OUT/'summary.json').write_text(json.dumps(meta,indent=2)+'\n')

    print('=== v0.83 Eq4.3 full-a no-new-events collider screen ===')
    print(S.to_string(index=False))
    print(f'robust screened candidate mass: {cand if cand is not None else "NONE"}')
    print('interpretation: global 0.03 dijet floor is sufficient-only; exact mass-specific dijet still required for final certification')
    return 0

if __name__=='__main__':
    raise SystemExit(main())
