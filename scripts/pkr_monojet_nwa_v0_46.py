#!/usr/bin/env python3
"""v0.46 public monojet narrow-width rescaling screen for the PKR 50-GeV benchmark.

This is deliberately NOT a detector recast.  It uses the public CMS EXO-20-004
vector-mediator coupling anchor at m_med=100 GeV, m_DM=m_med/3, g_chi=1:
observed g_q^95 = 0.018.  For a narrow on-shell spin-1 mediator, monojet yield
scales as g_q^2 BR_inv to leading NWA accuracy, while the hard-recoil kinematics
are controlled by pp -> Z'+j.  We therefore calculate the exact BR-rescaled
coupling threshold that a 50-GeV public limit must beat to exclude each PKR
width convention.

The calculation reads v0.45; no new arbitrary collider limit is invented.
"""
from __future__ import annotations
import json, math
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
INP = ROOT / "results" / "v0.45" / "width_br_dijet_gate.csv"
OUT = ROOT / "results" / "v0.46"
OUT.mkdir(parents=True, exist_ok=True)

# CMS EXO-20-004 published anchor (JHEP 11 (2021) 153):
CMS_M_GEV = 100.0
CMS_GQ95 = 0.018
CMS_GCHI = 1.0
CMS_MDM_OVER_M = 1.0 / 3.0

# Standard vector mediator widths used only to compute the reference BR_inv.
def gamma_vector_dirac_over_M(gchi: float, r: float) -> float:
    if 4.0*r >= 1.0:
        return 0.0
    return gchi*gchi/(12.0*math.pi) * (1.0 + 2.0*r) * math.sqrt(1.0 - 4.0*r)

def gamma_vector_quarks_over_M(gq: float, nflav: int=5) -> float:
    # massless u,d,s,c,b; top closed at 100 GeV
    return nflav * 3.0 * gq*gq/(12.0*math.pi)

def main() -> None:
    if not INP.exists():
        raise RuntimeError(f"Missing {INP}; run v0.45 first")
    d = pd.read_csv(INP)

    r = CMS_MDM_OVER_M**2
    ginv_ref = gamma_vector_dirac_over_M(CMS_GCHI, r)
    gqq_ref = gamma_vector_quarks_over_M(CMS_GQ95, 5)
    br_ref = ginv_ref / (ginv_ref + gqq_ref)

    rows=[]
    for _, x in d.iterrows():
        gq=float(x.gq); br=float(x.BR_invisible)
        # Same NWA yield as a CMS vector mediator with gchi=1:
        # gq_PKR^2 BR_PKR = gq_limit^2 BR_ref.
        critical = gq * math.sqrt(br / br_ref)
        mu_vs_100 = (gq*gq*br) / (CMS_GQ95*CMS_GQ95*br_ref)
        gamma_inv_over_M = float(x.Gamma_nu_GeV) / float(x.M_GeV)
        gchi_eq_massless = math.sqrt(max(0.0, 12.0*math.pi*gamma_inv_over_M))
        rows.append(dict(
            model=x.model,
            convention=x.convention,
            M_GeV=float(x.M_GeV),
            gq=gq,
            BR_invisible=br,
            Gamma_over_M=float(x.Gamma_over_M),
            gchi_equiv_massless_vector=gchi_eq_massless,
            CMS_anchor_M_GeV=CMS_M_GEV,
            CMS_anchor_gq95=CMS_GQ95,
            CMS_anchor_BRinv=br_ref,
            critical_gq95_at_50GeV=critical,
            required_weakening_vs_CMS100=critical/CMS_GQ95-1.0,
            NWA_signal_ratio_to_CMS100_boundary=mu_vs_100,
            decision_rule=(f"excluded if observed gq95(50 GeV) < {critical:.6f}"),
        ))

    R=pd.DataFrame(rows)
    R.to_csv(OUT/'monojet_nwa_mapping.csv',index=False)

    # Convention-independent summary: use the full bracket, not a cherry-picked width.
    summary=[]
    for model,q in R.groupby('model'):
        lo=float(q.critical_gq95_at_50GeV.min())
        hi=float(q.critical_gq95_at_50GeV.max())
        summary.append(dict(
            model=model,
            critical_gq95_min=lo,
            critical_gq95_max=hi,
            CMS100_anchor=CMS_GQ95,
            minimum_weakening_needed_to_survive=lo/CMS_GQ95-1.0,
            maximum_weakening_needed_to_survive=hi/CMS_GQ95-1.0,
            status='EXACT_50GEV_LIMIT_NEEDED',
        ))
    S=pd.DataFrame(summary)
    S.to_csv(OUT/'monojet_gate_summary.csv',index=False)

    meta={
      'version':'v0.46',
      'scope':'public monojet narrow-width rescaling screen; not a detector-level recast',
      'input':'results/v0.45/width_br_dijet_gate.csv',
      'cms_reference':{
        'analysis':'CMS EXO-20-004 / JHEP 11 (2021) 153',
        'luminosity_fb_inv':137,
        'mmed_GeV':100,
        'mDM_over_mmed':1/3,
        'gchi':1.0,
        'observed_gq95':0.018,
        'BRinv_computed':br_ref,
      },
      'nwa_mapping':'sigma(pp->Z+j)*BRinv proportional to gq^2*BRinv for a narrow on-shell mediator',
      'important_caveat':'The published exact 0.018 anchor is at 100 GeV. v0.46 does not invent a 50-GeV CMS point; it outputs the exact 50-GeV gq95 threshold needed to decide each PKR convention.',
      'next_if_marginal':'Obtain/digitize the exact low-mass coupling curve or perform the validated CMS EXO-20-004 MadAnalysis/simplified-likelihood recast with generated pp->Z+j events.',
    }
    (OUT/'metadata.json').write_text(json.dumps(meta,indent=2))

    print('\n=== v0.46 PKR monojet NWA mapping ===')
    print(R[['model','convention','BR_invisible','gchi_equiv_massless_vector',
             'critical_gq95_at_50GeV','required_weakening_vs_CMS100',
             'NWA_signal_ratio_to_CMS100_boundary']].to_string(index=False))
    print('\n=== v0.46 monojet gate summary ===')
    print(S.to_string(index=False))

if __name__=='__main__':
    main()
