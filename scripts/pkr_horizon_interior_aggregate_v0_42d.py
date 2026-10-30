#!/usr/bin/env python3
from __future__ import annotations
import argparse,json
from pathlib import Path
import pandas as pd

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--indir',required=True); ap.add_argument('--outdir',required=True); a=ap.parse_args()
    src=Path(a.indir); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    fs=sorted(src.rglob('summary_*.csv'))
    if not fs: raise RuntimeError('no summaries')
    r=pd.concat([pd.read_csv(f) for f in fs],ignore_index=True); r.to_csv(out/'interior_convergence.csv',index=False)
    med=r[r.setting=='medium'].set_index(['model','topology']); fine=r[r.setting=='fine'].set_index(['model','topology'])
    rows=[]
    for idx in fine.index:
        m=med.loc[idx]; f=fine.loc[idx]
        rows.append(dict(model=idx[0],topology=idx[1],
            fine_Erec_GeV=f.target_Erec_GeV,fine_tau_SM_eff=f.target_tau_SM_eff,
            fine_reco_frac=f.target_reco_frac,
            medium_to_fine_reco=abs(f.target_reco_frac-m.target_reco_frac)/max(abs(f.target_reco_frac),1e-12),
            medium_to_fine_weighted=abs(f.weighted_abs_horizon_reco-m.weighted_abs_horizon_reco)/max(abs(f.weighted_abs_horizon_reco),1e-12),
            fine_true_tau_frac_3PeV=f.true_tau_frac_3PeV,fine_true_tau_frac_5PeV=f.true_tau_frac_5PeV,
            fine_sec_numu_frac_3PeV=f.sec_numu_frac_3PeV,fine_sec_numu_frac_5PeV=f.sec_numu_frac_5PeV))
    c=pd.DataFrame(rows); c.to_csv(out/'interior_certification.csv',index=False)
    (out/'metadata.json').write_text(json.dumps(dict(version='v0.42d',mu=.35,target_Erec_GeV=10**6.7,Emin_GeV=1e6,Emax_GeV=3e8,sigma_log10E=.30,criterion='medium-to-fine <=5%',boundary_note='target is 0.7 dex above Emin; two-sigma lower migration scale is approximately 1.26 PeV'),indent=2))
    print(c.to_string(index=False),flush=True)
if __name__=='__main__': main()
