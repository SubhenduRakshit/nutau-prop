#!/usr/bin/env python3
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np,pandas as pd

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--indir',required=True);ap.add_argument('--outdir',required=True);a=ap.parse_args();src=Path(a.indir);out=Path(a.outdir);out.mkdir(parents=True,exist_ok=True)
 fs=sorted(src.rglob('summary_*.csv'))
 if not fs:raise RuntimeError('no convergence summaries')
 r=pd.concat([pd.read_csv(f) for f in fs],ignore_index=True);r.to_csv(out/'convergence.csv',index=False)
 fine=r[r.setting=='fine'].set_index(['model','topology']);med=r[r.setting=='medium'].set_index(['model','topology']);rows=[]
 for idx in fine.index:
  f=fine.loc[idx];m=med.loc[idx]
  rows.append(dict(model=idx[0],topology=idx[1],
   rel_change_reco_1p259=abs(f.reco_frac_1p259PeV-m.reco_frac_1p259PeV)/max(abs(f.reco_frac_1p259PeV),1e-12),
   rel_change_weighted_horizon=abs(f.SM_event_weighted_abs_horizon_reco-m.SM_event_weighted_abs_horizon_reco)/max(abs(f.SM_event_weighted_abs_horizon_reco),1e-12),
   rel_change_true_1PeV=abs(f.tau_frac_1PeV-m.tau_frac_1PeV)/max(abs(f.tau_frac_1PeV),1e-12),
   rel_change_true_1p5PeV=abs(f.tau_frac_1p5PeV-m.tau_frac_1p5PeV)/max(abs(f.tau_frac_1p5PeV),1e-12)))
 c=pd.DataFrame(rows);c.to_csv(out/'medium_vs_fine.csv',index=False)
 for f in src.rglob('true_*.csv'):
  (out/f.name).write_bytes(f.read_bytes())
 (out/'metadata.json').write_text(json.dumps(dict(version='v0.42cp',mu=.70,energy_min_GeV=1e6,energy_max_GeV=3e8,criterion='medium-to-fine relative change <=5% preferred'),indent=2))
 print(r.to_string(index=False));print('\nMEDIUM vs FINE');print(c.to_string(index=False))
if __name__=='__main__':main()
