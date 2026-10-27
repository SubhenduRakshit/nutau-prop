#!/usr/bin/env python3
import argparse, json
from pathlib import Path
import pandas as pd
from nutau_prop.config import load_config

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--config",default="config/default.yaml")
    ap.add_argument("--scan",default="results/bsm_scan.csv")
    args=ap.parse_args()
    c=load_config(args.config); g=c["gate"]
    df=pd.read_csv(args.scan)
    sel=df[
      (df.E_GeV>=float(g["energy_min_gev"])) &
      (df.E_GeV<=float(g["energy_max_gev"])) &
      (df.sm_transmission>=float(g["min_sm_transmission"]))
    ].copy()
    sel["absfrac"]=sel.frac_nu.abs()
    best=sel.sort_values("absfrac",ascending=False).head(20)
    threshold=float(g["minimum_distortion"])
    passed=bool((sel.absfrac>=threshold).any())
    out={
      "passed":passed,
      "threshold":threshold,
      "max_abs_distortion":float(sel.absfrac.max()) if len(sel) else None,
      "best_rows":best.to_dict(orient="records")
    }
    Path("results").mkdir(exist_ok=True)
    Path("results/gate.json").write_text(json.dumps(out,indent=2))
    print("PASS" if passed else "FAIL")
    if len(best):
        print(best[["model","cosz_abs","E_GeV","frac_nu","sm_transmission"]].to_string(index=False))
    print("Wrote results/gate.json")

if __name__=="__main__":
    main()
