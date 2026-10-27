#!/usr/bin/env python3
import argparse, json, subprocess, sys
from pathlib import Path
import numpy as np
import pandas as pd
from nutau_prop.config import load_config
from nutau_prop.tables import load_table
from nutau_prop.transport import TauTransport
from nutau_prop.util import powerlaw_per_logE

def ensure_table(path, args):
    if path.exists(): return
    cmd=[sys.executable,"scripts/precompute_xsec.py","--out",str(path)]+args
    subprocess.run(cmd,check=True)

def transport(table,c,mu):
    t=load_table(table)
    tr=TauTransport(t,float(c["physics"]["tau_beta_cm2_per_g"]),
                    float(c["physics"]["tau_polarization"]))
    p0=powerlaw_per_logE(t["energies"],float(c["physics"]["spectral_index"]))
    r=tr.propagate(mu,p0,max_step_km=float(c["transport"]["max_step_km"]),
                  max_p=float(c["transport"]["max_interaction_probability"]))
    return t["energies"],p0,r

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--config",default="config/default.yaml")
    ap.add_argument("--outdir",default="results")
    args=ap.parse_args()
    c=load_config(args.config); out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)
    sm=Path("cache/sm.npz"); ensure_table(sm,["--kind","sm"])
    rows=[]
    models=[]
    for r2 in c["benchmarks"]["charge_radius_cm2"]:
        if float(r2)==0: continue
        key=f"cr_{float(r2):+.1e}".replace("+","p").replace("-","m")
        p=Path("cache")/(key+".npz")
        ensure_table(p,["--kind","charge_radius",f"--r2={float(r2)}"])
        models.append((key,p,dict(kind="charge_radius",r2_cm2=float(r2))))
    for eps in c["benchmarks"]["derivative"]["eps_em_inf"]:
        for mx in c["benchmarks"]["derivative"]["mediator_mass_gev"]:
            key=f"deriv_eps{eps}_mx{mx}"
            p=Path("cache")/(key+".npz")
            ensure_table(p,["--kind","derivative_em","--eps-em-inf",str(eps),"--mx",str(mx)])
            models.append((key,p,dict(kind="derivative_em",eps_em_inf=eps,mx_gev=mx)))
    for mu in c["transport"]["cos_zenith_abs"]:
        E,p0,rsm=transport(sm,c,float(mu))
        for key,p,meta in models:
            _,_,rb=transport(p,c,float(mu))
            frac=np.divide(rb["nu"]-rsm["nu"],rsm["nu"],
                           out=np.full_like(E,np.nan),where=rsm["nu"]>0)
            trans=np.divide(rsm["nu"],p0,out=np.zeros_like(E),where=p0>0)
            for ee,ff,tt in zip(E,frac,trans):
                rows.append(dict(model=key,cosz_abs=float(mu),E_GeV=ee,
                                 frac_nu=ff,sm_transmission=tt,**meta))
    df=pd.DataFrame(rows)
    df.to_csv(out/"bsm_scan.csv",index=False)
    print(out/"bsm_scan.csv")

if __name__=="__main__":
    main()
