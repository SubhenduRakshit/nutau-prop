#!/usr/bin/env python3
import argparse, subprocess, sys
from pathlib import Path
import numpy as np, pandas as pd
from nutau_prop.config import load_config
from nutau_prop.tables import load_table
from nutau_prop.transport import TauTransport
from nutau_prop.util import powerlaw_per_logE

def ensure(path,args,config="config/default.yaml",workers=8):
    if path.exists(): return
    cmd=[sys.executable,"scripts/precompute_xsec.py","--config",config,
         "--workers",str(workers),"--out",str(path)]+list(args)
    subprocess.run(cmd,check=True)

def propagate(path,c,mu,P):
    t=load_table(path)
    tr=TauTransport(t,float(c["physics"]["tau_beta_cm2_per_g"]),P)
    p0=powerlaw_per_logE(t["energies"],float(c["physics"]["spectral_index"]))
    r=tr.propagate(mu,p0,max_step_km=float(c["transport"]["max_step_km"]),
                   max_p=float(c["transport"]["max_interaction_probability"]))
    return t,p0,r["nu"]

def make_models(c):
    models=[]
    for r2 in c["benchmarks"]["charge_radius_cm2"]:
        if float(r2)==0: continue
        key=f"cr_{float(r2):+.1e}".replace("+","p").replace("-","m")
        models.append((key,["--kind","charge_radius",f"--r2={float(r2)}"],
                       dict(kind="charge_radius",r2_cm2=float(r2))))
    for eps in c["benchmarks"]["derivative"]["eps_em_inf"]:
        for mx in c["benchmarks"]["derivative"]["mediator_mass_gev"]:
            key=f"deriv_eps{eps}_mx{mx}"
            models.append((key,["--kind","derivative_em","--eps-em-inf",str(eps),"--mx",str(mx)],
                           dict(kind="derivative_em",eps_em_inf=float(eps),mx_gev=float(mx))))
    return models

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--config",default="config/default.yaml")
    ap.add_argument("--out",default="results/bsm_scan_combined.csv")
    ap.add_argument("--workers",type=int,default=8)
    a=ap.parse_args()
    c=load_config(a.config); Path("cache").mkdir(exist_ok=True); Path("results").mkdir(exist_ok=True)

    smn=Path("cache/sm.npz")
    sma=Path("cache/sm_antinu.npz")
    ensure(smn,["--kind","sm"],a.config,a.workers)
    ensure(sma,["--kind","sm","--antinu"],a.config,a.workers)

    model_paths=[]
    for key,args,meta in make_models(c):
        pn=Path("cache")/(key+".npz")
        pa=Path("cache")/(key+"_antinu.npz")
        ensure(pn,args,a.config,a.workers)
        ensure(pa,args+["--antinu"],a.config,a.workers)
        model_paths.append((key,pn,pa,meta))

    rows=[]
    for mu in map(float,c["transport"]["cos_zenith_abs"]):
        tsn,p0,fsn=propagate(smn,c,mu,-1.0)
        tsa,_,fsa=propagate(sma,c,mu,+1.0)
        for key,pn,pa,meta in model_paths:
            tbn,_,fbn=propagate(pn,c,mu,-1.0)
            tba,_,fba=propagate(pa,c,mu,+1.0)

            # Equal incident nu and antinu fluxes.
            prop_ratio=(fbn+fba)/np.maximum(fsn+fsa,1e-300)

            # Detector interaction proxy.  CC is unchanged; NC is modified.
            rn=(tbn["sigma_cc"]+tbn["sigma_nc"])/np.maximum(tsn["sigma_cc"]+tsn["sigma_nc"],1e-300)
            ra=(tba["sigma_cc"]+tba["sigma_nc"])/np.maximum(tsa["sigma_cc"]+tsa["sigma_nc"],1e-300)
            det_ratio=0.5*(rn+ra)
            nc_ratio=0.5*(tbn["sigma_nc"]/np.maximum(tsn["sigma_nc"],1e-300)
                         +tba["sigma_nc"]/np.maximum(tsa["sigma_nc"],1e-300))

            for i,E in enumerate(tsn["energies"]):
                rows.append(dict(
                    model=key,cosz_abs=mu,E_GeV=float(E),
                    prop_ratio=float(prop_ratio[i]),
                    prop_frac=float(prop_ratio[i]-1.0),
                    detector_total_xsec_ratio=float(det_ratio[i]),
                    detector_nc_ratio=float(nc_ratio[i]),
                    combined_proxy_ratio=float(prop_ratio[i]*det_ratio[i]),
                    sm_trans_nu=float(fsn[i]/p0[i]),
                    sm_trans_antinu=float(fsa[i]/p0[i]),
                    **meta))
    pd.DataFrame(rows).to_csv(a.out,index=False)
    print(a.out)

if __name__=="__main__":
    main()
