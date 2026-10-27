#!/usr/bin/env python3
import argparse,json
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np,pandas as pd
import matplotlib.pyplot as plt

from nutau_prop.config import load_config
from nutau_prop.tables import load_table
from nutau_prop.transport import TauTransport
from nutau_prop.historical import RRZIterator,historical_flux_per_logE,nadir_to_mu,relative_max
from nutau_prop.util import nproc_default

def _one_species(table,phi0,mu,P,c):
    tr=TauTransport(table,beta_tau=0.0,tau_polarization=P,
                    secondary_branching=float(c["physics"]["secondary_numu_branching"]))
    direct=tr.propagate(mu,phi0,
        max_step_km=float(c["transport"]["max_step_km"]),
        max_p=float(c["transport"]["max_interaction_probability"]))["nu"]
    it=RRZIterator(table,P).solve(mu,phi0,
        depth_steps=int(c["historical"]["z_depth_steps"]),
        n_iter=int(c["historical"]["z_iterations"]))
    return direct,it

def _case(payload):
    gamma,theta,nu_path,anu_path,c=payload
    tn=load_table(nu_path); ta=load_table(anu_path)
    E=tn["energies"]
    phi0=historical_flux_per_logE(E,gamma,c["historical"]["cutoff_gev"])
    mu=nadir_to_mu(theta)

    dn,zn=_one_species(tn,phi0,mu,-1.0,c)  # nu_tau -> tau-
    da,za=_one_species(ta,phi0,mu,+1.0,c)  # anti-nu_tau -> tau+

    direct=(dn+da)/(2*phi0)
    Ftot=[(zn["F"][n][-1]+za["F"][n][-1])/2 for n in range(4)]
    ratioZ=[f/phi0 for f in Ftot]

    # "Relevant" bins: published energy range, with non-negligible final flux.
    mask=(E>=1e3)&(E<=3e7)&(ratioZ[2]>1e-5)
    f32=relative_max(Ftot[3],Ftot[2],mask)
    f21=relative_max(Ftot[2],Ftot[1],mask)
    d2=relative_max((dn+da)/2,Ftot[2],mask)

    return dict(gamma=float(gamma),theta=float(theta),E=E,direct=direct,
                ratioZ=ratioZ,f32=f32,f21=f21,direct_z2=d2)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--config",default="config/historical.yaml")
    ap.add_argument("--nu-table",default="cache/historical_nu.npz")
    ap.add_argument("--antinu-table",default="cache/historical_antinu.npz")
    ap.add_argument("--outdir",default="results/historical")
    ap.add_argument("--workers",type=int,default=nproc_default())
    a=ap.parse_args(); c=load_config(a.config)
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)

    cases=[(g,t,a.nu_table,a.antinu_table,c)
           for g in (1.0,2.0) for t in (0.0,30.0,60.0)]
    with ProcessPoolExecutor(max_workers=min(int(a.workers),len(cases))) as ex:
        res=list(ex.map(_case,cases))

    rows=[]; zrows=[]
    fig,axes=plt.subplots(2,3,figsize=(12,7),sharex=True)
    for i,r in enumerate(res):
        E=r["E"]; ax=axes[i//3,i%3]
        ax.loglog(E,r["direct"],label="direct")
        ax.loglog(E,r["ratioZ"][1],"--",label="RR-like Z1")
        ax.loglog(E,r["ratioZ"][2],"-.",label="RR-like Z2")
        ax.loglog(E,r["ratioZ"][3],":",label="RR-like Z3")
        ax.axhline(1,lw=.6)
        ax.set_title(f"E^-{r['gamma']:.0f}, theta={r['theta']:.0f} deg")
        ax.grid(True,which="both",alpha=.2)
        for e,d,z1,z2,z3 in zip(E,r["direct"],r["ratioZ"][1],r["ratioZ"][2],r["ratioZ"][3]):
            rows.append(dict(gamma=r["gamma"],nadir_deg=r["theta"],E_GeV=e,
                             direct_ratio=d,Z1_ratio=z1,Z2_ratio=z2,Z3_ratio=z3))
    axes[0,0].legend(fontsize=7)
    for ax in axes[-1,:]: ax.set_xlabel("E [GeV]")
    for ax in axes[:,0]: ax.set_ylabel("(nu_tau + antinu_tau) / incident")
    fig.tight_layout(); fig.savefig(out/"rr_reproduction.png",dpi=180)
    pd.DataFrame(rows).to_csv(out/"rr_reproduction.csv",index=False)

    hard0=next(r for r in res if r["gamma"]==1 and r["theta"]==0)
    soft0=next(r for r in res if r["gamma"]==2 and r["theta"]==0)
    bump=(hard0["E"]>=1e4)&(hard0["E"]<=1e5)
    hard_bump=float(np.max(hard0["direct"][bump]))
    soft_bump=float(np.max(soft0["direct"][bump]))

    d={
      "max_abs_F3_over_F2_minus1":max(r["f32"] for r in res),
      "hard_Eminus1_vertical_max_abs_F2_over_F1_minus1":hard0["f21"],
      "max_direct_vs_Z2_fraction":max(r["direct_z2"] for r in res),
      "hard_Eminus1_vertical_bump_1e4_1e5":hard_bump,
      "soft_Eminus2_vertical_max_1e4_1e5":soft_bump,
      "criteria":{
        "F3_over_F2_lt":float(c["historical"]["rr_F3_over_F2_max_fraction"]),
        "F2_over_F1_gt":float(c["historical"]["rr_F2_over_F1_min_effect"]),
        "F2_over_F1_lt":float(c["historical"]["rr_F2_over_F1_max_effect"]),
        "direct_vs_Z2_lt":float(c["historical"]["rr_direct_vs_Z2_max_fraction"]),
        "hard_bump_gt":float(c["historical"]["rr_hard_bump_min"]),
        "soft_bump_lt":float(c["historical"]["rr_soft_bump_max"])
      }
    }
    d["pass_convergence"]=d["max_abs_F3_over_F2_minus1"]<d["criteria"]["F3_over_F2_lt"]
    d["pass_second_order"]=(
        d["hard_Eminus1_vertical_max_abs_F2_over_F1_minus1"]>d["criteria"]["F2_over_F1_gt"]
        and d["hard_Eminus1_vertical_max_abs_F2_over_F1_minus1"]<d["criteria"]["F2_over_F1_lt"])
    d["pass_independent_solver"]=d["max_direct_vs_Z2_fraction"]<d["criteria"]["direct_vs_Z2_lt"]
    d["pass_hard_bump"]=d["hard_Eminus1_vertical_bump_1e4_1e5"]>d["criteria"]["hard_bump_gt"]
    d["pass_soft_no_bump"]=d["soft_Eminus2_vertical_max_1e4_1e5"]<d["criteria"]["soft_bump_lt"]
    d["PASS"]=all(d[k] for k in
        ("pass_convergence","pass_second_order","pass_independent_solver",
         "pass_hard_bump","pass_soft_no_bump"))
    (out/"rr_gate.json").write_text(json.dumps(d,indent=2))
    print(json.dumps(d,indent=2))
    if not d["PASS"]:
        raise SystemExit("RR historical gate FAILED")
if __name__=="__main__": main()
