#!/usr/bin/env python3
import argparse,json
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np,pandas as pd
import matplotlib.pyplot as plt

from nutau_prop.config import load_config
from nutau_prop.tables import load_table
from nutau_prop.transport import TauTransport,NuOnlyTransport
from nutau_prop.historical import historical_flux_per_logE,nadir_to_mu,interp_logE
from nutau_prop.util import nproc_default

def _species(table,phi0,mu,P,c):
    t=TauTransport(table,beta_tau=0.0,tau_polarization=P,
                   secondary_branching=float(c["physics"]["secondary_numu_branching"]))
    full=t.propagate(mu,phi0,
        max_step_km=float(c["transport"]["max_step_km"]),
        max_p=float(c["transport"]["max_interaction_probability"]),
        track_secondary_numu=True)
    primary=NuOnlyTransport(table).propagate(mu,phi0,
        max_step_km=float(c["transport"]["max_step_km"]),
        max_p=float(c["transport"]["max_interaction_probability"]))["nu"]
    return full,primary

def _case(payload):
    gamma,theta,nu_path,anu_path,c=payload
    tn=load_table(nu_path); ta=load_table(anu_path); E=tn["energies"]
    p=historical_flux_per_logE(E,gamma,c["historical"]["cutoff_gev"])
    mu=nadir_to_mu(theta)

    # nu_tau creates tau- -> anti-nu_mu; anti-nu_tau creates tau+ -> nu_mu.
    fn,pn=_species(tn,p,mu,-1.0,c)
    fa,pa=_species(ta,p,mu,+1.0,c)

    nutau=(fn["nu"]+fa["nu"])/(2*p)
    primary=(pn+pa)/(2*p)
    secondary=(fn["numu"]+fa["numu"])/(2*p)

    # DRS Eq. (5), applied separately then combined.  It intentionally
    # overestimates slightly because the secondary spectrum is softer.
    B=float(c["physics"]["secondary_numu_branching"])
    eq5=B*np.maximum(((fn["nu"]-pn)+(fa["nu"]-pa))/(2*p),0.0)

    mask=(eq5>1e-5)&(E<=3e7)
    frac=np.abs(secondary-eq5)/np.maximum(eq5,1e-30)
    med=float(np.median(frac[mask])) if np.any(mask) else np.inf
    return dict(gamma=float(gamma),theta=float(theta),E=E,nutau=nutau,
                primary=primary,secondary=secondary,eq5=eq5,eq5_med=med)

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

    rows=[]; fig,axes=plt.subplots(2,3,figsize=(12,7),sharex=True)
    for i,r in enumerate(res):
        ax=axes[i//3,i%3]; E=r["E"]
        ax.loglog(E,r["nutau"],"--",label="nu_tau total")
        ax.loglog(E,r["primary"],":",label="primary nu_mu/e")
        ax.loglog(E,r["secondary"],"-",label="secondary direct")
        ax.loglog(E,r["eq5"],"-.",label="DRS Eq.(5)")
        ax.set_title(f"E^-{r['gamma']:.0f}, theta={r['theta']:.0f} deg")
        ax.grid(True,which="both",alpha=.2)
        for e,x,y,z,q in zip(E,r["nutau"],r["primary"],r["secondary"],r["eq5"]):
            rows.append(dict(gamma=r["gamma"],nadir_deg=r["theta"],E_GeV=e,
                             nutau_ratio=x,primary_numu_ratio=y,
                             secondary_ratio=z,DRS_eq5_ratio=q))
    axes[0,0].legend(fontsize=7)
    for ax in axes[-1,:]: ax.set_xlabel("E [GeV]")
    for ax in axes[:,0]: ax.set_ylabel("Flux / incident flux")
    fig.tight_layout(); fig.savefig(out/"drs_fig1_reproduction.png",dpi=180)
    pd.DataFrame(rows).to_csv(out/"drs_fig1_reproduction.csv",index=False)

    hard0=next(r for r in res if r["gamma"]==1 and r["theta"]==0)
    soft0=next(r for r in res if r["gamma"]==2 and r["theta"]==0)
    nutau_1p=interp_logE(hard0["E"],hard0["nutau"],1e6)
    numu_1p=interp_logE(hard0["E"],hard0["primary"],1e6)

    d={
      "Eminus1_theta0_nutau_ratio_at_1PeV":nutau_1p,
      "Eminus1_theta0_primary_numu_ratio_at_1PeV":numu_1p,
      "max_median_fraction_direct_secondary_vs_DRS_Eq5":max(r["eq5_med"] for r in res),
      "Eminus1_theta0_peak_secondary_ratio":float(np.max(hard0["secondary"])),
      "Eminus2_theta0_peak_secondary_ratio":float(np.max(soft0["secondary"])),
      "criteria":{
        "nutau_1PeV_range":[float(c["historical"]["drs_nutau_1pev_min"]),
                            float(c["historical"]["drs_nutau_1pev_max"])],
        "numu_1PeV_range":[float(c["historical"]["drs_numu_1pev_min"]),
                           float(c["historical"]["drs_numu_1pev_max"])],
        "Eq5_median_fraction_lt":float(c["historical"]["drs_eq5_median_fractional_tolerance"])
      }
    }
    lo,hi=d["criteria"]["nutau_1PeV_range"]
    d["pass_nutau_1PeV"]=lo<nutau_1p<hi
    lo,hi=d["criteria"]["numu_1PeV_range"]
    d["pass_numu_1PeV"]=lo<numu_1p<hi
    d["pass_eq5"]=d["max_median_fraction_direct_secondary_vs_DRS_Eq5"]<d["criteria"]["Eq5_median_fraction_lt"]
    d["pass_secondary_spectral_hierarchy"]=(
        d["Eminus2_theta0_peak_secondary_ratio"] <
        d["Eminus1_theta0_peak_secondary_ratio"])
    d["PASS"]=all(d[k] for k in
        ("pass_nutau_1PeV","pass_numu_1PeV","pass_eq5","pass_secondary_spectral_hierarchy"))
    (out/"drs_gate.json").write_text(json.dumps(d,indent=2))
    print(json.dumps(d,indent=2))
    if not d["PASS"]:
        raise SystemExit("DRS historical gate FAILED")
if __name__=="__main__": main()
