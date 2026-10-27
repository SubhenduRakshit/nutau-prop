#!/usr/bin/env python3
import argparse,json,math
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar

from nutau_prop.config import load_config
from nutau_prop.icecube import (
    fetch_icecube_2015,parse_total_effective_area,powerlaw_integral,
    solid_angle_from_cosz,poisson_asimov_q
)

def interp_model(scan,model,E,mu,col):
    d=scan[scan.model==model]
    mus=np.sort(d.cosz_abs.unique())
    # interpolate in logE at each available mu, then across mu.
    vals=[]
    for m in mus:
        q=d[d.cosz_abs==m].sort_values("E_GeV")
        vals.append(np.interp(np.log10(E),np.log10(q.E_GeV),q[col],
                              left=q[col].iloc[0],right=q[col].iloc[-1]))
    return float(np.interp(mu,mus,vals,left=vals[0],right=vals[-1]))

def bin_sm_counts(row,cfg,gamma=None):
    ice=cfg["icecube"]; flux=cfg["flux"]
    g=float(flux["gamma"] if gamma is None else gamma)
    elo=max(float(row.E_lo_GeV),float(flux["energy_min_gev"]))
    ehi=min(float(row.E_hi_GeV),float(flux["energy_max_gev"]))
    if ehi<=elo: return 0.0
    fint=powerlaw_integral(float(flux["phi0_per_gev_cm2_s_sr"]),
                           float(flux["pivot_gev"]),g,elo,ehi)
    domega=solid_angle_from_cosz(row.cosz_lo,row.cosz_hi)
    T=float(ice["exposure_days"])*86400.0
    return fint*domega*T*float(row.Aeff_m2)*1e4

def mu_for_earth(row):
    c=0.5*(float(row.cosz_lo)+float(row.cosz_hi))
    return max(0.0,-c)

def energy_center(row):
    return math.sqrt(float(row.E_lo_GeV)*float(row.E_hi_GeV))

def profiled_norm_q(true,test):
    st=float(np.sum(true)); ss=float(np.sum(test))
    if ss<=0: return np.inf,1.0
    alpha=st/ss
    return poisson_asimov_q(true,alpha*test),alpha

def profile_norm_tilt(true,aeff,cfg):
    amin=float(cfg["analysis"]["profile_gamma_min"])
    amax=float(cfg["analysis"]["profile_gamma_max"])
    def objective(gamma):
        sm=np.array([bin_sm_counts(r,cfg,gamma) for _,r in aeff.iterrows()])
        q,_=profiled_norm_q(true,sm)
        return q
    res=minimize_scalar(objective,bounds=(amin,amax),method="bounded",
                        options={"xatol":1e-4})
    sm=np.array([bin_sm_counts(r,cfg,res.x) for _,r in aeff.iterrows()])
    q,alpha=profiled_norm_q(true,sm)
    return float(q),float(alpha),float(res.x)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--config",default="config/detector.yaml")
    ap.add_argument("--scan",default="results/bsm_scan_combined.csv")
    ap.add_argument("--data-dir",default="data/icecube_2015")
    ap.add_argument("--outdir",default="results")
    ap.add_argument("--no-download",action="store_true")
    a=ap.parse_args()
    cfg=load_config(a.config); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)

    eff=Path(a.data_dir)/cfg["icecube"]["effective_area_file"]
    if not eff.exists():
        hits=list(Path(a.data_dir).rglob(cfg["icecube"]["effective_area_file"]))
        if hits: eff=hits[0]
    if not eff.exists():
        if a.no_download:
            raise FileNotFoundError(f"Missing {eff}")
        eff=fetch_icecube_2015(a.data_dir,cfg["icecube"]["release_url"])

    A=parse_total_effective_area(eff)
    # Restrict to bins overlapping the published astrophysical-fit interval.
    emin=float(cfg["flux"]["energy_min_gev"]); emax=float(cfg["flux"]["energy_max_gev"])
    A=A[(A.E_hi_GeV>emin)&(A.E_lo_GeV<emax)&(A.Aeff_m2>0)].copy().reset_index(drop=True)
    A["E_center_GeV"]=[energy_center(r) for _,r in A.iterrows()]
    A["cosz_center"]=0.5*(A.cosz_lo+A.cosz_hi)
    A["mu_earth"]=[mu_for_earth(r) for _,r in A.iterrows()]
    A["N_SM"]=[bin_sm_counts(r,cfg) for _,r in A.iterrows()]

    scan=pd.read_csv(a.scan)
    models=sorted(scan.model.unique())
    detail=[]; summary=[]
    projection=float(cfg["analysis"]["projection_exposure_years"])
    scale_proj=projection*365.25/float(cfg["icecube"]["exposure_days"])

    for model in models:
        prop=[]; det=[]; comb=[]
        for _,r in A.iterrows():
            E=float(r.E_center_GeV); mu=float(r.mu_earth)
            if mu<=0:
                rp=1.0
            else:
                rp=interp_model(scan,model,E,mu,"prop_ratio")
            rd=interp_model(scan,model,E,max(mu,0.2),"detector_total_xsec_ratio")
            prop.append(rp); det.append(rd); comb.append(rp*rd)
        nsm=A.N_SM.to_numpy(float)
        nb_prop=nsm*np.array(prop)
        nb=nsm*np.array(comb)

        qfix=poisson_asimov_q(nb,nsm)
        qnorm,alpha=profiled_norm_q(nb,nsm)
        qtilt,alpha_t,gamma_t=profile_norm_tilt(nb,A,cfg)

        # Separate propagation-only diagnostic.
        qprop=poisson_asimov_q(nb_prop,nsm)
        qprop_norm,_=profiled_norm_q(nb_prop,nsm)

        for i,r in A.iterrows():
            detail.append(dict(
                model=model,
                E_lo_GeV=r.E_lo_GeV,E_hi_GeV=r.E_hi_GeV,
                cosz_lo=r.cosz_lo,cosz_hi=r.cosz_hi,
                mu_earth=r.mu_earth,Aeff_m2=r.Aeff_m2,
                N_SM=nsm[i],
                propagation_ratio=prop[i],
                detector_total_xsec_proxy=det[i],
                total_BSM_ratio=comb[i],
                N_BSM_prop_only=nb_prop[i],
                N_BSM=nb[i],
                deltaN=nb[i]-nsm[i],
                frac=(nb[i]/nsm[i]-1) if nsm[i]>0 else np.nan
            ))

        z=math.sqrt(max(qfix,0)); zn=math.sqrt(max(qnorm,0)); zt=math.sqrt(max(qtilt,0))
        summary.append(dict(
            model=model,
            exposure_days=float(cfg["icecube"]["exposure_days"]),
            N_SM_total=float(nsm.sum()),
            N_BSM_total=float(nb.sum()),
            total_fractional_shift=float(nb.sum()/nsm.sum()-1) if nsm.sum()>0 else np.nan,
            Z_fixed_flux=z,
            Z_profile_norm=zn,
            Z_profile_norm_tilt=zt,
            bestfit_norm_profile=float(alpha),
            bestfit_norm_tilt=float(alpha_t),
            bestfit_gamma=float(gamma_t),
            Z_prop_only_fixed=math.sqrt(max(qprop,0)),
            Z_prop_only_profile_norm=math.sqrt(max(qprop_norm,0)),
            projected_years=projection,
            projected_Z_profile_norm_tilt=zt*math.sqrt(scale_proj),
            years_for_2sigma=(float(cfg["icecube"]["exposure_days"])/365.25)*(2/zt)**2 if zt>0 else np.inf,
            years_for_3sigma=(float(cfg["icecube"]["exposure_days"])/365.25)*(3/zt)**2 if zt>0 else np.inf,
            max_abs_bin_fraction=float(np.nanmax(np.abs(nb/np.maximum(nsm,1e-300)-1))),
            max_abs_delta_events=float(np.max(np.abs(nb-nsm))),
            bins_with_SM_ge_0p01=int(np.sum(nsm>=0.01)),
            note="Detector factor is a total-cross-section proxy; published total Aeff is nu/antinu averaged."
        ))

    detail_df=pd.DataFrame(detail)
    summ=pd.DataFrame(summary).sort_values("Z_profile_norm_tilt",ascending=False)
    detail_df.to_csv(out/"event_weighted_counts.csv",index=False)
    summ.to_csv(out/"event_weighted_summary.csv",index=False)

    payload={
      "icecube_release":cfg["icecube"]["release_name"],
      "effective_area_source":cfg["icecube"]["release_url"],
      "effective_area_energy_format":A.attrs.get("energy_format"),
      "effective_area_angle_format":A.attrs.get("angle_format"),
      "flux":cfg["flux"],
      "important_method_note":"IceCube Aeff already contains SM Earth absorption; propagation enters only through T_BSM/T_SM.",
      "detector_reweight_note":"BSM detector response is approximated by the inclusive CC+NC cross-section ratio. Exact CC/NC/topology reweighting is a later refinement.",
      "models":summ.to_dict(orient="records")
    }
    (out/"event_weighted_summary.json").write_text(json.dumps(payload,indent=2))

    # Plot the most sensitive model: SM vs BSM counts summed over zenith per true-E bin.
    if len(summ):
        best=summ.iloc[0].model
        d=detail_df[detail_df.model==best].copy()
        grp=d.groupby(["E_lo_GeV","E_hi_GeV"],as_index=False)[["N_SM","N_BSM"]].sum()
        Ec=np.sqrt(grp.E_lo_GeV*grp.E_hi_GeV)
        fig,ax=plt.subplots()
        ax.step(Ec,grp.N_SM,where="mid",label="SM")
        ax.step(Ec,grp.N_BSM,where="mid",label=f"BSM: {best}")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("true neutrino energy [GeV]")
        ax.set_ylabel("expected tau-flavor events / bin")
        ax.grid(True,which="both",alpha=.25); ax.legend()
        fig.tight_layout(); fig.savefig(out/"event_weighted_best_model.png",dpi=180)

    print(summ[[
        "model","N_SM_total","N_BSM_total","total_fractional_shift",
        "Z_fixed_flux","Z_profile_norm","Z_profile_norm_tilt",
        "projected_Z_profile_norm_tilt","years_for_2sigma"
    ]].to_string(index=False))
    print(out/"event_weighted_summary.json")

if __name__=="__main__":
    main()
