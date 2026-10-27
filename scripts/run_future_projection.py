#!/usr/bin/env python3
import argparse, json, math
from pathlib import Path
import numpy as np, pandas as pd
from nutau_prop.config import load_config
from nutau_prop.icecube import powerlaw_integral, solid_angle_from_cosz

BANDS=[
    ("25-100_TeV",2.5e4,1.0e5),
    ("100-300_TeV",1.0e5,3.0e5),
    ("300-1000_TeV",3.0e5,1.0e6),
    ("1-1.4_PeV",1.0e6,1.4e6),
]

def sm_bin_at_gamma(row,cfg,gamma):
    flux=cfg["flux"]; ice=cfg["icecube"]
    elo=max(float(row.E_lo_GeV),float(flux["energy_min_gev"]))
    ehi=min(float(row.E_hi_GeV),float(flux["energy_max_gev"]))
    if ehi<=elo: return 0.0
    fint=powerlaw_integral(float(flux["phi0_per_gev_cm2_s_sr"]),
                           float(flux["pivot_gev"]),float(gamma),elo,ehi)
    domega=solid_angle_from_cosz(row.cosz_lo,row.cosz_hi)
    return fint*domega*float(ice["exposure_days"])*86400.0*float(row.Aeff_m2)*1e4

def qbin(true,test):
    if true>0 and test>0:
        return 2*(test-true+true*math.log(true/test))
    if true<=0 and test>0:
        return 2*test
    return 0.0

def first_meta(scan,model):
    q=scan[scan.model==model].iloc[0]
    d={"kind":q.get("kind","unknown")}
    for k in ("r2_cm2","eps_em_inf","mx_gev"):
        if k in q.index:
            try: d[k]=float(q[k])
            except Exception: d[k]=None
    return d

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--summary",default="results/future_current/event_weighted_summary.csv")
    ap.add_argument("--counts",default="results/future_current/event_weighted_counts.csv")
    ap.add_argument("--scan",default="results/bsm_scan_future.csv")
    ap.add_argument("--detector-config",default="config/detector.yaml")
    ap.add_argument("--future-config",default="config/future_detectors.yaml")
    ap.add_argument("--outdir",default="results/future")
    a=ap.parse_args()

    det=load_config(a.detector_config)
    fut=load_config(a.future_config)
    S=pd.read_csv(a.summary)
    C=pd.read_csv(a.counts)
    X=pd.read_csv(a.scan)
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)

    current_years=float(det["icecube"]["exposure_days"])/365.25
    gen2_years=float(fut["icecube_gen2"]["exposure_years"])
    gen2_base=float(fut["icecube_gen2"]["baseline_area_multiplier"])
    gen2_opt=float(fut["icecube_gen2"]["optimistic_area_multiplier"])
    tau_tag=float(fut["icecube_gen2"]["tau_tag_rate_multiplier"])
    cr_bound=float(fut["constraints"]["tau_charge_radius_abs_1sigma_cm2"])

    rows=[]; bandrows=[]
    for _,s in S.iterrows():
        model=s["model"]; z=float(s["Z_profile_norm_tilt"])
        meta=first_meta(X,model)
        kind=meta.get("kind","unknown")
        r2=meta.get("r2_cm2",np.nan)
        eps=meta.get("eps_em_inf",np.nan)

        scale_time=gen2_years/current_years
        z5=z*math.sqrt(scale_time*gen2_base)
        z10=z*math.sqrt(scale_time*gen2_opt)
        ztau=z*math.sqrt(scale_time*tau_tag)

        req2=((2.0/z)**2/scale_time) if z>0 else np.inf
        req3=((3.0/z)**2/scale_time) if z>0 else np.inf

        if kind=="charge_radius":
            allowed=bool(abs(float(r2))<=cr_bound*(1+1e-12))
            status="within_2025_global_1sigma_tau_NCR" if allowed else "outside_2025_global_1sigma_tau_NCR"
        elif kind=="derivative_em":
            allowed=False
            status="exploratory_forward_blind_no_global_bound_applied"
        else:
            allowed=False; status="unclassified"

        # Rebuild the nuisance-profiled SM prediction so q can be localized in energy.
        d=C[C.model==model].copy()
        gamma=float(s["bestfit_gamma"])
        alpha=float(s["bestfit_norm_tilt"])
        d["N_SM_profile"]=[alpha*sm_bin_at_gamma(r,det,gamma) for _,r in d.iterrows()]
        d["q_profile"]=[qbin(float(r.N_BSM),float(r.N_SM_profile)) for _,r in d.iterrows()]
        qtot=float(d.q_profile.sum())
        for label,lo,hi in BANDS:
            # Assign by geometric energy center.
            ec=np.sqrt(d.E_lo_GeV*d.E_hi_GeV)
            m=(ec>=lo)&(ec<hi)
            qb=float(d.loc[m,"q_profile"].sum())
            bandrows.append(dict(model=model,band=label,q=qb,
                                 q_fraction=qb/qtot if qtot>0 else 0.0,
                                 Z_equiv=math.sqrt(max(qb,0))))
        # High-statistic bins: where current IceCube expectation is appreciable.
        hi_stat=d[d.N_SM>=0.1]
        maxfrac=float(np.max(np.abs(hi_stat.frac))) if len(hi_stat) else np.nan
        medE=float(d.loc[d.q_profile.idxmax(),"E_lo_GeV"]*
                   d.loc[d.q_profile.idxmax(),"E_hi_GeV"])**0.5 if qtot>0 else np.nan

        rows.append(dict(
            model=model,kind=kind,r2_cm2=r2,eps_em_inf=eps,
            external_constraint_status=status,
            current_Z_profile_norm_tilt=z,
            current_total_fractional_shift=float(s["total_fractional_shift"]),
            Gen2_10yr_Z_area5=z5,
            Gen2_10yr_Z_area10=z10,
            Gen2_10yr_Z_tauTagScale5=ztau,
            required_area_multiplier_10yr_for_2sigma=req2,
            required_area_multiplier_10yr_for_3sigma=req3,
            max_abs_fraction_in_bins_NSM_ge_0p1=maxfrac,
            dominant_profiled_q_energy_GeV=medE
        ))

    R=pd.DataFrame(rows).sort_values("Gen2_10yr_Z_area10",ascending=False)
    B=pd.DataFrame(bandrows)
    R.to_csv(out/"future_sensitivity_summary.csv",index=False)
    B.to_csv(out/"future_sensitivity_energy_bands.csv",index=False)

    # Generic detector scaling table.
    grows=[]
    for _,r in R.iterrows():
        z=float(r.current_Z_profile_norm_tilt)
        for mult in fut["generic_future_tau"]["area_multipliers"]:
            years=float(fut["generic_future_tau"]["exposure_years"])
            zg=z*math.sqrt((years/current_years)*float(mult))
            grows.append(dict(model=r.model,area_multiplier=float(mult),
                              exposure_years=years,Z_profile_norm_tilt=zg))
    G=pd.DataFrame(grows)
    G.to_csv(out/"generic_tau_telescope_scaling.csv",index=False)

    # Decision restricted to charge-radius points with an external current bound.
    A=R[R.external_constraint_status=="within_2025_global_1sigma_tau_NCR"].copy()
    best=A.iloc[0].to_dict() if len(A) else None
    decision={
      "constraint_used":fut["constraints"],
      "gen2_scaling_basis":{
        "baseline_area_multiplier":gen2_base,
        "optimistic_area_multiplier":gen2_opt,
        "exposure_years":gen2_years,
        "tau_tag_rate_multiplier":tau_tag,
        "note":"Transparent scaling of the validated IceCube response; not a full Gen2 response matrix."
      },
      "best_externally_bounded_charge_radius_point":best,
      "any_bounded_point_reaches_2sigma_Gen2_baseline":
          bool(len(A) and np.any(A.Gen2_10yr_Z_area5>=2)),
      "any_bounded_point_reaches_2sigma_Gen2_optimistic":
          bool(len(A) and np.any(A.Gen2_10yr_Z_area10>=2)),
      "any_bounded_point_reaches_3sigma_Gen2_optimistic":
          bool(len(A) and np.any(A.Gen2_10yr_Z_area10>=3)),
      "derivative_points_are_exploratory":
          "No global laboratory bound is imposed on eps_em_inf; current-IceCube detectability alone is not called an allowed-region proof."
    }
    (out/"future_publishability_gate.json").write_text(json.dumps(decision,indent=2,default=str))

    cols=["model","external_constraint_status","current_Z_profile_norm_tilt",
          "Gen2_10yr_Z_area5","Gen2_10yr_Z_area10",
          "required_area_multiplier_10yr_for_2sigma",
          "dominant_profiled_q_energy_GeV"]
    print(R[cols].to_string(index=False))
    print(out/"future_publishability_gate.json")

if __name__=="__main__":
    main()
