#!/usr/bin/env python3
import argparse,json,math
from pathlib import Path
import numpy as np,pandas as pd
from scipy.optimize import minimize_scalar

from nutau_prop.config import load_config
from nutau_prop.arca import (
    fetch_zenodo_record,response_inventory,fig7_response,calibrate_fig7_response,
    powerlaw,integrate_loggrid
)
from nutau_prop.icecube import poisson_asimov_q

def interp_scan(scan,model,E,mu,col):
    d=scan[scan.model==model]
    mus=np.sort(d.cosz_abs.unique())
    vals=[]
    for m in mus:
        q=d[d.cosz_abs==m].sort_values("E_GeV")
        vals.append(np.interp(np.log10(E),np.log10(q.E_GeV),q[col],
                              left=q[col].iloc[0],right=q[col].iloc[-1]))
    # horizon -> no Earth propagation for prop ratio
    if col=="prop_ratio":
        mus=np.r_[0.0,mus]; vals=np.r_[1.0,vals]
    return float(np.interp(mu,mus,vals,left=vals[0],right=vals[-1]))

def avg_upgoing_prop(scan,model,E,col="prop_ratio"):
    # half-sky average in cos(theta), using 8-point midpoint quadrature
    mus=(np.arange(8)+0.5)/8
    return float(np.mean([interp_scan(scan,model,E,m,col) for m in mus]))

def detector_nc_ratio(scan,model,E):
    # high-Q detector interaction does not depend on Earth chord; use a stable
    # available scan slice.
    return interp_scan(scan,model,E,0.5,"detector_nc_ratio")

def bin_edges():
    # true-energy bins dense enough to exploit spectral + up/down information
    return np.logspace(np.log10(2.5e4),np.log10(1.4e6),17)

def model_counts(scan,model,cfg,gamma,bsmmode=True):
    expo=float(cfg["arca"]["exposure_years"])*365.25*86400
    phi0=float(cfg["analysis_flux"]["phi0_per_gev_cm2_s_sr"])
    piv=float(cfg["analysis_flux"]["pivot_gev"])
    sf_sh,sf_tr=calibrate_fig7_response(cfg)
    edges=bin_edges()
    rows=[]

    # 18-point quadrature inside each energy bin.
    for lo,hi in zip(edges[:-1],edges[1:]):
        E=np.geomspace(lo,hi,18)
        phi=powerlaw(phi0,piv,gamma,E)
        R=fig7_response(E)
        sh=R["shower"]

        # Downgoing showers: no Earth propagation.
        if bsmmode:
            rnc=np.array([detector_nc_ratio(scan,model,e) for e in E])
        else:
            rnc=np.ones_like(E)
        # Only one of the three equal flavor NC components is tau and BSM-sensitive.
        nc_down_factor=(2.0+rnc)/3.0
        A_down=(sh["down_ecc"]+sh["down_tcc"]+sh["down_nc"]*nc_down_factor)*sf_sh
        Nd=expo*2*math.pi*integrate_loggrid(E,phi*A_down)*1e4
        rows.append(dict(channel="shower_down",E_lo=lo,E_hi=hi,N=Nd))

        # Upgoing showers: tau CC carries propagation; only tau third of NC
        # carries propagation + BSM detector NC shift.
        if bsmmode:
            rp=np.array([avg_upgoing_prop(scan,model,e) for e in E])
            rnc=np.array([detector_nc_ratio(scan,model,e) for e in E])
        else:
            rp=np.ones_like(E); rnc=np.ones_like(E)
        A_up=(sh["up_ecc"]
              + sh["up_tcc"]*rp
              + sh["up_nc"]*((2.0+rp*rnc)/3.0))*sf_sh
        Nu=expo*2*math.pi*integrate_loggrid(E,phi*A_up)*1e4
        rows.append(dict(channel="shower_up",E_lo=lo,E_hi=hi,N=Nu))

        # nu_mu CC track control bins. Tau-only BSM leaves these unchanged.
        for key,(clo,chi) in R["track_cos_bins"].items():
            A=R["track"][key]*sf_tr
            domega=2*math.pi*(chi-clo)
            Nt=expo*domega*integrate_loggrid(E,phi*A)*1e4
            rows.append(dict(channel="track_"+key,E_lo=lo,E_hi=hi,N=Nt))
    return pd.DataFrame(rows)

def profiled_q(scan,model,cfg,true_bsm,include_track=True):
    gamma0=float(cfg["analysis_flux"]["gamma"])
    true=model_counts(scan,model,cfg,gamma0,bsmmode=true_bsm)
    if not include_track:
        true=true[true.channel.str.startswith("shower")].reset_index(drop=True)
    t=true.N.to_numpy(float)

    gmin=float(cfg["profile"]["gamma_min"]); gmax=float(cfg["profile"]["gamma_max"])
    def objective(g):
        sm=model_counts(scan,model,cfg,g,bsmmode=False)
        if not include_track:
            sm=sm[sm.channel.str.startswith("shower")].reset_index(drop=True)
        s=sm.N.to_numpy(float)
        alpha=t.sum()/s.sum()
        return poisson_asimov_q(t,alpha*s)

    opt=minimize_scalar(objective,bounds=(gmin,gmax),method="bounded",
                        options={"xatol":2e-4})
    sm=model_counts(scan,model,cfg,opt.x,bsmmode=False)
    if not include_track:
        sm=sm[sm.channel.str.startswith("shower")].reset_index(drop=True)
    s=sm.N.to_numpy(float); alpha=t.sum()/s.sum()
    q=poisson_asimov_q(t,alpha*s)
    return math.sqrt(max(q,0)),float(alpha),float(opt.x),true,sm

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--config",default="config/arca.yaml")
    ap.add_argument("--scan",default="results/bsm_scan_future.csv")
    ap.add_argument("--outdir",default="results/arca")
    ap.add_argument("--data-dir",default="data/arca_zenodo")
    ap.add_argument("--skip-zenodo",action="store_true")
    a=ap.parse_args()
    cfg=load_config(a.config); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)

    # Acquire official record for provenance/inventory.  Current v0.6 still uses
    # the calibrated Fig.7 response unless a later parser explicitly promotes
    # an exact Zenodo table.  This avoids silently misreading unknown schemas.
    zenodo_ok=False; zerr=None
    if not a.skip_zenodo:
        try:
            fetch_zenodo_record(a.data_dir,int(cfg["arca"]["zenodo_record"]))
            zenodo_ok=True
        except Exception as e:
            zerr=repr(e)
    inv=response_inventory(a.data_dir)
    (out/"zenodo_inventory.json").write_text(json.dumps(inv,indent=2))

    scan=pd.read_csv(a.scan)
    models=sorted(scan.model.unique())
    rows=[]; details=[]
    for model in models:
        zboth,ab,gb,true,sm=profiled_q(scan,model,cfg,True,include_track=True)
        zshow,as_,gs,_,_=profiled_q(scan,model,cfg,True,include_track=False)

        # fixed-flux/no-nuisance upper-level diagnostic
        sm0=model_counts(scan,model,cfg,float(cfg["analysis_flux"]["gamma"]),False)
        t=true.N.to_numpy(float); s=sm0.N.to_numpy(float)
        zfix=math.sqrt(max(poisson_asimov_q(t,s),0))

        sh=true[true.channel.str.startswith("shower")]
        shsm=sm0[sm0.channel.str.startswith("shower")]
        frac=(sh.N.sum()/shsm.N.sum()-1) if shsm.N.sum()>0 else np.nan

        # up/down ratio is the detector-specific Earth-shape observable
        up=true[true.channel=="shower_up"].N.sum()
        dn=true[true.channel=="shower_down"].N.sum()
        upsm=sm0[sm0.channel=="shower_up"].N.sum()
        dnsm=sm0[sm0.channel=="shower_down"].N.sum()
        rd=(up/dn)/(upsm/dnsm)-1 if dn>0 and dnsm>0 and upsm>0 else np.nan

        rows.append(dict(
            model=model,response_source="Fig7_digitisation_calibrated_to_Table1_2",
            zenodo_download_ok=zenodo_ok,
            Z_fixed_flux_all_channels=zfix,
            Z_profile_norm_tilt_showers_only=zshow,
            Z_profile_norm_tilt_tracks_plus_showers=zboth,
            bestfit_norm_tracks_plus_showers=ab,
            bestfit_gamma_tracks_plus_showers=gb,
            shower_total_fractional_shift=frac,
            up_down_ratio_fractional_shift=rd,
            N_BSM_showers=float(sh.N.sum()),
            N_SM_showers=float(shsm.N.sum()),
            N_BSM_tracks=float(true[true.channel.str.startswith("track")].N.sum()),
            N_SM_tracks=float(sm0[sm0.channel.str.startswith("track")].N.sum())
        ))
        q=true.copy()
        q["model"]=model
        q["N_SM_same_flux"]=sm0.N.to_numpy(float)
        q["frac_same_flux"]=q.N/q.N_SM_same_flux-1
        details.append(q)

    R=pd.DataFrame(rows).sort_values("Z_profile_norm_tilt_tracks_plus_showers",ascending=False)
    D=pd.concat(details,ignore_index=True)
    R.to_csv(out/"arca_sensitivity_summary.csv",index=False)
    D.to_csv(out/"arca_event_bins.csv",index=False)

    decision={
      "response_source":"coarse published Fig.7 digitisation, separately calibrated to exact published 220/year track and 60/year shower cosmic rates",
      "zenodo_record_downloaded":zenodo_ok,
      "zenodo_error":zerr,
      "important_limit":"This is a signal-only optimistic upper-bound likelihood; atmospheric backgrounds are not yet included bin-by-bin.",
      "interpretation":"If Z<2 even here, ARCA cannot rescue the benchmark. If Z>=2, exact Zenodo response/background parsing is mandatory before a positive claim.",
      "best_model":R.iloc[0].to_dict() if len(R) else None,
      "any_model_ge_2sigma_upper_bound":bool(len(R) and np.any(R.Z_profile_norm_tilt_tracks_plus_showers>=2)),
      "any_charge_radius_ge_2sigma_upper_bound":bool(len(R) and np.any(
          R[R.model.str.startswith("cr_")].Z_profile_norm_tilt_tracks_plus_showers>=2))
    }
    (out/"arca_publishability_gate.json").write_text(json.dumps(decision,indent=2,default=str))

    print(R[[
      "model","Z_profile_norm_tilt_showers_only",
      "Z_profile_norm_tilt_tracks_plus_showers",
      "shower_total_fractional_shift","up_down_ratio_fractional_shift",
      "N_BSM_showers","N_BSM_tracks"
    ]].to_string(index=False))
    print(out/"arca_publishability_gate.json")

if __name__=="__main__":
    main()
