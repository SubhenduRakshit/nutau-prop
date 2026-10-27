#!/usr/bin/env python3
"""v0.11 detector-NSI stage for Case 6, before modified Earth transport.

The new interaction is added only to tau-neutrino NC interactions.  ARCA's
published shower decomposition is therefore reweighted channel by channel:
  e CC unchanged; tau CC unchanged; the tau third of all-flavor NC is changed.
Muon-track control channels are unchanged because Case 6 is NC-like.

This is the most detector-consistent public-response calculation possible with
currently available channel tables.  The residual limitation is that published
Aeff tables have already integrated over the SM inelasticity distribution, so
changed-y selection efficiency cannot be reconstructed exactly without MC.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

from nutau_prop.arca import fig7_response, calibrate_fig7_response, integrate_loggrid, powerlaw
from nutau_prop.config import load_config

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "v0.11" / "tables"
OUT = ROOT / "results" / "v0.11"
UHE = ROOT / "data" / "uhe_official"
YEAR = 365.25 * 86400.0


def tag(lam):
    return str(lam).replace(".", "p")


def load_case6(lam, species="nu"):
    z = np.load(TABLES / f"case6_L{tag(lam)}_{species}.npz")
    return {k: z[k] for k in z.files}


def logcurve(E, grid, values, left=None, right=None):
    E=np.asarray(E,float); grid=np.asarray(grid,float); values=np.asarray(values,float)
    lo = values[0] if left is None else left; hi = values[-1] if right is None else right
    return 10.0**np.interp(np.log10(E),np.log10(grid),np.log10(np.maximum(values,1e-300)),
                           left=np.log10(max(lo,1e-300)),right=np.log10(max(hi,1e-300)))


def ratios(table, E):
    scc=logcurve(E,table["energies"],table["sigma_cc"])
    snc0=logcurve(E,table["energies"],table["sigma_nc_sm"])
    snp=logcurve(E,table["energies"],table["sigma_np"],left=0.0,right=0.0)
    r_nc=(snc0+snp)/np.maximum(snc0,1e-300)
    r_incl=(scc+snc0+snp)/np.maximum(scc+snc0,1e-300)
    return r_nc,r_incl


def ehlert_flux(E):
    q=pd.read_json(UHE/"flux.json")
    ef=q["Energy [GeV]"].to_numpy(float)
    sed=q["log10(E^2 F(E) [GeV.cm^-2.s^-1.sr^-1])"].to_numpy(float)
    s=interp1d(ef,sed,bounds_error=False,fill_value=-np.inf)(E)
    return 10.0**s/E**2


def arca_counts(lam,cfg,flux_name):
    sf_sh,sf_tr=calibrate_fig7_response(cfg)
    expo=float(cfg["arca"]["exposure_years"])*YEAR
    if flux_name=="legacy_powerlaw":
        E=np.logspace(np.log10(float(cfg["analysis_flux"]["energy_min_gev"])),
                      np.log10(float(cfg["analysis_flux"]["energy_max_gev"])),1201)
        phi=powerlaw(float(cfg["analysis_flux"]["phi0_per_gev_cm2_s_sr"]),
                     float(cfg["analysis_flux"]["pivot_gev"]),
                     float(cfg["analysis_flux"]["gamma"]),E)
    elif flux_name=="Ehlert_UHE_to_100PeV":
        E=np.logspace(5,8,1801)
        phi=ehlert_flux(E)
    else:
        raise ValueError(flux_name)
    R=fig7_response(E); sh=R["shower"]
    nu=load_case6(lam,"nu"); anu=load_case6(lam,"antinu")
    rnu,_=ratios(nu,E); ranu,_=ratios(anu,E); rnc=0.5*(rnu+ranu)
    # Equal incident flavors. Only the tau one-third of the published all-flavor
    # NC response is modified at the detector in this stage.
    nc_factor=(2.0+rnc)/3.0
    asm_down=(sh["down_ecc"]+sh["down_tcc"]+sh["down_nc"])*sf_sh
    absm_down=(sh["down_ecc"]+sh["down_tcc"]+sh["down_nc"]*nc_factor)*sf_sh
    asm_up=(sh["up_ecc"]+sh["up_tcc"]+sh["up_nc"])*sf_sh
    absm_up=(sh["up_ecc"]+sh["up_tcc"]+sh["up_nc"]*nc_factor)*sf_sh
    factor=expo*2.0*math.pi*1e4
    nsm_down=factor*integrate_loggrid(E,phi*asm_down)
    nb_down=factor*integrate_loggrid(E,phi*absm_down)
    nsm_up=factor*integrate_loggrid(E,phi*asm_up)
    nb_up=factor*integrate_loggrid(E,phi*absm_up)
    # Track control remains exactly unchanged at detector level.
    ntrack=0.0
    for key,(clo,chi) in R["track_cos_bins"].items():
        domega=2.0*math.pi*(chi-clo)
        ntrack += expo*domega*integrate_loggrid(E,phi*R["track"][key]*sf_tr)*1e4
    return dict(flux=flux_name,Lambda_TeV=lam,
                E_min_GeV=float(E.min()),E_max_GeV=float(E.max()),
                N_shower_SM=nsm_down+nsm_up,N_shower_detectorNSI=nb_down+nb_up,
                shower_fractional_shift=(nb_down+nb_up)/(nsm_down+nsm_up)-1.0,
                up_down_SM=nsm_up/nsm_down,up_down_detectorNSI=nb_up/nb_down,
                N_track_SM=ntrack,N_track_detectorNSI=ntrack,
                direct_track_fractional_shift=0.0)


def icecube_tau_response(lam):
    # This is retained as an explicitly labelled inclusive tau-response proxy;
    # IceCube's public EHE table does not separate tau CC and tau NC acceptance.
    values=np.loadtxt(UHE/"icecube_ehe_effective_area.csv",delimiter=",",skiprows=1)
    Eic=values[:,0]; A=values[:,1:]*1e4
    E=np.logspace(5,11,6001); phi=ehlert_flux(E)
    atot=interp1d(Eic,A[:,0],bounds_error=False,fill_value=0.0)(E)
    atau=interp1d(Eic,A[:,3],bounds_error=False,fill_value=0.0)(E)
    nu=load_case6(lam,"nu"); anu=load_case6(lam,"antinu")
    _,rinu=ratios(nu,E); _,rianu=ratios(anu,E); rincl=0.5*(rinu+rianu)
    # Only tau-initial area is eligible. This still cannot resolve CC/NC-specific
    # selection or altered Earth propagation, hence 'proxy' in every field name.
    absm=atot+atau*(rincl-1.0)
    T=12.6*YEAR
    fold=lambda area:4.0*math.pi*T*np.trapezoid(area*phi,E)
    return dict(Lambda_TeV=lam,IceCube_EHE_SM=fold(atot),
                IceCube_tau_detectorNSI_proxy=fold(absm),
                IceCube_tau_detectorNSI_proxy_fractional_shift=fold(absm)/fold(atot)-1.0)


def y_diagnostics(lam):
    q=load_case6(lam,"nu")
    out={"Lambda_TeV":lam}
    for epev in (1.0,10.0,100.0,220.0,300.0):
        my=float(np.interp(np.log10(epev*1e6),np.log10(q["energies"]),q["mean_y_np"]))
        out[f"mean_y_NP_{epev:g}PeV"]=my
    return out


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--config",default="config/arca.yaml")
    ap.add_argument("--lambdas",nargs="+",type=float,default=[1.5,2.0,2.44,3.0])
    a=ap.parse_args(); cfg=load_config(a.config); OUT.mkdir(parents=True,exist_ok=True)
    required=[TABLES/f"case6_L{tag(l)}_nu.npz" for l in a.lambdas]
    if not all(p.exists() for p in required):
        raise SystemExit("Missing v0.11 differential tables; run build_case6_tables_v0_11.py first")
    rows=[]
    for lam in a.lambdas:
        for flux in ("legacy_powerlaw","Ehlert_UHE_to_100PeV"):
            r=arca_counts(lam,cfg,flux); r.update(y_diagnostics(lam)); rows.append(r)
    D=pd.DataFrame(rows); D.to_csv(OUT/"case6_detector_channel_fold_v0_11.csv",index=False)
    IC=pd.DataFrame([icecube_tau_response(l) for l in a.lambdas])
    IC.to_csv(OUT/"case6_icecube_tau_detector_proxy_v0_11.csv",index=False)
    report={
      "version":"v0.11",
      "detector_NSI_status":"channel-resolved for public ARCA shower components; direct muon-track increment forced to zero",
      "remaining_public_response_limit":"ARCA/IceCube Aeff tables integrate over SM y; exact changed-y selection requires detector MC",
      "important":"No modified Earth propagation is included in v0.11. That is deliberately the next v0.12 stage.",
    }
    (OUT/"detector_NSI_status_v0_11.json").write_text(json.dumps(report,indent=2))
    print("\nARCA channel-resolved detector-only Case-6 fold")
    print(D[["flux","Lambda_TeV","N_shower_SM","N_shower_detectorNSI","shower_fractional_shift","N_track_SM","N_track_detectorNSI"]].to_string(index=False,float_format=lambda x:f"{x:.5g}"))
    print("\nIceCube EHE tau-initial detector-response proxy (no modified Earth transport)")
    print(IC.to_string(index=False,float_format=lambda x:f"{x:.5g}"))
    print("\nResidual response systematic: published Aeff does not expose y-dependent selection. Case-6 direct track shift = 0 by topology.")
    print("WROTE",OUT)


if __name__=="__main__":
    main()
