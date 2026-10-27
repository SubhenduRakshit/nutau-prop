#!/usr/bin/env python3
"""v0.8 UHE Case-6 detector fold.

Downloads official IceCube EHE and KM3NeT/KM3-230213A response inputs, validates
KM3NeT against the collaboration's Ehlert expectation, then folds Case-6 r(E).
Primary EFT prediction: sqrt(shat)<Lambda.  Q<Lambda is an aggressive envelope.

Important topology caveat: Case 6 is NC-like; KM3NeT's published UHE response here
is a bright-track response.  Therefore NC-equivalent KM3 values are acceptance
proxies, not claimed physical bright-track enhancements.
"""
import json, math
from pathlib import Path
from urllib.request import Request, urlopen
from urllib.parse import urlencode, quote
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

D=Path("data/uhe_official"); D.mkdir(parents=True,exist_ok=True)
O=Path("results/v0.8"); O.mkdir(parents=True,exist_ok=True)

def GET(url):
    req=Request(url,headers={"User-Agent":"nutau-prop-v0.8"})
    with urlopen(req,timeout=90) as r: return r.read()

def get_icecube():
    meta_url="https://dataverse.harvard.edu/api/datasets/:persistentId/?"+urlencode({"persistentId":"doi:10.7910/DVN/JHK49D"})
    meta=json.loads(GET(meta_url)); files=meta["data"]["latestVersion"]["files"]
    cand=[]
    for f in files:
        df=f["dataFile"]
        name=" ".join(str(df.get(k,"")) for k in ("filename","description")) + " " + str(f.get("label",""))
        if "effective" in name.lower() and "area" in name.lower(): cand.append(f)
    if not cand: raise RuntimeError("IceCube EHE effective-area file not found in Dataverse metadata")
    # Prefer CSV and the smallest exact effective_area-named product.
    cand.sort(key=lambda f: ("csv" not in f["dataFile"].get("filename","").lower(), f["dataFile"].get("filesize",10**18)))
    ff=cand[0]; p=D/"icecube_ehe_effective_area.csv"
    p.write_bytes(GET(f"https://dataverse.harvard.edu/api/access/datafile/{ff['dataFile']['id']}?format=original"))
    return p

def numeric_response(path):
    df=pd.read_csv(path,sep=None,engine="python",comment="#")
    num=df.apply(pd.to_numeric,errors="coerce").dropna(how="all")
    usable=[c for c in num.columns if num[c].notna().sum()>=3]
    if len(usable)<2: raise RuntimeError(f"Could not identify two numeric columns in {path}")
    # Prefer names containing energy and effective area.
    ec=next((c for c in usable if "ener" in str(c).lower()),usable[0])
    ac=next((c for c in usable if ("area" in str(c).lower() or "aeff" in str(c).lower()) and c!=ec),None)
    if ac is None: ac=next(c for c in usable if c!=ec)
    q=num[[ec,ac]].dropna(); E=q[ec].to_numpy(float); A=q[ac].to_numpy(float).copy()*1e4  # IceCube official area is m^2
    en=str(ec).lower(); an=str(ac).lower()
    # Energy normalisation.
    if "ev" in en and "gev" not in en: E=E/1e9
    elif "tev" in en: E=E*1e3
    # Area normalisation.  Official EHE tables commonly use m^2; retain cm^2 if explicit.
    if False: A=A*1e4
    return E,A,ec,ac

def download_km3():
    base="https://raw.githubusercontent.com/KM3NeT/KM3-230213A-data/main/"
    paths={"km3":"data/supplementary/simulations/effective_area_brighttrackselection_allflavour_skyavg.json",
           "flux":"data/external/flux_models/cosmogenic_Ehlert+[2304.07321].json"}
    for key,path in paths.items(): (D/f"{key}.json").write_bytes(GET(base+quote(path,safe="/+[]")))

def rcurve(cs,lam,col,E):
    q=cs[np.isclose(cs.Lambda_TeV,lam)].sort_values("E_PeV")
    x=q.E_PeV.to_numpy(float)*1e6; y=q[col].to_numpy(float)
    good=np.isfinite(x)&np.isfinite(y)&(x>0)&(y>0)
    x=x[good]; y=y[good]; out=np.zeros_like(E)
    m=(E>=x.min())&(E<=x.max())
    out[m]=10**np.interp(np.log10(E[m]),np.log10(x),np.log10(y))
    return out

def main():
    icfile=get_icecube(); download_km3()
    Eic,Aic,ec,ac=numeric_response(icfile)
    km=pd.read_json(D/"km3.json"); Ekm=km["Energy [GeV]"].to_numpy(float); Akm=km["Aeff [cm^2]"].to_numpy(float)
    fl=pd.read_json(D/"flux.json"); Ef=fl["Energy [GeV]"].to_numpy(float); Sf=fl["log10(E^2 F(E) [GeV.cm^-2.s^-1.sr^-1])"].to_numpy(float)
    sed=interp1d(Ef,Sf,bounds_error=False,fill_value=-np.inf)
    E=np.logspace(5,11,6001); phi=10.0**sed(E)/E**2
    aeff_ic=interp1d(Eic,Aic,bounds_error=False,fill_value=0.0)(E)
    aeff_km=interp1d(Ekm,Akm,bounds_error=False,fill_value=0.0)(E)
    TIC=12.6*365.25*86400.; TKM=335.*86400.
    def fold(aeff,T,w=None,elo=None,ehi=None):
        if w is None: w=np.ones_like(E)
        mask=np.ones(E.size,dtype=bool)
        if elo is not None: mask&=E>=elo
        if ehi is not None: mask&=E<=ehi
        return 4*np.pi*T*np.trapezoid((aeff*phi*w)[mask],E[mask])
    NIC0=fold(aeff_ic,TIC); NKM0=fold(aeff_km,TKM); NKM90=fold(aeff_km,TKM,elo=7.24e7,ehi=2.57e9)
    print(f"ICECUBE columns: {ec!r}, {ac!r}")
    print(f"VALIDATION KM3 Ehlert SM 335d = {NKM0:.6f}  target=0.242767")
    print(f"VALIDATION KM3 central90      = {NKM90:.6f}  target=0.211117")
    print(f"IceCube EHE Ehlert SM 12.6y  = {NIC0:.6f}")
    if abs(NKM0/0.242767-1)>0.02: raise RuntimeError("KM3 validation failed by >2%")
    csfile=O/"case6_pdf_cut_v0_8.csv"
    if not csfile.exists(): raise RuntimeError(f"Missing {csfile}; run scripts/run_case6_pdf_cut_v0_8.py first")
    cs=pd.read_csv(csfile); rows=[]
    for lam in (1.5,2.0,2.44,3.0):
        q=cs[np.isclose(cs.Lambda_TeV,lam)].iloc[0]; rs=rcurve(cs,lam,"r_shatcut",E); rq=rcurve(cs,lam,"r_Qcut",E)
        NICs=fold(aeff_ic,TIC,1+rs); NICq=fold(aeff_ic,TIC,1+rq)
        NKMs=fold(aeff_km,TKM,1+rs); NKMq=fold(aeff_km,TKM,1+rq); NKMs90=fold(aeff_km,TKM,1+rs,7.24e7,2.57e9)
        scale90=2.302585093/NICs if NICs>0 else math.inf
        rows.append(dict(Lambda_TeV=lam,c4_over_4pi=float(q.c4_over_4pi),IceCube_SM=NIC0,
            IceCube_shat_proxy=NICs,IceCube_Q_proxy=NICq,IceCube_delta_shat=NICs-NIC0,
            IceCube_frac_shat=NICs/NIC0-1,IceCube_P0_shat=math.exp(-NICs),naive_flux_scale_90_if_zero=scale90,
            KM3_brighttrack_SM=NKM0,KM3_brighttrack_direct_NP=NKM0,
            KM3_NC_equiv_shat_proxy=NKMs,KM3_NC_equiv_Q_proxy=NKMq,
            KM3_NC_equiv_shat_central90=NKMs90,KM3_Pge1_NC_equiv_shat=1-math.exp(-NKMs),
            KM3_NC_equiv_at_naive_IC90=NKMs*scale90))
    R=pd.DataFrame(rows); R.to_csv(O/"uhe_event_counts_v0_8.csv",index=False)
    print("\nPRIMARY: conservative sqrt(shat)<Lambda")
    print(R[["Lambda_TeV","c4_over_4pi","IceCube_SM","IceCube_shat_proxy","IceCube_frac_shat","IceCube_P0_shat","KM3_brighttrack_SM","KM3_NC_equiv_shat_proxy","KM3_NC_equiv_shat_central90","naive_flux_scale_90_if_zero"]].to_string(index=False,float_format=lambda x:f"{x:.4g}"))
    print("\nAGGRESSIVE Q<Lambda envelope")
    print(R[["Lambda_TeV","IceCube_Q_proxy","KM3_NC_equiv_Q_proxy"]].to_string(index=False,float_format=lambda x:f"{x:.4g}"))
    print("\nIMPORTANT: Case 6 is NC-like; KM3 NC-equivalent values are acceptance proxies, not new bright-track counts.")
    print("WROTE",O/"uhe_event_counts_v0_8.csv")
if __name__=="__main__": main()
