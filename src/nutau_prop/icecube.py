from __future__ import annotations
from pathlib import Path
import io, math, urllib.request, zipfile
import numpy as np
import pandas as pd

ICECUBE_2015_URL = (
    "https://icecube.wisc.edu/data-releases/"
    "20150219_Search_for_contained_neutrino_events_at_energies_greater_than_1_TeV_in_2_years_of_data.zip"
)

def fetch_icecube_2015(dest="data/icecube_2015", url=ICECUBE_2015_URL):
    """
    Download and extract the official IceCube 2015 contained-event data release.
    Existing extracted files are reused.  No project files are deleted.
    """
    dest=Path(dest)
    target=dest/"effective_area.nu_tau.txt"
    if target.exists():
        return target
    dest.mkdir(parents=True,exist_ok=True)
    archive=dest/"icecube_2015_contained.zip"
    if not archive.exists():
        req=urllib.request.Request(url,headers={"User-Agent":"nutau-prop/0.4"})
        with urllib.request.urlopen(req,timeout=120) as r, open(archive,"wb") as f:
            shutil_copyfileobj(r,f)
    with zipfile.ZipFile(archive) as z:
        z.extractall(dest)
    hits=list(dest.rglob("effective_area.nu_tau.txt"))
    if not hits:
        raise FileNotFoundError(
            "Downloaded IceCube archive but effective_area.nu_tau.txt was not found."
        )
    # Keep original extraction hierarchy; return actual discovered file.
    return hits[0]

def shutil_copyfileobj(src,dst,length=1024*1024):
    while True:
        b=src.read(length)
        if not b: break
        dst.write(b)

def _numeric_rows(path):
    rows=[]
    with open(path,"r",errors="replace") as f:
        for line in f:
            q=line.strip()
            if not q or q.startswith(("#","%",";")):
                continue
            parts=q.replace(","," ").split()
            vals=[]
            try:
                vals=[float(x) for x in parts]
            except ValueError:
                continue
            if len(vals)>=5:
                rows.append(vals)
    if not rows:
        raise ValueError(f"No >=5-column numeric rows found in {path}")
    n=min(len(r) for r in rows)
    return np.asarray([r[:n] for r in rows],float)

def parse_total_effective_area(path):
    """
    Parse IceCube effective_area.nu_tau.txt.

    The official release describes bins in true neutrino energy and zenith.
    This parser deliberately supports the common representations:
      [E_lo,E_hi,cosz_lo,cosz_hi,Aeff_m2],
      [log10E_lo,log10E_hi,cosz_lo,cosz_hi,Aeff_m2],
      or zenith bounds in radians/degrees.

    It returns canonical GeV energy bounds and cos(zenith) bounds.
    IceCube zenith convention is used: cosz<0 is upgoing through Earth.
    """
    a=_numeric_rows(path)
    if a.shape[1] < 5:
        raise ValueError("Expected at least 5 numeric columns")
    e0,e1,z0,z1,aeff=[a[:,i].copy() for i in range(5)]

    # Energy format: log10(GeV) values are O(3--9), while GeV values are >=1e3.
    if np.nanmax(np.r_[e0,e1]) < 50.0:
        elo=10.0**e0; ehi=10.0**e1
        energy_format="log10_GeV"
    else:
        elo=e0; ehi=e1
        energy_format="GeV"

    zmax=float(np.nanmax(np.abs(np.r_[z0,z1])))
    zmin=float(np.nanmin(np.r_[z0,z1]))
    if zmax <= 1.000001 and zmin >= -1.000001:
        c0,c1=z0,z1; angle_format="coszenith"
    elif zmax <= math.pi+1e-6 and zmin >= -1e-9:
        c0,c1=np.cos(z0),np.cos(z1); angle_format="radian_zenith"
    elif zmax <= 180.0001 and zmin >= -1e-9:
        c0,c1=np.cos(np.deg2rad(z0)),np.cos(np.deg2rad(z1)); angle_format="degree_zenith"
    else:
        raise ValueError(f"Unrecognized zenith representation: range [{zmin},{zmax}]")

    clo=np.minimum(c0,c1); chi=np.maximum(c0,c1)
    out=pd.DataFrame(dict(
        E_lo_GeV=elo,E_hi_GeV=ehi,
        cosz_lo=clo,cosz_hi=chi,
        Aeff_m2=aeff
    ))
    out.attrs["energy_format"]=energy_format
    out.attrs["angle_format"]=angle_format
    return out

def powerlaw_integral(phi0,pivot,gamma,elo,ehi):
    """Integral dE phi0*(E/pivot)^(-gamma), units cm^-2 s^-1 sr^-1."""
    elo=float(elo); ehi=float(ehi)
    if ehi<=elo: return 0.0
    pref=float(phi0)*float(pivot)**float(gamma)
    if abs(gamma-1.0)<1e-12:
        return pref*math.log(ehi/elo)
    return pref*(ehi**(1-gamma)-elo**(1-gamma))/(1-gamma)

def solid_angle_from_cosz(clo,chi):
    return 2.0*math.pi*abs(float(chi)-float(clo))

def poisson_asimov_q(mu_true,mu_test):
    """2*Poisson KL divergence for Asimov counts."""
    a=np.asarray(mu_true,float); b=np.asarray(mu_test,float)
    m=(a>0)&(b>0)
    q=np.sum(2*(b[m]-a[m]+a[m]*np.log(a[m]/b[m])))
    q += np.sum(2*b[(a<=0)&(b>0)])
    return float(max(q,0.0))
