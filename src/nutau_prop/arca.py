from __future__ import annotations
from pathlib import Path
import json, math, urllib.request, zipfile, tarfile, shutil
import numpy as np
import pandas as pd

ZENODO_RECORD = 10949703
ZENODO_API = f"https://zenodo.org/api/records/{ZENODO_RECORD}"

def fetch_zenodo_record(dest="data/arca_zenodo", record=ZENODO_RECORD):
    """
    Download every public file attached to the official KM3NeT/ARCA Zenodo
    record and recursively extract zip/tar archives.  Returns an inventory.
    If Zenodo is temporarily unavailable the caller may use the documented
    Fig.7 fallback.
    """
    dest=Path(dest); dest.mkdir(parents=True,exist_ok=True)
    meta_path=dest/"record.json"
    url=f"https://zenodo.org/api/records/{int(record)}"
    req=urllib.request.Request(url,headers={"User-Agent":"nutau-prop/0.6"})
    with urllib.request.urlopen(req,timeout=120) as r:
        meta=json.loads(r.read().decode("utf-8"))
    meta_path.write_text(json.dumps(meta,indent=2))
    inv=[]
    for f in meta.get("files",[]):
        key=f.get("key") or f.get("filename") or "unknown"
        link=(f.get("links") or {}).get("self") or (f.get("links") or {}).get("content")
        if not link:
            continue
        out=dest/key
        out.parent.mkdir(parents=True,exist_ok=True)
        if not out.exists():
            req=urllib.request.Request(link,headers={"User-Agent":"nutau-prop/0.6"})
            with urllib.request.urlopen(req,timeout=180) as r, open(out,"wb") as w:
                shutil.copyfileobj(r,w)
        inv.append(dict(name=key,size=out.stat().st_size,path=str(out)))
        _extract_if_archive(out,dest/"extracted")
    # inventory extracted tree too
    ex=dest/"extracted"
    if ex.exists():
        for q in ex.rglob("*"):
            if q.is_file():
                inv.append(dict(name=q.name,size=q.stat().st_size,path=str(q)))
    (dest/"inventory.json").write_text(json.dumps(inv,indent=2))
    return inv

def _extract_if_archive(path,dest):
    try:
        if zipfile.is_zipfile(path):
            d=dest/path.stem
            d.mkdir(parents=True,exist_ok=True)
            with zipfile.ZipFile(path) as z: z.extractall(d)
        elif tarfile.is_tarfile(path):
            d=dest/path.stem
            d.mkdir(parents=True,exist_ok=True)
            with tarfile.open(path) as t: t.extractall(d)
    except Exception:
        pass

def response_inventory(dest="data/arca_zenodo"):
    dest=Path(dest)
    inv=dest/"inventory.json"
    if inv.exists():
        return json.loads(inv.read_text())
    rows=[]
    for q in dest.rglob("*"):
        if q.is_file():
            rows.append(dict(name=q.name,size=q.stat().st_size,path=str(q)))
    return rows

# ---------- Published Fig. 7 fallback ----------
# Coarse log-log digitisation of Fig.7(b,c), calibrated below to the exact
# Table-1/2 cosmic rates.  It is not advertised as the Zenodo response.
_E=np.array([1e3,3e3,1e4,3e4,1e5,3e5,1e6,3e6,1e7,3e7,1e8],float)

# Shower selected: eCC, tauCC, all-flavour NC; down/up hemispheres.
# Units m^2 before rate calibration.
_SHOWER={
 "down_ecc": np.array([.12,.35,1.0,3.2,9,22,50,95,160,210,250]),
 "down_tcc": np.array([.05,.16,.55,1.8,5.5,14,35,70,125,205,320]),
 "down_nc":  np.array([.035,.12,.42,1.4,4.5,12,30,62,115,195,290]),
 "up_ecc":   np.array([.18,.55,1.3,3.8,8.5,14,20,24,25,23,20]),
 "up_tcc":   np.array([.06,.20,.58,1.7,4.8,9,15,20,25,30,34]),
 "up_nc":    np.array([.035,.12,.38,1.2,3.8,7.5,12,16,20,23,25]),
}

# Track selected nu_mu CC in four published cos(theta) bins.
_TRACK={
 "c008_m025": np.array([.25,1.2,7,28,100,270,650,1250,2200,3300,4500]),
 "cm025_m050":np.array([.25,1.1,6.5,24,82,215,480,760,680,400,210]),
 "cm050_m075":np.array([.25,1.1,6.3,22,72,165,205,145,58,17,2.5]),
 "cm075_m100":np.array([.25,1.05,6.0,20,62,115,82,28,3.5,.25,.008]),
}
_TRACK_COS_BINS={
 "c008_m025":(-0.25,0.08),
 "cm025_m050":(-0.50,-0.25),
 "cm050_m075":(-0.75,-0.50),
 "cm075_m100":(-1.00,-0.75),
}

def _loginterp(vals,E):
    E=np.asarray(E,float)
    return 10**np.interp(np.log10(E),np.log10(_E),
                          np.log10(np.maximum(vals,1e-12)))

def fig7_response(E):
    E=np.asarray(E,float)
    return {
      "shower":{k:_loginterp(v,E) for k,v in _SHOWER.items()},
      "track":{k:_loginterp(v,E) for k,v in _TRACK.items()},
      "track_cos_bins":dict(_TRACK_COS_BINS),
      "source":"published_Fig7_coarse_digitisation"
    }

def powerlaw(phi0,pivot,gamma,E):
    return float(phi0)*(np.asarray(E,float)/float(pivot))**(-float(gamma))

def integrate_loggrid(E,y):
    return float(np.trapezoid(y,E))

def calibrate_fig7_response(cfg, ngrid=800):
    """
    Calibrate the digitised shapes to the exact published BDT-selected
    cosmic-neutrino rates: 220 tracks/year and 60 showers/year for the
    paper's Eq.(6) E^-2 reference flux.  This makes absolute normalization
    anchored to Tables 1/2 even when Zenodo data cannot be machine-read.
    """
    E=np.logspace(3,8,ngrid)
    R=fig7_response(E)
    phi0=float(cfg["arca"]["reference_flux_phi0_at_1GeV"])
    phi=phi0*E**(-float(cfg["arca"]["reference_gamma"]))
    sec=365.25*86400.0

    # panel c: each curve is convolved with the corresponding per-flavor flux;
    # all-flavour NC curve is already the sum over flavors.
    sh=R["shower"]
    # one hemisphere = 2pi sr
    raw_sh=sec*2*math.pi*integrate_loggrid(
        E,phi*(sh["down_ecc"]+sh["down_tcc"]+sh["down_nc"]
              +sh["up_ecc"]+sh["up_tcc"]+sh["up_nc"])
    )*1e4
    sf_sh=float(cfg["arca"]["reference_shower_rate_per_year"])/raw_sh

    tr=R["track"]; raw_tr=0.0
    for k,(clo,chi) in R["track_cos_bins"].items():
        domega=2*math.pi*(chi-clo)
        raw_tr += sec*domega*integrate_loggrid(E,phi*tr[k])*1e4
    sf_tr=float(cfg["arca"]["reference_track_rate_per_year"])/raw_tr
    return sf_sh,sf_tr

def published_rate_crosscheck(cfg):
    sf_sh,sf_tr=calibrate_fig7_response(cfg)
    return dict(shower_scale=sf_sh,track_scale=sf_tr,
                shower_target=float(cfg["arca"]["reference_shower_rate_per_year"]),
                track_target=float(cfg["arca"]["reference_track_rate_per_year"]))
