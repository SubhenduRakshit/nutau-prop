#!/usr/bin/env python3
"""v0.47: exact public CMS EXO-20-004 HEPData low-mass monojet gate.

Downloads the HEPData record for CMS EXO-20-004, locates the vector-mediator
Figure 20-b table (observed g_q^95 versus mediator mass), extracts the 50 GeV
point if present, and compares it to the v0.46 PKR NWA exclusion thresholds.

No interpolation is used when an exact 50 GeV row exists. If the public table
does not contain 50 GeV, the script reports the nearest bracketing points and
leaves the gate OPEN rather than inventing an interpolation.
"""
from __future__ import annotations

import csv
import io
import json
import math
import re
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "v0.47"
OUT.mkdir(parents=True, exist_ok=True)
IN46 = ROOT / "results" / "v0.46" / "monojet_gate_summary.csv"
RECORD_URL = "https://www.hepdata.net/record/ins1894408?format=json"
UA = "nutau-prop/0.47 (+public CMS HEPData reinterpretation)"


def fetch_bytes(url: str) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "*/*"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return r.read()


def fetch_json(url: str):
    return json.loads(fetch_bytes(url).decode("utf-8"))


def text_of(x) -> str:
    if isinstance(x, dict):
        return " ".join(text_of(v) for v in x.values())
    if isinstance(x, list):
        return " ".join(text_of(v) for v in x)
    return str(x)


def candidate_tables(record: dict) -> list[dict]:
    tabs = record.get("data_tables") or record.get("tables") or []
    scored = []
    for t in tabs:
        s = text_of(t).lower()
        score = 0
        if "figure 20" in s or "fig. 20" in s or "fig 20" in s:
            score += 4
        if "20-b" in s or "20b" in s:
            score += 4
        if "vector" in s:
            score += 2
        if "g_q" in s or "gq" in s or "quark coupling" in s or "coupling" in s:
            score += 2
        if score:
            scored.append((score, t))
    scored.sort(key=lambda z: z[0], reverse=True)
    return [t for _, t in scored]


def data_links(t: dict) -> list[str]:
    out = []
    d = t.get("data")
    if isinstance(d, dict):
        # Prefer machine-readable table formats.
        for k in ("csv", "json", "yaml", "yml"):
            v = d.get(k)
            if isinstance(v, str):
                out.append(v)
        for v in d.values():
            if isinstance(v, str) and v not in out:
                out.append(v)
    elif isinstance(d, str):
        out.append(d)
    for k in ("csv", "json", "yaml"):
        v = t.get(k)
        if isinstance(v, str) and v not in out:
            out.append(v)
    # HEPData sometimes provides a canonical table URL only. Try format query.
    base = t.get("url") or t.get("record") or t.get("table")
    if isinstance(base, str):
        for fmt in ("csv", "json"):
            sep = "&" if "?" in base else "?"
            out.append(f"{base}{sep}format={fmt}")
    return out


def read_csv_relaxed(raw: bytes) -> pd.DataFrame:
    txt = raw.decode("utf-8", errors="replace")
    lines = [ln for ln in txt.splitlines() if ln.strip() and not ln.lstrip().startswith("#")]
    return pd.read_csv(io.StringIO("\n".join(lines)))


def numeric_series(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s.astype(str).str.replace(r"[^0-9eE+\-.]", "", regex=True), errors="coerce")


def extract_50_from_df(df: pd.DataFrame):
    cols = list(df.columns)
    num = {c: numeric_series(df[c]) for c in cols}

    # Find a mediator-mass column containing exactly 50 GeV.
    mass_candidates = []
    for c in cols:
        name = str(c).lower()
        vals = num[c]
        has50 = bool(((vals - 50.0).abs() < 1e-8).any())
        if has50:
            score = (3 if "med" in name else 0) + (2 if "mass" in name else 0) + (1 if "gev" in name else 0)
            mass_candidates.append((score, c))
    if not mass_candidates:
        return None
    mass_candidates.sort(reverse=True)
    mcol = mass_candidates[0][1]
    mvals = num[mcol]
    rows = df[(mvals - 50.0).abs() < 1e-8].copy()
    if rows.empty:
        return None

    # Identify observed gq/limit column. Prefer explicit observed + gq/coupling naming.
    value_candidates = []
    for c in cols:
        if c == mcol:
            continue
        name = str(c).lower()
        vals = numeric_series(rows[c]).dropna()
        if vals.empty:
            continue
        v = float(vals.iloc[0])
        if not (0.0 < v < 5.0):
            continue
        score = 0
        if "observ" in name: score += 5
        if "g_q" in name or "gq" in name: score += 4
        if "coupl" in name: score += 2
        if "limit" in name or "95" in name: score += 1
        # expected/uncertainty columns are disfavoured unless nothing else exists
        if "expect" in name: score -= 3
        if "error" in name or "unc" in name: score -= 3
        value_candidates.append((score, c, v))
    if not value_candidates:
        return None
    value_candidates.sort(key=lambda z: z[0], reverse=True)
    score, vcol, val = value_candidates[0]
    return dict(mass_column=str(mcol), value_column=str(vcol), observed_gq95=float(val), row=rows.iloc[0].to_dict(), score=int(score))


def nearest_mass_rows(df: pd.DataFrame):
    out=[]
    for c in df.columns:
        name=str(c).lower()
        if "mass" not in name and "med" not in name:
            continue
        vals=numeric_series(df[c])
        q=df.assign(_m=vals).dropna(subset=["_m"])
        q=q.iloc[(q._m-50.0).abs().argsort()[:4]]
        if len(q):
            out.append(dict(column=str(c), rows=q.to_dict(orient="records")))
    return out


def main():
    if not IN46.exists():
        raise RuntimeError(f"Missing {IN46}; run v0.46 first")

    rec = fetch_json(RECORD_URL)
    (OUT / "hepdata_record.json").write_text(json.dumps(rec, indent=2))
    tabs = candidate_tables(rec)
    inventory=[]
    exact=[]

    for i,t in enumerate(tabs):
        meta={k:t.get(k) for k in ("name","location","description","doi") if k in t}
        links=data_links(t)
        entry={"rank":i,"meta":meta,"links":links,"attempts":[]}
        for j,u in enumerate(links):
            try:
                raw=fetch_bytes(u)
                suf="csv" if ("csv" in u.lower() or raw[:200].count(b",")>2) else "dat"
                path=OUT/f"candidate_{i:02d}_{j:02d}.{suf}"
                path.write_bytes(raw)
                try:
                    df=read_csv_relaxed(raw)
                    hit=extract_50_from_df(df)
                    entry["attempts"].append({"url":u,"saved":path.name,"columns":[str(c) for c in df.columns],"nrows":len(df),"exact50":hit is not None})
                    if hit is not None:
                        exact.append(dict(rank=i,table_meta=meta,url=u,file=path.name,**hit))
                    else:
                        entry["nearest"]=nearest_mass_rows(df)
                except Exception as e:
                    entry["attempts"].append({"url":u,"saved":path.name,"parse_error":repr(e)})
            except Exception as e:
                entry["attempts"].append({"url":u,"fetch_error":repr(e)})
        inventory.append(entry)

    (OUT/"candidate_inventory.json").write_text(json.dumps(inventory,indent=2,default=str))

    v46=pd.read_csv(IN46)
    rows=[]
    if exact:
        # Prefer the highest-ranked table and strongest semantic column score.
        exact.sort(key=lambda z:(-z["score"], z["rank"]))
        best=exact[0]
        g95=float(best["observed_gq95"])
        for _,r in v46.iterrows():
            lo=float(r.critical_gq95_min); hi=float(r.critical_gq95_max)
            if g95 < lo:
                verdict="EXCLUDED_FOR_BOTH_WIDTH_CONVENTIONS"
            elif g95 > hi:
                verdict="SURVIVES_FOR_BOTH_WIDTH_CONVENTIONS"
            else:
                verdict="WIDTH_CONVENTION_DECIDES"
            rows.append(dict(model=r.model,CMS50_observed_gq95=g95,critical_gq95_min=lo,critical_gq95_max=hi,verdict=verdict,source_table=best["table_meta"].get("name"),mass_column=best["mass_column"],value_column=best["value_column"]))
        pd.DataFrame(rows).to_csv(OUT/"monojet_exact50_gate.csv",index=False)
        (OUT/"exact50_extraction.json").write_text(json.dumps(best,indent=2,default=str))
    else:
        pd.DataFrame(columns=["model","CMS50_observed_gq95","critical_gq95_min","critical_gq95_max","verdict"]).to_csv(OUT/"monojet_exact50_gate.csv",index=False)

    meta={
      "version":"v0.47",
      "analysis":"CMS EXO-20-004 / JHEP 11 (2021) 153",
      "hepdata_record":"ins1894408",
      "target":"Figure 20-b vector mediator observed gq^95 at mmed=50 GeV",
      "policy":"use exact 50 GeV row only; do not interpolate if absent",
      "exact_50_found":bool(exact),
      "next_if_absent":"validated MadAnalysis5 CMS-EXO-20-004 simplified-likelihood recast",
    }
    (OUT/"metadata.json").write_text(json.dumps(meta,indent=2))

    print("\n=== v0.47 HEPData discovery ===")
    print(f"candidate tables: {len(tabs)}; exact 50 GeV extraction: {bool(exact)}")
    if exact:
        best=exact[0]
        print(f"table={best['table_meta'].get('name')}  observed gq95(50 GeV)={best['observed_gq95']:.8g}")
        print("\n=== v0.47 exact 50 GeV monojet gate ===")
        print(pd.DataFrame(rows).to_string(index=False))
    else:
        print("No exact 50 GeV numeric row was found. See candidate_inventory.json; monojet gate remains OPEN and v0.48 should use the validated MadAnalysis5 implementation.")

if __name__=="__main__":
    main()
