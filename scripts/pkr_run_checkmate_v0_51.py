#!/usr/bin/env python3
from __future__ import annotations

import gzip, json, os, re, shutil, subprocess, sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "cache" / "v0.51"
OUT = ROOT / "results" / "v0.51"
MANIFEST = WORK / "benchmark_manifest.csv"
PREP = ROOT / "scripts" / "pkr_prepare_checkmate_v0_51.py"
SHOWER = WORK / "bin" / "pythia_lhe_to_hepmc_v0_51"


def run(cmd, cwd=None, log=None, env=None):
    p = subprocess.run(cmd, cwd=cwd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=env)
    if log:
        Path(log).parent.mkdir(parents=True, exist_ok=True)
        Path(log).write_text(p.stdout)
    if p.returncode:
        tail = "\n".join(p.stdout.splitlines()[-30:])
        raise RuntimeError(f"command failed ({p.returncode}): {' '.join(map(str,cmd))}\n{tail}")
    return p.stdout


def executable_from_env(name, candidates):
    e = os.environ.get(name)
    if e and Path(e).exists(): return Path(e)
    for c in candidates:
        c = Path(c).expanduser()
        if c.exists(): return c
    return None


def find_lhe(procdir: Path) -> Path:
    cand = list(procdir.glob("Events/run_*/unweighted_events.lhe.gz")) + list(procdir.glob("Events/run_*/unweighted_events.lhe"))
    if not cand:
        raise RuntimeError(f"No LHE found below {procdir}")
    return sorted(cand)[-1]


def ensure_plain_lhe(src: Path) -> Path:
    if src.suffix != ".gz": return src
    dst = src.with_suffix("")
    if not dst.exists() or dst.stat().st_mtime < src.stat().st_mtime:
        with gzip.open(src, "rb") as fi, open(dst, "wb") as fo:
            shutil.copyfileobj(fi, fo)
    return dst


def parse_lhe_xsec_pb(path: Path) -> float:
    op = gzip.open if path.suffix == ".gz" else open
    with op(path, "rt", errors="replace") as f:
        in_init = False
        nprup = None
        rows = []
        for line in f:
            s = line.strip()
            if s == "<init>":
                in_init = True; continue
            if in_init and nprup is None:
                parts = s.split()
                if len(parts) >= 10:
                    nprup = int(parts[-1])
                continue
            if in_init and nprup is not None:
                if s.startswith("</init>"): break
                parts = s.split()
                if len(parts) >= 4:
                    rows.append(float(parts[0]))
                    if len(rows) >= nprup: break
        if not rows:
            raise RuntimeError(f"Could not parse LHE <init> cross section from {path}")
        return float(sum(rows))


def write_checkmate_card(row, xsec_pb: float, hepmc: Path) -> Path:
    tag = row.tag
    card = Path(row.checkmate_card)
    cmout = WORK / "checkmate" / tag
    text = f"""[Parameters]
Name: pkr_v051_{tag}
Analyses: cms_2107_13021
OutputDirectory: {cmout}
OutputExists: Overwrite
SkipParamCheck: True
EffTab: True
Multibin: detailed
Expected: True
Uplim: True
Model: simple
InvisiblePIDs: 5000521
RandomSeed: 5101
TotalResultFileColumns: analysis, sr, o, b, db, s, ds, s95obs, s95exp, robscons, rexpcons
BestPerAnalysisResultFileColumns: analysis, sr, o, b, db, s, ds, s95obs, s95exp, robscons, rexpcons

[signal]
Events: {hepmc}
XSect: {xsec_pb:.12g} PB
XSectErr: 10 %
"""
    card.write_text(text)
    return card


def parse_multibin(root: Path):
    files = list(root.rglob("multibin_limits/results.dat"))
    if not files: return None, None, None
    text = files[-1].read_text(errors="replace")
    blocks = re.split(r"=+", text)
    best = None
    for b in blocks:
        if "cms_2107_13021" not in b: continue
        mo = re.search(r"Observed upper limit:\s*mu\s*=\s*([0-9.eE+\-]+)", b)
        mc = re.search(r"Observed CLs for mu = 1:\s*([0-9.eE+\-]+)", b)
        if mo:
            mu = float(mo.group(1)); cls = float(mc.group(1)) if mc else None
            best = (mu, 1.0/mu if mu > 0 else None, cls)
    return best if best else (None, None, None)


def parse_singlebin(root: Path):
    files = list(root.rglob("evaluation/total_results.txt"))
    if not files: return None
    mx = None
    for line in files[-1].read_text(errors="replace").splitlines()[1:]:
        p = line.split()
        if not p or "cms_2107_13021" not in p[0]: continue
        # Requested column order puts robscons at index 9 (analysis,sr,o,b,db,s,ds,s95obs,s95exp,robscons,...)
        try: r = float(p[9])
        except Exception: continue
        mx = r if mx is None else max(mx, r)
    return mx


def verdict(r):
    if r is None: return "NO_MULTIBIN_RESULT"
    if r > 1.20: return "ROBUSTLY_EXCLUDED_CHECKMATE"
    if r < 0.80: return "ROBUSTLY_SURVIVES_CHECKMATE"
    return "MARGINAL_OPERATOR_LEVEL_RECAST_REQUIRED"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)
    run([sys.executable, str(PREP)], log=OUT/"prepare.log")
    d = pd.read_csv(MANIFEST)

    mg5 = executable_from_env("PKR_MG5", [WORK/"MG5_aMC"/"bin"/"mg5_aMC", Path.home()/"MG5_aMC"/"bin"/"mg5_aMC"])
    checkmate = executable_from_env("PKR_CHECKMATE", [WORK/"checkmate-install"/"bin"/"CheckMATE", WORK/"checkmate2"/"bin"/"CheckMATE"])
    if mg5 is None: raise RuntimeError("MadGraph not found. Run bash setup_pkr_checkmate_v0_51.sh")
    if checkmate is None: raise RuntimeError("CheckMATE not found. Run bash setup_pkr_checkmate_v0_51.sh")
    if not SHOWER.exists(): raise RuntimeError("Pythia shower bridge missing. Run bash setup_pkr_checkmate_v0_51.sh")

    rows = []
    force = os.environ.get("PKR_CM_FORCE", "0") == "1"
    for i, row in d.iterrows():
        tag = row.tag
        procdir = Path(row.mg5_dir)
        mglog = OUT / f"mg5_{tag}.log"
        if force or not list(procdir.glob("Events/run_*/unweighted_events.lhe*")):
            run([str(mg5), str(row.mg5_card)], log=mglog)
        lhe_src = find_lhe(procdir)
        xsec_pb = parse_lhe_xsec_pb(lhe_src)
        lhe = ensure_plain_lhe(lhe_src)
        hepmc = Path(row.hepmc)
        if force or not hepmc.exists() or hepmc.stat().st_size < 1000:
            run([str(SHOWER), str(lhe), str(hepmc), str(5101+i)], log=OUT/f"pythia_{tag}.log")
        card = write_checkmate_card(row, xsec_pb, hepmc)
        cmroot = WORK / "checkmate" / tag
        if force and cmroot.exists(): shutil.rmtree(cmroot)
        run([str(checkmate), str(card)], log=OUT/f"checkmate_{tag}.log")
        mu95, rmb, cls = parse_multibin(cmroot)
        rsingle = parse_singlebin(cmroot)
        rows.append(dict(
            tag=tag, model=row.model, convention=row.convention,
            xsec_pb=xsec_pb, BR_invisible=float(row.BR_invisible),
            gchi_width_matched=float(row.gchi_width_matched),
            multibin_mu95=mu95, multibin_r=rmb, multibin_CLs_mu1=cls,
            max_singlebin_robscons=rsingle, verdict=verdict(rmb),
        ))
        print(f"{tag}: xsec={xsec_pb:.6g} pb  r_MB={rmb}  r_single={rsingle}  {verdict(rmb)}")

    R = pd.DataFrame(rows)
    R.to_csv(OUT/"checkmate_gate.csv", index=False)
    summary=[]
    for model,q in R.groupby("model"):
        vals=[x for x in q.multibin_r.tolist() if pd.notna(x)]
        if len(vals)!=len(q): st="INCOMPLETE"
        elif max(vals)<0.80: st="COLLIDER_MONOJET_PASS_BOTH_WIDTH_CONVENTIONS"
        elif min(vals)>1.20: st="COLLIDER_MONOJET_FAIL_BOTH_WIDTH_CONVENTIONS"
        else: st="MARGINAL_OR_WIDTH_CONVENTION_DEPENDENT"
        summary.append(dict(model=model,min_multibin_r=min(vals) if vals else None,max_multibin_r=max(vals) if vals else None,status=st))
    S=pd.DataFrame(summary); S.to_csv(OUT/"summary.csv",index=False)
    meta={
        "version":"v0.51",
        "analysis":"CheckMATE2 cms_2107_13021 = CMS EXO-20-004 jets+MET",
        "signal":"MG5 LO pp -> Xd Xd~ j, Pythia8 shower/hadronization, CheckMATE Delphes+analysis",
        "matching":"DMsimp spin-1 vector quark production; gVXd chosen to reproduce each PKR invisible partial width; total width fixed to v0.45 PKR total width",
        "nwa_justification":"all v0.45 Gamma/M < 9e-4; invisible decay tensor does not affect visible recoil in the on-shell limit",
        "primary_gate":"multibin observed upper limit mu95; r=1/mu95",
        "robust_thresholds":{"survive":"r<0.8","exclude":"r>1.2","marginal":"0.8<=r<=1.2"},
        "limitations":["LO one-hard-jet matrix element without MLM merging","10% signal cross-section systematic","not an off-shell custom PKR UFO; dedicated operator UFO required if marginal"],
    }
    (OUT/"metadata.json").write_text(json.dumps(meta,indent=2))
    print("\n=== v0.51 CheckMATE PKR monojet gate ==="); print(R.to_string(index=False))
    print("\n=== v0.51 summary ==="); print(S.to_string(index=False))

if __name__ == "__main__": main()
