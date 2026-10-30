#!/usr/bin/env python3
"""Prepare v0.51 PKR-matched monojet recast cards.

Physics logic
-------------
The PKR mediator is extremely narrow (Gamma/M < 9e-4). For an on-shell invisible
Z' produced through the same vector quark current, the detector-level monojet
recoil depends on production and BR(Z'->invisible), not on the unobserved
spinor structure of the decay. We therefore use the standard DMsimp_s_spin1
UFO with a nearly massless Dirac invisible fermion Xd and choose gVXd such that
its partial width exactly equals the PKR invisible width from v0.45. The total
mediator width is set to the PKR total width, so production x BR is matched.

This is a narrow-width-equivalent collider recast, not an off-shell validation
of the charge-conjugate derivative/dipole UFO vertex itself. If CheckMATE lands
near r=1, a dedicated operator UFO is then required.
"""
from __future__ import annotations

import csv, math, os
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
IN45 = ROOT / "results" / "v0.45" / "width_br_dijet_gate.csv"
WORK = ROOT / "cache" / "v0.51"
OUT = ROOT / "results" / "v0.51"
MODEL = WORK / "models" / "DMsimp_s_spin1"
NEVENTS = int(os.environ.get("PKR_CM_NEVENTS", "100000"))
MCHI = float(os.environ.get("PKR_CM_MCHI_GEV", "0.001"))
PTJMIN = float(os.environ.get("PKR_CM_PTJMIN_GEV", "100.0"))
PDF_LHAID = int(os.environ.get("PKR_CM_LHAID", "14400"))  # CT18NLO central


def gchi_for_width(M: float, mchi: float, gamma_inv: float) -> float:
    r = (mchi / M) ** 2
    if 4*r >= 1:
        raise ValueError("invisible state above threshold")
    phase = (1.0 + 2.0*r) * math.sqrt(1.0 - 4.0*r)
    return math.sqrt(gamma_inv * 12.0 * math.pi / (M * phase))


def slug(model: str, conv: str) -> str:
    return f"{model}__{conv}".replace("-", "_")


def mg5_card(tag: str, row: pd.Series, gchi: float) -> str:
    M = float(row.M_GeV); gq = float(row.gq); width = float(row.Gamma_total_GeV)
    proc = WORK / "mg5" / tag
    # Explicitly zero all non-PKR visible couplings, including FCNC terms that
    # exist in some historical DMsimp UFO mirrors.
    sets = [
        f"set MY1 {M:.16g}",
        f"set MXd {MCHI:.16g}",
        f"set WY1 {width:.16g}",
        f"set gVXd {gchi:.16g}",
        "set gAXd 0",
        "set gVXc 0",
    ]
    for p in ["gVd11","gVu11","gVd22","gVu22","gVd33","gVu33"]:
        sets.append(f"set {p} {gq:.16g}")
    for p in [
        "gAd11","gAu11","gAd22","gAu22","gAd33","gAu33",
        "gVl11","gVl22","gVl33","gAl11","gAl22","gAl33",
        "gnu11","gnu22","gnu33","gVu31","gAu31","gVd31","gAd31","gVh",
    ]:
        sets.append(f"set {p} 0")
    runsets = [
        f"set nevents {NEVENTS}",
        "set ebeam1 6500",
        "set ebeam2 6500",
        "set pdlabel lhapdf",
        f"set lhaid {PDF_LHAID}",
        f"set ptj {PTJMIN:.16g}",
        "set etaj 5.0",
        "set maxjetflavor 5",
        "set ickkw 0",
    ]
    return "\n".join([
        "set automatic_html_opening False",
        f"import model {MODEL}",
        "generate p p > xd xd~ j",
        f"output {proc} -f",
        f"launch {proc} -f",
        *runsets,
        *sets,
        "done",
        "",
    ])


def checkmate_card(tag: str, hepmc: Path, xsec_pb: float) -> str:
    outdir = WORK / "checkmate" / tag
    return f"""[Parameters]
Name: pkr_v051_{tag}
Analyses: cms_2107_13021
OutputDirectory: {outdir}
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


def main() -> None:
    if not IN45.exists():
        raise RuntimeError(f"missing {IN45}")
    WORK.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    (WORK / "cards").mkdir(parents=True, exist_ok=True)
    (WORK / "mg5").mkdir(parents=True, exist_ok=True)
    (WORK / "hepmc").mkdir(parents=True, exist_ok=True)
    (WORK / "checkmate").mkdir(parents=True, exist_ok=True)

    d = pd.read_csv(IN45)
    manifest = []
    for _, r in d.iterrows():
        tag = slug(str(r.model), str(r.convention))
        gchi = gchi_for_width(float(r.M_GeV), MCHI, float(r.Gamma_nu_GeV))
        card = WORK / "cards" / f"mg5_{tag}.txt"
        card.write_text(mg5_card(tag, r, gchi))
        manifest.append({
            "tag": tag,
            "model": r.model,
            "convention": r.convention,
            "M_GeV": float(r.M_GeV),
            "gq": float(r.gq),
            "kappa_GeV_inv": float(r.kappa_GeV_inv),
            "Gamma_nu_GeV": float(r.Gamma_nu_GeV),
            "Gamma_total_GeV": float(r.Gamma_total_GeV),
            "Gamma_over_M": float(r.Gamma_over_M),
            "BR_invisible": float(r.BR_invisible),
            "gchi_width_matched": gchi,
            "mchi_GeV": MCHI,
            "nevents": NEVENTS,
            "ptjmin_GeV": PTJMIN,
            "lhaid": PDF_LHAID,
            "mg5_card": str(card),
            "mg5_dir": str(WORK / "mg5" / tag),
            "hepmc": str(WORK / "hepmc" / f"{tag}.hepmc"),
            "checkmate_card": str(WORK / "cards" / f"checkmate_{tag}.in"),
        })

    M = pd.DataFrame(manifest)
    M.to_csv(WORK / "benchmark_manifest.csv", index=False)
    M.to_csv(OUT / "benchmark_manifest.csv", index=False)
    print("\n=== v0.51 prepared benchmarks ===")
    print(M[["model","convention","gchi_width_matched","BR_invisible","Gamma_over_M"]].to_string(index=False))
    print(f"\nPrepared {len(M)} MG5 cards; nevents={NEVENTS}, pTj(min)={PTJMIN:g} GeV, LHAID={PDF_LHAID}.")


if __name__ == "__main__":
    main()
