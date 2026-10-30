#!/usr/bin/env python3
"""v0.68: exactly one conservative CT18NLO CheckMATE benchmark.

Runs only the PKR Eq. 4.4 derivative benchmark with the `one_Weyl_pair`
width convention.  This is the weaker of the two retained invisible-width
conventions, so it is the conservative detector-level test.

Workflow
--------
1. Build a fresh DMsimp spin-1 surrogate process for pp -> xd xd~ j.
2. Read the exact Eq. 4.4/one-Weyl M=50 GeV benchmark from v0.45.
3. Persist the exact mass, gq, invisible-width-matched gchi and total width.
4. Use CT18NLO (LHAID 14400) and the validated v0.67 LHAPDF compatibility patch.
5. Generate a fresh LHE sample and require xsec/BR consistency with the
   independent v0.67 CT18 smoke sample to 2%.
6. Shower/hadronize with Pythia8 and run exactly one CheckMATE
   cms_2107_13201 (CMS EXO-20-004) analysis.

No scan is performed.  A multibin result is primary.  If CheckMATE exits
nonzero only during multibin statistics after producing the ordinary analysis,
the single-bin number is preserved as diagnostic evidence only, never as a
final exclusion/survival conclusion.
"""
from __future__ import annotations

import gzip
import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "cache" / "v0.68_checkmate_single"
OUT = ROOT / "results" / "v0.68_checkmate_single"
IN45 = ROOT / "results" / "v0.45" / "width_br_dijet_gate.csv"
SMOKE67 = ROOT / "results" / "v0.67_ct18_smoke" / "ct18_smoke.json"
TAG = "eq4p4_derivative__one_Weyl_pair"
ANALYSIS = os.environ.get("PKR_CM_ANALYSIS", "cms_2107_13201")
XSEC_BR_TOL = 0.02

# Reuse the now-certified CT18 bridge patch and deterministic MG5 machinery.
P67 = ROOT / "scripts" / "pkr_ct18_smoke_v0_67.py"
spec67 = importlib.util.spec_from_file_location("pkr_ct18_v067", P67)
if spec67 is None or spec67.loader is None:
    raise RuntimeError(f"cannot import {P67}")
v067 = importlib.util.module_from_spec(spec67)
spec67.loader.exec_module(v067)
v065 = v067.v065
base = v065.base
cfg = v065.cfg

# Reuse only parsing helpers from the old CheckMATE driver, not its invalid v0.51 samples.
P51 = ROOT / "scripts" / "pkr_run_checkmate_v0_51.py"
spec51 = importlib.util.spec_from_file_location("pkr_cm_v051", P51)
if spec51 is None or spec51.loader is None:
    raise RuntimeError(f"cannot import {P51}")
cm51 = importlib.util.module_from_spec(spec51)
spec51.loader.exec_module(cm51)


def run(cmd, *, cwd=None, log: Path | None = None, allow_nonzero=False):
    p = subprocess.run([str(x) for x in cmd], cwd=cwd, text=True,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if log is not None:
        log.parent.mkdir(parents=True, exist_ok=True)
        log.write_text(p.stdout)
    if p.returncode and not allow_nonzero:
        tail = "\n".join(p.stdout.splitlines()[-60:])
        raise RuntimeError(f"command failed rc={p.returncode}: {' '.join(map(str, cmd))}\n{tail}")
    return p.returncode, p.stdout


def ensure_plain_lhe(src: Path, dst: Path) -> Path:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.suffix == ".gz":
        with gzip.open(src, "rb") as fi, dst.open("wb") as fo:
            shutil.copyfileobj(fi, fo)
    else:
        shutil.copy2(src, dst)
    return dst


def exact_benchmark():
    d = pd.read_csv(IN45)
    q = d[(d["model"] == "eq4p4_derivative") & (d["convention"] == "one_Weyl_pair")]
    if len(q) != 1:
        raise RuntimeError(f"expected exactly one Eq4.4 one_Weyl benchmark in {IN45}, found {len(q)}")
    return q.iloc[0]


def write_checkmate_card(hepmc: Path, xsec_pb: float) -> tuple[Path, Path]:
    card = WORK / "cards" / f"checkmate_{TAG}.in"
    cmout = WORK / "checkmate" / TAG
    card.parent.mkdir(parents=True, exist_ok=True)
    text = f"""[Parameters]
Name: pkr_v068_{TAG}
Analyses: {ANALYSIS}
OutputDirectory: {cmout}
OutputExists: Overwrite
SkipParamCheck: True
EffTab: True
Multibin: detailed
Expected: True
Uplim: True
Model: simple
InvisiblePIDs: 5000521
RandomSeed: 6801
TotalResultFileColumns: analysis, sr, o, b, db, s, ds, s95obs, s95exp, robscons, rexpcons
BestPerAnalysisResultFileColumns: analysis, sr, o, b, db, s, ds, s95obs, s95exp, robscons, rexpcons

[signal]
Events: {hepmc}
XSect: {xsec_pb:.12g} PB
XSectErr: 10 %%
"""
    card.write_text(text)
    return card, cmout


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)
    if not IN45.exists() or not SMOKE67.exists():
        raise RuntimeError("missing v0.45 benchmark table or v0.67 CT18 smoke result")

    row = exact_benchmark()
    M = float(row.M_GeV)
    gq = float(row.gq)
    gamma_inv = float(row.Gamma_nu_GeV)
    gamma_total = float(row.Gamma_total_GeV)
    br = float(row.BR_invisible)
    gchi = cfg.gchi_for_width(M, cfg.MCHI, gamma_inv)

    procdir = WORK / "mg5" / TAG
    if procdir.exists():
        shutil.rmtree(procdir)
    (WORK / "cards").mkdir(parents=True, exist_ok=True)
    build_card = WORK / "cards" / f"build_{TAG}.txt"
    base.build_process(procdir, build_card, OUT / "mg5_build.log")

    param = procdir / "Cards" / "param_card.dat"
    run_card = procdir / "Cards" / "run_card.dat"
    base.edit_param_card(param, M, cfg.MCHI, gq, gchi, gamma_total)
    base.configure_run_card(run_card, seed=6801)
    ok_before, *_ = base.certify_persisted_param_card(param, M, cfg.MCHI, gq, gchi, gamma_total)
    if not ok_before:
        raise RuntimeError("Eq4.4 parameter card failed before generation")

    # Apply the validated v0.67 CT18 bridge patch, explicitly rebuild libpdf,
    # and retain v0.66 process-local Safari suppression through base.generate_events.
    v067.patch_lhapdf_lambda_fallback_v067(procdir)
    base.generate_events(procdir, OUT / "mg5_generate.log")

    ok_after, actual_M, actual_mchi, actual_gchi, actual_width, actual_gqs = \
        base.certify_persisted_param_card(param, M, cfg.MCHI, gq, gchi, gamma_total)
    if not ok_after:
        raise RuntimeError("Eq4.4 parameter card changed during generation")

    lhe_src = cfg.find_lhe(procdir)
    banner = cfg.find_banner(procdir)
    xsec_lhe = cfg.parse_lhe_init_xsec_pb(lhe_src)
    xsec_banner = cfg.parse_banner_xsec_pb(banner)
    xsec_rel = abs(xsec_banner - xsec_lhe) / max(abs(xsec_lhe), 1e-30)
    banner_text = banner.read_text(errors="replace")
    has_14400 = bool(re.search(r"(?im)^\s*14400\s*=\s*lhaid\b", banner_text))

    smoke = json.loads(SMOKE67.read_text())
    smoke_xsec = float(smoke["metrics"]["lhe_init_xsec_pb"])
    smoke_br = float(smoke["benchmark"]["BR_invisible"])
    ratio_here = xsec_lhe / br
    ratio_smoke = smoke_xsec / smoke_br
    ratio_rel = abs(ratio_here - ratio_smoke) / ratio_smoke
    if ratio_rel > XSEC_BR_TOL:
        raise RuntimeError(
            f"CT18 xsec/BR consistency failed: here={ratio_here:.8g} pb, "
            f"v0.67={ratio_smoke:.8g} pb, rel={ratio_rel:.3%}"
        )

    lhe = ensure_plain_lhe(lhe_src, WORK / "lhe" / f"{TAG}.lhe")
    shower = Path(os.environ.get("PKR_SHOWER", ROOT / "cache" / "v0.51" / "bin" / "pythia_lhe_to_hepmc_v0_51"))
    checkmate = Path(os.environ.get("PKR_CHECKMATE", ROOT / "cache" / "v0.51" / "checkmate-install" / "bin" / "CheckMATE"))
    if not shower.exists():
        raise RuntimeError(f"Pythia LHE->HepMC bridge missing: {shower}")
    if not checkmate.exists():
        raise RuntimeError(f"CheckMATE missing: {checkmate}")

    hepmc = WORK / "hepmc" / f"{TAG}.hepmc"
    hepmc.parent.mkdir(parents=True, exist_ok=True)
    run([shower, lhe, hepmc, "6801"], log=OUT / "pythia.log")
    if not hepmc.exists() or hepmc.stat().st_size < 1000:
        raise RuntimeError("Pythia produced no usable HepMC sample")

    card, cmout = write_checkmate_card(hepmc, xsec_lhe)
    if cmout.exists():
        shutil.rmtree(cmout)
    cm_rc, cm_stdout = run([checkmate, card], log=OUT / "checkmate.log", allow_nonzero=True)

    mu95, rmb, cls = cm51.parse_multibin(cmout)
    rsingle = cm51.parse_singlebin(cmout)
    if cm_rc != 0 and rmb is None and rsingle is None:
        tail = "\n".join(cm_stdout.splitlines()[-60:])
        raise RuntimeError(f"CheckMATE failed before any usable analysis output (rc={cm_rc})\n{tail}")

    if rmb is not None:
        if rmb > 1.2:
            detector_status = "ROBUST_MULTIBIN_EXCLUSION_OF_CONSERVATIVE_WIDTH_CONVENTION"
        elif rmb < 0.8:
            detector_status = "ROBUST_MULTIBIN_SURVIVAL_OF_CONSERVATIVE_WIDTH_CONVENTION"
        else:
            detector_status = "MARGINAL_MULTIBIN_DEDICATED_OPERATOR_RECAST_REQUIRED"
        finality = "multibin_available"
    else:
        detector_status = "MULTIBIN_UNRESOLVED_SINGLEBIN_DIAGNOSTIC_ONLY"
        finality = "not_final"

    result = {
        "version": "v0.68_checkmate_single",
        "status": "PASS" if (rmb is not None or rsingle is not None) else "FAIL",
        "benchmark": {
            "tag": TAG,
            "model": "eq4p4_derivative",
            "convention": "one_Weyl_pair",
            "why_this_one": "lower BR_inv than operator_plus_hc; conservative NWA monojet test",
            "M_GeV": actual_M,
            "mchi_GeV": actual_mchi,
            "gq": gq,
            "gchi_width_matched": actual_gchi,
            "Gamma_total_GeV": actual_width,
            "Gamma_over_M": actual_width / actual_M,
            "BR_invisible": br,
            "lhaid": 14400,
        },
        "generator_checks": {
            "param_card_before": bool(ok_before),
            "param_card_after": bool(ok_after),
            "banner_records_lhaid_14400": has_14400,
            "banner_vs_lhe_rel_diff": xsec_rel,
            "xsec_over_BR_consistent_with_v067_within_2pct": bool(ratio_rel <= XSEC_BR_TOL),
            "xsec_over_BR_relative_difference_vs_v067": ratio_rel,
        },
        "cross_section_pb": xsec_lhe,
        "checkmate": {
            "analysis": ANALYSIS,
            "returncode": cm_rc,
            "multibin_mu95": mu95,
            "multibin_r": rmb,
            "multibin_CLs_mu1": cls,
            "max_singlebin_robscons": rsingle,
            "detector_status": detector_status,
            "finality": finality,
        },
        "interpretation_rule": (
            "Only a valid multibin result can close the collider gate. If multibin is unavailable, "
            "single-bin robscons is diagnostic only. Because this is the lower-BR width convention, "
            "a robust multibin exclusion here also excludes the higher-BR operator_plus_hc convention "
            "within the validated NWA mapping."
        ),
    }
    (OUT / "single_benchmark.json").write_text(json.dumps(result, indent=2) + "\n")
    pd.DataFrame([{
        "tag": TAG, "BR_invisible": br, "xsec_pb": xsec_lhe,
        "xsec_over_BR_pb": ratio_here, "xsec_over_BR_rel_vs_v067": ratio_rel,
        "checkmate_returncode": cm_rc, "multibin_mu95": mu95,
        "multibin_r": rmb, "multibin_CLs_mu1": cls,
        "max_singlebin_robscons": rsingle, "detector_status": detector_status,
    }]).to_csv(OUT / "single_benchmark.csv", index=False)

    print("=== v0.68 single conservative CheckMATE benchmark ===")
    print(f"benchmark: Eq4.4 derivative / one_Weyl_pair, BRinv={br:.9f}")
    print(f"CT18 xsec: {xsec_lhe:.8g} pb")
    print(f"CT18 xsec/BR consistency vs v0.67: {ratio_rel:.3%}")
    print(f"CheckMATE rc={cm_rc}, multibin_r={rmb}, singlebin_robscons={rsingle}")
    print(f"detector status: {detector_status}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
