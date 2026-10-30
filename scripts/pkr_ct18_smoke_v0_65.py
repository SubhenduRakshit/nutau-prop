#!/usr/bin/env python3
"""v0.65: one-sample CT18NLO/LHAPDF compatibility smoke test.

Purpose
-------
v0.64 certified the generator normalization/NWA BR scaling independently of
LHAPDF.  The remaining collider-environment blocker is CT18NLO (LHAID 14400):
LHAPDF 6.5.6 raises MetadataError when MG5's legacy LHAGLUE bridge requests
AlphaS_Lambda4/5, although MG5's own pdf_lhapdf62.cc explicitly supplies -1.0
as the intended fallback for missing metadata.

This test patches only those two legacy metadata getters to catch the exception
and return the same -1.0 fallback already encoded by MG5.  It then generates
ONE BR=0.45 sample with CT18NLO and checks parameter persistence plus
banner-vs-LHE normalization.  No CheckMATE run occurs.
"""
from __future__ import annotations

import importlib.util
import json
import re
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = ROOT / "scripts" / "pkr_certify_mg5_normalization_v0_62.py"
spec = importlib.util.spec_from_file_location("pkr_norm_v062", BASE_PATH)
if spec is None or spec.loader is None:
    raise RuntimeError(f"cannot import {BASE_PATH}")
v062 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v062)

base = v062.base          # deterministic v0.55 implementation
cfg = base.base           # v0.52 configuration/constants

WORK = ROOT / "cache" / "v0.65_ct18_smoke"
OUT = ROOT / "results" / "v0.65_ct18_smoke"
BR = 0.45


def patch_lhapdf_lambda_fallback(procdir: Path) -> None:
    src = procdir / "Source" / "PDF" / "pdf_lhapdf62.cc"
    if not src.exists():
        raise RuntimeError(f"expected MG5 LHAPDF bridge missing: {src}")
    text = src.read_text(errors="replace")

    old4 = 'qcdl4 = ACTIVESETS[nset].activeMember()->info().get_entry_as<double>("AlphaS_Lambda4", -1.0);'
    new4 = ('try { qcdl4 = ACTIVESETS[nset].activeMember()->info().get_entry_as<double>("AlphaS_Lambda4"); } '
            'catch (...) { qcdl4 = -1.0; }')
    old5 = 'qcdl5 = ACTIVESETS[nset].activeMember()->info().get_entry_as<double>("AlphaS_Lambda5", -1.0);'
    new5 = ('try { qcdl5 = ACTIVESETS[nset].activeMember()->info().get_entry_as<double>("AlphaS_Lambda5"); } '
            'catch (...) { qcdl5 = -1.0; }')

    if text.count(old4) != 1 or text.count(old5) != 1:
        raise RuntimeError(
            f"refusing blind LHAPDF patch: AlphaS_Lambda4 matches={text.count(old4)}, "
            f"AlphaS_Lambda5 matches={text.count(old5)}"
        )
    text = text.replace(old4, new4, 1).replace(old5, new5, 1)
    src.write_text(text)
    print(f"{procdir.name}: patched missing AlphaS_Lambda4/5 metadata to MG5's intended -1 fallback")


def main() -> int:
    if not cfg.MG5.exists():
        raise RuntimeError(f"MG5 not found at {cfg.MG5}")
    if not cfg.MODEL.exists():
        raise RuntimeError(f"model not found at {cfg.MODEL}")

    OUT.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)
    procdir = WORK / "mg5" / "br45"
    if procdir.exists():
        shutil.rmtree(procdir)
    (WORK / "cards").mkdir(parents=True, exist_ok=True)

    build_card = WORK / "cards" / "build_br45.txt"
    base.build_process(procdir, build_card, OUT / "build_br45.log")

    param = procdir / "Cards" / "param_card.dat"
    run_card = procdir / "Cards" / "run_card.dat"
    if not param.exists() or not run_card.exists():
        raise RuntimeError(f"missing generated cards below {procdir}/Cards")

    gamma_vis = cfg.visible_width_reference()
    gamma_inv = gamma_vis * BR / (1.0 - BR)
    gamma_total = gamma_vis + gamma_inv
    gchi = cfg.gchi_for_width(cfg.M, cfg.MCHI, gamma_inv)

    base.edit_param_card(param, cfg.M, cfg.MCHI, cfg.GQ, gchi, gamma_total)
    base.configure_run_card(run_card, seed=6501)  # inherited CT18NLO/LHAID 14400 setup
    ok_before, *_ = base.certify_persisted_param_card(
        param, cfg.M, cfg.MCHI, cfg.GQ, gchi, gamma_total
    )
    if not ok_before:
        raise RuntimeError("persisted parameter card failed before CT18 generation")

    patch_lhapdf_lambda_fallback(procdir)
    base.generate_events(procdir, OUT / "generate_br45.log")

    ok_after, actual_M, actual_mchi, actual_gchi, actual_width, actual_gq = base.certify_persisted_param_card(
        param, cfg.M, cfg.MCHI, cfg.GQ, gchi, gamma_total
    )
    lhe = cfg.find_lhe(procdir)
    banner = cfg.find_banner(procdir)
    xb = cfg.parse_banner_xsec_pb(banner)
    xi = cfg.parse_lhe_init_xsec_pb(lhe)
    rel = abs(xb - xi) / max(abs(xb), 1.0e-30)

    banner_text = banner.read_text(errors="replace")
    has_lhaid_14400 = bool(re.search(r"(?im)^\s*14400\s*=\s*lhaid\b", banner_text))
    has_lhapdf = "lhapdf" in banner_text.lower()
    event_hash = cfg.event_payload_hash(lhe)

    checks = {
        "param_card_matches_target_before": bool(ok_before),
        "param_card_matches_target_after": bool(ok_after),
        "banner_vs_lhe_init": bool(rel <= cfg.XSEC_MATCH_TOL),
        "banner_records_lhaid_14400": has_lhaid_14400,
        "banner_records_lhapdf": has_lhapdf,
        "event_payload_nonempty": bool(event_hash),
    }
    status = "PASS" if all(checks.values()) else "FAIL"
    result = {
        "version": "v0.65_ct18_smoke",
        "status": status,
        "purpose": "CT18NLO/LHAPDF compatibility smoke test only; no CheckMATE conclusion",
        "benchmark": {
            "BR_invisible": BR,
            "M_GeV": actual_M,
            "mchi_GeV": actual_mchi,
            "gchi": actual_gchi,
            "width_GeV": actual_width,
            "gq_min": min(actual_gq),
            "gq_max": max(actual_gq),
            "lhaid": 14400,
        },
        "checks": checks,
        "metrics": {
            "mg5_banner_xsec_pb": xb,
            "lhe_init_xsec_pb": xi,
            "banner_vs_init_rel_diff": rel,
            "event_payload_sha256": event_hash,
        },
        "patch_note": (
            "pdf_lhapdf62.cc catches missing AlphaS_Lambda4/5 metadata and returns -1.0, "
            "the fallback already specified by upstream MG5."
        ),
        "next_step": "If PASS, run exactly one CT18NLO CheckMATE benchmark; no scan.",
    }
    (OUT / "ct18_smoke.json").write_text(json.dumps(result, indent=2) + "\n")

    print("=== v0.65 CT18NLO smoke test ===")
    print(f"status: {status}")
    print(f"xsec: {xi:.8g} pb")
    print(f"banner-vs-LHE relative difference: {rel:.3e}")
    print(f"LHAID 14400 in banner: {has_lhaid_14400}")
    print(f"param card after generation: {ok_after}")
    return 0 if status == "PASS" else 2


if __name__ == "__main__":
    sys.exit(main())
