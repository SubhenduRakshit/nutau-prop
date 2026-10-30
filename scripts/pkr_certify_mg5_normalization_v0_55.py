#!/usr/bin/env python3
"""v0.55 deterministic MG5 normalization certification.

v0.52 proved that `launch ... -f` ran with UFO defaults.
v0.53-v0.54 showed that scripting MG5's interactive launch state machine is
not reliable enough for a certification gate.

v0.55 removes that ambiguity completely:
  1. MG5 is used only to build a fresh process directory.
  2. Python edits the persisted Cards/param_card.dat and Cards/run_card.dat.
  3. The edited cards are read back and certified before event generation.
  4. The process-local bin/generate_events executable is then run with -f.
  5. Banner/LHE normalization, unique-event hashes, and sigma/BR scaling are
     checked exactly as before.

No CheckMATE run is performed here.
"""
from __future__ import annotations

import importlib.util
import json
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = ROOT / "scripts" / "pkr_certify_mg5_normalization_v0_52.py"
spec = importlib.util.spec_from_file_location("pkr_norm_v052", BASE_PATH)
if spec is None or spec.loader is None:
    raise RuntimeError(f"cannot import {BASE_PATH}")
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)

WORK = ROOT / "cache" / "v0.55_mg5_norm"
OUT = ROOT / "results" / "v0.55_mg5_norm"
PASS_MARKER = WORK / "NORMALIZATION_CERTIFIED.json"


def run(cmd, log: Path, cwd: Path | None = None) -> None:
    p = subprocess.run([str(x) for x in cmd], cwd=cwd, text=True,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    log.parent.mkdir(parents=True, exist_ok=True)
    log.write_text(p.stdout)
    if p.returncode:
        tail = "\n".join(p.stdout.splitlines()[-60:])
        raise RuntimeError(f"command failed ({p.returncode}): {' '.join(map(str, cmd))}\n{tail}")


def build_process(procdir: Path, card: Path, log: Path) -> None:
    card.write_text("\n".join([
        "set automatic_html_opening False",
        f"import model {base.MODEL}",
        "generate p p > xd xd~ j",
        f"output {procdir} -f",
        "quit",
        "",
    ]))
    run([base.MG5, card], log)


def _replace_numeric_line(raw: str, value: float) -> str:
    comment = ""
    if "#" in raw:
        _, c = raw.split("#", 1)
        comment = " #" + c.rstrip("\n")
    return f"    {value:.16e}{comment}\n"


def edit_param_card(path: Path, M: float, mchi: float, gq: float,
                    gchi: float, width: float) -> None:
    dm = {1: 0.0, 2: gchi, 3: 0.0}
    for k in range(4, 10):
        dm[k] = gq
    for k in range(10, 30):
        dm[k] = 0.0
    mass = {5000001: M, 5000521: mchi}

    out = []
    block = None
    seen = {"M": False, "mchi": False, "gchi": False, "width": False,
            "gq": set()}
    for raw in path.read_text(errors="replace").splitlines(keepends=True):
        clean = raw.split("#", 1)[0].strip()
        parts = clean.split()
        if parts and parts[0].upper() == "BLOCK" and len(parts) >= 2:
            block = parts[1].upper()
            out.append(raw)
            continue
        if parts and parts[0].upper() == "DECAY" and len(parts) >= 3:
            pid = int(parts[1])
            if pid == 5000001:
                comment = ""
                if "#" in raw:
                    comment = " #" + raw.split("#", 1)[1].rstrip("\n")
                out.append(f"DECAY  5000001  {width:.16e}{comment}\n")
                seen["width"] = True
            else:
                out.append(raw)
            block = None
            continue
        if block == "MASS" and len(parts) >= 2:
            try:
                idx = int(parts[0])
            except ValueError:
                out.append(raw); continue
            if idx in mass:
                comment = ""
                if "#" in raw:
                    comment = " #" + raw.split("#", 1)[1].rstrip("\n")
                out.append(f"  {idx:9d}  {mass[idx]:.16e}{comment}\n")
                seen["M" if idx == 5000001 else "mchi"] = True
                continue
        if block == "DMINPUTS" and len(parts) >= 2:
            try:
                idx = int(parts[0])
            except ValueError:
                out.append(raw); continue
            if idx in dm:
                comment = ""
                if "#" in raw:
                    comment = " #" + raw.split("#", 1)[1].rstrip("\n")
                out.append(f"  {idx:4d}  {dm[idx]:.16e}{comment}\n")
                if idx == 2:
                    seen["gchi"] = True
                if 4 <= idx <= 9:
                    seen["gq"].add(idx)
                continue
        out.append(raw)

    if not (seen["M"] and seen["mchi"] and seen["gchi"] and seen["width"] and len(seen["gq"]) == 6):
        raise RuntimeError(f"failed to locate required param-card entries in {path}: {seen}")
    path.write_text("".join(out))


def set_run_card_value(path: Path, key: str, value: str) -> None:
    lines = path.read_text(errors="replace").splitlines(keepends=True)
    pat = re.compile(rf"^(?P<lhs>.*?)=(?P<rhs>\s*{re.escape(key)}\b.*)$", re.I)
    found = False
    out = []
    for raw in lines:
        m = pat.match(raw)
        if m and not raw.lstrip().startswith("#"):
            suffix = m.group("rhs")
            nl = "\n" if raw.endswith("\n") else ""
            suffix = suffix.rstrip("\n")
            out.append(f" {value} ={suffix}{nl}")
            found = True
        else:
            out.append(raw)
    if not found:
        raise RuntimeError(f"run-card key {key!r} not found in {path}")
    path.write_text("".join(out))


def configure_run_card(path: Path, seed: int) -> None:
    vals = {
        "nevents": str(base.NEVENTS),
        "iseed": str(seed),
        "ebeam1": "6500",
        "ebeam2": "6500",
        "pdlabel": "'lhapdf'",
        "lhaid": str(base.PDF_LHAID),
        "ptj": f"{base.PTJMIN:.16g}",
        "etaj": "5.0",
        "maxjetflavor": "5",
        "ickkw": "0",
    }
    for k, v in vals.items():
        set_run_card_value(path, k, v)


def generate_events(procdir: Path, log: Path) -> None:
    gen = procdir / "bin" / "generate_events"
    madevent = procdir / "bin" / "madevent"
    if gen.exists():
        run([gen, "-f"], log, cwd=procdir)
    elif madevent.exists():
        run([madevent, "generate_events", "-f"], log, cwd=procdir)
    else:
        raise RuntimeError(f"no process-local event generator found below {procdir}/bin")


def certify_persisted_param_card(param: Path, M: float, mchi: float, gq: float,
                                  gchi: float, width: float):
    blocks, decays = base.parse_slha(param)
    actual_M = blocks.get("MASS", {}).get((5000001,))
    actual_mchi = blocks.get("MASS", {}).get((5000521,))
    actual_gchi = blocks.get("DMINPUTS", {}).get((2,))
    actual_gq = [blocks.get("DMINPUTS", {}).get((k,)) for k in range(4, 10)]
    actual_width = decays.get(5000001)
    required = [actual_M, actual_mchi, actual_gchi, actual_width, *actual_gq]
    if any(x is None for x in required):
        raise RuntimeError(f"missing required persisted parameters in {param}: {required}")
    ok = (base.rel_close(actual_M, M) and base.rel_close(actual_mchi, mchi)
          and base.rel_close(actual_gchi, gchi) and base.rel_close(actual_width, width)
          and all(base.rel_close(x, gq) for x in actual_gq))
    return ok, actual_M, actual_mchi, actual_gchi, actual_width, actual_gq


def main() -> int:
    if not base.MG5.exists():
        raise RuntimeError(f"MG5 not found at {base.MG5}")
    if not base.MODEL.exists():
        raise RuntimeError(f"model not found at {base.MODEL}")

    OUT.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)
    (WORK / "cards").mkdir(exist_ok=True)
    (WORK / "mg5").mkdir(exist_ok=True)
    PASS_MARKER.unlink(missing_ok=True)

    gamma_vis = base.visible_width_reference()
    rows = []
    for i, br in enumerate(base.BR_TARGETS, start=1):
        gamma_inv = gamma_vis * br / (1.0 - br)
        gamma_total = gamma_vis + gamma_inv
        gchi = base.gchi_for_width(base.M, base.MCHI, gamma_inv)
        tag = f"br{int(round(100*br)):02d}"
        procdir = WORK / "mg5" / tag
        if procdir.exists():
            shutil.rmtree(procdir)
        build_card = WORK / "cards" / f"build_{tag}.txt"
        build_process(procdir, build_card, OUT / f"build_{tag}.log")

        param = procdir / "Cards" / "param_card.dat"
        run_card = procdir / "Cards" / "run_card.dat"
        if not param.exists() or not run_card.exists():
            raise RuntimeError(f"missing generated cards below {procdir}/Cards")

        edit_param_card(param, base.M, base.MCHI, base.GQ, gchi, gamma_total)
        configure_run_card(run_card, seed=5500 + i)
        param_ok, actual_M, actual_mchi, actual_gchi, actual_width, actual_gq = certify_persisted_param_card(
            param, base.M, base.MCHI, base.GQ, gchi, gamma_total)
        if not param_ok:
            raise RuntimeError(f"persisted param card failed readback before generation for {tag}")
        param_hash_before = base.sha256_bytes(param.read_bytes())
        run_hash_before = base.sha256_bytes(run_card.read_bytes())

        generate_events(procdir, OUT / f"generate_{tag}.log")

        # Re-read after generation: the actual generated sample must still be tied to these cards.
        param_ok_after, *_ = certify_persisted_param_card(param, base.M, base.MCHI, base.GQ, gchi, gamma_total)
        param_hash_after = base.sha256_bytes(param.read_bytes())
        run_hash_after = base.sha256_bytes(run_card.read_bytes())
        cards_stable = (param_hash_before == param_hash_after and run_hash_before == run_hash_after)

        lhe = base.find_lhe(procdir)
        banner = base.find_banner(procdir)
        x_banner = base.parse_banner_xsec_pb(banner)
        x_init = base.parse_lhe_init_xsec_pb(lhe)
        xsec_rel_diff = abs(x_banner - x_init) / max(abs(x_banner), 1.0e-30)
        event_hash = base.event_payload_hash(lhe)

        rows.append({
            "tag": tag,
            "target_BR_invisible": br,
            "gamma_visible_GeV": gamma_vis,
            "target_gamma_invisible_GeV": gamma_inv,
            "target_gamma_total_GeV": gamma_total,
            "target_Gamma_over_M": gamma_total / base.M,
            "target_gchi": gchi,
            "actual_M_GeV": actual_M,
            "actual_mchi_GeV": actual_mchi,
            "actual_gchi": actual_gchi,
            "actual_width_GeV": actual_width,
            "actual_gq_min": min(actual_gq),
            "actual_gq_max": max(actual_gq),
            "param_card_ok_before_generation": param_ok,
            "param_card_ok_after_generation": param_ok_after,
            "cards_stable_through_generation": cards_stable,
            "param_card_sha256": param_hash_after,
            "run_card_sha256": run_hash_after,
            "event_payload_sha256": event_hash,
            "mg5_banner_xsec_pb": x_banner,
            "lhe_init_xsec_pb": x_init,
            "banner_vs_init_rel_diff": xsec_rel_diff,
            "sigma_over_BR_pb": x_init / br,
        })
        print(f"{tag}: BR={br:.3f} xsec={x_init:.8g} pb xsec/BR={x_init/br:.8g} pb param_ok={param_ok_after} cards_stable={cards_stable}")

    R = pd.DataFrame(rows)
    ratios = R["sigma_over_BR_pb"].astype(float)
    ratio_mean = float(ratios.mean())
    ratio_max_dev = float(((ratios / ratio_mean) - 1.0).abs().max())
    banner_init_max = float(R["banner_vs_init_rel_diff"].max())
    params_ok = bool(R["param_card_ok_before_generation"].all() and R["param_card_ok_after_generation"].all())
    cards_stable = bool(R["cards_stable_through_generation"].all())
    event_hashes_unique = R["event_payload_sha256"].nunique() == len(R)
    param_hashes_unique = R["param_card_sha256"].nunique() == len(R)
    run_hashes_unique = R["run_card_sha256"].nunique() == len(R)
    xsec_sources_ok = banner_init_max <= base.XSEC_MATCH_TOL
    nwa_scaling_ok = ratio_max_dev <= base.NWA_TOL
    xs = R["lhe_init_xsec_pb"].astype(float).tolist()
    monotonic_ok = all(b > a for a, b in zip(xs, xs[1:]))

    passed = all([params_ok, cards_stable, event_hashes_unique, param_hashes_unique,
                  run_hashes_unique, xsec_sources_ok, nwa_scaling_ok, monotonic_ok])
    summary = {
        "version": "v0.55_mg5_norm",
        "status": "PASS" if passed else "FAIL",
        "purpose": "deterministic generator-level normalization certification only; no CheckMATE conclusion",
        "execution_mode": "build process; edit persisted cards directly; verify; process-local generate_events -f",
        "checks": {
            "persisted_param_cards_match_targets": params_ok,
            "cards_stable_through_generation": cards_stable,
            "event_payload_hashes_unique": bool(event_hashes_unique),
            "param_card_hashes_unique": bool(param_hashes_unique),
            "run_card_hashes_unique": bool(run_hashes_unique),
            "mg5_banner_vs_lhe_init": xsec_sources_ok,
            "sigma_over_BR_constant_within_2pct": nwa_scaling_ok,
            "cross_section_monotonic_with_BR": monotonic_ok,
        },
        "metrics": {
            "max_banner_vs_init_relative_difference": banner_init_max,
            "sigma_over_BR_mean_pb": ratio_mean,
            "max_relative_deviation_sigma_over_BR": ratio_max_dev,
        },
        "tolerances": {
            "sigma_over_BR": base.NWA_TOL,
            "banner_vs_init": base.XSEC_MATCH_TOL,
            "param_card_relative": base.PARAM_RTOL,
        },
        "history": {
            "v0.52": "failed: launch -f ran with UFO defaults",
            "v0.53": "failed: interactive launch did not produce LHE",
            "v0.54": "failed: interactive launch state sequencing still did not produce LHE",
        },
        "next_step": "If PASS, run exactly one CheckMATE benchmark. If FAIL, stop and diagnose generator normalization.",
    }
    R.to_csv(OUT / "normalization_certification.csv", index=False)
    (OUT / "normalization_certification.json").write_text(json.dumps(summary, indent=2) + "\n")
    if passed:
        PASS_MARKER.write_text(json.dumps(summary, indent=2) + "\n")

    print("\n=== v0.55 MG5 normalization certification ===")
    print(R[["tag", "target_BR_invisible", "actual_M_GeV", "actual_gchi", "actual_width_GeV",
             "lhe_init_xsec_pb", "sigma_over_BR_pb"]].to_string(index=False))
    print(json.dumps(summary, indent=2))
    return 0 if passed else 2


if __name__ == "__main__":
    sys.exit(main())
