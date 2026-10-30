#!/usr/bin/env python3
"""v0.52: certify MG5 normalization before any new CheckMATE run.

This is deliberately a generator-level gate only.  It makes three fresh,
independent pp -> Xd Xd~ j samples at fixed M and g_q while varying only the
invisible branching ratio.  It then requires:

  1. the persisted MG5 param_card.dat to contain the intended M, g_q, g_chi
     and total width;
  2. the generated event payloads to have distinct hashes (no sample reuse);
  3. MG5's persisted banner cross section to agree with the independent LHE
     <init> extraction;
  4. sigma/BR_inv to be constant to 2% across the three narrow-width points.

No detector simulation is run here.  A PASS marker is written only if every
check succeeds.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
IN45 = ROOT / "results" / "v0.45" / "width_br_dijet_gate.csv"
WORK = ROOT / "cache" / "v0.52_mg5_norm"
OUT = ROOT / "results" / "v0.52_mg5_norm"
PASS_MARKER = WORK / "NORMALIZATION_CERTIFIED.json"

MG5 = Path(os.environ.get("PKR_MG5", ROOT / "cache" / "v0.51" / "MG5_aMC" / "bin" / "mg5_aMC")).expanduser()
MODEL = Path(os.environ.get("PKR_CM_MODEL", ROOT / "cache" / "v0.51" / "models" / "DMsimp_s_spin1")).expanduser()
NEVENTS = int(os.environ.get("PKR_NORM_NEVENTS", "5000"))
M = float(os.environ.get("PKR_NORM_M_GEV", "50.0"))
GQ = float(os.environ.get("PKR_NORM_GQ", "0.03"))
MCHI = float(os.environ.get("PKR_CM_MCHI_GEV", "0.001"))
PTJMIN = float(os.environ.get("PKR_CM_PTJMIN_GEV", "100.0"))
PDF_LHAID = int(os.environ.get("PKR_CM_LHAID", "14400"))
BR_TARGETS = [0.25, 0.45, 0.60]
NWA_TOL = 0.02
PARAM_RTOL = 1.0e-7
XSEC_MATCH_TOL = 1.0e-6


def ffloat(s: str) -> float:
    return float(s.replace("D", "E").replace("d", "e"))


def gchi_for_width(mass: float, mchi: float, gamma_inv: float) -> float:
    r = (mchi / mass) ** 2
    phase = (1.0 + 2.0 * r) * math.sqrt(1.0 - 4.0 * r)
    return math.sqrt(gamma_inv * 12.0 * math.pi / (mass * phase))


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def read_bytes(path: Path) -> bytes:
    if path.suffix == ".gz":
        with gzip.open(path, "rb") as f:
            return f.read()
    return path.read_bytes()


def find_lhe(procdir: Path) -> Path:
    cand = list(procdir.glob("Events/run_*/unweighted_events.lhe.gz")) + list(procdir.glob("Events/run_*/unweighted_events.lhe"))
    if len(cand) != 1:
        raise RuntimeError(f"expected exactly one LHE below {procdir}, found {len(cand)}: {cand}")
    return cand[0]


def find_banner(procdir: Path) -> Path:
    cand = list(procdir.glob("Events/run_*/*_banner.txt"))
    if len(cand) != 1:
        raise RuntimeError(f"expected exactly one MG5 banner below {procdir}, found {len(cand)}: {cand}")
    return cand[0]


def parse_banner_xsec_pb(path: Path) -> float:
    text = path.read_text(errors="replace")
    m = re.search(r"#\s*Integrated weight \(pb\)\s*:\s*([0-9.eEdD+\-]+)", text)
    if not m:
        raise RuntimeError(f"could not find Integrated weight (pb) in {path}")
    return ffloat(m.group(1))


def parse_lhe_init_xsec_pb(path: Path) -> float:
    text = read_bytes(path).decode("utf-8", errors="replace")
    m = re.search(r"<init>\s*(.*?)\s*</init>", text, flags=re.S)
    if not m:
        raise RuntimeError(f"no <init> block in {path}")
    lines = [x.strip() for x in m.group(1).splitlines() if x.strip() and not x.lstrip().startswith("#")]
    if not lines:
        raise RuntimeError(f"empty <init> block in {path}")
    head = lines[0].split()
    if len(head) < 10:
        raise RuntimeError(f"malformed first <init> line in {path}: {lines[0]}")
    nprup = int(head[-1])
    rows = lines[1:1 + nprup]
    if len(rows) != nprup:
        raise RuntimeError(f"expected {nprup} subprocess rows in <init>, found {len(rows)}")
    return sum(ffloat(row.split()[0]) for row in rows)


def event_payload_hash(path: Path) -> str:
    data = read_bytes(path)
    marker = b"</init>"
    i = data.find(marker)
    if i < 0:
        raise RuntimeError(f"cannot isolate event payload in {path}")
    payload = data[i + len(marker):]
    return sha256_bytes(payload)


def parse_slha(path: Path):
    blocks: dict[str, dict[tuple[int, ...], float]] = {}
    decays: dict[int, float] = {}
    block = None
    for raw in path.read_text(errors="replace").splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        parts = line.split()
        key = parts[0].upper()
        if key == "BLOCK" and len(parts) >= 2:
            block = parts[1].upper()
            blocks.setdefault(block, {})
            continue
        if key == "DECAY" and len(parts) >= 3:
            decays[int(parts[1])] = ffloat(parts[2])
            block = None
            continue
        if block is None or len(parts) < 2:
            continue
        try:
            value = ffloat(parts[-1])
            idx = tuple(int(x) for x in parts[:-1])
        except ValueError:
            continue
        blocks[block][idx] = value
    return blocks, decays


def rel_close(actual: float, expected: float, rtol: float = PARAM_RTOL) -> bool:
    scale = max(abs(expected), 1.0e-15)
    return abs(actual - expected) / scale <= rtol


def make_card(procdir: Path, gamma_total: float, gchi: float) -> str:
    sets = [
        f"set MY1 {M:.16g}",
        f"set MXd {MCHI:.16g}",
        f"set WY1 {gamma_total:.16g}",
        f"set gVXd {gchi:.16g}",
        "set gAXd 0",
        "set gVXc 0",
    ]
    for p in ["gVd11", "gVu11", "gVd22", "gVu22", "gVd33", "gVu33"]:
        sets.append(f"set {p} {GQ:.16g}")
    for p in [
        "gAd11", "gAu11", "gAd22", "gAu22", "gAd33", "gAu33",
        "gVl11", "gVl22", "gVl33", "gAl11", "gAl22", "gAl33",
        "gnu11", "gnu22", "gnu33", "gVu31", "gAu31", "gVd31", "gAd31", "gVh",
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
        f"output {procdir} -f",
        f"launch {procdir} -f",
        *runsets,
        *sets,
        "done",
        "",
    ])


def run_mg5(card: Path, log: Path) -> None:
    p = subprocess.run([str(MG5), str(card)], text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    log.write_text(p.stdout)
    if p.returncode:
        tail = "\n".join(p.stdout.splitlines()[-50:])
        raise RuntimeError(f"MG5 failed for {card} (rc={p.returncode})\n{tail}")


def visible_width_reference() -> float:
    if not IN45.exists():
        raise RuntimeError(f"missing benchmark-width input {IN45}")
    d = pd.read_csv(IN45)
    q = d[(d["M_GeV"] - M).abs() < 1e-9 * max(M, 1.0)]
    q = q[(q["gq"] - GQ).abs() < 1e-9 * max(abs(GQ), 1.0)]
    if q.empty:
        raise RuntimeError(f"no v0.45 rows at M={M}, gq={GQ}")
    vis = (q["Gamma_total_GeV"] - q["Gamma_nu_GeV"]).astype(float)
    spread = (vis.max() - vis.min()) / vis.mean()
    if spread > 1.0e-8:
        raise RuntimeError(f"v0.45 visible widths are inconsistent: relative spread={spread:.3e}")
    return float(vis.mean())


def main() -> int:
    if not MG5.exists():
        raise RuntimeError(f"MG5 not found at {MG5}; run the v0.51 setup first or export PKR_MG5")
    if not MODEL.exists():
        raise RuntimeError(f"DMsimp model not found at {MODEL}; run the v0.51 setup first or export PKR_CM_MODEL")

    OUT.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)
    (WORK / "cards").mkdir(exist_ok=True)
    (WORK / "mg5").mkdir(exist_ok=True)
    PASS_MARKER.unlink(missing_ok=True)

    gamma_vis = visible_width_reference()
    rows = []

    for i, br in enumerate(BR_TARGETS, start=1):
        gamma_inv = gamma_vis * br / (1.0 - br)
        gamma_total = gamma_vis + gamma_inv
        gchi = gchi_for_width(M, MCHI, gamma_inv)
        tag = f"br{int(round(100*br)):02d}"
        procdir = WORK / "mg5" / tag
        card = WORK / "cards" / f"mg5_norm_{tag}.txt"
        log = OUT / f"mg5_norm_{tag}.log"

        # A certification sample must be fresh.  Never reuse a previous process tree.
        if procdir.exists():
            shutil.rmtree(procdir)
        card.write_text(make_card(procdir, gamma_total, gchi))
        run_mg5(card, log)

        param = procdir / "Cards" / "param_card.dat"
        if not param.exists():
            raise RuntimeError(f"persisted param_card.dat missing for {tag}")
        blocks, decays = parse_slha(param)

        actual_M = blocks.get("MASS", {}).get((5000001,))
        actual_mchi = blocks.get("MASS", {}).get((5000521,))
        actual_gchi = blocks.get("DMINPUTS", {}).get((2,))
        actual_gq = [blocks.get("DMINPUTS", {}).get((k,)) for k in range(4, 10)]
        actual_width = decays.get(5000001)
        required = [actual_M, actual_mchi, actual_gchi, actual_width, *actual_gq]
        if any(x is None for x in required):
            raise RuntimeError(f"missing required persisted parameters in {param}: {required}")

        param_ok = (
            rel_close(actual_M, M)
            and rel_close(actual_mchi, MCHI)
            and rel_close(actual_gchi, gchi)
            and rel_close(actual_width, gamma_total)
            and all(rel_close(x, GQ) for x in actual_gq)
        )

        lhe = find_lhe(procdir)
        banner = find_banner(procdir)
        x_banner = parse_banner_xsec_pb(banner)
        x_init = parse_lhe_init_xsec_pb(lhe)
        xsec_rel_diff = abs(x_banner - x_init) / max(abs(x_banner), 1.0e-30)
        event_hash = event_payload_hash(lhe)
        param_hash = sha256_bytes(param.read_bytes())

        rows.append({
            "tag": tag,
            "target_BR_invisible": br,
            "gamma_visible_GeV": gamma_vis,
            "target_gamma_invisible_GeV": gamma_inv,
            "target_gamma_total_GeV": gamma_total,
            "target_Gamma_over_M": gamma_total / M,
            "target_gchi": gchi,
            "actual_M_GeV": actual_M,
            "actual_mchi_GeV": actual_mchi,
            "actual_gchi": actual_gchi,
            "actual_width_GeV": actual_width,
            "actual_gq_min": min(actual_gq),
            "actual_gq_max": max(actual_gq),
            "param_card_ok": param_ok,
            "param_card_sha256": param_hash,
            "event_payload_sha256": event_hash,
            "mg5_banner_xsec_pb": x_banner,
            "lhe_init_xsec_pb": x_init,
            "banner_vs_init_rel_diff": xsec_rel_diff,
            "sigma_over_BR_pb": x_init / br,
            "lhe_path": str(lhe),
            "param_card_path": str(param),
        })
        print(f"{tag}: BR={br:.3f}  xsec(init)={x_init:.8g} pb  xsec/BR={x_init/br:.8g} pb  param_ok={param_ok}")

    R = pd.DataFrame(rows)
    ratios = R["sigma_over_BR_pb"].astype(float)
    ratio_mean = float(ratios.mean())
    ratio_max_dev = float(((ratios / ratio_mean) - 1.0).abs().max())
    banner_init_max = float(R["banner_vs_init_rel_diff"].max())
    params_ok = bool(R["param_card_ok"].all())
    event_hashes_unique = R["event_payload_sha256"].nunique() == len(R)
    param_hashes_unique = R["param_card_sha256"].nunique() == len(R)
    xsec_sources_ok = banner_init_max <= XSEC_MATCH_TOL
    nwa_scaling_ok = ratio_max_dev <= NWA_TOL
    monotonic_ok = all(b > a for a, b in zip(R["lhe_init_xsec_pb"].tolist(), R["lhe_init_xsec_pb"].tolist()[1:]))

    passed = all([
        params_ok,
        event_hashes_unique,
        param_hashes_unique,
        xsec_sources_ok,
        nwa_scaling_ok,
        monotonic_ok,
    ])

    summary = {
        "version": "v0.52_mg5_norm",
        "status": "PASS" if passed else "FAIL",
        "purpose": "generator-level normalization certification only; no CheckMATE conclusion",
        "fixed_parameters": {"M_GeV": M, "gq": GQ, "mchi_GeV": MCHI, "ptjmin_GeV": PTJMIN, "lhaid": PDF_LHAID},
        "BR_targets": BR_TARGETS,
        "nevents_per_point": NEVENTS,
        "checks": {
            "persisted_param_cards_match_targets": params_ok,
            "event_payload_hashes_unique": bool(event_hashes_unique),
            "param_card_hashes_unique": bool(param_hashes_unique),
            "mg5_banner_vs_lhe_init": xsec_sources_ok,
            "sigma_over_BR_constant_within_2pct": nwa_scaling_ok,
            "cross_section_monotonic_with_BR": monotonic_ok,
        },
        "metrics": {
            "max_banner_vs_init_relative_difference": banner_init_max,
            "sigma_over_BR_mean_pb": ratio_mean,
            "max_relative_deviation_sigma_over_BR": ratio_max_dev,
        },
        "tolerances": {"sigma_over_BR": NWA_TOL, "banner_vs_init": XSEC_MATCH_TOL, "param_card_relative": PARAM_RTOL},
        "next_step": "If PASS, run exactly one CheckMATE benchmark. If FAIL, stop and diagnose MG5 normalization; do not interpret collider r values.",
    }

    R.to_csv(OUT / "normalization_certification.csv", index=False)
    (OUT / "normalization_certification.json").write_text(json.dumps(summary, indent=2) + "\n")
    if passed:
        PASS_MARKER.write_text(json.dumps(summary, indent=2) + "\n")

    print("\n=== v0.52 MG5 normalization certification ===")
    print(R[["tag", "target_BR_invisible", "param_card_ok", "mg5_banner_xsec_pb", "lhe_init_xsec_pb", "sigma_over_BR_pb"]].to_string(index=False))
    print(json.dumps(summary, indent=2))
    return 0 if passed else 2


if __name__ == "__main__":
    sys.exit(main())
