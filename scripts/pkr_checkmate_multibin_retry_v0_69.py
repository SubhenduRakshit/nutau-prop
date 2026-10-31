#!/usr/bin/env python3
"""v0.69: retry only the CheckMATE multibin stage from the existing v0.68 HepMC sample.

No MG5 and no Pythia are run.  CheckMATE's spey_wrapper.py imports jax
unconditionally although its own documentation treats jax as optional.  On
this installation that makes the optional statistics wrapper fail to import.
We patch only that import to be optional and force Backend: numpy, then rerun
exactly the same CMS EXO-20-004 CheckMATE benchmark from the saved v0.68 HepMC.
"""
from __future__ import annotations

import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
V68 = ROOT / "results" / "v0.68_checkmate_single" / "single_benchmark.json"
WORK68 = ROOT / "cache" / "v0.68_checkmate_single"
WORK = ROOT / "cache" / "v0.69_multibin_retry"
OUT = ROOT / "results" / "v0.69_multibin_retry"
TAG = "eq4p4_derivative__one_Weyl_pair"

P51 = ROOT / "scripts" / "pkr_run_checkmate_v0_51.py"
spec = importlib.util.spec_from_file_location("cm51", P51)
if spec is None or spec.loader is None:
    raise RuntimeError(f"cannot import {P51}")
cm51 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cm51)


def patch_optional_jax() -> list[str]:
    roots = [ROOT / "cache" / "v0.51" / "checkmate-install", ROOT / "cache" / "v0.51" / "checkmate2"]
    patched = []
    for root in roots:
        if not root.exists():
            continue
        for p in root.rglob("spey_wrapper.py"):
            text = p.read_text(errors="replace")
            marker = "# v0.69: CheckMATE docs treat jax as optional"
            if marker in text:
                patched.append(str(p)); continue
            if re.search(r"(?m)^import jax\s*$", text):
                repl = marker + "\ntry:\n    import jax\nexcept ImportError:\n    jax = None"
                text = re.sub(r"(?m)^import jax\s*$", repl, text, count=1)
                p.write_text(text)
                patched.append(str(p))
    if not patched:
        raise RuntimeError("could not locate/patch CheckMATE spey_wrapper.py")
    return patched


def write_card(v68: dict) -> tuple[Path, Path]:
    src = WORK68 / "cards" / f"checkmate_{TAG}.in"
    hepmc = WORK68 / "hepmc" / f"{TAG}.hepmc"
    if not src.exists() or not hepmc.exists():
        raise RuntimeError("missing saved v0.68 CheckMATE card or HepMC sample")
    cmout = WORK / "checkmate" / TAG
    card = WORK / "cards" / f"checkmate_{TAG}.in"
    card.parent.mkdir(parents=True, exist_ok=True)
    text = src.read_text(errors="replace")
    text = re.sub(r"(?m)^Name:.*$", f"Name: pkr_v069_{TAG}", text, count=1)
    text = re.sub(r"(?m)^OutputDirectory:.*$", f"OutputDirectory: {cmout}", text, count=1)
    if re.search(r"(?im)^Backend\s*:", text):
        text = re.sub(r"(?im)^Backend\s*:.*$", "Backend: numpy", text, count=1)
    else:
        text = text.replace("[Parameters]\n", "[Parameters]\nBackend: numpy\n", 1)
    card.write_text(text)
    return card, cmout


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)
    if not V68.exists():
        raise RuntimeError(f"missing {V68}")
    v68 = json.loads(V68.read_text())

    # Verify the true required Spey stack before running CheckMATE.
    for mod in ("spey", "spey_pyhf", "pyhf"):
        __import__(mod)
    patched = patch_optional_jax()

    card, cmout = write_card(v68)
    if cmout.exists():
        shutil.rmtree(cmout)
    checkmate = Path(os.environ.get("PKR_CHECKMATE", ROOT / "cache" / "v0.51" / "checkmate-install" / "bin" / "CheckMATE"))
    if not checkmate.exists():
        raise RuntimeError(f"CheckMATE missing: {checkmate}")

    p = subprocess.run([str(checkmate), str(card)], text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (OUT / "checkmate.log").write_text(p.stdout)

    mu95, rmb, cls = cm51.parse_multibin(cmout)
    rsingle = cm51.parse_singlebin(cmout)
    if rmb is None:
        status = "MULTIBIN_STILL_UNRESOLVED"
    elif rmb > 1.2:
        status = "ROBUST_MULTIBIN_EXCLUSION"
    elif rmb < 0.8:
        status = "ROBUST_MULTIBIN_SURVIVAL"
    else:
        status = "MARGINAL_MULTIBIN"

    result = {
        "version": "v0.69_multibin_retry",
        "source": "existing v0.68 HepMC; no MG5/Pythia rerun",
        "backend": "numpy",
        "patched_spey_wrapper_paths": patched,
        "checkmate_returncode": p.returncode,
        "multibin_mu95": mu95,
        "multibin_r": rmb,
        "multibin_CLs_mu1": cls,
        "max_singlebin_robscons": rsingle,
        "status": status,
    }
    (OUT / "multibin_retry.json").write_text(json.dumps(result, indent=2) + "\n")
    print("=== v0.69 CheckMATE multibin retry ===")
    print("MG5 rerun: no; Pythia rerun: no")
    print(f"CheckMATE rc={p.returncode}")
    print(f"multibin_r={rmb}, multibin_mu95={mu95}, CLs(mu=1)={cls}")
    print(f"singlebin_robscons={rsingle}")
    print(f"status: {status}")
    return 0 if rmb is not None else 2


if __name__ == "__main__":
    sys.exit(main())
