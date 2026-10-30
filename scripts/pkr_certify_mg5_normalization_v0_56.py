#!/usr/bin/env python3
"""v0.56 deterministic MG5 normalization certification.

v0.55 successfully built a fresh process directory and reached the persisted
run card, but its editor assumed the historical MadGraph syntax

    value = variable

for every setting.  At least one generated-card setting (`pdlabel`) was not
written in that orientation, so v0.55 stopped before event generation.

v0.56 preserves the v0.55 deterministic workflow and changes only the run-card
editor: it recognizes either `value = variable` or `variable = value`, keeps
whatever orientation is already present, and still refuses to invent a missing
key.  No detector simulation is run here.
"""
from __future__ import annotations

import importlib.util
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = ROOT / "scripts" / "pkr_certify_mg5_normalization_v0_55.py"
spec = importlib.util.spec_from_file_location("pkr_norm_v055", BASE_PATH)
if spec is None or spec.loader is None:
    raise RuntimeError(f"cannot import {BASE_PATH}")
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)

base.WORK = ROOT / "cache" / "v0.56_mg5_norm"
base.OUT = ROOT / "results" / "v0.56_mg5_norm"
base.PASS_MARKER = base.WORK / "NORMALIZATION_CERTIFIED.json"


def set_run_card_value_v056(path: Path, key: str, value: str) -> None:
    lines = path.read_text(errors="replace").splitlines(keepends=True)
    key_re = re.compile(rf"\b{re.escape(key)}\b", re.I)
    found = False
    out: list[str] = []

    for raw in lines:
        if raw.lstrip().startswith("#") or "=" not in raw:
            out.append(raw)
            continue

        body = raw.rstrip("\n")
        nl = "\n" if raw.endswith("\n") else ""
        before_comment, sep_comment, comment = body.partition("!")
        lhs, sep, rhs = before_comment.partition("=")
        if not sep:
            out.append(raw)
            continue

        lhs_s = lhs.strip()
        rhs_s = rhs.strip()
        lhs_has = bool(key_re.fullmatch(lhs_s))
        rhs_token = rhs_s.split()[0] if rhs_s else ""
        rhs_has = bool(key_re.fullmatch(rhs_token))

        if lhs_has:
            new_core = f"{key} = {value}"
            suffix = f" !{comment}" if sep_comment else ""
            out.append(new_core + suffix + nl)
            found = True
        elif rhs_has:
            # Preserve canonical MadGraph value = variable orientation.
            trailing = rhs_s[len(rhs_token):]
            new_core = f" {value} = {rhs_token}{trailing}"
            suffix = f" !{comment}" if sep_comment else ""
            out.append(new_core + suffix + nl)
            found = True
        else:
            out.append(raw)

    if not found:
        nearby = [ln.strip() for ln in lines if key.lower() in ln.lower()][:5]
        raise RuntimeError(f"run-card key {key!r} not found in {path}; nearby={nearby}")
    path.write_text("".join(out))


base.set_run_card_value = set_run_card_value_v056


if __name__ == "__main__":
    sys.exit(base.main())
