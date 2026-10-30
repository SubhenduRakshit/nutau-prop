#!/usr/bin/env python3
"""v0.66: corrected one-sample CT18NLO/LHAPDF compatibility smoke test.

v0.65 established two remaining implementation issues:

1. Safari still opened because process-local MadEvent reads
   Cards/me5_configuration.txt; setting automatic_html_opening only in the
   top-level MG5 build command is not sufficient.
2. The CT18NLO crash persisted because pdf_lhapdf62.cc contains multiple
   legacy AlphaS_Lambda4/5 metadata reads.  v0.65 patched only getlam4/5,
   while the earlier PDF-initialisation reads could still throw.

v0.66 changes only those compatibility details.  It patches every expected
AlphaS_Lambda4/5 read in the freshly generated process source to preserve the
fallback already encoded by MG5, removes stale libpdf build products so the
patched source is necessarily recompiled, and disables process-local automatic
HTML opening.  The physics benchmark and CT18NLO LHAID 14400 are unchanged.
No CheckMATE run occurs.
"""
from __future__ import annotations

import importlib.util
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = ROOT / "scripts" / "pkr_ct18_smoke_v0_65.py"
spec = importlib.util.spec_from_file_location("pkr_ct18_v065", BASE_PATH)
if spec is None or spec.loader is None:
    raise RuntimeError(f"cannot import {BASE_PATH}")
v065 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v065)

v065.WORK = ROOT / "cache" / "v0.66_ct18_smoke"
v065.OUT = ROOT / "results" / "v0.66_ct18_smoke"


def disable_process_html_opening(procdir: Path) -> None:
    cfg = procdir / "Cards" / "me5_configuration.txt"
    if not cfg.exists():
        raise RuntimeError(f"process-local MG5 configuration missing: {cfg}")
    text = cfg.read_text(errors="replace")
    pat = re.compile(r"(?m)^\s*#?\s*automatic_html_opening\s*=.*$")
    if pat.search(text):
        text = pat.sub("automatic_html_opening = False", text, count=1)
    else:
        text += "\nautomatic_html_opening = False\n"
    npat = re.compile(r"(?m)^\s*#?\s*notification_center\s*=.*$")
    if npat.search(text):
        text = npat.sub("notification_center = False", text, count=1)
    else:
        text += "notification_center = False\n"
    cfg.write_text(text)
    print(f"{procdir.name}: disabled process-local HTML opening and notifications")


def patch_lhapdf_lambda_fallback_v066(procdir: Path) -> None:
    src = procdir / "Source" / "PDF" / "pdf_lhapdf62.cc"
    if not src.exists():
        raise RuntimeError(f"expected MG5 LHAPDF bridge missing: {src}")
    text = src.read_text(errors="replace")

    replacements = [
        (
            'w50512_.qcdl4 = pdf->info().get_entry_as<double>("AlphaS_Lambda4", 0.0);',
            'try { w50512_.qcdl4 = pdf->info().get_entry_as<double>("AlphaS_Lambda4"); } catch (...) { w50512_.qcdl4 = 0.0; }',
        ),
        (
            'w50512_.qcdl5 = pdf->info().get_entry_as<double>("AlphaS_Lambda5", 0.0);',
            'try { w50512_.qcdl5 = pdf->info().get_entry_as<double>("AlphaS_Lambda5"); } catch (...) { w50512_.qcdl5 = 0.0; }',
        ),
        (
            'lhapdfr_.qcdlha4 = pdf->info().get_entry_as<double>("AlphaS_Lambda4", 0.0);',
            'try { lhapdfr_.qcdlha4 = pdf->info().get_entry_as<double>("AlphaS_Lambda4"); } catch (...) { lhapdfr_.qcdlha4 = 0.0; }',
        ),
        (
            'lhapdfr_.qcdlha5 = pdf->info().get_entry_as<double>("AlphaS_Lambda5", 0.0);',
            'try { lhapdfr_.qcdlha5 = pdf->info().get_entry_as<double>("AlphaS_Lambda5"); } catch (...) { lhapdfr_.qcdlha5 = 0.0; }',
        ),
        (
            'qcdl4 = ACTIVESETS[nset].activeMember()->info().get_entry_as<double>("AlphaS_Lambda4", -1.0);',
            'try { qcdl4 = ACTIVESETS[nset].activeMember()->info().get_entry_as<double>("AlphaS_Lambda4"); } catch (...) { qcdl4 = -1.0; }',
        ),
        (
            'qcdl5 = ACTIVESETS[nset].activeMember()->info().get_entry_as<double>("AlphaS_Lambda5", -1.0);',
            'try { qcdl5 = ACTIVESETS[nset].activeMember()->info().get_entry_as<double>("AlphaS_Lambda5"); } catch (...) { qcdl5 = -1.0; }',
        ),
    ]

    for old, new in replacements:
        n = text.count(old)
        if n != 1:
            raise RuntimeError(f"refusing blind LHAPDF patch: expected exactly one occurrence, found {n}: {old}")
        text = text.replace(old, new, 1)

    src.write_text(text)

    # Force the patched bridge to be recompiled even if MG5 copied stale build products.
    for stale in [
        procdir / "Source" / "PDF" / "pdf_lhapdf62.o",
        procdir / "lib" / "libpdf.a",
        procdir / "lib" / "libpdf.dylib",
    ]:
        stale.unlink(missing_ok=True)

    print(f"{procdir.name}: patched all six legacy AlphaS_Lambda4/5 metadata reads and forced libpdf rebuild")


v065.patch_lhapdf_lambda_fallback = patch_lhapdf_lambda_fallback_v066

# Wrap process generation so process-local MadEvent cannot open Safari.
previous_generate_events = v065.base.generate_events


def generate_events_v066(procdir: Path, log: Path) -> None:
    disable_process_html_opening(procdir)
    previous_generate_events(procdir, log)


v065.base.generate_events = generate_events_v066


if __name__ == "__main__":
    sys.exit(v065.main())
