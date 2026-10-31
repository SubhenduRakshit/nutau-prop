#!/usr/bin/env python3
"""v0.71: audit existing v0.69 CheckMATE multibin output only.

No MG5/Pythia/Delphes/CheckMATE rerun. Locate every multibin result-like file
under the existing v0.69 output tree and print the relevant Spey/covariance
lines from the saved CheckMATE log. This distinguishes a parser/path issue from
an actual missing likelihood result.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = "eq4p4_derivative__one_Weyl_pair"
CMROOT = ROOT / "cache" / "v0.69_multibin_retry" / "checkmate" / TAG
LOG = ROOT / "results" / "v0.69_multibin_retry" / "checkmate.log"


def main() -> int:
    print("=== v0.71 multibin output audit ===")
    print(f"cmroot: {CMROOT}")
    files = []
    if CMROOT.exists():
        for p in CMROOT.rglob("*"):
            if p.is_file() and ("multibin" in str(p).lower() or p.name == "results.dat"):
                files.append(p)
    print(f"multibin/result-like files: {len(files)}")
    for p in files:
        print(f"FILE {p}")
        if p.name == "results.dat":
            text = p.read_text(errors="replace")
            for line in text.splitlines():
                if any(k in line.lower() for k in ("analysis:", "upper limit", "cls", "likelihood")):
                    print("  " + line)

    if LOG.exists():
        print("--- relevant checkmate.log lines ---")
        for i, line in enumerate(LOG.read_text(errors="replace").splitlines(), 1):
            low = line.lower()
            if any(k in low for k in ("spey", "covariance", "upper limit", "multibin", "signal events below", "no signal events", "error is dominated")):
                print(f"{i}: {line}")
    else:
        print(f"missing log: {LOG}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
