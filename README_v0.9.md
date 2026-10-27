# nutau-prop v0.9: public-data Case-6 audit

This version preserves the v0.8 baseline in Git history and adds an analytic
topology and flux-normalization audit. It folds the frozen Case-6 cross-section
grid with IceCube's released differential-limit windows. The result is a
**public-data screen, not a reproduction of IceCube's binned likelihood**.

## Reproduce

From the repository root, with `numpy`, `scipy`, and `pandas` installed:

```sh
./run_v0_9.sh
python -m pytest tests/test_public_ehe_v0_9.py
latexmk -pdf -interaction=nonstopmode -halt-on-error -cd manuscript/v0.9/main.tex
```

The v0.9 fold is one-dimensional and uses one CPU thread. The preserved
LHAPDF grid was the expensive part; rerunning it requires CT10nlo and may
benefit from at most 14 workers if a future parallel implementation is added.
Do not set 14 BLAS threads while also using 14 worker processes.

## Versioned outputs

- `results/v0.9/icecube_public_windows_v0_9.csv`: every published IceCube
  sliding-window true-energy screen.
- `results/v0.9/case6_public_audit_v0_9.csv`: benchmark summary, including
  all-flavor and tau-initial response proxies and the unchanged direct track.
- `manuscript/v0.9/main.tex`, `references.bib`, and `main.pdf`: JHEP-style
  draft by Subhendu Rakshit.
- `STATUS_AUDIT_v0.9.md` and `PROGRESS.log`: decisions, failures, and limits.

The official IceCube release is DOI `10.7910/DVN/JHK49D`; its public files do
not include reconstructed-bin observations, migration response, or background
templates. The exact collaboration likelihood therefore remains outstanding.
