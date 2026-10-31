# v0.86 status — 31 October 2026

- **Resolved:** exact decay/DIS normalization for explicitly specified
  distinct-state, massless chiral Eq. 4.3/4.4 vertices. Relative rate is four;
  legacy absolute DIS rates were too large by eight/sixteen respectively.
- **Calculated on Mac Studio:** fresh CT18 folds with 14 worker processes;
  old 100-PeV 5% points become approximately 0.62%/0.31%.
- **Validated:** general-angle symbolic traces, Ward/flavour identities,
  decay crossing, 20 tests, quadrature and envelope interpolation checks.
- **Collider scope:** conditional NWA screens only. No nominal 5% survivor
  in the tested mass/coefficient domains with the stated inputs. This is
  not a robust exclusion; relaxed M50 Eq43 screening can exceed 5%.
- **Not resolved:** original PKR UFO/field mapping, exact mass-specific
  dijet recast, full UV validity, tau-specific IceCube, updated final
  propagation/detector likelihood. No claim that PKR itself is disproved.
- **Deliverables:** `PKR_NORMALIZATION_RESOLUTION_v0.86.md`,
  `results/v0.86_normalization/`, `manuscript/v0.86/normalization.tex` and PDF.
  Historical failures and unrelated files preserved. Runner has no automatic
  Git side effects.
