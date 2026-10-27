# nutau-prop v0.18 — corrected production transport after finite-beta_tau repair

Date: 2026-10-28

## Preconditions

The finite-beta_tau transport bug found in the audit was repaired in v0.17 by evolving continuous tau energy loss and decay together in each constant-density segment and composing tau evolution with neutrino interactions using Strang-style splitting. The v0.17 validation passed all substantive analytic, positivity, max-p and matter-step convergence gates. A single auxiliary monotonic-refinement heuristic failed only because both successive errors were already below 1e-6; v0.17.1 therefore classifies points with both successive differences below 1e-4 as already converged.

The v0.14 hard-kernel/PDF results and the v0.16 resolved-mediator cross-section tables remain valid. Earlier finite-beta_tau transport outputs are retained as historical/provisional and are not used as physics decision gates.

## Ordered production rerun

1. Re-run Case 5 with the repaired transport for both the Dirac-sink and Majorana tau-mu brackets.
2. Re-run Case 6 for both the Dirac-sink and Majorana tau-antitau flip brackets.
3. Re-run Case 7 for both the Dirac-sink and Majorana helicity-flip brackets.
4. Re-run only the v0.16 resolved-mediator transport, reusing its saved cross-section tables.
5. Only after 1-4 are inspected, perform detector-response weighting.

## Numerical settings

Production uses `max_p = 0.02` and `max_step_km = 12`. The repaired v0.17 solver showed typical `0.02 -> 0.01` transport differences of `~1e-6--2e-5` and matter-step differences below `~6e-4` in the representative Case-6 tests, so these settings are comfortably inside the numerical accuracy required for the phenomenology.

Independent parameter points are parallelized with `NUTAU_NPROC` (default 12).

## Zenith coverage

The corrected rerun explicitly includes shallow and long Earth chords:

`|cos z| = 0.05, 0.1, 0.2, 0.5, 0.8, 1.0`.

This prevents the interpretation from being dominated by vertical bins in which the SM flux is already effectively zero.

## Flux normalization

The official Ehlert cosmogenic input in `data/uhe_official/flux.json` is treated as the all-flavor `nu+antinu` diffuse flux used with the public all-flavor KM3NeT response. At Earth it is split equally among the six active species. For each species,

`dPhi_species/dlnE = (1/6) * 10^[log10(E^2 F_all)] / E`.

The full CSV outputs therefore retain physical units `cm^-2 s^-1 sr^-1 per dlnE`, rather than the unit-normalized spectrum used in the earlier transport-screen scripts.

## Required outputs

For every operator/mass, zenith and external-state bracket v0.18 records:

- all-flavor incident spectral weight;
- relevant incident tau or tau+mu spectral weight;
- matched SM transmitted spectral weight;
- BSM transmitted spectral weight;
- transmission relative to the incident spectrum;
- raw BSM/SM ratios;
- interpretation ratios masked whenever the matched-SM transmission is below `1e-12`.

The raw ratio is preserved for reproducibility. The masked ratio is the quantity intended for physics inspection; this prevents enormous ratios in bins with essentially zero absolute surviving flux from being mistaken for experimental sensitivity.

## Scope of the v0.18 decision

v0.18 is a propagation result only. It may identify promising or irrelevant energy/zenith regions, but it does **not** by itself establish exclusion, discovery reach, or a no-go statement. A detector-response fold over energy and zenith is required after the corrected transport is inspected.
