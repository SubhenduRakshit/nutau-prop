# nutau-prop v0.15.1 — audit correction before interpreting contact-operator transport

Date: 2026-10-27

## Status

The qualitative and quantitative conclusions recorded from the first v0.15 chirality-aware transport run are **provisional and must not be used in the manuscript or as a decision gate** until the corrected transport is rerun.

The exact-PDF Case-5/Case-7 hard-kernel prescreen in v0.14 remains valid.  The issue is in the new v0.15 transport layer and in how its output was interpreted.

## Audit finding 1: secondary muon neutrinos were accumulated but not propagated

In `TauHelicityFlipTransport`, tau-decay secondary `nu_mu`/`antinu_mu` flux was added to arrays but was not subsequently evolved through SM CC/NC interactions.  Consequently the saved `secondary_mu_ratio` columns are not reliable transport observables.

Correction: existing secondary muon-neutrino flux is now propagated through SM CC/NC on every substep before new tau-decay daughters are added.

## Audit finding 2: SM and BSM runs used different numerical substep schedules

The original v0.15 solver selected the number of substeps from the interaction rate actually enabled in that run.  Thus

- the matched SM run used a step size based on `sigma_SM`,
- the BSM run used a smaller step size based on `sigma_SM + sigma_NP`.

Because the algorithm is operator split and finite-step, comparing solutions at different numerical resolutions can generate spurious percent-level differences.  This is the likely origin of saved `Dirac sink` entries above unity.

For a pure additional sink, the discrete positive transport problem should not produce a component larger than the matched SM solution when both are evaluated at the same numerical resolution.  Therefore any robust Dirac-sink ratio greater than one is a hard consistency failure, not a physical signal.

Correction: SM and BSM runs now use the same substep schedule, determined conservatively from the full SM+NP rate even when the NP term is disabled in the baseline.

## Audit finding 3: opacity ratios alone are not an experimental conclusion

Several earlier statements interpreted very small BSM/SM ratios on long Earth chords as evidence that the contact operator is self-defeating.  This is insufficient by itself.

A ratio in a bin where the SM Earth-transmitted flux is already extremely small can have negligible detector weight.  Conversely, the scientifically relevant signal may live near the horizon or at lower column depth where both the absolute surviving flux and detector acceptance are substantial.

The corrected analysis must therefore show, at minimum:

1. matched SM and BSM propagated spectra;
2. transmission relative to the incident spectrum, not only BSM/SM ratios;
3. zenith dependence including shallow chords;
4. detector-weighted or response-weighted sensitivity before declaring an operator excluded or unpromising.

## What remains valid

The following v0.14 results are unaffected because they are hard-kernel/PDF calculations rather than transport conclusions:

- Case 5 has the same `y(1-y)` shape as Case 6 and the expected benchmark normalization ratio;
- at its nominal `c/(4 pi)=1` endpoint, Case 5 reaches `r_NP/SM ~ 0.055, 0.068, 0.073` at `100, 220, 300 PeV` with `<y> ~ 0.5`;
- at its endpoint, Case 7 reaches `r_NP/SM ~ 0.031, 0.054, 0.064` with `<y> ~ 0.8`;
- the chirality/external-state distinction (Dirac sink versus Majorana active transition/flip) remains physically necessary;
- the hard `sqrt(shat)<Lambda_UV` restriction is an EFT-validity boundary, not a physical turnover.

## What is withdrawn pending rerun

Pending v0.15.1 validation, do **not** quote as established:

- the numerical Case-5/Case-7 transport ratios from the first v0.15 run;
- the statement that Cases 5 and 7 fail the positive-signal gate;
- the combined statement that the tested contact-operator route is generically self-defeating through Earth opacity;
- secondary-muon conclusions from the original Case-7 helicity run.

The older Case-6 v0.12/v0.13 calculation has separate validated machinery, but the v0.15 chirality reinterpretation of Case 6 must also be rerun with the matched-step correction before it is used as a robustness argument.

## Required validation gate

Before any new full scan or v0.16 interpretation:

1. rerun only the transport stage using the existing v0.15 tables;
2. require `Dirac sink <= matched SM` component-by-component within numerical tolerance;
3. repeat representative points with smaller `max_p` and verify convergence;
4. propagate tau-decay secondary muon neutrinos through SM interactions;
5. report incident, SM-surviving and BSM-surviving spectral weights/transmissions;
6. only then perform detector weighting and decide whether opacity or detector enhancement wins.

The resolved-mediator v0.16 cross-section tables may be retained, but v0.16 transport/phenomenological interpretation is paused until this correction gate passes.
