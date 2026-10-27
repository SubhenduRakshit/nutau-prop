# nutau-prop v0.15.1 — audit correction before interpreting contact-operator transport

Date: 2026-10-27

> **v0.17 addendum.** The v0.15.1 matched-step/secondary-muon corrections were necessary but not sufficient. The v0.15.1 convergence gate subsequently failed systematically and exposed a deeper finite-`beta_tau` operator-splitting problem: the old solver applied tau energy loss over an entire neutrino substep before allowing the tau to decay. Both historical Rakshit-Reya and Dutta-Reno-Sarcevic validators explicitly used `beta_tau=0.0`, so they did not validate this finite-loss evolution. Therefore the earlier statement below that v0.12/v0.13 Case-6 finite-loss machinery was independently protected by the historical validation is withdrawn. See `FINITE_BETA_TAU_REPAIR_v0.17.md` for the repaired solver and new validation gate.

## Status

The qualitative and quantitative conclusions recorded from the first v0.15 chirality-aware transport run are **provisional and must not be used in the manuscript or as a decision gate** until the corrected transport is rerun.

The exact-PDF Case-5/Case-7 hard-kernel prescreen in v0.14 remains valid. The issue is in the transport layer and in how its output was interpreted.

## Audit finding 1: secondary muon neutrinos were accumulated but not propagated

In `TauHelicityFlipTransport`, tau-decay secondary `nu_mu`/`antinu_mu` flux was added to arrays but was not subsequently evolved through SM CC/NC interactions. Consequently the saved `secondary_mu_ratio` columns are not reliable transport observables.

Correction: existing secondary muon-neutrino flux is propagated through SM CC/NC on every substep before new tau-decay daughters are added.

## Audit finding 2: SM and BSM runs used different numerical substep schedules

The original v0.15 solver selected the number of substeps from the interaction rate actually enabled in that run. Thus the matched SM run used a step size based on `sigma_SM`, while the BSM run used a smaller step size based on `sigma_SM + sigma_NP`.

Because the algorithm is operator split and finite-step, comparing solutions at different numerical resolutions can generate spurious percent-level differences. This is the likely origin of saved `Dirac sink` entries above unity.

Correction: SM and BSM runs use the same substep schedule, determined conservatively from the full SM+NP rate even when the NP term is disabled in the baseline.

## Audit finding 3: opacity ratios alone are not an experimental conclusion

Several earlier statements interpreted very small BSM/SM ratios on long Earth chords as evidence that the contact operator is self-defeating. This is insufficient by itself.

A ratio in a bin where the SM Earth-transmitted flux is already extremely small can have negligible detector weight. Conversely, the scientifically relevant signal may live near the horizon or at lower column depth where both the absolute surviving flux and detector acceptance are substantial.

The corrected analysis must therefore show, at minimum:

1. matched SM and BSM propagated spectra;
2. transmission relative to the incident spectrum, not only BSM/SM ratios;
3. zenith dependence including shallow chords;
4. detector-weighted or response-weighted sensitivity before declaring an operator excluded or unpromising.

## Audit finding 4: finite-beta_tau evolution was not historically validated

The v0.15.1 convergence failure exposed a deeper problem inherited by finite-loss transport runs. The old algorithm applied continuous tau energy loss over a whole neutrino-interaction substep and only then applied decay. At PeV energies the tau may physically decay after tens of metres to kilometres, so this ordering can over-propagate and over-cool a tau before allowing regeneration.

Both `scripts/validate_rakshit_reya.py` and `scripts/validate_dutta_reno_sarcevic.py` instantiate `TauTransport(..., beta_tau=0.0, ...)`. Their PASS status validates the zero-loss decay/regeneration implementation, not the finite-`beta_tau` UHE solver.

v0.17 replaces this ordering by exact constant-density coupled tau loss+decay evolution and Strang-style composition. See `FINITE_BETA_TAU_REPAIR_v0.17.md`.

## What remains valid

The following hard-kernel/PDF results are unaffected:

- Case 5 has the same `y(1-y)` shape as Case 6 and the expected benchmark normalization ratio;
- at its nominal `c/(4 pi)=1` endpoint, Case 5 reaches `r_NP/SM ~ 0.055, 0.068, 0.073` at `100, 220, 300 PeV` with `<y> ~ 0.5`;
- at its endpoint, Case 7 reaches `r_NP/SM ~ 0.031, 0.054, 0.064` with `<y> ~ 0.8`;
- the chirality/external-state distinction remains physically necessary;
- the hard `sqrt(shat)<Lambda_UV` restriction is an EFT-validity boundary, not a physical turnover;
- v0.16 resolved-mediator cross-section/shape tables remain hard-kernel results.

## What is withdrawn pending v0.17 validation

Do **not** quote as established until the v0.17 gate passes and transport is rerun:

- v0.12 finite-loss Case-6 Earth-propagation numbers;
- transport-dependent v0.13 propagation+detector interpretation;
- v0.15 Case-5/Case-7 transport ratios and the claim that they fail;
- v0.15 Case-6 chirality robustness transport;
- v0.16 resolved-mediator transport results;
- the broad statement that the tested contact-operator route is self-defeating through Earth opacity.

## Required validation gate

Use `run_tau_transport_validation_v0_17.sh`. It requires:

1. analytic constant-density single-tau survival/energy-loss agreement;
2. decay-probability conservation and decay-energy-moment agreement;
3. beta_tau=0 regression to historical RR/DRS gates;
4. one-step/two-half-step refinement consistency;
5. representative full finite-loss transport convergence in `max_p`;
6. PREM matter-step refinement;
7. positivity/finite-number checks.

Only after this passes should any Case-5/6/7 or resolved-mediator transport be interpreted phenomenologically.
