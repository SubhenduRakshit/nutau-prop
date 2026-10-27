# nutau-prop v0.17 — finite-beta_tau transport repair and revalidation

Date: 2026-10-28

## Why v0.17 was required

An external audit of the v0.15 chirality-aware transport exposed two immediate issues:

1. tau-decay secondary muon neutrinos were accumulated in one solver but not subsequently propagated through SM interactions;
2. matched SM and BSM runs used different substep schedules, allowing a pure additional Dirac sink to appear spuriously above the matched SM solution.

The v0.15.1 correction fixed those two issues, but its convergence gate then failed systematically.  Inspection exposed a deeper finite-`beta_tau` weakness inherited from the older transport implementation.

The old algorithm evolved a tau over a full neutrino-interaction substep by:

1. applying continuous energy loss over the entire substep;
2. only afterwards applying tau decay over that substep.

At PeV energies the tau decay length can be tens of metres to kilometres, while a neutrino interaction/matter substep can be much longer.  A newly created tau could therefore be forced to lose energy over a path much longer than its physical lifetime before it was allowed to decay.  This is a first-order operator-splitting error and explains the strong sensitivity to `max_p` found by v0.15.1.

## Historical-validation limitation discovered

The repository's Rakshit-Reya and Dutta-Reno-Sarcevic historical validators both explicitly instantiate the transport with

`beta_tau = 0.0`.

Those historical PASS gates remain useful validations of the interaction kernels, tau-decay/regeneration spectra and secondary-channel structure in the zero-loss limit, but they **did not validate finite tau energy loss**.  Therefore they cannot be cited as validation of the finite-`beta_tau` UHE transport used in v0.12-v0.16.

## Exact constant-density tau evolution

For continuous tau energy loss

`dE/ds = - beta_tau rho E`,

write

`b = beta_tau rho`,

so that

`E(s) = E0 exp(-b s)`.

The instantaneous decay hazard is

`1/L_tau(E) = m_tau / (c tau_tau E)`.

Defining

`a = m_tau / (c tau_tau E0)`,

the accumulated decay optical depth over a constant-density segment is

`U(s) = (a/b) [exp(b s)-1]` for `b>0`,

with the smooth zero-loss limit `U=a s`.

Hence

`P_surv(s) = exp[-U(s)]`.

For a decay occurring at optical depth `u`, the tau energy at the decay point is

`E_decay(u) = E0 / [1 + (b/a) u]`.

The implementation integrates the decay-energy distribution using

`t = 1-exp(-u)`,

for which the physical decay measure is simply

`exp(-u) du = dt`.

Gauss-Legendre quadrature in `t` therefore integrates the daughter-energy distribution directly, without requiring metre-scale spatial stepping.

Implementation: `src/nutau_prop/tau_evolution.py`.

## Transport composition

All repaired tau-neutrino transport paths now use Strang-style splitting inside each constant-density neutrino substep:

`A(h/2) -> B(h) -> A(h/2)`,

where

- `A` is the exact coupled tau energy-loss + decay evolution;
- `B` is neutrino CC/NC/NP interaction and redistribution.

This means an existing tau is evolved for half a step, neutrino interactions occur at the midpoint, newly created taus are added there, and all taus are evolved for the remaining half step.

Updated modules:

- `src/nutau_prop/transport.py`;
- `src/nutau_prop/transport_case6.py`;
- `src/nutau_prop/chirality_transport.py`.

The v0.15.1 matched SM/BSM step schedule is retained, and secondary muon neutrinos are propagated after production.

## Independent v0.17 validation gate

New script: `scripts/validate_tau_transport_v0_17.py`.

It performs no PDF integration and reuses existing v0.11 tables.  The gate contains:

1. direct comparison of single-tau survival against the analytic constant-density formula;
2. probability conservation between surviving and decayed tau populations;
3. comparison of the decay-energy logarithmic moment against an independent `scipy.integrate.quad` calculation;
4. exact `beta_tau=0` regression tests;
5. one-step versus two-half-step semigroup/refinement tests;
6. representative full Case-6 transport convergence under `max_p = 0.04, 0.02, 0.01`;
7. independent PREM matter-step refinement from 12 km to 6 km;
8. positivity and finite-number checks;
9. separate reruns of the historical Rakshit-Reya and Dutta-Reno-Sarcevic `beta_tau=0` gates.

Runner: `run_tau_transport_validation_v0_17.sh`.

Independent parameter points are multiprocessing-parallelized through `NUTAU_NPROC`.

## Status of previous results

### Still valid

The following do not depend on finite-`beta_tau` Earth transport and remain valid unless separately audited:

- Case-5/6/7 hard matrix-element derivations;
- exact-PDF cross-section tables;
- `d sigma/dy` and mean-inelasticity results;
- EFT-domain/perturbativity diagnostics;
- v0.14 Case-5/7 prescreen;
- v0.16 resolved-mediator cross-section/shape tables.

### Provisional pending v0.17 PASS

Do not use as established physics conclusions until the v0.17 finite-`beta_tau` gate passes and the relevant transport is rerun:

- v0.12 Case-6 finite-loss Earth-propagation numbers;
- transport-dependent parts of v0.13 detector+propagation conclusions;
- v0.15 Case-5/7 transport ratios and the claim that those operators fail;
- v0.15 Case-6 chirality robustness transport;
- v0.16 resolved-mediator transport results.

This does **not** mean those physical conclusions are necessarily wrong.  It means their finite-`beta_tau` numerical support must be regenerated with a validated solver.

## Interpretation discipline after the repair

Even after the transport solver passes, a small BSM/SM ratio in a bin where the SM transmitted flux is already negligible is not by itself a detector-level exclusion or a useful negative result.

Subsequent phenomenology must report:

- incident spectral weight;
- SM transmitted weight;
- BSM transmitted weight;
- zenith dependence including shallow chords;
- detector-response weighting before deciding whether detector enhancement or Earth opacity wins.

No general no-go theorem is claimed.
