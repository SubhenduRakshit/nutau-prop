# nutau-prop — canonical UHE BSM / propagation gap-closure ToDo v0.90

Date: 2026-10-31
Status: LIVING CANONICAL TODO LIST
Companion record: `UHE_BSM_EXPLORATION_MASTER_v0.90.md`

Purpose: maintain one place for every unresolved physics, numerical, experimental, novelty and publication gap. Update this file whenever a gap is closed, a new blocker is found, or a branch is demoted/promoted. Preserve failed attempts in the historical files; do not delete them.

## Status vocabulary

- `TODO` — not yet started / not yet certified
- `IN PROGRESS` — active work exists but gate not closed
- `BLOCKED` — requires unavailable detector response, external data, or model specification
- `PASS` — required gate closed positively
- `FAIL` — required gate closed negatively; preserve details
- `DEFER` — not worth current resources unless higher-priority branches fail

## Priority order

### P0 — charged-scalar / charged-Higgs loophole (highest priority)

Goal: determine whether a modern 2026 gauge-complete `H^+ (cs, tau nu_tau)` model can still produce the large UHE CC effect already seen in our screen.

- [ ] `TODO` Reconstruct Huang et al. (2112.09476) scalar model exactly: field content, gauge quantum numbers, Yukawa texture, chirality and benchmark definitions.
- [ ] `TODO` Derive the exact `nu_tau s -> tau c` and conjugate differential cross sections with finite `m_tau`, `m_c`, full propagator and consistent spin/color averages.
- [ ] `TODO` Reproduce Huang et al.'s published benchmark cross sections / event trends as an external validation.
- [ ] `TODO` Map our broken-phase `(g_tau,g_L,g_R,M_H)` parameterization onto a gauge-complete scalar sector; identify all neutral/charged partners and CKM-related couplings.
- [ ] `TODO` Apply current `D_s -> tau nu`, `K/tau` decay, flavour, electroweak and neutrino-mass constraints relevant to the actual texture.
- [ ] `TODO` Apply current direct LHC constraints, including modern `H^+ -> cs`, `H^+ -> tau nu`, top-decay/top-associated searches where applicable; do not copy bounds from incompatible production modes.
- [ ] `TODO` Check perturbativity, total width, vacuum stability/unitarity of the scalar sector and whether the `~90 GeV` benchmark is internally consistent.
- [ ] `TODO` Determine the allowed envelope in `(M_H, y_cs, y_tau)` and compare directly with the coupling needed for 5%, 10%, 20% UHE effects.
- [ ] `TODO` If an allowed benchmark survives, implement exact CC kernel in the Earth transport solver.
- [ ] `TODO` Compute emerging-tau probability `P(nu_tau -> tau_exit)(E_nu,E_tau,beta)`, not only outgoing neutrino flux.
- [ ] `TODO` Fold with GRAND/POEMMA/Trinity acceptance and compare with Huang et al. using current detector projections.
- [ ] `TODO` Novelty gate: identify what is genuinely new relative to Huang et al.; likely candidates are updated 2026 constraints, full all-flavour secondaries, and/or a distinct viable flavour texture.

Completion criterion: one explicit gauge-complete benchmark that simultaneously passes present constraints and yields a detector-level effect large enough for a next-generation tau-neutrino telescope, OR a documented exclusion of the entire targeted texture.

### P1 — resonant models

#### P1a strange-aligned `R2-tilde` leptoquark

- [ ] `TODO` Re-derive exact resonance cross section and finite-width result for the stored `M~1.5 TeV, lambda~2` benchmark.
- [ ] `TODO` Include both weak-doublet components and all gauge-related channels.
- [ ] `TODO` Apply current pair-production, single-production, `tau+jet`, MET, ditau and flavour constraints for the exact strange/tau texture.
- [ ] `TODO` Compare with Kirk–Okawa–Wu and identify whether our flavour choice genuinely avoids their dominant exclusions.
- [ ] `TODO` If viable, propagate with full `nu_tau-tau-nu_tau` transport and compute emerging-tau detector observables.
- [ ] `TODO` Determine whether its signal can be distinguished from the charged-scalar scenario and from small-x QCD uncertainty.

Completion criterion: allowed resonance with a detector-level signature, or robust closure of the benchmark/texture.

#### P1b tau-neutrinogluon / leptogluon

- [ ] `TODO` Specify a complete gauge-consistent colored multiplet and transition operator normalization.
- [ ] `TODO` Justify or exclude the large stored matching coefficient `a~4` from a plausible UV mechanism.
- [ ] `TODO` Compute exact width and finite-width resonance fold.
- [ ] `TODO` Apply current colored-partner pair-production, `tau+jet`, jet+MET and relevant contact-tail constraints.
- [ ] `TODO` If viable, run full Earth propagation and detector fold.

Completion criterion: a UV-credible allowed benchmark with >= useful detector distortion, or closure.

### P2 — CGC / small-x QCD control and discriminator

Current v0.88b result is a total-cross-section surrogate with SM conditional `y` shapes.

- [x] `PASS` Full PREM all-flavour propagation with CC removal, NC feed-down, `nu_tau-tau-nu_tau` regeneration and tau-decay secondary e/mu channels for the surrogate.
- [x] `PASS` Store `gamma=2.0,2.5`, fixed chords and upgoing-angle averages.
- [ ] `TODO` Produce an `E-cos(theta)` map together with absolute survival / event-weight proxy; mask extremely opaque bins so huge ratios are not misinterpreted.
- [ ] `TODO` Verify the 100-PeV feed-down feature quantitatively by decomposing contributions from parent energies above 100 PeV.
- [ ] `TODO` Replace the surrogate with a published/digitized or directly implemented nonlinear-BK/CGC `d sigma/dy` (or structure-function) calculation.
- [ ] `TODO` Add modern NNLO + small-x resummation + PDF/nuclear uncertainty as the SM/QCD baseline before claiming a CGC discrepancy.
- [ ] `TODO` Investigate composition-dependent nuclear saturation along PREM: use `R_CGC(E,A(r))` rather than one nucleon-level suppression everywhere.
- [ ] `TODO` Check whether a core/mantle feature survives realistic detector angular smearing.
- [ ] `TODO` Literature novelty gate: exhaustively verify whether PREM-composition-dependent nuclear CGC + full all-flavour regeneration has already been done.
- [ ] `TODO` Compare CGC against any surviving BSM model using the same `(E, cos theta, flavour, topology)` likelihood.

Completion criterion: either a precision QCD-control model suitable for profiling as a nuisance, or a defensible new nuclear-CGC propagation result.

### P3 — all-flavour regeneration / secondary-neutrino novelty

- [x] `PASS` Tau decay regenerates `nu_tau` and generates secondary `nu_mu`; corresponding e secondary included in the CGC all-flavour runner.
- [ ] `TODO` Audit exact charge-conjugation bookkeeping for all `tau-` and `tau+` decay secondaries in every transport branch.
- [ ] `TODO` Quantify how much of the final `nu_e` and `nu_mu` flux comes specifically from tau-decay secondaries as a function of `E` and zenith.
- [ ] `TODO` Define robust observables such as `D_tau/mu = (Phi_tau/Phi_mu)_model/(Phi_tau/Phi_mu)_SM` and secondary-flavour fractions.
- [ ] `TODO` Test whether these observables distinguish tau-specific BSM from flavour-universal CGC/PDF changes after nuisance profiling.
- [ ] `TODO` Check literature specifically for using regenerated secondary `nu_e/nu_mu` as a discriminator of BSM in UHE Earth propagation.
- [ ] `TODO` Translate flux differences into detector topology changes (tracks/showers/tau tags), rather than only flavour flux ratios.

Completion criterion: identify a genuinely new and experimentally useful all-flavour observable, or document that existing literature already covers it.

### P4 — GRAND / POEMMA / Trinity observability

- [ ] `TODO` Implement emerging-tau kernel / exit probability using the existing tau propagation machinery.
- [ ] `TODO` Validate against a published SM emerging-tau benchmark.
- [ ] `TODO` Obtain current public acceptances / effective areas for GRAND, POEMMA and Trinity (or use clearly labelled published projections when machine-readable responses are unavailable).
- [ ] `TODO` Fold each surviving BSM/CGC model with detector acceptance and exposure.
- [ ] `TODO` Use energy + emergence-angle likelihoods and profile astrophysical flux normalization/slope.
- [ ] `TODO` Report absolute expected event counts; do not interpret huge relative flux ratios in opaque bins as observability.
- [ ] `TODO` Compare our forecast directly with Huang et al. and later tau-telescope studies.

Completion criterion: detector-level significance / exclusion forecast for each surviving model.

### P5 — PKR Eq. 4.3 / 4.4 residual low-mass loopholes

- [x] `PASS` Correct project normalization and preserve superseded results.
- [x] `PASS` Exact CMS dijet inputs and tested 50/100/220-GeV collider-compatible maxima.
- [ ] `TODO` Resolve 5-GeV and ~40-GeV parameter space with appropriate low-mass dijet/quarkonium/meson constraints and operator-specific collider searches.
- [ ] `TODO` Replace pilot monojet shapes at higher masses with a full validated Run-2 recast if this branch remains worth pursuing.
- [ ] `BLOCKED` Exact IceCube arbitrary-y detector reweight requires channel-resolved response / MC truth not present in the public HESE schema inspected so far.
- [ ] `TODO` If a stable sterile `N` is used, implement bidirectional `nu_tau <-> N` transport rather than an irreversible sink and include possible N-initiated detector interactions.
- [ ] `TODO` Reassess whether any remaining low-mass point can reach detector relevance after corrected normalization; if not, close this branch explicitly.

Completion criterion: close or certify the remaining low-mass window without overclaiming IceCube constraints.

### P6 — strong/endpoint EFT branches

Includes active gluon d=8, sterile derivative-tensor and sterile Rayleigh/gluon.

- [ ] `DEFER` Do not spend long transport runs unless an explicit strong/composite UV model is specified.
- [ ] `TODO` For active gluon d=8, if revisited: finish crossed/coupled-channel unitarity and operator-specific collider constraints.
- [ ] `TODO` Quantify regulator/cutoff dependence as a theory uncertainty, not as alternate predictions.
- [ ] `TODO` Seek explicit composite/resolved mediator realizations only if they provide a controlled rate larger than the contact EFT.

Completion criterion: explicit UV completion or permanent demotion to illustrative endpoint benchmarks.

### P7 — modern-QCD / detector-systematics baseline

- [ ] `TODO` Build one consistent SM reference using modern neutrino DIS, small-x resummation and nuclear corrections over the full energy range used in transport.
- [ ] `TODO` Propagate PDF/QCD uncertainties through the Earth solver.
- [ ] `TODO` Include astrophysical flux normalization, spectral index and source-composition uncertainties in detector forecasts.
- [ ] `TODO` Keep absolute event counts and detector backgrounds wherever public information allows.
- [ ] `TODO` Distinguish theory uncertainty from BSM signal before any discovery claim.

### P8 — publication / novelty gates

- [ ] `TODO` Before drafting a paper, formulate one sharp question rather than a broad operator survey.
- [ ] `TODO` Candidate paper A: modern 2026 viability + UHE tau-telescope consequences of the `cs-tau nu` charged-scalar loophole.
- [ ] `TODO` Candidate paper B: distinguish nonlinear/nuclear small-x QCD from a surviving tau-specific BSM interaction using full all-flavour transport.
- [ ] `TODO` Candidate paper C: secondary-flavour regeneration as a discriminator, only if literature check confirms novelty and detector effect survives.
- [ ] `TODO` For the chosen paper, perform a fresh literature search immediately before claiming novelty.
- [ ] `TODO` Maintain a claims table: calculated / screened / constrained / detector-observable / novel.
- [ ] `TODO` Do not call a broad “no-go theorem” unless the exact model space and assumptions are stated.

## Closed branches retained for reference, not active TODO

- sterile d=6 scalar: controlled effect ~0.006%; CLOSED.
- sterile d=6 tensor: corrected ~0.5–0.7%; CLOSED for several-percent target.
- active cR/sR two-derivative d=8 at tested benchmark: <0.25% total; CLOSED for target.
- sterile electroweak dipole as hard-DIS enhancement: CLOSED for target.
- universal `B-3L_tau` vector under tested assumptions: sub-percent; CLOSED on size.

These should only be reopened if a genuinely different UV realization or experimental assumption changes the relevant bound or mechanism.

## Update policy

Whenever work is completed:

1. change the checkbox/status here;
2. add the exact versioned artifact/code/result path;
3. record PASS/FAIL and the numerical reason;
4. update `UHE_BSM_EXPLORATION_MASTER_v0.90.md` if the scientific conclusion changes;
5. preserve superseded/failed calculations rather than deleting them;
6. increment version number when the decision tree changes materially.

## Immediate next action

Highest-value next task:

**Re-certify the 90-GeV `cs-tau nu_tau` charged-scalar model against current 2026 flavour + collider constraints, because it is the only already-computed nonresonant weak-scale mechanism in our portfolio that naturally produced O(10%)–O(20%) UHE charged-current effects without relying on an EFT endpoint.**
