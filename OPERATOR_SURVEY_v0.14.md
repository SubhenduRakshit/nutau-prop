# nutau-prop v0.14 — momentum-growing operator prescreen

Date: 2026-10-27

## Purpose

Continue the agreed operator survey after the completed Case-6 v0.13 joint propagation+detector fold. Priority remains: Case 5 first, then Case 7, then a resolved mediator/physical-turnover model.

## Case 5: derivative tensor/vector quark operator

Pandey–Karmakar–Rakshit Case 5:

L5 = C5 * d^nu( bar(nu_i^c) sigma_{mu nu} P_L nu_j ) (bar q gamma^mu q),  i != j,

with the 2019 monojet benchmark bound C5 = 1.8e-10 GeV^-3.

The complete low-energy neutrino EFT of Altmannshofer–Tammaro–Zupan contains the corresponding derivative tensor/vector DIS structure Q10. Its hard-scattering kernel has the same y(1-y) dependence as the derivative vector Case-6-like Q8 structure, but half the matrix-element normalization:

|M5|^2 = 32 C5^2 shat^3 y(1-y),
|M6|^2 = 64 C6^2 shat^3 y(1-y).

Therefore, for the same PDF and the same sqrt(shat)<Lambda_UV domain restriction,

sigma5/sigma6 = (1/2) (C5/C6)^2 = 0.0219037317

when C5=1.8e-10 and C6=8.6e-10 GeV^-3.

This does NOT automatically make Case 5 irrelevant: Case 6 produced exponentially large opacity on long chords, so a few percent of its extra optical depth can still matter. Also, because C5 is smaller, the nominal c=4pi diagnostic is reached only at Lambda_UV = 4.118 TeV, compared with 2.445 TeV for Case 6. The exact-PDF v0.14 prescreen therefore extends Case 5 to its own perturbative-domain edge before deciding whether full transport is warranted.

### External-state caution — mandatory before transport

The previous shorthand that Case 5 simply produces an active nu_mu daughter was too quick. The charge-conjugated tensor bilinear requires an explicit Dirac/Majorana interpretation. In the standard low-energy EFT classification, derivative tensor operators are chirality flipping for Dirac neutrinos, while the Majorana tensor analogue vanishes flavor-diagonally and only mediates flavor-changing transitions. Before any Case-5 Earth transport we must derive the exact external state of the Pandey operator and decide whether the daughter is an active antineutrino/Majorana helicity state, an active different flavor, or a right-handed sterile state in the chosen UV interpretation. No track claim is allowed before this audit.

This observation also triggers a robustness audit of Case 6, whose Pandey operator likewise contains bar(nu^c). The completed Case-6 opacity result remains a valid calculation for the diagonal NC-like kernel assumption used in v0.11-v0.13, but publication wording must not hide the charge-conjugation/external-state assumption. If the operator implies nu<->antinu conversion, a coupled nu/antinu transport correction is required and should be quantified.

## Case 7: neutrino-gluon Rayleigh operator

Pandey–Karmakar–Rakshit Case 7:

L7 = C7 * (bar(nu^c) P_L nu) G_{mu nu} G^{mu nu},

with the 2019 monojet benchmark bound C7 = 1.6e-10 GeV^-3.

Using the standard low-energy EFT Rayleigh normalization and absorbing alpha_s/(12 pi) into the coefficient convention used by Pandey gives

|M7|^2 = 16 C7^2 shat^3 y^3.

Thus the Case-7 kernel is much harder in inelasticity than Case 6: before PDF/cut effects, y^3 alone has <y>=4/5, whereas y(1-y) has <y>=1/2. This suggests stronger energy degradation per BSM scattering. On the other hand, the coefficient is smaller and the gluon PDF can partially compensate. The competition cannot be settled without the explicit CT10nlo gluon fold, which is implemented in v0.14.

The nominal c=4pi scale for C7=1.6e-10 GeV^-3 is Lambda_UV = 4.282 TeV.

### UV-completion caution

Pandey et al. discussed a Type-II-seesaw realization of Case 7 and found the effective coefficient proportional to the neutrino Yukawa times scalar mixing. Neutrino-mass and scalar-mixing constraints make that particular renormalizable completion far too small to generate a significant NC cross-section deviation. Therefore Case 7 is presently a useful kinematic/PDF counterexample, not yet a compelling UV-complete positive benchmark.

## Gauge-invariance caution

The complete low-energy EFT analysis shows that derivative dimension-7 neutrino-current operators can arise from electroweak-gauge-invariant dimension-8 structures after the Higgs gets a vev. Consequently, a low-energy coefficient C/Lambda^3 and a hard monojet contact bound should not automatically be interpreted as a complete gauge-invariant UV theory. This strengthens the motivation for the later resolved-mediator stage.

## v0.14 numerical gate

Files:
- src/nutau_prop/momentum_operators.py
- scripts/run_momentum_operator_prescreen_v0_14.py
- run_momentum_operator_prescreen_v0_14.sh

The prescreen computes exact CT10nlo PDF folds for Cases 5 and 7 at E = 0.1, 1, 10, 100, 220, 300 PeV, using the same Qmin and primary sqrt(shat)<Lambda_UV prescription as Case 6. It reports sigma_NP/sigma_SM, <y>, c/(4pi), and the ratio to the same-cut Case-6 kernel. Case 5 carries an exact regression gate against the analytic factor 0.0219037317.

Decision after this gate:
1. If Case 5 reaches a propagation-relevant optical depth before its perturbativity edge, resolve its external-state structure and run coupled flavor/helicity transport.
2. If Case 7 reaches a propagation-relevant cross section, run its high-y gluon transport; otherwise retain it as a documented negative counterexample.
3. Only after these contact cases, move to the resolved-mediator physical-turnover calculation.
