# nutau-prop v0.15 — momentum-growing contact operators after exact-PDF prescreen

Date: 2026-10-27

## What v0.14 established

The exact CT10nlo prescreen for Pandey–Karmakar–Rakshit Cases 5 and 7 passed its kernel gate.

At the nominal perturbative endpoints c/(4pi)=1:

- Case 5: Lambda_UV = 4.1176 TeV, r_NP/SM = 0.0552, 0.0678, 0.0732 at 100, 220, 300 PeV, with <y> ~ 0.50.
- Case 7: Lambda_UV = 4.2825 TeV, r_NP/SM = 0.0311, 0.0538, 0.0642 at 100, 220, 300 PeV, with <y> ~ 0.80.

At equal EFT cut, Case 5 is exactly 0.0219037 times the Case-6 NP cross section at the respective 2019 LHC benchmark coefficients, because the Case-5 hard kernel has the same y(1-y) shape and half the matrix-element normalization.

The prescreen therefore does not justify a broad parameter scan.  Only the 3 TeV point and the c/(4pi)=1 endpoint are carried into transport.

## External-state correction

A crucial operator-basis point must be kept explicit.

Altmannshofer–Tammaro–Zupan (arXiv:1812.02778) show that the relevant dimension-5 and dimension-7 neutrino bilinears are chirality flipping for Dirac neutrinos: an incoming left-handed active neutrino is converted into a right-handed neutrino.  For Majorana neutrinos the corresponding operators violate lepton number by two units; the tensor derivative operators vanish for flavor diagonal alpha=beta and are therefore flavor changing.

The original Pandey–Karmakar–Rakshit Case 5 explicitly has i != j, following the same hermiticity argument as their transition dipole operator.  Therefore the earlier shorthand 'tau -> mu escape' is not sufficient.

### Brackets used in v0.15

**Dirac bracket:** NP scattering is an active-flux sink because the daughter is right handed and sterile to ordinary weak detection to leading order.

**Majorana Case 5 bracket:** choose the tau-mu transition as the most track-favorable flavor assignment, but propagate both directions required by the off-diagonal Hermitian coupling:

- nu_tau -> antinu_mu,
- antinu_tau -> nu_mu,
- nu_mu -> antinu_tau,
- antinu_mu -> nu_tau.

Thus Case 5 is *not* treated as a one-way escape from further BSM scattering.

**Majorana Case 7 bracket:** flavor-diagonal nu_tau <-> antinu_tau helicity flip with the Case-7 y kernel.

## Why Case 5 remains the first contact-operator transport test

Case 5 has the larger allowed NP fraction at 100–300 PeV and a much less destructive inelasticity distribution than Case 7.  Its special possible observable is flavor redistribution, especially the muon-track fraction, rather than a large total event-rate enhancement.

## Why Case 7 remains second

Case 7 benefits from the gluon PDF but its y^3 kernel gives <y> ~ 0.8.  Each scattering therefore removes most of the parent neutrino energy from the UHE bin.  In addition, the explicit Type-II-seesaw completion discussed in JHEP 11 (2019) 046 suppresses its coefficient strongly through the neutrino Yukawa and scalar mixing.  It is retained as an EFT counterexample, not presently as the preferred UV benchmark.

## Case-6 robustness audit added in v0.15

The same chirality issue affects the physical interpretation of the Case-6 daughter.  The published v0.12 calculation treated the NP daughter as the same active tau-neutrino species.  v0.15 therefore reuses the exact Case-6 tables and brackets that assumption by:

1. Dirac active-to-RH sink;
2. Majorana nu_tau <-> antinu_tau helicity flip.

This is a robustness audit, not a retraction of the v0.12 numerical kernel.  If both physically motivated brackets are at least as suppressive at UHE, the conclusion that contact Case 6 fails the positive-signal gate is strengthened.

## v0.15 decision rule

Run full differential Earth transport only for:

- Case 5: Lambda_UV = 3 TeV and 4.1176 TeV;
- Case 7: Lambda_UV = 3 TeV and 4.2825 TeV.

If neither produces a distinctive surviving energy/zenith/flavor signal at the perturbative endpoint, stop broadening the contact-operator survey and move to the already planned resolved-mediator calculation.  The resolved mediator addresses the actual Case-6 failure mechanism by providing a physical high-Q turnover rather than a hard EFT-domain cut.
