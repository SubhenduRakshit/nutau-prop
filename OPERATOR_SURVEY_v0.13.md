# nutau-prop v0.13 — momentum-growing operator survey

Date: 2026-10-27

## Status entering this survey

The v0.11 detector-NSI stage and v0.12 flavor-aware Earth-transport stage have completed successfully on the Mac Studio. The ordered runner uses `set -euo pipefail`; therefore completion implies the historical Rakshit–Reya / Dutta–Reno–Sarcevic gate and the v0.11 differential-kernel regression gate both passed before v0.12 ran.

The new v0.12 result changes the interpretation of the earlier v0.8/v0.9 effective-area proxy. Case 6 produces very strong zenith-dependent opacity at UHE. At 220 PeV, the propagated tau-neutrino ratio is approximately 0.152, 8.97e-3, 1.45e-4, 5.43e-8 for |cos z|=0.5 and Lambda=1.5,2.0,2.44,3.0 TeV respectively, and is much smaller for a vertical chord. Hence a detector-only `(1+r)` rescaling is not a faithful UHE prediction.

## Momentum-growing operators in Pandey–Karmakar–Rakshit, JHEP 11 (2019) 046

### Case 5: derivative tensor neutrino–quark contact operator

L_eff = (c3/Lambda^3) d_nu (bar(nu_i^c) sigma^{mu nu} P_L nu_j) (bar q gamma_mu q),  i != j.

Published bounds in that paper:
- CHARM: c3/Lambda^3 <~ 2.9e-7 GeV^-3
- LHC monojet: c3/Lambda^3 <~ 1.8e-10 GeV^-3
- IceCube cascade-era bound: c3/Lambda^3 <~ 5.3e-8 GeV^-3

Physics priority: VERY HIGH. It has genuine momentum enhancement but a different Lorentz / inelasticity structure from Case 6. Because i != j, the hard interaction is intrinsically flavor converting. For a tau-initiated coupling this means the daughter can emerge in another neutrino flavor and need not remain subject to the same tau-specific opacity. A tau->mu assignment can also feed a muon-track channel after propagation even though the hard vertex itself remains neutral-current-like. This is precisely the qualitative loophole that Case 6 does not possess, so Case 5 is the first alternative operator that should be computed fully.

### Case 6: derivative vector neutrino–quark contact operator

L_eff = (c4/Lambda^3) (bar(nu^c) i<->d_mu nu) (bar q gamma^mu q).

Published LHC benchmark bound used in v0.7–v0.12:
- c4/Lambda^3 <~ 8.6e-10 GeV^-3.

Physics status: FULL differential Case-6 kernel, detector-only channel fold, and flavor-aware Earth transport now implemented. Remaining detector task is the combined propagation + detector response, because the v0.12 opacity is too large for the old factorized A_eff proxy.

### Case 7: neutrino–gluon contact operator

L_eff = (c5/Lambda^3) (bar(nu^c) P_L nu) G_{mu nu} G^{mu nu}.

Published bounds in that paper:
- CHARM: c5/Lambda^3 <~ 1.6e-6 GeV^-3
- LHC monojet: c5/Lambda^3 <~ 1.6e-10 GeV^-3
- IceCube cascade-era bound: c5/Lambda^3 <~ 5.5e-8 GeV^-3

Physics priority: MEDIUM/LOW for a UV-complete benchmark. The 2019 paper already discusses a Type-II-seesaw realization and finds the coefficient strongly suppressed by the neutrino Yukawa and scalar-mixing constraints, too small there to give a significant NC cross-section deviation. A generic EFT calculation is still useful as a kinematic counterexample because the gluon PDF weights a very different small-x region.

## Resolved-mediator momentum-dependent structures

The same 2019 paper also studied two dimension-5 neutrino-Z' vertices:

1. dipole-like: (c1/Lambda) (bar nu_i^c sigma_{mu nu} P_L nu_j) Z'^{mu nu}, i != j;
2. derivative-vector: (c2/Lambda) (bar nu_i^c <->d_mu nu) Z'^mu.

Both have additional momentum dependence, but unlike the contact operators their high-Q behavior is softened by the mediator propagator. The paper explicitly notes that the contact interactions rise faster with energy while the Z'-mediated cases are relaxed by propagator suppression.

These resolved-mediator structures are therefore not merely additional operator variants: they are direct candidates for the later physical-turnover step requested after the contact-operator survey.

## Ordered next calculation

1. Derive and independently validate the partonic |M|^2 and d^2 sigma/(dx dy) for Case 5. Do not infer it by rescaling Case 6; the tensor structure can materially alter y weighting.
2. Implement Case 5 in the same v0.11/v0.12 machinery, including explicit flavor conversion i != j and tau->mu / tau->e assignments.
3. Implement Case 7 only after the Case-5 result unless an analytic estimate shows a compelling gluon-enhanced window despite the stronger LHC bound and UV-completion suppression.
4. Compare Cases 5,6,7 at fixed currently allowed coefficient using: total r(E), <y>(E), zenith-dependent propagated spectra by flavor, secondary muon-neutrino spectrum, and channel-resolved detector yields.
5. Then move to the requested physical-turnover stage using an explicit resolved mediator rather than an ad hoc hard cutoff.

## Submission-gate implication after v0.12

The project is scientifically closer to submission because the historical transport, detector NSI topology, and full Case-6 Earth transport are now internally consistent. However, it is not yet submission-ready. The decisive remaining calculation is a joint propagation+detector fold with the same kernel. The v0.12 opacity is large enough that the old IceCube/KM3NeT common-flux screen may be quantitatively altered, so the earlier UHE 'clash' must not be quoted as final until that joint fold is performed.
