# nutau-prop v0.13 — momentum-growing operator survey

Date: 2026-10-27

## Status entering this survey

The v0.11 detector-NSI stage and v0.12 flavor-aware Earth-transport stage have completed successfully on the Mac Studio. The ordered runner uses `set -euo pipefail`; therefore completion implies the historical Rakshit–Reya / Dutta–Reno–Sarcevic gate and the v0.11 differential-kernel regression gate both passed before v0.12 ran.

The v0.13 joint propagation+detector fold has also passed its numerical bookkeeping gate. The full calculation confirms that Earth opacity does not rescue the contact Case-6 benchmark from the public IceCube UHE screen. Even the most favorable public-response envelope allows only about 0.30–0.31 of the Ehlert benchmark normalization, while the corresponding KM3NeT bright-track expectation remains only about 0.07 events. Thus the old v0.8/v0.9 tension survives the self-consistent propagation calculation, although the physical reason is now clearer: detector-vertex enhancement competes with strong zenith-dependent depletion of the incident tau-neutrino flux.

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

Published LHC benchmark bound used in v0.7–v0.13:
- c4/Lambda^3 <~ 8.6e-10 GeV^-3.

Physics status: COMPLETE for the present public-response study. The differential Case-6 kernel, detector-only topology treatment, flavor-aware Earth transport, secondary-nu_mu production, and the joint propagation+detector fold have all been implemented. The joint v0.13 result shows that the contact interaction gives a modest regeneration excess near the PeV scale but develops severe zenith-dependent opacity at tens to hundreds of PeV. The old detector-only `(1+r)` picture is therefore not physically adequate, but the corrected treatment still leaves a strong UHE constraint.

#### What the cutoff scale means in the Case-6 calculation

There are two distinct quantities that should not be conflated.

First, the interaction strength used in the numerical scan is the dimensionful coefficient

C4 = 8.6e-10 GeV^-3,

which is held fixed. With that choice the partonic matrix element grows as

|M|^2 = 64 C4^2 shat^3 y(1-y),

so the new contribution rises rapidly with partonic energy. This momentum growth is the physical feature being tested.

Second, the scan parameter Lambda_UV is an EFT-validity boundary, imposed primarily as

sqrt(shat) < Lambda_UV.

It does not represent a physical suppression mechanism and it does not make the interaction turn over. Phase-space points with sqrt(shat) above Lambda_UV are simply discarded because the contact description is no longer trusted there. Since shat = x s ~ 2 m_N E_nu x, the cut corresponds at fixed neutrino energy to

x < Lambda_UV^2 / (2 m_N E_nu).

As E_nu increases, more of the high-x, high-shat phase space lies outside the EFT domain. Conversely, increasing Lambda_UV admits more of that rapidly growing phase space and therefore produces a larger cross section and stronger Earth opacity. This explains why the 1.5, 2.0, 2.44 and 3.0 TeV scans become progressively more absorptive even though C4 itself is unchanged.

The notation C4 = c4 / Lambda^3 is useful for judging the size of an underlying dimensionless coefficient, but it must not be misread as the way the numerical scan is performed. In the present scan C4 is fixed and

c4 = C4 Lambda_UV^3

therefore changes with Lambda_UV. The corresponding c4/(4 pi) values are 0.231, 0.547, 0.994 and 1.848 for Lambda_UV = 1.5, 2.0, 2.44 and 3.0 TeV. Thus 2.44 TeV is already near the nominal 4 pi coefficient diagnostic and 3 TeV is beyond it. This is only a perturbativity diagnostic, not a proof of unitarity or UV consistency.

The conceptual distinction is therefore:

- contact Case 6 + hard sqrt(shat) cut = rapidly growing EFT interaction with an imposed boundary of calculational validity;
- resolved mediator = a physical amplitude whose propagator changes the high-momentum behavior and can produce a genuine turnover or saturation.

The hard EFT cut by itself cannot be advertised as the mechanism that keeps the interaction large in the PeV propagation window while making it safe at hundreds of PeV or EeV. That requires an explicit UV completion or resolved mediator.

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

These resolved-mediator structures are therefore not merely additional operator variants: they are direct candidates for the later physical-turnover step requested after the contact-operator survey. In particular, they provide the correct language for replacing the present hard EFT-validity cut by a physical high-momentum modification of the amplitude.

## Ordered next calculation

1. Derive and independently validate the partonic |M|^2 and d^2 sigma/(dx dy) for Case 5. Do not infer it by rescaling Case 6; the tensor structure can materially alter y weighting.
2. Implement Case 5 in the same v0.11/v0.12 machinery, including explicit flavor conversion i != j and tau->mu / tau->e assignments.
3. Implement Case 7 only after the Case-5 result unless an analytic estimate shows a compelling gluon-enhanced window despite the stronger LHC bound and UV-completion suppression.
4. Compare Cases 5,6,7 at fixed currently allowed coefficient using: total r(E), <y>(E), zenith-dependent propagated spectra by flavor, secondary muon-neutrino spectrum, and channel-resolved detector yields.
5. Then move to the requested physical-turnover stage using an explicit resolved mediator rather than an ad hoc hard cutoff.

## Submission-gate implication after v0.13

The Case-6 calculation itself is now essentially complete at the level possible with public detector responses: historical transport validation, differential cross section, topology-aware detector treatment, Earth propagation and regeneration, secondary-flavor production, and the joint propagation+detector fold are all in place. The conclusion is that contact Case 6 can generate a distinctive regeneration-to-opacity transition but does not provide the desired viable UHE enhancement once the IceCube public constraint is imposed.

This is not a general no-go theorem for momentum-growing tau-neutrino interactions. Case 5 is qualitatively different because it is flavor changing, and a resolved mediator is qualitatively different because it can physically soften the high-momentum growth. Those two mechanisms are therefore the decisive remaining tests before the project-level submission decision.
