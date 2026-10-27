# nutau-prop v0.16 — resolved derivative-Z' turnover test

Date: 2026-10-27

## Why this stage exists

Cases 5, 6 and 7 established a structural problem for momentum-growing contact interactions: once the new interaction is large enough to modify 10--300 PeV propagation, Earth opacity becomes self-defeating at higher energy. The next physics question is therefore not whether the interaction can grow, but whether that growth can be physically softened before the UHE regime.

## Model

We use the derivative-vector neutrino-Z' structure already studied in Pandey, Karmakar and Rakshit, JHEP 11 (2019) 046:

L ⊃ (c2/Lambda) (bar(nu^c) i<->partial_mu nu) Z'^mu
    + g_q (bar q gamma_mu q) Z'^mu.

Define

kappa = c2 g_q / Lambda  [GeV^-1].

For Q^2 << M_Z'^2, integrating out the mediator gives the Case-6 contact coefficient

C4_lowQ = kappa / M_Z'^2.

The v0.16 primary scan fixes C4_lowQ to the Case-6 benchmark,

C4_lowQ = 8.6e-10 GeV^-3,

and varies M_Z'. This is a shape test: it asks whether the same low-Q interaction can avoid runaway high-Q growth once the mediator is resolved.

The exact differential kernel is the Case-6 kernel multiplied by

[M_Z'^2/(M_Z'^2 + Q^2)]^2.

Equivalently,

d^2 sigma_res/dxdy = d^2 sigma_C6(contact)/dxdy * [M_Z'^2/(M_Z'^2 + Q^2)]^2.

At Q^2 << M_Z'^2 the contact result is recovered. At Q^2 >> M_Z'^2 the explicit high-energy growth is softened by the propagator. This is a physical turnover mechanism, unlike the hard sqrt(shat)<Lambda_UV restriction used only as an EFT-validity diagnostic for the contact operator.

## Scan

Primary masses:

M_Z' = 0.1, 0.2, 0.5, 1, 2, 3 TeV.

The required low-Q-matching coupling product is printed for each mass:

kappa = C4_lowQ M_Z'^2.

For reference this gives approximately

- 0.1 TeV: 8.6e-6 GeV^-1
- 0.2 TeV: 3.44e-5 GeV^-1
- 0.5 TeV: 2.15e-4 GeV^-1
- 1.0 TeV: 8.6e-4 GeV^-1
- 2.0 TeV: 3.44e-3 GeV^-1
- 3.0 TeV: 7.74e-3 GeV^-1

These are not asserted to be currently allowed. Coupling viability is a separate gate after the shape/transport gate.

## Historical constraint context

Pandey et al. found for this derivative-Z' interaction at M_Z'=5 GeV an LHC monojet benchmark kappa <~ 3.3e-3 GeV^-1 and an old IceCube cross-section benchmark kappa <~ 2.5e-4 GeV^-1. Their mass scan showed that the LHC bound becomes more important than IceCube in approximately the 40--220 GeV range. These 2019 numbers are historical benchmarks, not a 2026 recast.

A targeted literature search performed at v0.16 did not identify a later paper that directly supersedes the exact derivative operator and its mass-dependent monojet/IceCube comparison. This should not be interpreted as proof that no stronger indirect or model-dependent constraint exists.

## Important EFT caveat

Resolving the Z' propagator removes the contact-mediator approximation, but the neutrino-Z' vertex itself is dimension five. Therefore the model is not automatically a complete UV theory. The primary v0.16 scan deliberately uses no hard sqrt(shat) cutoff so that the propagator effect is isolated. The code exposes an optional vertex cutoff only as a diagnostic.

## Transport treatment

The resolved derivative-vector interaction inherits the Case-6 external-state/chirality issue. v0.16 therefore keeps the same two brackets used in the Case-6 robustness audit:

1. Dirac sink: the chirality-flipped daughter is sterile to SM interactions.
2. Majorana flip: nu_tau <-> antinu_tau with the exact NP energy-degradation kernel.

SM CC/NC propagation, tau loss/decay, regeneration, PREM chords and secondary nu_mu follow the already validated machinery.

## Decision gate

A useful mediator mass must satisfy both qualitative requirements:

1. retain a non-negligible 10--300 PeV interaction/propagation distortion;
2. make the 1--10 EeV behavior materially less catastrophic than contact Case 6.

Only masses that pass this shape/transport gate proceed to:

- a joint propagation+detector fold;
- a current collider/low-energy/UHE coupling audit;
- a more explicit UV completion of the dimension-5 neutrino-Z' vertex if needed.

If no mass passes at the Case-6 low-Q normalization, the next controlled test is to vary the low-Q normalization. If no normalization yields a useful window, then this derivative-Z' parent does not solve the structural opacity problem.
