# PKR width / branching-ratio / dijet prescreen — v0.45

## Motivation
The v0.44 detector result is a practical diffuse-sensitivity failure, not an official detector exclusion. Before any further detector refinement, the conditional M_Z'=50 GeV, g_q=0.03 PKR benchmark must pass basic model certification.

The first cheap certification step is the explicit mediator width and branching ratio. This also corrects the interpretation of the CMS low-mass dijet bound, which assumes BR(Z'->qq)=1.

## Width formulas
For the universal vector quark interaction

`g_q Z'_mu qbar gamma^mu q`,

per open quark flavor:

Gamma(Z'->q qbar) = Nc g_q^2 M/(12 pi) (1+2 r_q) sqrt(1-4 r_q).

For one massless Weyl neutrino pair, with `a=c_i/Lambda=kappa_i/g_q`, the PKR dim-5 vertices give

- Eq. (4.4) derivative: Gamma_nu = a^2 M^3/(48 pi),
- Eq. (4.3) transition dipole: Gamma_nu = a^2 M^3/(12 pi).

Because charge-conjugation / h.c. bookkeeping can differ by a factor two depending on how the invisible channel is counted, v0.45 reports both one-pair and operator+h.c. widths. A conclusion stable under this bracket is treated as robust for the prescreen.

## Dijet logic
CMS EXO-24-007 (2026) quotes observed limits on universal `g_q` ranging from 0.03 to 0.13 over 50--250 GeV for a spin-1 resonance with BR(Z'->qq)=1.

For a narrow resonance with invisible decays, the resonance yield scales approximately as

sigma(pp->Z') BR_qq ~ g_q^2 BR_qq.

Therefore define the BR-corrected equivalent coupling

`g_q,eff = g_q sqrt(BR_qq)`.

v0.45 compares this with 0.03, the strongest observed limit anywhere in the quoted mass interval. If `g_q,eff < 0.03`, the benchmark passes a deliberately conservative dijet prescreen without requiring the exact 50-GeV HEPData point.

## What this does not certify
This step does **not** replace:

- an operator-specific Run-2 monojet recast;
- kinetic-mixing / transition-moment constraints;
- gauge completion and anomaly cancellation;
- flavor multiplicity choices beyond the explicit width bracket.

If v0.45 passes, the next decisive calculation is the modern operator-specific monojet recast.
