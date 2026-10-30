# nutau-prop v0.35 — self-consistent resolved derivative-Z' gate

Date: 2026-10-30
Status: quantitative size gate + historical-constraint audit; modern likelihood recast not yet certified.

## Executive decision

The resolved derivative active-to-sterile mediator remains the leading **non-resonant** candidate after enforcing the requested coupling-dependent validity domain inside the PDF integral. The controlled CT18 scan still finds a 5% solution for mediator masses through roughly 500 GeV, but not at 1–2 TeV within the scanned perturbative domain.

This does **not** yet certify an allowed model: the historical PKR constraint curve must be replaced by a modern operator-specific ATLAS/IceCube recast, and the product coefficient kappa must be unpacked into c2, gq and Lambda before mediator-width/direct-search statements are made.

## Interaction and conventions

The historical PKR Eq. (4.4) structure is

\[
\mathcal L\supset {c_2\over\Lambda}
(\bar\nu^c i\overleftrightarrow{\partial_\mu}\nu)Z'^\mu
+g_q(\bar q\gamma_\mu q)Z'^\mu,
\qquad
\kappa\equiv{c_2g_q\over\Lambda}.
\]

For the project phenomenology the useful resolved variant is interpreted as an active-to-sterile transition, \(\nu_\tau q\to N q\). This is not silently identified with a diagonal active neutral current; the external-state/gauge-completion issue is kept explicit.

The resolved DIS kernel used here is

\[
{d^2\sigma\over dxdy}={4\over\pi}\,
\kappa^2 x s^2 y(1-y)
{\sum_q x f_q(x,Q^2)\over(M_{Z'}^2+Q^2)^2}.
\]

## Crucial v0.35 change: moving validity cut

For **every** trial kappa the partonic domain is recomputed with

\[
\boxed{\sqrt{\hat s}<1/\kappa}.
\]

No fixed-cut quadratic rescaling is used. Increasing kappa therefore both strengthens the matrix element and shrinks its allowed partonic phase space. The code searches directly for the first crossing

\[
\sigma_{\rm NP}/\sigma_{\rm SM}=0.05
\]

using a logarithmic kappa scan followed by log-space bisection, while reevaluating the PDF integral at every trial point.

Code: `scripts/resolved_derivative_selfconsistent_v0_35.py`.

## CT18NLO decision scan

The first completed self-consistent scan gives:

| M_Z' [GeV] | kappa_5%(100 PeV) [GeV^-1] | kappa_5%(300 PeV) [GeV^-1] |
|---:|---:|---:|
| 5 | 3.0943e-5 | 3.1251e-5 |
| 40 | 5.4254e-5 | 5.1974e-5 |
| 100 | 7.9706e-5 | 7.4021e-5 |
| 220 | 1.2008e-4 | 1.0898e-4 |
| 500 | 2.1643e-4 | 1.9065e-4 |
| 1000 | no first crossing | no first crossing |
| 2000 | no first crossing | no first crossing |

The `none` entries mean no first 5% crossing was found within the scanned kappa range while the coupling-dependent validity condition was imposed. They are not zero-cross-section statements.

Numerical result commit: `eff8d54542359012e7863a9725f8baa0cbaf0038`.

## Comparison with the historical PKR baseline

Verified directly from PKR text:

- at M_Z'=5 GeV, kappa_LHC < 3.3e-3 GeV^-1;
- at M_Z'=5 GeV, kappa_IceCube < 2.5e-4 GeV^-1;
- Fig. 2(c) states that the LHC bound overtakes IceCube approximately in the 40–220 GeV interval;
- Fig. 4 gives the explicit 5-GeV neutrino-nucleon cross-section comparison.

At 5 GeV the new self-consistent 5% requirement,

\[
\kappa_{5\%}\simeq3.1\times10^{-5}\ {\rm GeV}^{-1},
\]

is about a factor eight below the historical IceCube upper bound and over two orders of magnitude below the historical LHC upper bound. Therefore the historical 5-GeV constraints do **not** close the 5% window.

A full numerical digitization of every point of Fig. 2(c) is not claimed here. The text anchors and crossover interval are verified, while the complete curve must be reconstructed from the figure before using it as quantitative data at intermediate masses. Do not invent interpolated PKR limits.

## Modern constraint update: what is and is not established

### LHC

The appropriate current public monojet baseline is the ATLAS full Run-2 13-TeV, 139-fb^-1 jets+missing-momentum search. It is substantially more powerful than the dataset used in PKR. However, no published bound on this exact derivative active-to-sterile operator supplies a ready-made kappa_max(M) curve.

Therefore any simple luminosity rescaling is only a **screening estimate**. A certification-quality update needs generation of pp -> j + missing energy with the resolved propagator and derivative vertex, ATLAS signal-region acceptance, detector/reconstruction treatment or validated efficiency maps, and the published likelihood/bin information.

### IceCube

Modern IceCube all-flavour neutrino-nucleon cross-section measurements extend through the multi-PeV region (notably the 7.5-year HESE analysis up to about 10 PeV), but likewise do not publish an operator-specific kappa_max(M) curve. A modern update requires folding this interaction through the energy/zenith-dependent event likelihood rather than treating the measurement as a single normalization number.

Hence the modern constraint curve is **identified but not yet reproduced**. The historical PKR curve remains a baseline, not a current certificate.

## Keep the model parameters separate

The one-parameter kappa scan is a cross-section diagnostic only:

\[
\kappa={c_2g_q\over\Lambda}.
\]

A physical model must retain independently:

1. c2/Lambda — derivative neutrino-transition strength and its UV validity;
2. g_q — quark coupling and hadron-collider/direct-Z' constraints;
3. Lambda — UV scale;
4. M_Z' and the total width Gamma_Z';
5. branching fractions and any additional Z' couplings;
6. direct resonance/dijet/low-mass mediator constraints.

In particular, Gamma_Z' cannot be inferred from kappa alone. The requested sqrt(shat)<1/kappa cutoff is retained as the conservative product-level portfolio test, not presented as a complete UV theorem.

## Current portfolio decision

The v0.35 size calculation strengthens rather than removes the low/intermediate-mass window:

\[
\boxed{M_{Z'}\sim5\text{--}500\ {\rm GeV}\quad\text{can reach a 5% UHE effect in the self-consistent size fold.}}
\]

But present status remains

\[
\boxed{\text{SIZE PASS / MODERN-CONSTRAINT CONDITIONAL / UV CONDITIONAL}.}
\]

The most promising candidate therefore remains the resolved derivative active-to-sterile mediator **for the non-resonant programme**. If a modern monojet+IceCube recast drives kappa_max(M) below the tabled kappa_5%(M) throughout 5–500 GeV, close this branch and promote the strange-aligned tilde-R2 leptoquark resonance as the lead model.

No Earth-transport calculation is justified at this stage.
