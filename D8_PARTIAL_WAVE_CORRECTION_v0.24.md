# nutau-prop v0.24 — d=8 partial-wave correction and certification reorder

Date: 2026-10-29
Status: analytic correction to v0.23; preserves v0.23 as history

## Executive correction

The active d=8 two-derivative semileptonic cross sections in `D8_UHE_SCREEN_v0.23.md` are retained, but the leading partial-wave statement is corrected.

For

\[
\mathcal L_{\rm NP}=\frac{C}{\Lambda^4}
(\bar\nu_{\tau L}\gamma_\mu i\overleftrightarrow D_\nu\nu_{\tau L})
(\bar q_R\gamma^\mu i\overleftrightarrow D^\nu q_R),
\]

with

\[
\frac{d\hat\sigma_q}{dy}=\frac{C^2\hat s^3}{8\pi\Lambda^8}(2-y)^2(1-y)^2,
\qquad
\frac{d\hat\sigma_{\bar q}}{dy}=\frac{C^2\hat s^3}{8\pi\Lambda^8}(2-y)^2,
\]

the helicity amplitude for the quark channel is

\[
\mathcal M_q=\sqrt2\,C\frac{\hat s^2}{\Lambda^4}
\frac{(3+\cos\theta)(1+\cos\theta)}{4}.
\]

Using

\[
d^1_{11}=\frac{1+\cos\theta}{2},\qquad
 d^2_{11}=\frac{(1+\cos\theta)(2\cos\theta-1)}{2},
\]

one has the exact decomposition

\[
\frac{(3+c)(1+c)}4=\frac74 d^1_{11}+\frac14 d^2_{11}.
\]

With the convention

\[
\mathcal M=16\pi\sum_J(2J+1)a_J d^J_{11},
\]

this gives

\[
a_1=\frac{7\sqrt2}{192\pi}C\left(\frac{\hat s}{\Lambda^2}\right)^2
=0.016412\,C\left(\frac{\hat s}{\Lambda^2}\right)^2,
\]

\[
a_2=\frac{\sqrt2}{320\pi}C\left(\frac{\hat s}{\Lambda^2}\right)^2
=0.0014067\,C\left(\frac{\hat s}{\Lambda^2}\right)^2.
\]

Therefore **J=1, not J=2, is the binding partial wave**, by a factor 11.7 in amplitude.

At \(\hat s=\Lambda^2\):

- strict elastic criterion \(|\mathrm{Re}\,a_1|\le1/2\): \(C\le30.46\);
- loose criterion \(|a_1|\le1\): \(C\le60.93\).

Alioli et al.'s statement that this operator class contains a new l=2 amplitude component remains correct; it does not imply that the l=2 coefficient dominates the unitarity bound after crossing to \(\nu q\to\nu q\).

## Antiquark enhancement

The integrated y weights are

\[
I_q=\int_0^1dy\,(2-y)^2(1-y)^2=\frac{31}{30},
\]

\[
I_{\bar q}=\int_0^1dy\,(2-y)^2=\frac73.
\]

Hence

\[
\frac{I_{\bar q}}{I_q}=\frac{70}{31}=2.258.
\]

For a symmetric sea q≈qbar, the total q+qbar contribution is larger than a quark-only estimate by

\[
\frac{I_q+I_{\bar q}}{I_q}=\frac{101}{31}=3.258.
\]

This strengthens the motivation for charm and strange sea operators.

## Important subtlety: unitarity versus xi truncation

The earlier statement that xi=1/4 is `excluded outright` is only true if Lambda is treated as the physical UV cutoff and one demands unitarity all the way to shat=Lambda^2 even though the numerical fold is truncated earlier.

If the retained EFT domain itself ends at

\[
\hat s_{\max}=\xi\Lambda^2,
\]

then the strict partial-wave condition is

\[
C\,\xi^2\le30.46,
\]

or

\[
C_{\max}(\xi)=30.46/\xi^2.
\]

Thus

- xi=1: Cmax=30.5;
- xi=1/2: Cmax=121.8;
- xi=1/4: Cmax=487.4.

The cleanest parameterization is therefore not separate (C,Lambda,xi), but the dimensionful coefficient

\[
G_8=C/\Lambda^4
\]

and an independently stated maximum trusted invariant mass \(M_{\rm cut}\), with

\[
G_8 M_{\rm cut}^4\le30.46
\]

for the strict partial-wave criterion.

For UV interpretations in which Lambda itself is the new-physics scale, it is still reasonable to require C<=30.46 at shat=Lambda^2 and use xi only as a truncation diagnostic. Both conventions must be stated explicitly; they answer slightly different questions.

## UHE consequence under the stronger `unitary up to Lambda` convention

Using the v0.23 CT18 pure-NP screen at Lambda=4 TeV and C=30:

### xi=1

Approximate pure-NP modification:

- 100 PeV c_R: 4.3% of SM NC, 1.35% of total SM;
- 100 PeV s_R: 5.9% of SM NC, 1.84% of total SM;
- 300 PeV c_R: 9.0% of SM NC, 2.85% of total SM;
- 300 PeV s_R: 11.3% of SM NC, 3.57% of total SM.

Thus the 5–10% **total-cross-section** headline does not survive strict unitarity-to-Lambda, but the NC kernel can still change at the 5–11% level. Since Earth propagation is sensitive to the NC redistribution kernel rather than only to the inclusive total cross section, the original >=5% spectral/zenith distortion goal is not yet ruled out.

### xi=1/2 with the same C<=30 convention

The pure-NP contribution is only about 1.2–2.3% of NC (roughly 0.4–0.7% of total).

### xi=1/4 with the same C<=30 convention

The pure-NP contribution is at the few-per-mille level of NC and negligible in the total cross section.

Therefore, if Lambda is required to be a genuine UV scale with the EFT unitary up to Lambda, the observable effect is endpoint-dominated and loses robustness under conservative truncation.

## Endpoint dominance

The x^3 weighting of the d=8 operator makes the cutoff-edge problem more severe than in the old d=7 branch. The external audit reports approximately:

- 97.5% of sigma_NP from shat>0.25 Lambda^2;
- 84% from shat>0.5 Lambda^2;
- fold-weighted <(shat/Lambda^2)^2>≈0.57.

These numbers should be independently reproduced with the exact dynamic-Q CT18/NNPDF integration before publication, but the qualitative conclusion is accepted: truncation uncertainty is potentially O(1) near the parameter region that gives the largest UHE signal.

## Collider ordering

The exact tau-tau recast remains necessary, but it is no longer the first gate.

The certification order is now:

1. **J=1 partial-wave unitarity** in the exact operator convention;
2. **exact normalization/crossing check** tying UHE and collider amplitudes to one convention;
3. **HighPT tau+tau- recast** for [O8,lu partial2]3322 and [O8,ld partial2]3322;
4. monojet/dineutrino screen;
5. positivity-sign mapping;
6. IceCube 60 TeV–10 PeV energy-resolved tau-specific NC-like reweight;
7. dynamic-Q CT18/NNPDF UHE fold and endpoint-support diagnostics;
8. Earth propagation only for survivors.

The external audit quotes approximate collider-equivalent C limits that can be stronger for strange and weaker for charm, but those values are not promoted to certified project bounds until their exact likelihood/operator normalization is reproduced.

## Project verdict after this correction

The active d=8 charm/strange twist-2 branch is **not dead**, but its most attractive region is much narrower than v0.23 suggested if Lambda is interpreted as a physical cutoff and the EFT must remain unitary to Lambda.

The clean surviving possibility is:

- NC cross-section modifications of order 5–10% near 100–300 PeV at xi≈1;
- total-cross-section modifications only at the few-percent level;
- possible >5% spectral/zenith distortions after Earth propagation still need to be tested because NC energy redistribution can amplify a modest inclusive change.

The endpoint sensitivity means a resolved spin-2 mediator is now a particularly well-motivated UV-completion benchmark, not an optional afterthought.
