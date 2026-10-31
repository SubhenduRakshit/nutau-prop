# PKR Eq. (4.4) invisible-width normalization — v0.73

## Scope

This note resolves the factor-of-two width ambiguity left open in v0.45–v0.72 for the resolved derivative interaction of Pandey–Karmakar–Rakshit (PKR), JHEP 11 (2019) 046, arXiv:1907.07700.

PKR Eq. (4.4) is

\[
\mathcal L \supset a\,(\bar\nu^c i\overleftrightarrow{\partial^\mu}\nu)Z'_\mu,
\qquad a\equiv c^{(2)}/\Lambda,
\]

plus the quark current. The printed equation does not explicitly show `+ h.c.`.

The present project interprets the neutrino vertex as an active-to-sterile/off-diagonal transition with two distinct effectively massless Weyl states. The question is therefore whether the physical invisible width is the previously tabulated `one_Weyl_pair` value or twice that value after Hermitian completion.

## Hermiticity

The operator

\[
O=(\bar\nu^c i\overleftrightarrow{\partial^\mu}\nu)Z'_\mu
\]

is not, by itself, a complete Hermitian interaction for the active-to-sterile Weyl realization. The physical Lagrangian is

\[
\mathcal L_{\rm int}=a O+a^* O^\dagger.
\]

This is also the standard convention in the general neutrino EFT basis: non-Hermitian operators enter the effective Lagrangian together with their Hermitian conjugates.

For distinct Weyl states, `O` and `O^dagger` generate conjugate invisible final states. In the massless limit they are orthogonal channels and have equal partial widths, so their rates add.

## One-channel width

For

\[
Z'(P)\to \nu_i(p_1)\,N(p_2),
\]

the derivative vertex gives, up to an irrelevant phase,

\[
\mathcal M=a\,\epsilon_\mu(P)(p_1-p_2)^\mu\langle12\rangle.
\]

For massless final states,

\[
|\langle12\rangle|^2=2p_1\!\cdot p_2=M^2,
\]

and the polarization average gives

\[
\frac13\sum_\lambda
|\epsilon^{(\lambda)}\!\cdot(p_1-p_2)|^2
=\frac{M^2}{3}.
\]

Hence

\[
\overline{|\mathcal M|^2}=\frac{a^2M^4}{3}.
\]

For two distinct massless particles,

\[
\Gamma_{\rm one}
=\frac{\overline{|\mathcal M|^2}}{16\pi M}
=\boxed{\frac{a^2M^3}{48\pi}}.
\]

This is exactly the old `one_Weyl_pair` entry in v0.45.

## Hermitian-completed invisible width

The Hermitian conjugate gives the conjugate helicity/charge channel with the same rate. Therefore, for the active-to-sterile/off-diagonal realization used throughout the present project,

\[
\boxed{
\Gamma_{\rm inv}
=2\Gamma_{\rm one}
=\frac{a^2M^3}{24\pi}.
}
\]

Thus the physically appropriate branch is the old `operator_plus_hc` branch.

This conclusion is scoped to the distinct-state active-to-sterile realization. A separately normalized diagonal self-conjugate Majorana operator would require its own identical-particle/operator-normalization convention and should not be silently substituted for this model.

## Numerical benchmark at M = 50 GeV

For the certified benchmark

- M = 50 GeV,
- g_q = 0.03,
- kappa_2 = 1.068919824203e-4 GeV^-1,
- a = kappa_2/g_q,

v0.45 gives

\[
\Gamma_q=0.0179038627\ {\rm GeV},
\]

\[
\Gamma_{\rm inv}^{\rm one}=0.0105236563\ {\rm GeV},
\]

while the Hermitian-completed value is

\[
\boxed{\Gamma_{\rm inv}=0.0210473126\ {\rm GeV}}.
\]

Therefore

\[
\Gamma_{\rm tot}=0.0389511754\ {\rm GeV},
\qquad
\Gamma_{\rm tot}/M=7.7902\times10^{-4},
\]

and

\[
\boxed{BR_{\rm inv}=0.5403511563}.
\]

The mediator remains safely narrow.

## Collider consequence from v0.72

The validated CMS EXO-20-004 covariance fit for the lower-width sample gave

\[
\mu_{95}=1.1784367171,
\qquad r_{\rm one}=1/\mu_{95}=0.848582.
\]

In the narrow-width limit the monojet signal shape is fixed by on-shell Z' production and the total invisible momentum equals the mediator momentum. Changing only the invisible width convention therefore rescales the signal normalization by the branching-ratio ratio.

Using the Hermitian-completed branch,

\[
r_{O+O^\dagger}
=r_{\rm one}\,
\frac{0.5403511563}{0.3701925697}
=1.238632.
\]

Under the pre-declared v0.51/v0.72 robustness thresholds (`r>1.2` exclusion, `r<0.8` survival), the benchmark crosses the exclusion threshold:

\[
\boxed{
\text{Eq. (4.4), M=50 GeV, g_q=0.03 benchmark: collider-monojet FAIL}
}
\]

for the active-to-sterile/off-diagonal Hermitian completion.

The margin above the robustness threshold is modest. This verdict should therefore be quoted specifically for this certified benchmark and NWA model, not as a model-independent exclusion of every derivative-neutrino interaction.

## Important operator caveat

PKR Eq. (4.4) is not the same object as the standard lepton-number-conserving derivative operator Q8 in the Altmannshofer–Tammaro–Zupan low-energy EFT basis, which is written as

\[
(\bar\nu_\beta i\overleftrightarrow{\partial}_\mu P_L\nu_\alpha)
(\bar f\gamma^\mu f).
\]

PKR instead prints a charge-conjugated neutrino bilinear. Consequently, one must not import Q8 flavor/hermiticity statements into Eq. (4.4) without specifying the neutrino-field realization. The present conclusion uses only the explicit active-to-sterile/off-diagonal realization adopted in this repository.

## Decision

1. The v0.45 width ambiguity is closed for the project model.
2. Use `operator_plus_hc` as the physical Eq. (4.4) width convention.
3. Do not generate more collider events for this benchmark.
4. The M=50 GeV, g_q=0.03 Eq. (4.4) benchmark fails the certified monojet gate at the pre-declared robustness threshold.
5. The broader momentum-growing interaction survey remains open; this does not establish a no-go theorem for other operators or mediator structures.
