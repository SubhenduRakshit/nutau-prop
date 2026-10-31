# PKR Eq. (4.3) transition-dipole width certification — v0.82

> **v0.86 correction:** the widths below survive for the specified distinct-state
> convention, but the cross-section factor-two claim is superseded. Exact
> scattering traces give factor four, and the inherited absolute DIS kernel
> also mismatched these widths. See `PKR_NORMALIZATION_RESOLUTION_v0.86.md`.
> The original text is retained as history, not an accepted DIS normalization.

Date: 2026-10-31
Status: analytic normalization certification

## Operator and project realization

We take the resolved transition-dipole interaction

\[
\mathcal L_{4.3}=a\,(\bar\nu_i^c\sigma_{\mu\nu}P_L\nu_j) Z'^{\mu\nu}+\mathrm{h.c.},
\qquad i\neq j,
\]

with
\[
a\equiv c^{(1)}/\Lambda,\qquad \kappa_1=g_q a.
\]

The Earth-attenuation realization is active-to-sterile/off-diagonal, e.g. \(\nu_\tau\to N\), not a diagonal active NC interaction.

## Exact on-shell relation to Eq. (4.4)

For massless external fermions and \(\sigma^{\mu\nu}=i[\gamma^\mu,\gamma^\nu]/2\),

\[
q_\mu\,\bar u(p_1)\sigma^{\mu\nu}P_Lv(p_2)
=-i(p_1-p_2)^\nu\,\bar u(p_1)P_Lv(p_2),
\qquad q=p_1+p_2.
\]

Since
\[
Z'_{\mu\nu}\to -i(q_\mu\epsilon_\nu-q_\nu\epsilon_\mu)
\]
and \(\sigma^{\mu\nu}\) is antisymmetric, the Eq. (4.3) decay amplitude is, up to an irrelevant phase,

\[
\mathcal M_{4.3}=2a\,\epsilon\!\cdot\!(p_1-p_2)\,\bar uP_Lv.
\]

The certified Eq. (4.4) one-channel amplitude has the same structure without the factor 2. Therefore

\[
\overline{|\mathcal M_{4.3}|^2}=4\,\overline{|\mathcal M_{4.4}|^2}
=\frac{4a^2M^4}{3}.
\]

Using the massless two-body phase space,

\[
\boxed{\Gamma_{4.3}^{\rm one}=\frac{a^2M^3}{12\pi}}.
\]

For the physical project convention \(aO+a^*O^\dagger\), the h.c. creates the distinct conjugate invisible final state, so

\[
\boxed{\Gamma_{4.3}^{O+O^\dagger}=\frac{a^2M^3}{6\pi}}.
\]

Thus the old v0.45 `operator_plus_hc` numerical width is the physical convention for the current active-to-sterile/off-diagonal realization.

## Cross-section relation

The already-derived resolved DIS relation remains

\[
\frac{d^2\sigma_{4.3}}{dxdy}=2\frac{d^2\sigma_{4.4}}{dxdy}
\]

at equal numerical \(M\), \(g_q\), and \(a\), with the same validity boundary

\[
\sqrt{\hat s}<\frac{1}{a}=\frac{g_q}{\kappa_1}.
\]

Therefore the corrected Eq. (4.4) base envelopes can be reused exactly for Eq. (4.3) by multiplying the UHE ratio by 2; no new CT18 integration is required.

## Certification decision

- Exact one-channel dipole width: PASS.
- Physical O+h.c. invisible width: PASS, \(a^2M^3/(6\pi)\).
- Coupling separation: retain \((M,g_q,a)\), never only \(\kappa_1\).
- Old v0.45 one-Weyl row: historical bracket only.
- Old v0.45 operator+h.c. row: physical normalization for the present project realization.

Next: build the complete Eq. (4.3) >=5%-at-100-PeV full-a envelope, apply width/vertex/dijet screens pointwise, then reuse the existing mass-dependent monojet production/acceptance vectors for a no-new-events collider screen.
