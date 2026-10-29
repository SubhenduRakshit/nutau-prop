# v0.27 — Active gluon dimension-8 momentum-growing screen

## Operator

The next candidate is the class-14 SMEFT operator

\[
Q_{l^2G^2D}^{33}=i(\bar L_3\gamma^\mu\overleftrightarrow D^\nu L_3)
G^A_{\mu\rho}G^{A\rho}{}_{\nu},
\]

with

\[
\mathcal L_{\rm NP}=\frac{C_G}{\Lambda^4}Q_{l^2G^2D}^{33}.
\]

This is a genuine dimension-8 SMEFT basis element (Murphy, JHEP 10 (2020) 174, arXiv:2005.00059). For p=r=3 it contains both the tau-neutrino and charged-tau pieces required by SU(2)L.

## Partonic process

The relevant UHE process is

\[
\nu_\tau g\to\nu_\tau g.
\]

For a fixed left-handed neutrino and an unpolarized initial gluon, direct contraction of the two field strengths gives

\[
\overline{|\mathcal M|^2}
=2\frac{C_G^2\hat s^4}{\Lambda^8}
(1-y)\left[1+(1-y)^2\right].
\]

Hence

\[
\boxed{
\frac{d\hat\sigma}{dy}
=\frac{C_G^2\hat s^3}{8\pi\Lambda^8}
(1-y)\left[1+(1-y)^2\right]
}
\]

and

\[
\int_0^1dy\,(1-y)[1+(1-y)^2]=\frac34.
\]

There is no tree-level SM elastic `nu g -> nu g` amplitude, so the first UHE screen contains only the positive pure-BSM contribution. Loop-level elastic SM interference is a later correction if this branch survives.

## Helicity partial waves

For incoming neutrino helicity h=-1/2, the nonzero gluon-helicity amplitudes are

\[
\mathcal M_{\lambda_g=-1\to-1}
=-2C_G\frac{\hat s^2}{\Lambda^4}\cos\frac\theta2,
\]

\[
\mathcal M_{\lambda_g=+1\to+1}
=-2C_G\frac{\hat s^2}{\Lambda^4}\cos^3\frac\theta2.
\]

They are pure Jacob-Wick partial waves:

\[
J=\frac12:\quad
a_{1/2}=-\frac{C_G}{16\pi}
\left(\frac{\hat s}{\Lambda^2}\right)^2,
\]

\[
J=\frac32:\quad
a_{3/2}=-\frac{C_G}{32\pi}
\left(\frac{\hat s}{\Lambda^2}\right)^2.
\]

Thus the strict condition \(|\mathrm{Re}\,a_J|\le1/2\) implies

\[
\boxed{
|C_G|\left(\frac{\hat s}{\Lambda^2}\right)^2\le8\pi
}.
\]

If unitarity is required through \(\hat s=\Lambda^2\),

\[
|C_G|\le8\pi=25.1327.
\]

## Positivity

The published forward-limit positivity condition for this Murphy-basis coefficient is

\[
-c_{l^2G^2D}\ge0,
\]

therefore

\[
\boxed{C_G\le0}.
\]

Reference: M. Chala and X. Li, Phys. Rev. D 109, 065015 (2024), Appendix B.

The strict benchmark is consequently

\[
\boxed{C_G=-8\pi},\qquad \Lambda=4~\mathrm{TeV}.
\]

## Nucleon fold

Let \(s=2m_NE_\nu\), \(\hat s=xs\), and \(Q^2=xys\). With LHAPDF returning \(xg(x,Q^2)\),

\[
\frac{d^2\sigma_G}{dx\,dy}
=\frac{C_G^2s^3x^2}{8\pi\Lambda^8}
[xg(x,Q^2)]
(1-y)[1+(1-y)^2],
\]

subject to

\[
Q>Q_{\min},\qquad \hat s<\xi\Lambda^2.
\]

The numerical screen must evaluate \(\xi=1,0.5,0.25\) at 1, 10, 100 and 300 PeV using the same CT18NLO and dynamic-Q conventions as the closed quark branch.

## Gate

1. Positivity: PASS, requires C_G<=0.
2. Exact normalization: derived above.
3. Partial-wave unitarity: PASS, strict ceiling |C_G|<=8pi for `unitary-to-Lambda`.
4. Exact CT18NLO UHE fold: NEXT COMPUTATION.
5. If max allowed Delta sigma_NC/sigma_NC remains below 5%, close the branch without Earth propagation.
6. If it exceeds 5%, then perform the collider gate, dominated by the SU(2)L-related `gg -> tau+ tau-` contact interaction and missing-energy/monojet channels, before any Earth propagation.

The gluon branch must therefore not be prejudged from the quark result: the coefficient ceiling is three times larger and the gluon PDF can be substantially larger than the charm/strange sea in the relevant x range.
