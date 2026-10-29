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

For a fixed left-handed neutrino and an unpolarized initial gluon,

\[
\overline{|\mathcal M|^2}
=2\frac{C_G^2\hat s^4}{\Lambda^8}
(1-y)\left[1+(1-y)^2\right],
\]

hence

\[
\boxed{
\frac{d\hat\sigma}{dy}
=\frac{C_G^2\hat s^3}{8\pi\Lambda^8}
(1-y)\left[1+(1-y)^2\right]
}.
\]

The angular integral is 3/4. There is no tree-level SM elastic `nu g -> nu g` amplitude, so the first UHE screen contains only the positive pure-BSM contribution.

## Helicity partial waves

The nonzero gluon-helicity amplitudes are

\[
\mathcal M_{-\to-}=-2C_G\frac{\hat s^2}{\Lambda^4}\cos\frac\theta2,
\qquad
\mathcal M_{+\to+}=-2C_G\frac{\hat s^2}{\Lambda^4}\cos^3\frac\theta2.
\]

They give

\[
a_{1/2}=-\frac{C_G}{16\pi}\rho^2,
\qquad
a_{3/2}=-\frac{C_G}{32\pi}\rho^2,
\qquad \rho=\frac{\hat s}{\Lambda^2}.
\]

Thus strict \(|\mathrm{Re}\,a_J|\le1/2\) implies

\[
\boxed{|C_G|\rho^2\le8\pi}.
\]

Requiring unitarity through \(\hat s=\Lambda^2\) gives \(|C_G|\le8\pi\).

## Positivity

Published forward-limit positivity in the Murphy basis gives

\[
-c_{l^2G^2D}\ge0,
\qquad\Rightarrow\qquad
\boxed{C_G\le0}.
\]

Reference: M. Chala and X. Li, Phys. Rev. D 109, 065015 (2024), Appendix B.

The strict benchmark is

\[
\boxed{C_G=-8\pi},\qquad \Lambda=4~\mathrm{TeV}.
\]

## Exact CT18NLO dynamic-Q UHE fold

The exact fold was executed in GitHub Actions using CT18NLO, dynamic \(Q^2=xys\), \(Q_{\min}=1.3\) GeV, and the same SM reference as the closed quark branch. Results are in `results/v0.27/d8_gluon_uhe_v0_27.csv`.

For `unitary-to-Lambda` and \(C_G=-8\pi\):

| xi | E_nu | Delta sigma_NC / sigma_NC | Delta sigma_tot / sigma_tot | support rho>0.25 | support rho>0.5 |
|---:|---:|---:|---:|---:|---:|
| 1.00 | 10 PeV | 0.255% | 0.070% | 49.3% | 9.6% |
| 1.00 | 100 PeV | **12.34%** | 3.46% | 93.1% | 71.8% |
| 1.00 | 300 PeV | **25.96%** | **7.32%** | 95.1% | 76.9% |
| 0.50 | 100 PeV | 3.11% | 0.871% | 75.6% | 0 |
| 0.50 | 300 PeV | **5.25%** | 1.48% | 78.7% | 0 |
| 0.25 | 100 PeV | 0.755% | 0.211% | 0 | 0 |
| 0.25 | 300 PeV | 1.11% | 0.314% | 0 | 0 |

The 5% NC target is therefore passed for xi=1 at 100 and 300 PeV, and marginally passed for xi=0.5 at 300 PeV. It fails for xi=0.25.

## Interpretation

Unlike the `cR/sR` contact branch, this operator cannot be closed at the UHE screen. The larger gluon luminosity and the weaker strict coefficient ceiling make an O(10%) NC modification possible.

However, the xi=1 signal is strongly endpoint dominated: 93--95% of the BSM cross section at 100--300 PeV comes from rho>0.25, and 72--77% from rho>0.5. Consequently the contact-EFT interpretation is sensitive to the UV completion/cutoff treatment.

## Gate status

1. Exact normalization: **PASS**.
2. Positivity: **PASS**, C_G<=0.
3. Partial-wave unitarity: **PASS**, strict |C_G|<=8pi for `unitary-to-Lambda`.
4. Exact CT18NLO UHE fold: **PASS phenomenological 5% gate for xi=1; marginal at xi=0.5, 300 PeV**.
5. Collider gate: **NEXT**. The same SU(2)L operator contains the charged-tau interaction and generates `gg -> tau+ tau-`; high-mass ditau data are therefore the primary constraint to establish. Missing-energy/monojet channels are secondary checks.
6. Earth propagation: **DO NOT RUN YET**. It is justified only if the collider/UV-validity gate leaves a coefficient large enough to retain >=5% UHE NC modification.

This branch is therefore **alive but not yet certified**.
