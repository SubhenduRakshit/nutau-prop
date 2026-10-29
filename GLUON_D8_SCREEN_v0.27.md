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

### QCD-coupling convention

In the canonical Murphy basis written above there is **no explicit factor of** `g_s` or `alpha_s` multiplying the operator. For the two-external-gluon contact vertex, both field strengths contribute their linear pieces,

\[
G^A_{\mu\nu}=\partial_\mu G^A_\nu-\partial_\nu G^A_\mu
+g_s f^{ABC}G^B_\mu G^C_\nu,
\]

so the `nu g -> nu g` and crossed `gg -> tau tau` two-gluon amplitudes from this operator contain no explicit `g_s` at tree level. The nonlinear term is relevant to vertices with additional gluons.

One may instead define a rescaled operator containing an explicit `g_s^2`; then the corresponding Wilson coefficient is inversely rescaled. In that convention a statement such as "alpha_s cancels under unitarity saturation" is physically equivalent. The project will quote all numerical bounds in the canonical Murphy convention to avoid convention-dependent alpha_s factors.

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

The bar includes the average over the two initial gluon helicities. The color average cancels the final-color sum because the two-gluon vertex is proportional to `delta^{AB}`. The angular integral is 3/4. There is no tree-level SM elastic `nu g -> nu g` amplitude, so the first UHE screen contains only the positive pure-BSM contribution.

## Helicity partial waves

The nonzero gluon-helicity amplitudes are

\[
\mathcal M_{-\to-}=-2C_G\frac{\hat s^2}{\Lambda^4}\cos\frac\theta2,
\qquad
\mathcal M_{+\to+}=-2C_G\frac{\hat s^2}{\Lambda^4}\cos^3\frac\theta2.
\]

Averaging their squared magnitudes over the two initial gluon helicities reproduces the unpolarized matrix element above, providing the direct spin-average normalization check analogous to the antiquark check in the closed quark branch.

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

The comparison with the quark operator does **not** reduce to gluon-versus-charm luminosity alone. The binding partial wave is different (`J=1/2` here versus antiquark `J=0` there), so the unitarity-saturated coefficient ceilings and angular kernels differ. In the common canonical normalization `|C_G|_max/|C_q|_max=3`, while the gluon angular integral is smaller. These effects must be retained in any scaling estimate; the exact PDF fold is preferred.

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

Unlike the `cR/sR` contact branch, this operator cannot be closed at the UHE screen. The gluon luminosity together with the different partial-wave ceiling makes an O(10%) NC modification possible in the most permissive contact-EFT treatment.

However, the xi=1 signal is strongly endpoint dominated: 93--95% of the BSM cross section at 100--300 PeV comes from rho>0.25, and 72--77% from rho>0.5. Consequently the contact-EFT interpretation is sensitive to the UV completion/cutoff treatment.

## Gate status

1. Exact normalization: **PASS internally**, including explicit initial-gluon polarization averaging; retain an independent Feynman-rule/crossing check before calling the branch fully certified.
2. Positivity: **PASS**, C_G<=0.
3. Partial-wave unitarity: **PASS**, strict |C_G|<=8pi for `unitary-to-Lambda`.
4. Exact CT18NLO UHE fold: **PASS phenomenological 5% gate for xi=1; marginal at xi=0.5, 300 PeV**.
5. Collider gate: **NEXT**. The same SU(2)L operator contains the charged-tau interaction and generates `gg -> tau+ tau-`; high-mass ditau data are therefore the primary constraint to establish. Missing-energy/monojet channels are secondary checks.
6. Earth propagation: **DO NOT RUN YET**. It is justified only if the collider/UV-validity gate leaves a coefficient large enough to retain >=5% UHE NC modification.

This branch is therefore **alive but not yet certified**.

## Resolved spin-2 follow-up if the contact gluon branch fails

A resolved spin-2 mediator should not be treated as simply another contact coefficient. For `nu q -> nu q` the mediator is exchanged in the **t channel**, so schematically

\[
\mathcal M\propto \frac{\text{momentum numerator}}{Q^2+M_T^2}.
\]

The contact-like momentum growth therefore turns over once `Q^2 ~ M_T^2`. This replaces the arbitrary hard `Theta(shat < xi Lambda^2)` endpoint regulator by physical propagator suppression and can convert a monotonic EFT rise into a localized distortion in energy and zenith.

The cost must be included from the start: crossing produces direct collider channels, and depending on the UV couplings an `s`-channel resonance in `c cbar -> tau+ tau-` or gluon-initiated production can be much more constraining than an EFT tail. Thus the resolved-mediator sequence is `choose gauge-consistent couplings -> direct-production limits -> exact propagator UHE fold -> Earth propagation`, not the reverse.
