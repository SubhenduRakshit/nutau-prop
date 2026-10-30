# nutau-prop v0.36 — PKR Eq. (4.3) resolved transition-dipole Z' screen

Date: 2026-10-30

## Scope

This checkpoint completes the remaining serious PKR non-resonant branch before making any broad no-go statement. It treats PKR Eq. (4.3),

\[
\mathcal L\supset {c^{(1)}\over\Lambda}(\bar\nu_i^c\sigma_{\mu\nu}P_L\nu_j)Z'^{\mu\nu}
+g_q\bar q\gamma_\mu q Z'^\mu,
\qquad i\ne j,
\]

with
\[
\kappa_1\equiv {c^{(1)}g_q\over\Lambda}.
\]

The interaction is transition/off-diagonal. For the active-to-sterile interpretation relevant to Earth attenuation, the hard process is \(\nu_\tau q\to Nq\). It must not be silently re-labelled as a diagonal active NC interaction.

## Exact kernel relation to PKR Eq. (4.4)

The validated contact limits in the repository are
\[
|\mathcal M_{4.6}|^2=32C_5^2\hat s^3y(1-y),\qquad
|\mathcal M_{4.7}|^2=64C_6^2\hat s^3y(1-y).
\]

Integrating out the resolved Eq. (4.3) mediator gives
\[
C_5^{\rm eff}={2\kappa_1\over M_{Z'}^2+Q^2},
\]
whereas Eq. (4.4) gives
\[
C_6^{\rm eff}={\kappa_2\over M_{Z'}^2+Q^2}.
\]
Therefore, at equal \(M_{Z'}\) and equal numerical \(\kappa\), for the same massless off-diagonal/sterile external states,
\[
\boxed{
{d^2\sigma_{4.3}\over dxdy}=2\,{d^2\sigma_{4.4}\over dxdy}.
}
\]
This fixes the normalization without a new fit or guess.

The v0.36 calculation therefore obtains the Eq. (4.3) 5% threshold by solving the already validated Eq. (4.4) fold at a 2.5% target, while recomputing the moving validity domain for every trial coupling.

## Validity prescription

For every trial \(\kappa_1\), the PDF integral imposes
\[
\sqrt{\hat s}<1/\kappa_1.
\]
No fixed-cut quadratic rescaling is used. This is a one-parameter screening validity prescription; a full UV model must separate \(c^{(1)}/\Lambda\) and \(g_q\).

## CT18 self-consistent 5% thresholds

| M_Z' [GeV] | kappa1_5%(100 PeV) [GeV^-1] | kappa1_5%(300 PeV) [GeV^-1] |
|---:|---:|---:|
| 5 | 2.1880e-5 | 2.2098e-5 |
| 15 | 2.7932e-5 | 2.7529e-5 |
| 40 | 3.8364e-5 | 3.6751e-5 |
| 100 | 5.6361e-5 | 5.2226e-5 |
| 220 | 8.4594e-5 | 7.5899e-5 |
| 500 | 1.4390e-4 | 1.2611e-4 |
| 1000 | 3.3001e-4 | 2.6618e-4 |
| 2000 | no crossing | no crossing |

Thus Eq. (4.3) is a stronger size candidate than Eq. (4.4): it needs roughly 1/sqrt(2) smaller effective coupling at low/intermediate mass and still admits a controlled 5% crossing near 1 TeV, where Eq. (4.4) did not.

## Historical PKR constraints

PKR quotes at M_Z'=5 GeV:
\[
\kappa_{1,\rm LHC}<1.3\times10^{-3}\ {\rm GeV}^{-1},\qquad
\kappa_{1,\rm IceCube}<4.8\times10^{-4}\ {\rm GeV}^{-1}.
\]
The v0.36 5% requirement, \(\kappa_{1,5\%}\simeq2.2\times10^{-5}\) GeV^-1, is therefore far below both historical 5-GeV bounds.

PKR Fig. 2(b) shows LHC becoming stronger than IceCube above roughly 15 GeV. PKR also notes that the momentum-enhanced dipole interaction can make the Z' width large at high mass; in their setup the collider sensitivity eventually saturates and ceases to give a useful constraint around the several-hundred-GeV regime.

These historical statements are baselines only; they are not modern recasts.

## Width and coupling separation

The cross section depends on the product
\[
\kappa_1=(c^{(1)}/\Lambda)g_q,
\]
but the neutrino partial width depends on the dipole vertex itself, parametrically
\[
{\Gamma_{Z'\to\nu N}\over M_{Z'}}\propto
\left({c^{(1)}\over\Lambda}\right)^2M_{Z'}^2
={\kappa_1^2M_{Z'}^2\over g_q^2}.
\]
Therefore a single \(\kappa_1\) curve is insufficient for certification. Small \(g_q\) at fixed \(\kappa_1\) forces a large transition dipole and can generate a broad mediator. Direct-Z' quark constraints pull in the opposite direction by limiting \(g_q\). A viable scan must retain at least \(M_{Z'}\), \(g_q\), \(c^{(1)}/\Lambda\), total width and branching fractions separately.

The exact numerical prefactor of the dipole width should be fixed in the eventual UV/likelihood implementation; only the scaling above is used here.

## Modern constraints to be recast

1. Full Run-2 ATLAS monojet/jets+MET (13 TeV, 139 fb^-1) is the appropriate current collider dataset. A trustworthy bound requires signal generation, acceptance/cuts and the binned likelihood; luminosity rescaling of PKR is only a screening estimate.
2. Current IceCube cross-section/event data require an operator-specific energy/zenith likelihood. No ready-made 2026 \(\kappa_1(M)\) curve exists for this exact transition interaction.
3. Modern active-sterile photon transition-moment constraints are relevant only after matching the Z' dipole to an induced photon dipole (e.g. via kinetic mixing). They do not directly constrain \(\kappa_1\) without specifying \(g_q\), the mixing and the sterile mass.
4. Direct low-mass Z' constraints on \(g_q\), meson decays and LLP searches must be applied separately.

## Recheck of PKR Eqs. (4.6)-(4.8)

- Eq. (4.6), derivative-tensor contact: LHC historically dominates; repository size passes only near the strong/4pi endpoint. Retain only as strong/composite benchmark.
- Eq. (4.7), derivative-vector contact: LHC historically dominates and the literal charge-conjugated external states require care. The resolved Eq. (4.4) is the controlled continuation; no clean contact loophole remains.
- Eq. (4.8), neutrino-gluon contact: LHC historically dominates by orders of magnitude over IceCube; the stored >5% result occurs only near the strong/endpoint regime. Retain as strong/composite Rayleigh benchmark.

No additional perturbative PKR contact candidate is promoted.

## Decision

Size alone gives
\[
\boxed{\text{Eq. (4.3) transition dipole is at least as important as Eq. (4.4), and quantitatively stronger.}}
\]

Model credibility does not automatically improve, because Eq. (4.3) has additional width and transition-moment/matching liabilities. The two resolved PKR mediators should therefore be treated as the final pair of serious non-resonant candidates:

1. Eq. (4.4) resolved derivative active-to-sterile Z' — cleaner phenomenologically, no dipole-width liability;
2. Eq. (4.3) resolved transition-dipole Z' — stronger UHE size reach, but width and dipole matching must be certified.

A non-resonant no-go claim is premature until modern monojet + IceCube recasts close both branches. No Earth transport is justified yet.
