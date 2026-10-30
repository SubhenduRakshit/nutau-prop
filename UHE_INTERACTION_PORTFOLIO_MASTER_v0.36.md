# nutau-prop — authoritative UHE interaction portfolio v0.36

Date: 2026-10-30
Status: authoritative portfolio through the PKR Eq. (4.3) resolved-dipole screen.

> **Audit policy:** `AUDIT_PORTFOLIO_v0.33_REVIEW.md` is **SUPERSEDED—DO NOT USE**. Preserve it and its failed reasoning; do not delete it. Earlier portfolio/checkpoint files remain historical records. This file is the current single decision record.

## Target

Find an interaction giving at least a 5% modification of the inclusive UHE neutrino–nucleon cross section around 100–300 PeV while maintaining the correct external states, gauge/UV consistency, partonic validity, and collider/flavour/neutrino constraints. No Earth transport is run until one benchmark clears the complete constraint gate.

## A. Serious non-resonant candidates still open

### A1. PKR Eq. (4.4): resolved derivative active-to-sterile Z'

\[
\mathcal L\supset {c_2\over\Lambda}(\bar N i\overleftrightarrow{\partial_\mu}P_L\nu_\tau)Z'^\mu
+g_q\bar q\gamma_\mu qZ'^\mu,
\qquad \kappa_2={c_2g_q\over\Lambda}.
\]

Topology: \(\nu_\tau q\to Nq\), a sterile sink.

Every trial coupling is folded with the moving domain
\[
\sqrt{\hat s}<1/\kappa_2.
\]

CT18 self-consistent 5% thresholds:

| M [GeV] | kappa2_5%(100 PeV) | kappa2_5%(300 PeV) |
|---:|---:|---:|
| 5 | 3.0943e-5 | 3.1251e-5 |
| 40 | 5.4254e-5 | 5.1974e-5 |
| 100 | 7.9706e-5 | 7.4021e-5 |
| 220 | 1.2008e-4 | 1.0898e-4 |
| 500 | 2.1643e-4 | 1.9065e-4 |
| 1000 | no controlled crossing | no crossing |

Historical PKR at 5 GeV: LHC <3.3e-3 GeV^-1, IceCube <2.5e-4 GeV^-1. Historical constraints therefore leave the low-mass 5% window open. Modern ATLAS/IceCube recasts and a credible UV origin of the derivative vertex are still required.

Verdict: **OPEN / SIZE PASS / best phenomenological non-resonant candidate / modern-constraint and UV conditional**.

### A2. PKR Eq. (4.3): resolved transition-dipole Z'

\[
\mathcal L\supset {c^{(1)}\over\Lambda}(\bar N^c\sigma_{\mu\nu}P_L\nu_\tau)Z'^{\mu\nu}
+g_q\bar q\gamma_\mu qZ'^\mu,
\qquad \kappa_1={c^{(1)}g_q\over\Lambda}.
\]

It is intrinsically transition/off-diagonal. It is not a diagonal active NC interaction.

Using the validated PKR contact kernels,
\[
\boxed{d^2\sigma_{4.3}/dxdy=2\,d^2\sigma_{4.4}/dxdy}
\]
at equal \(M\) and numerical \(\kappa\) for the same massless off-diagonal/sterile external states. The v0.36 scan therefore solves the validated Eq. (4.4) fold at 2.5% while imposing \(\sqrt{\hat s}<1/\kappa_1\) for every trial coupling.

CT18 self-consistent Eq. (4.3) 5% thresholds:

| M [GeV] | kappa1_5%(100 PeV) | kappa1_5%(300 PeV) |
|---:|---:|---:|
| 5 | 2.1880e-5 | 2.2098e-5 |
| 15 | 2.7932e-5 | 2.7529e-5 |
| 40 | 3.8364e-5 | 3.6751e-5 |
| 100 | 5.6361e-5 | 5.2226e-5 |
| 220 | 8.4594e-5 | 7.5899e-5 |
| 500 | 1.4390e-4 | 1.2611e-4 |
| 1000 | 3.3001e-4 | 2.6618e-4 |
| 2000 | no crossing | no crossing |

Historical PKR at 5 GeV: LHC <1.3e-3 GeV^-1, IceCube <4.8e-4 GeV^-1. LHC overtakes IceCube around 15 GeV in the historical analysis.

Eq. (4.3) is quantitatively stronger than Eq. (4.4) and extends the controlled size reach to about 1 TeV. However the width depends separately on the dipole vertex:
\[
\Gamma_{\nu N}/M\propto (c^{(1)}/\Lambda)^2M^2
=\kappa_1^2M^2/g_q^2.
\]
Thus \(\kappa_1\) alone cannot certify the model. The separate parameters \(c^{(1)}/\Lambda\), \(g_q\), total width, branching ratios, induced photon transition moment, and direct-Z' limits must be retained.

Verdict: **OPEN / SIZE PASS / strongest PKR non-resonant size candidate / width and matching conditional**.

### Decision between A1 and A2

- Eq. (4.3) wins on raw UHE size.
- Eq. (4.4) is presently cleaner phenomenologically because it avoids the transition-dipole width/matching liability.
- Both need the same modern ATLAS monojet and IceCube operator-specific recasts before either can be certified.

Therefore the project should carry them as a **two-candidate non-resonant frontier**, not prematurely choose one and discard the other.

## B. Resonant candidates

### B1. Strange-aligned tilde-R2 scalar leptoquark

Renormalizable and gauge complete. Stored M=1.5 TeV, lambda=2 benchmark gives about 7.17% at 100 PeV and 9.71% at 300 PeV. Finite width, both weak-doublet components, pair/single production and flavour constraints remain to be combined.

Verdict: **OPEN / best UV-clean resonant fallback**.

### B2. Tau-neutrinogluon / leptogluon resonance

Stored M=3 TeV, Lambda=9 TeV, a=4 benchmark gives 6.14% and 12.95% at 100/300 PeV. Strong size, but the large chromomagnetic transition coefficient and coloured-partner limits are unresolved.

Verdict: **OPEN / SIZE PASS / strong-matching conditional**.

## C. Strong/composite non-resonant benchmarks only

- PKR Eq. (4.6) / sterile derivative-tensor contact: 5.52% / 7.32%, but only near C Lambda^3 ~4pi; historical LHC bound dominates.
- PKR Eq. (4.8) / sterile gluon-Rayleigh contact: 3.11% / 6.42%, passing only near the strong/validity endpoint; historical monojet dominates and weakly coupled UV completions are loop suppressed.
- Active gluon twist-2 d=8: >5% only at the aggressive validity/strong-coupling edge; tightening the shat cut collapses the effect.

Verdict: **not closed as strong/composite physics, but not perturbative EFT survivors**.

## D. Other conditional branch

Charged scalar: stored tuned benchmark ~7.6% / 7.4%, but gauge completion, CKM-linked flavour and collider constraints make it lower priority.

Verdict: **conditional / likely tuned**.

## E. Failed / closed for the >5% perturbative goal

1. Gauge-invariant sterile d=6 scalar: ~0.005–0.007% at current controlled bounds.
2. Gauge-invariant sterile d=6 tensor: ~0.54–0.71%; even 4pi contact domain below 5%.
3. Preliminary d=6 tensor 8–11% claim: **withdrawn**, factor-16 normalization error from mixing low-energy epsilon and SMEFT conventions.
4. Active quark twist-2 d=8: <0.25% total after interference/positivity/unitarity/dynamic-Q treatment.
5. Sterile electroweak/photon dipole portal: low-Q dominated; fails the hard inclusive >5% objective. This is distinct from the PKR resolved Z' transition dipole above.
6. Universal B-3L_tau vector: sub-percent at allowed couplings, with u/d NC interference cancellation.
7. Literal PKR Eq. (4.7) derivative-vector contact: collider dominated and external-state interpretation nontrivial; its useful controlled continuation is resolved Eq. (4.4).
8. Generic weakly-coupled UV completions of the endpoint tensor/Rayleigh contacts: loop suppression prevents the required size without strong dynamics.

## F. Current no-go status

The survey now supports a strong **practical pattern**:

\[
\text{weakly coupled non-resonant interactions usually give }<5\%
\]
unless they use a resolved momentum-dependent mediator or approach a strong/validity endpoint.

But a no-go theorem is **not yet justified**, because two resolved PKR branches remain open:
\[
\boxed{\text{Eq. (4.4) derivative }Z'\quad\text{and}\quad\text{Eq. (4.3) transition-dipole }Z'.}
\]

If modern collider + IceCube recasts plus width/UV matching close both, the surveyed perturbative non-resonant space will be very close to a defensible no-go statement, leaving resonances or strong/composite dynamics as the natural mechanisms for >5% effects.

## G. Next justified work

1. Modern ATLAS 139 fb^-1 monojet/jets+MET recast for both Eq. (4.3) and Eq. (4.4), with actual signal acceptance and binned likelihood.
2. Modern IceCube energy/zenith likelihood recast for both resolved branches.
3. For Eq. (4.3), scan \((M,g_q,c^{(1)}/\Lambda)\) and impose mediator width, induced photon-transition-moment matching and direct-Z' constraints.
4. For Eq. (4.4), specify/assess the UV origin of the derivative neutrino vertex while keeping \(c_2/\Lambda\) and \(g_q\) separate.
5. Compare modern \(\kappa_{\max}(M)\) to the self-consistent \(\kappa_{5\%}(M)\) curves.
6. **No Earth propagation until at least one benchmark satisfies this complete gate.**

Detailed Eq. (4.3) record: `PKR_RESOLVED_DIPOLE_v0.36.md`.
