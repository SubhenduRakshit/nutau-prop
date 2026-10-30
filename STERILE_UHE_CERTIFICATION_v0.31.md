# nutau-prop v0.31 — sterile-transition UHE certification exploration

Date: 2026-10-30
Status: theory/constraint exploration checkpoint; not a final all-constraints certification.

## Goal

Explore the three leading non-resonant sterile-compatible candidates identified in v0.30:

1. derivative-tensor transition;
2. derivative-vector transition;
3. sterile gluon/Rayleigh transition.

The project criterion remains an additional UHE neutrino-nucleon contribution of at least 5% at 100–300 PeV, while keeping EFT/partial-wave control and avoiding obviously fatal collider/low-energy constraints. A sterile daughter is treated as active-flux disappearance for Earth attenuation unless it has appreciable SM-active mixing/interactions.

---

## 1. Derivative-tensor sterile transition

Repository hard kernel (PKR Case 5 / ATZ-Q10-like normalization):

\[
|M_5|^2=32 C_5^2\hat s^3 y(1-y),
\]

with benchmark

\[
C_5=1.8\times10^{-10}\;\mathrm{GeV}^{-3}.
\]

Exact-PDF repository result at its perturbative endpoint:

- Lambda_UV = 4.1176 TeV;
- 100 PeV: Delta sigma / sigma_SM = 5.52%;
- 220 PeV: 6.78%;
- 300 PeV: 7.32%;
- mean y approximately 0.50.

The dimensionless endpoint coefficient is

\[
c_5=C_5\Lambda^3\simeq 12.6\simeq4\pi.
\]

Therefore the >5% result sits at the strong-coupling/perturbative edge, not at a generic O(1) Wilson coefficient. For c=1 at fixed Lambda the rate would be suppressed by roughly 1/(4pi)^2, far below the 5% gate.

Gauge-invariant sterile-neutrino EFT classifications show that four-fermion derivative operators of the psi^4 D class are loop-generated in weakly-coupled renormalizable UV completions. Hence a weakly-coupled UV completion would naturally suppress the coefficient further. A strong/composite UV origin or a resolved mediator is required if this branch is to remain phenomenologically interesting.

Constraint status:

- PKR found LHC monojet constraints stronger than IceCube for contact-type two-neutrino/two-parton operators in their study.
- The sterile interpretation removes the problematic diagonal active-state identification and gives a genuine sink channel, but it does not remove the monojet/EFT constraint problem.
- A modern exact tau-flavor sterile monojet recast is still required.

Verdict: SIZE PASS, NATURALNESS/UV CONTROL WEAK. Retain as a serious non-resonant benchmark only for strong-coupling or resolved-UV scenarios.

---

## 2. Derivative-vector sterile transition

Low-energy structure:

\[
\mathcal L_D=C_D(\bar N i\overleftrightarrow{\partial_\mu}P_L\nu_\tau)(\bar q\gamma^\mu q),
\]

with repository Case-6 normalization

\[
|M_D|^2\propto64 C_D^2\hat s^3 y(1-y).
\]

Benchmark coefficient:

\[
C_D=8.6\times10^{-10}\;\mathrm{GeV}^{-3}.
\]

The previous exact PDF fold gave approximately

\[
\Delta\sigma/\sigma_{SM}\sim0.35-0.72
\]

at 100–300 PeV for representative EFT cutoffs around 2.0–2.44 TeV.

At Lambda=2.445 TeV,

\[
c_D=C_D\Lambda^3\simeq12.6\simeq4\pi.
\]

At Lambda=2.0 TeV the same coefficient corresponds to c approximately 6.9, still large but below 4pi. Thus this branch can remain above 5% before the absolute perturbative endpoint, unlike the tensor/Rayleigh cases.

Physics advantage:

\[
\nu_\tau q\to N q
\]

is true active-flux disappearance for a sufficiently sterile, long-lived N. The same raw cross-section addition therefore produces a stronger Earth attenuation effect than an active NC-like redistribution channel.

Main problems:

1. strong energy growth and endpoint sensitivity;
2. monojet/dineutrino constraints;
3. electroweak gauge completion of the low-energy derivative operator;
4. derivative psi^4 D structures are loop-generated in weakly-coupled renormalizable UV completions;
5. a fully resolved parent is needed before treating the highest-energy behavior as physical.

Resolved-Z' comparison: the repository v0.16 parent multiplies the contact kernel by

\[
\left[{M_{Z'}^2\over M_{Z'}^2+Q^2}\right]^2,
\]

which physically turns over the growth at high Q^2. This improves UV behavior but suppresses the desired UHE enhancement. A useful resolved model must therefore live in the intermediate window where the derivative numerator is important but propagator suppression has not yet erased the 5% effect.

Verdict: strongest SIZE candidate and still the highest-priority non-resonant sterile model, but only if a resolved strong/UV completion can be built. Pure contact interpretation should not be treated as a final model.

---

## 3. Sterile gluon / Rayleigh transition

Representative active-sterile gauge-invariant parent requires a Higgs insertion, schematically

\[
{C_G\over\Lambda^4}(\bar N\widetilde H^\dagger L_\tau)G^A_{\mu\nu}G^{A\mu\nu}
\]

or the CP-odd dual-gluon analogue.

Repository hard kernel:

\[
{d\hat\sigma\over dy}={|\kappa|^2\hat s^2y^3\over64\pi},
\qquad
\hat\sigma={|\kappa|^2\hat s^2\over256\pi}.
\]

Exact Case-7 PDF result at the perturbative endpoint:

- Lambda_UV = 4.2825 TeV;
- 100 PeV: 3.11%;
- 220 PeV: 5.38%;
- 300 PeV: 6.42%;
- mean y approximately 0.80.

Again,

\[
c_G^{eff}\simeq C_7\Lambda^3\simeq4\pi
\]

for the low-energy d=7 normalization used in the repository. The v0.29 common-scale weakly-coupled-style benchmark was much smaller: only about 0.7% and 1.4% of total SM at 100 and 300 PeV for the quoted 4-TeV example.

The hard y^3 kernel is nevertheless valuable for Earth attenuation because a sterile daughter plus high inelasticity removes active energy very efficiently.

Gauge/UV point: sterile-neutrino EFT classifications contain neutrino-gluon Rayleigh structures and classify the corresponding field-strength operators as loop-generated in weakly-coupled renormalizable UV theories. The active-sterile transition with a Higgs insertion is at still higher gauge-invariant dimension, so a large weak-scale coefficient is not naturally expected from a perturbative UV completion.

Verdict: SIZE PASS only near the strong-coupling endpoint; excellent attenuation shape but weak naturalness. Keep as a secondary strong/composite benchmark, not the leading weakly-coupled candidate.

---

## 4. Dimension-6 sterile semileptonic operators — why they are not promoted

Complete nuSMEFT/LNEFT bases contain scalar, vector, and tensor active-sterile quark operators. Literature single-coefficient analyses find representative NP-scale limits of about 0.5–1.5 TeV for neutral-current neutrino-quark operators, while charged-current operators can be far more strongly constrained by pseudoscalar-meson and tau decays, reaching multi-TeV to tens-of-TeV scales depending on flavor and chirality.

Their partonic growth is only

\[
\hat\sigma\propto \hat s/\Lambda^4,
\]

not the \hat s^2-type growth of the derivative/Rayleigh branches. If a light-quark neutral-current sterile operator is truly allowed near the low end of the 0.5–1.5 TeV scale range it could potentially approach or exceed the 5% UHE gate, but this requires an exact flavor-specific modern fit and PDF fold. Charged-current d=6 branches are much less promising because low-energy bounds typically drive Lambda too high.

Verdict: neutral-current d=6 sterile operators remain worth one dedicated exact fold; charged-current branches are demoted.

---

## 5. Dipole portal

Gauge-invariant operators O_NB and O_NW give photon/Z-mediated active-sterile scattering. The momentum in the dipole vertex is offset by the t-channel propagator, pushing the integral toward low Q^2 rather than creating sustained hard-DIS growth. Existing collider/low-energy/cosmology constraints on transition moments further suppress the allowed coefficient.

Verdict: reject as a route to a >5% inclusive UHE nu-N enhancement, though it can remain interesting for low-Q signatures.

---

## 6. Updated ranking

### Pure size ranking

1. derivative-vector sterile: tens of percent possible;
2. derivative-tensor sterile: 5.5–7.3%;
3. sterile gluon/Rayleigh: up to 6.4% at 300 PeV.

### Credibility ranking after UV/naturalness considerations

1. resolved derivative-vector parent, if a viable UV completion can keep >5%;
2. neutral-current d=6 sterile operator near the weakest allowed flavor directions — needs exact fold;
3. derivative-tensor sterile in a strong/composite UV scenario;
4. sterile gluon/Rayleigh in a strong/composite UV scenario;
5. dipole portal — reject for the present hard-UHE goal.

The key correction relative to v0.30 is that the tensor and Rayleigh >5% points occur essentially at c≈4pi and should not be presented as generic weakly-coupled EFT successes.

---

## 7. Immediate next quantitative calculation

The most informative next calculation is NOT Earth propagation yet. It is:

1. exact dynamic-Q CT18/NNPDF fold for the least-constrained neutral-current d=6 active-sterile quark operators using current one-coefficient bounds;
2. a resolved derivative-vector mediator scan with couplings constrained rather than fixed to the old Case-6 contact normalization;
3. only for survivors with Delta sigma/sigma_SM >=5%, run sterile-sink Earth attenuation with the repaired transport solver.

This sequence directly tests whether any non-resonant sterile model survives without requiring c≈4pi.