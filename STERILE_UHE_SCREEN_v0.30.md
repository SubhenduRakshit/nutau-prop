# nutau-prop v0.30 — comprehensive light-sterile / non-resonant UHE cross-section screen

Date: 2026-10-30
Status: analytic + repository-kernel screening checkpoint. This is **not** an all-constraints certification.

## Scope

Goal: identify every physically distinct operator/process class presently worth pursuing for an O(5%) or larger modification of the UHE neutrino-nucleon cross section at E_nu ~ 100--300 PeV, with emphasis on light sterile-neutrino daughters. Resonant models are kept only as comparators.

This checkpoint consolidates previously validated hard kernels in the repository and adds the sterile-transition interpretation. It distinguishes:

1. exact/reused repository PDF-fold results;
2. analytic scaling estimates where a new dedicated dynamic-Q LHAPDF fold has not yet been run;
3. constraint status.

The central physical distinction is whether the BSM daughter is active (NC-like redistribution) or sterile (true active-flux disappearance in Earth propagation).

## Numerical reference conventions

Where an existing repository result already quotes r_NP/SM, that ratio is used directly and is not renormalized to a different historical SM table. For v0.29 absolute cross sections, use the same repository SM reference implicit in that checkpoint. Do not mix old CT10/CT18 SM normalizations across versions.

Target research gate: Delta sigma / sigma_SM >= 5% at 100 or 300 PeV.

---

# A. Sterile-transition derivative quark operators

## A1. Derivative-vector transition (strongest non-resonant candidate)

Low-energy structure:

L_D = C_D (bar N i<->partial_mu P_L nu_tau)(bar q gamma^mu q),

with [C_D] = GeV^-3. A fully electroweak-gauge-invariant parent requires the appropriate Higgs insertion / nuSMEFT completion.

For m_N << sqrt(shat), the hard kernel is the same momentum-growing kernel as repository Case 6, now with an unambiguous sterile daughter:

|M|^2 ~ 64 C_D^2 shat^3 y(1-y)

for the normalization used in the Case-6 code.

The previous CT10nlo repository calculation with pointwise sqrt(shat)<Lambda gives

sigma_NP/sigma_SM ~ 0.35--0.72 at 100--300 PeV

for representative Lambda ~ 2.0--2.44 TeV and the Case-6 benchmark coefficient C_D = 8.6e-10 GeV^-3.

### Verdict

SIZE: PASS by a very large margin.

PHYSICS ADVANTAGE: for a Dirac sterile N, nu_tau q -> N q is a genuine disappearance channel rather than active NC regeneration.

MAIN PROBLEM: the same momentum growth drives the contact theory toward its EFT/partial-wave boundary. Monojet + gauge-completion + unitarity certification remain decisive. The old Case-6 diagonal-active interpretation should not be used as the physical sterile model; the hard kernel itself is reusable.

Priority: 1 among non-resonant sterile models.

---

## A2. Derivative-tensor transition (PKR Case-5 / ATZ-Q10-like branch)

Repository normalization:

L_5 = C_5 partial^nu(bar nu^c sigma_{mu nu} P_L nu)(bar q gamma^mu q),

with published benchmark C_5 = 1.8e-10 GeV^-3. In the Dirac interpretation the chirality-flipped daughter can be a sterile/right-handed state.

The exact repository relation at equal cutoff is

sigma_5/sigma_6 = 0.5 (C_5/C_6)^2 = 0.0219037317.

At its own perturbative endpoint Lambda_UV = 4.1176 TeV, the validated exact-PDF prescreen gives

- 100 PeV: r_NP/SM = 0.0552
- 220 PeV: r_NP/SM = 0.0678
- 300 PeV: r_NP/SM = 0.0732
- mean y ~ 0.50

### Verdict

SIZE: PASS at both 100 and 300 PeV.

This is currently the cleanest already-computed sterile-transition benchmark that crosses 5% without relying on a resonance.

Main caveats: exact electroweak gauge completion, modern monojet recast, and partial-wave/EFT interpretation of the endpoint.

Priority: 2.

---

# B. Sterile-transition gluon / Rayleigh operators

## B1. Scalar Rayleigh transition

Representative gauge-invariant parent:

(C_G/Lambda^4) (bar N Htilde^dagger L_tau) G^A_{mu nu} G^{A mu nu}

(or the corresponding CP-odd dual-G operator).

After EWSB this generates nu_tau g -> N g. In the physical-amplitude normalization used in v0.29,

M = kappa delta^{ab} B_nu [(k.k')(eps.eps') - ...]

and

|M|^2_bar = |kappa|^2 (-t)^3/4,

d sigmahat/dy = |kappa|^2 shat^2 y^3/(64 pi),

sigmahat = |kappa|^2 shat^2/(256 pi).

The repository Case-7 / Rayleigh exact-PDF result at its perturbative endpoint Lambda_UV = 4.2825 TeV is

- 100 PeV: r_NP/SM = 0.0311
- 220 PeV: r_NP/SM = 0.0538
- 300 PeV: r_NP/SM = 0.0642
- mean y ~ 0.80

### Verdict

SIZE: FAIL at 100 PeV but PASS at 300 PeV.

The high-y kernel makes this particularly effective as an attenuation mechanism once the final state is sterile: each scattering removes or strongly degrades the active parent.

A separate v0.29 normalization example with kappa=6.4e-10 GeV^-3 gave only 0.71%/1.42% of the old SM total at 100/300 PeV for the 4-TeV common-scale benchmark. Therefore the >5% result occurs only near the perturbative/contact endpoint, not for a generic weakly-coupled common-scale UV benchmark.

Priority: 3, but with serious UV-normalization scrutiny.

## B2. CP-odd dual-gluon transition

Replacing G G by G Gtilde gives the same unpolarized pure-BSM hard kernel for massless external states. Therefore the cross-section numbers above apply at equal physical amplitude coefficient.

SIZE: same as B1.

Do not identify their neutrino-mass/CP constraints; only the unpolarized UHE kernel is identical.

---

# C. Dimension-6 sterile semileptonic operators

The complete nuSMEFT contains gauge-invariant active-sterile four-fermion operators such as scalar/tensor structures built from L, N, Q, u, d. After electroweak decomposition they generate processes of the form

nu_tau q -> N q'

without a resonance.

For a generic contact coefficient C_6/Lambda^2,

M ~ C_6 shat/Lambda^2,

sigmahat ~ C_6^2 shat/(16--64 pi Lambda^4),

with Lorentz-structure-dependent angular factors.

A useful analytic scale estimate at 100 PeV is

s = 2 m_N E_nu ~ 1.88e8 GeV^2.

For Lambda=1 TeV and C_6=1, the prefactor before PDF momentum fractions is O(10^3 pb). Multiplication by the relevant parton momentum fraction gives typically:

- heavy sea (s,c): O(10--50 pb), much less than 5% of total SM;
- light valence/sea: O(100 pb), still generally below 5% unless angular/tensor factors or C_6>1 help;
- at 300 PeV the rate grows roughly linearly before PDF evolution/cuts and can approach the 5% gate for light-quark structures at Lambda~1 TeV.

However, the same light-quark operators are directly exposed to LHC MET/charged-current channels and low-energy flavor constraints. Raising the allowed scale from 1 TeV to 2 TeV suppresses the contact rate by 1/16, eliminating a 5% effect.

### Verdict

SIZE: only conditionally competitive for Lambda ~ 1 TeV, order-one/tensor-enhanced coefficients, especially light-quark channels.

PHENOMENOLOGICAL PROSPECT: weak compared with A1/A2 because the energy growth is only sigmahat~shat rather than shat^2.

No dimension-6 sterile four-fermion benchmark is presently promoted above the 5% certified-candidate line without a dedicated current collider/flavor fit.

Priority: secondary.

---

# D. Sterile electroweak dipole portal

Gauge-invariant operators:

O_NB = (bar L_tau sigma^{mu nu} N) Htilde B_{mu nu},
O_NW = (bar L_tau sigma^{mu nu} tau^I N) Htilde W^I_{mu nu}.

They generate nu_tau q -> N q through photon/Z exchange.

Although the dipole vertex contains momentum, the t-channel propagator produces a low-Q enhancement rather than sustained hard-DIS growth. The UHE integral is therefore dominated by comparatively small Q^2, and the cross section does not acquire the polynomial shat^2 growth of A/B.

Modern dipole-portal collider studies probe transition coefficients around d_gamma ~ 1e-7 GeV^-1 at future colliders, and direct/cosmological constraints can be stronger depending on m_N and lifetime. A coefficient at that scale is far below what would be needed for a robust O(5%) addition to the inclusive 100--300 PeV DIS cross section.

### Verdict

SIZE: FAIL for realistic allowed dipole coefficients as an inclusive UHE nu-N enhancement mechanism.

Potentially interesting for low-Q spectral signatures, not for this project's >5% hard-DIS gate.

Priority: reject for the present purpose.

---

# E. Resolved derivative Z' transition

Model:

(kappa) Z'_mu (bar N i<->partial^mu P_L nu_tau) + g_q Z'_mu bar q gamma^mu q

(or the corresponding charge-conjugated PKR structure).

The exact resolved kernel is the derivative contact kernel multiplied by

[M_Z'^2/(M_Z'^2+Q^2)]^2

when normalized to the same low-Q contact coefficient.

Repository v0.16 scanned M_Z' = 0.1--3 TeV at fixed low-Q C4 = 8.6e-10 GeV^-3. The propagator physically turns over the contact growth at Q^2 >> M_Z'^2.

### Verdict

The resolved mediator is theoretically healthier than an indefinitely growing contact operator, but the same propagator also removes much of the UHE enhancement. It remains a conditional candidate only if a coupling/mass point can simultaneously retain >5% at 100--300 PeV and pass collider/UV-vertex constraints.

No new >5% certified point is claimed here because the repository v0.16 hard tables were not stored in this checkpoint in a directly reusable summary form and the derivative vertex itself is still dimension five.

Priority: model-building follow-up, below A1/A2.

---

# F. Active-state non-resonant comparators already in the repo

These are retained because they establish how difficult the >5% gate is without a sterile sink.

## F1. Active two-derivative d=8 quark operators

O_8 ~ (bar L_3 gamma_mu <->D_nu L_3)(bar q_R gamma^mu <->D^nu q_R).

Strict partial-wave unitarity through shat=Lambda^2 gives

|C| <= 8 pi/3 = 8.3776.

At Lambda=4 TeV, the pure d=8 square then contributes only

- c_R: 0.11% / 0.22% of total SM at 100/300 PeV
- s_R: 0.14% / 0.28% of total SM at 100/300 PeV

The only possible rescue is positivity-fixed SM interference, whose exact dynamic-Q fold remains the decisive numerical gate.

Verdict: pure term FAILS decisively; interference not yet certified above 5%.

## F2. Active gluon twist-2 d=8

For C=-8 pi, M=4 TeV, repository v0.29 gives

xi=1:
- 100 PeV: 229.75 pb = 3.30% of old SM total
- 300 PeV: 733.984 pb = 7.05%

xi=0.5:
- 100 PeV: 62.828 pb = 0.90%
- 300 PeV: 164.157 pb = 1.58%

Verdict: PASS only in the aggressive full-domain xi=1 treatment at 300 PeV; not robust to the stricter cutoff.

## F3. Resolved charged scalar

Broken-phase benchmark A=|g_tau|^2(|g_L|^2+|g_R|^2).

For M_H=90 GeV, A=1:
- 100 PeV: 1658.4 pb, roughly 24% of the v0.29 SM reference
- 300 PeV: 2393.8 pb, roughly 23%

A tuned g_tau=1, g_L=g_R=0.4 gives A=0.32:
- 100 PeV: 530.7 pb, roughly 7.6%
- 300 PeV: 766.0 pb, roughly 7.4%

The simple Ds pseudoscalar constraint can cancel for g_L=g_R, but a gauge-complete two-doublet texture generates additional tau/K/neutral-scalar constraints. This is not an allowed benchmark.

Verdict: SIZE PASS, constraint viability unresolved and likely difficult.

---

# G. Resonant comparator: tau-neutrinogluon

Not part of the non-resonant ranking, but retained as the current benchmark to beat.

For M=3 TeV, Lambda=9 TeV, a=4:

- 100 PeV: 427.7 pb = 6.14%
- 300 PeV: 1347.5 pb = 12.95%

Gamma/M ~ 0.0345.

This obtains the enhancement from an on-shell color-octet leptogluon and increasing gluon luminosity rather than uncontrolled contact growth. Charged-partner collider limits and the origin of a=4 remain decisive.

---

# Ranked result

## Clear >5% non-resonant candidates

1. **Sterile derivative-vector quark transition** — ~35--72% at 100--300 PeV in the existing Case-6 benchmark; strongest but closest to EFT/unitarity danger.
2. **Sterile derivative-tensor quark transition** — 5.52% at 100 PeV, 7.32% at 300 PeV at its perturbative endpoint; best already-computed moderate-size sterile candidate.
3. **Sterile/Dirac gluon Rayleigh transition** — 3.11% at 100 PeV, 6.42% at 300 PeV; hard y~0.8 makes attenuation especially strong.
4. **Resolved charged scalar** — ~7--24% depending on coupling benchmark; gauge/flavor viability unresolved.
5. **Active gluon twist-2** — 7.05% at 300 PeV only for xi=1; falls below 2% for xi=0.5.

## Conditional / currently below gate

- dimension-6 sterile scalar/tensor semileptonic operators: can approach >5% only for ~1 TeV scales / favorable light-quark or tensor coefficients; current constraints likely remove much of that region;
- resolved derivative Z': needs exact coupling/mass optimization; propagator softens the growth;
- active d=8 two-derivative quark: pure term <0.3% of total under strict unitarity; only interference can rescue it.

## Reject for this purpose

- ordinary sterile electroweak dipole portal as a hard-DIS enhancement mechanism: low-Q dominated and too small at realistic coefficients.

---

# Important propagation consequence

For sterile final states, Delta sigma/sigma_SM is not the whole story. A process

nu_tau N -> N X

removes the active tau neutrino instead of regenerating a lower-energy active neutrino. Therefore a 5% sterile-conversion cross section can produce a larger change in Earth attenuation than a 5% NC-like active redistribution term.

This does **not** justify transport before certification. It means the sterile candidates should be ranked by an eventual attenuation optical depth, not only by inclusive cross-section percentage.

---

# Next certification order

Do not run Earth propagation yet. For the three sterile survivors, use:

1. exact external-state/gauge-invariant matching;
2. partial-wave unitarity with the correct fixed-helicity channels;
3. exact modern monojet / tau+MET / dijet or gluon constraints as applicable;
4. cosmology/BBN only for genuinely light, long-lived sterile daughters;
5. exact dynamic-Q CT18/NNPDF fold using one common SM reference;
6. only then Earth attenuation/transport.

Recommended priority: A2 first (moderate, already >5% and less violent than A1), then A1, then B1/B2.

---

# Literature anchors

- complete nuSMEFT/nuLEFT operator bases with sterile neutrinos: Li et al., arXiv:2105.09329; Liao & Ma, arXiv:1612.04527; earlier d=7 basis arXiv:1505.05264.
- PKR momentum-growing neutrino-parton interactions: Pandey, Karmakar, Rakshit, arXiv:1907.07700.
- sterile dipole portal collider context: arXiv:2404.09609; arXiv:2502.07024.
- repository source checkpoints: OPERATOR_SURVEY_v0.14.md, CONTACT_OPERATOR_RESULTS_v0.15.md, RESOLVED_MEDIATOR_v0.16.md, D8_CERTIFICATION_v0.25.md, CHAT_UHE_CROSS_SECTIONS_v0.29_c1_2026-10-30.md.

## Bottom line

The statement 'only resonances can give >5%' is false. The repository already contains two non-resonant sterile-compatible hard kernels that clear the gate: the derivative-tensor transition at 100 and 300 PeV and the Rayleigh/gluon transition at 300 PeV; the stronger derivative-vector transition can be tens of percent. The question is no longer size, but whether any of these survive gauge-completion, unitarity and collider/cosmology certification.
