# Broader PKR/ATZ and mediator assessment — v0.29 supplement
Date: 30 October 2026. Analytical/literature review only; no simulations or code changes.
This supplements AUDIT_UHE_OPERATORS_v0.29_2026-10-30.md. It corrects that audit's overly narrow candidate selection. No new parameter point is certified as satisfying all constraints.

## What the earlier audit missed

The failed 4-TeV active quark benchmarks do not exclude other Lorentz structures, lepton-number assignments, mediator masses, or production thresholds. The main search should compare physically distinct mechanisms, not only operators with the highest power of energy.

## Operator-family triage

Labels below follow [ATZ](https://arxiv.org/pdf/1812.02778); equations follow [PKR](https://arxiv.org/pdf/1907.07700). These papers supply the operator definitions, not current universal bounds.

| Family | Necessary interpretation/correction | Decision |
|---|---|---|
| Vector/axial neutrino currents, ATZ dimension-six Q1,2; PKR 4.1/4.5 | Active NC; include SM interference, electroweak partners and mediator dependence | Retain a consistent vector mediator as a control, not a promised sharp UHE turn-on |
| Scalar/pseudoscalar gluons, ATZ dimension-seven Q3,4; PKR 4.8 | Dirac chirality-flip and active LNV are distinct theories; alpha_s normalization differs between conventions | High-priority missing comparison |
| Derivatives, ATZ Q8–11; PKR 4.6/4.7 | Fix external states, flavour symmetry, on-shell identities and collider coefficient mapping first | One analytical reconciliation before further scans |
| Scalar/pseudoscalar/tensor quark operators, ATZ Q5–7 | Explicit quark-mass factors are a coefficient convention; physical Yukawa suppression requires a flavour assumption | Screen electroweak scalar/leptoquark parents rather than dismiss by notation |
| Photon dipole/Rayleigh, ATZ dimension-five Q1 and dimension-seven Q1,2 | Low-recoil/coherent physics, electromagnetic constraints and photon-initiated channels differ from gluon DIS | Secondary; do not count a forward/coherent rate as equivalent to hard UHE attenuation |
| Dark-vector dipole/derivative vertices, PKR 4.3/4.4 | A resolved vector plus dimension-five vertex is still an EFT | Consider only with a parent completion and controlled width/cutoff |

### A useful exact derivative check

Set sigma^{mu nu}=i[gamma^mu,gamma^nu]/2, q=p'-p, P=p'+p and neglect external neutrino masses. The free on-shell Dirac equations give
\[
q_\mu\bar u(p')\sigma^{\mu\nu}P_Lu(p)
=-iP^\nu\bar u(p')P_Lu(p).
\]
The derivative acting on the bilinear contributes iq_mu. Thus the tensor-divergence and scalar-overleftrightarrow-derivative structures have the same massless on-shell neutrino factor with these conventions. This applies to the corresponding axial quark current as well.

This is a restricted amplitude identity, not a claim that arbitrary massive/interacting off-shell EFT bases are identical. Additional interactions and mass terms require their own treatment. Nevertheless, independently assigning different massless DIS shapes or unrelated normalizations to ATZ Q8 versus Q10 (or Q9 versus Q11) must be reconciled with this identity. For PKR, compare identical chiral/LNV external states before using it; an unprojected interaction can contain additional channels. Neither the old factor 64 nor a proposed factor 16 is a substitute for this check.

## The overlooked active LNV scalar-gluon branch

Let chi_i=tilde H^dagger L_i. A manifestly electroweak-invariant completion is
\[
\mathcal L\supset \frac{c^{(+)}_{ij}}{M^5}
(\overline{\chi_i^c}\chi_j)G^A_{\mu\nu}G^{A\mu\nu}
+\frac{c^{(-)}_{ij}}{M^5}
(\overline{\chi_i^c}\chi_j)G^A_{\mu\nu}\widetilde G^{A\mu\nu}
+\mathrm{h.c.}
\]
Coefficients here omit conventional alpha_s factors and assume an explicitly fixed independent-flavour normalization. After electroweak breaking, g7=c v^2/(2M^5) multiplies the dimension-seven interaction.

Crucially, the scalar two-left-handed-neutrino bilinear is flavour symmetric. The tau-tau term does NOT vanish:
\[
\nu_\tau g\to\bar\nu_\tau g
\]
is possible. It retains an active daughter. This differs both from the Dirac ATZ active-to-sterile interaction and from the off-diagonal-only LNV derivative operator. A Dirac right-handed-singlet parent instead occurs at dimension eight and produces a sterile daughter.

Dimensional analysis gives M_parton ~ g7 shat^{3/2} and sigma_parton ~ |g7|^2 shat^2, before angular factors and validity restrictions. This is genuine momentum growth even though it is slower than active dimension-eight twist-two scattering.

The important restrictions are collider invisible production, neutrino scattering, electroweak completion, perturbative validity and neutrino-mass generation. Tau-only LNV does not inherit an electron-flavour neutrinoless-double-beta bound unchanged. Conversely, it is not exempt from mass/mixing constraints.

An analytic discriminator: the CP-even gluon condensate can contribute a Majorana mass proportional to g7 <G^2>. In a CP-conserving QCD vacuum <G tilde G>=0, so the dual-gluon operator lacks that particular condensate contribution. This does NOT eliminate UV threshold corrections or other mass constraints. It motivates checking the two structures separately rather than excluding or accepting them together.

A scalar/pseudoscalar mediator with neutrino coupling and a colored-particle loop generating S G^2 or a G tilde G is a matching route, not yet a complete benchmark. Resolve both the mediator propagator and the colored-loop form factor when their thresholds are approached. [Farzan et al.](https://arxiv.org/abs/1802.05171) discuss UV constructions for scalar neutrino/quark interactions, but do not certify this particular gluon model.

## Does a UV model help?

Yes, principally by fixing validity, correlations and thresholds—not by guaranteeing larger rates.

For a single ordinary t-channel mediator matched at zero momentum,
\[
G_{\rm eff}=\frac{g_\nu g_q}{M^2},\qquad
\frac{\mathcal M_{\rm full}}{\mathcal M_{\rm contact}}
=\frac{M^2}{M^2+Q^2}
\]
for the same numerator and couplings. The pure contribution is suppressed by the square of this ratio. For this propagator the local expansion parameter is Q^2/M^2, not automatically shat/M^2. A blanket shat cut is a conservative contact-screen prescription, not a necessary-and-sufficient validity condition for every mediator topology. Other heavy vertices, loop thresholds and perturbative unitarity remain independent restrictions.

For s-channel production,
\[
\hat s=2xm_NE_\nu\simeq M_R^2
\]
can instead access a real resonance. Its width must be computed from the same couplings. This genuinely accesses physics omitted by the contact theory. Comparing LHC and neutrino processes then requires their different parton luminosities and branching fractions; transferring a contact bound unchanged is invalid.

A very broad resonance or a nonrenormalizable vertex has not become UV safe merely because a propagator is included. Spin-two models with dimension-five matter couplings likewise need a specified higher cutoff/completion.

## Model shortlist beyond the existing contact calculation

### 1. Flavour-structured charged/neutral scalar sector

A charged scalar with couplings to tau-neutrino and charm-strange currents can produce
\[
\nu_\tau+s\to\tau^-+c
\]
through t-channel exchange (with the appropriate chirality chosen explicitly). This adds a charged-current channel and real taus, unlike an NC-only operator.

[Huang et al., 2112.09476](https://arxiv.org/abs/2112.09476) identify scalar and leptoquark scenarios with interesting tau-neutrino signals under the constraints used in that work. Their scalar example motivates revisiting a gauge-complete flavour texture rather than an arbitrary isolated H^+ vertex. This is not evidence that those historical benchmark points survive 2026 data.

Required first checks: gauge-related neutral scalars, flavour/CKM effects, D_s -> tau nu, tau/meson decays, neutrino masses if relevant, and collider tau+MET/dijet/pair production. [ATLAS 2412.17584](https://arxiv.org/abs/2412.17584) targets top-associated/top-decay production; [ATLAS 2407.10096](https://arxiv.org/abs/2407.10096) targets top-decay H^+ -> cs. Their limits must be applied according to the chosen production couplings, not ignored or copied indiscriminately. This is a promising *model-audit priority*, not a newly established allowed region.

### 2. Tau-coupled scalar leptoquarks R2/S1

They permit resonance and associated-production mechanisms, including gluon-initiated channels, and supply renormalizable gauge-invariant interactions. An R2 coupling to c_R and L3 also generates dimension-six interactions: do not pretend it realizes only the isolated dimension-eight derivative benchmark.

Important contrary evidence: [Kirk, Okawa and Wu, 2307.11152](https://arxiv.org/abs/2307.11152) revisit R2/S1 including collider and previously missed flavour constraints and find collider bounds more restrictive than projected neutrino reach for their scenarios. Appendix A also considers removing the coupling used for flavour anomalies; collider constraints remain important. Therefore leptoquarks are a useful controlled resonance comparison, not the default loophole. An analysis must identify exactly which different flavour texture or channel changes that conclusion, and include current direct searches.

### 3. Consistent vector mediator

Useful for active NC enhancement and for exposing the excessive conservatism of a universal shat cut in some t-channel applications. A large small-y cross section can yield much less energy degradation than its total rate suggests. Derive both sigma and integral y d sigma/dy.

Gauge/anomaly cancellation and charged-lepton couplings cannot be discarded independently. In particular a pure L_mu-L_tau gauge boson has no tree-level quark coupling without mixing or additional matter. Keep neutrino trident, meson, dilepton/dijet and relevant low-energy constraints according to the actual model.

### 4. Tau dipole portal to a massive neutral lepton

Threshold production nu_tau N -> N_heavy X can bypass constraints that assume a light final neutrino and provides a distinct experimental topology. It is not ordinary active NC scattering; include decay, missing energy and regeneration. The electroweak dipole still requires a UV completion. [Barducci et al., PRD 108, 115009](https://doi.org/10.1103/PhysRevD.108.115009) find remaining tau-portal sensitivity, but this does not establish a large allowed UHE inclusive rate. Retain as a secondary threshold candidate.

## Revised research decision

Compare three primary mechanisms analytically before any long computation:
1. existing active gluon twist-two;
2. scalar/dual-gluon LNV, including its tau-diagonal active daughter;
3. a gauge-complete flavour-structured scalar model with a charged-current contribution.

Use leptoquarks and a consistent vector model as controlled resonance/propagator comparisons. Collapse redundant derivative benchmarks using explicit on-shell checks. Do not repeatedly scan the same unverified normalization under different names.

For each, specify states, mediator spectrum, matching, coefficient conventions, invariant validity conditions, widths and the actually relevant laboratory constraints. Only after this analytic envelope leaves a potentially measurable effect should a small PDF fold be authorized. Cross-section enhancement and detector-event enhancement remain separate objectives. This review identifies better questions and candidates; it does not claim a fully allowed new benchmark.
