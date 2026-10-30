# nutau-prop v0.33 — consolidated UHE-interaction portfolio

Date: 2026-10-30  
Status: analytic and evidence audit; no new numerical cross-section or transport run in this checkpoint.

## Purpose and decision standard

This note reconciles `STERILE_UHE_SCREEN_v0.30.md`,
`STERILE_UHE_CERTIFICATION_v0.31.md`, the corrected
`STERILE_D6_RESOLVED_CERTIFICATION_v0.32.md`, the `chat-nutau3`
discussion, and the other operator/UV-model audits through repository commit
`8ed593c691d38e6f0a932acbb554a10152338ea4`.

The target is an interaction that changes the inclusive UHE neutrino--nucleon
cross section by at least 5% around 100--300 PeV while satisfying electroweak
gauge invariance, partial-wave/EFT control, applicable collider/flavour/neutrino
constraints, and the correct daughter topology in transport and detector
folding.

Four labels are kept separate:

- **SIZE PASS:** a calculated benchmark exceeds 5%;
- **CONDITIONAL:** a binding constraint or UV completion remains unevaluated;
- **STRONG/ENDPOINT:** the result needs a coefficient near `4 pi` or is
  dominated by the edge of the validity cut;
- **FAILED:** the best controlled/current-bound result is below 5%, or the
  proposed enhancement disappears after the correct treatment.

No interaction is yet certified against every requirement.

## Consolidated ranking

### Priority 1 — resolved derivative active-to-sterile mediator (non-resonant)

\[
\nu_\tau q\to N q,
\qquad
\mathcal L\supset
\kappa Z'_\mu(\bar N i\overleftrightarrow{\partial^\mu}P_L\nu_\tau)
+g_q Z'_\mu\bar q\gamma^\mu q.
\]

The exact v0.32 dynamic-\(Q\) fold retains
\((M_{Z'}^2+Q^2)^{-2}\) and gives, at the historical PKR normalization
\(\kappa=2.5\times10^{-4}\,\mathrm{GeV}^{-1}\):

| \(M_{Z'}\) | 100 PeV | 300 PeV |
|---:|---:|---:|
| 100 GeV | 45.7% | 50.6% |
| 220 GeV | 19.3% | 22.4% |
| 500 GeV | 6.26% | 7.60% |
| 1 TeV | 1.86% | 2.33% |

The coupling needed for 5% is approximately
\(\kappa=(0.8--2.2)\times10^{-4}\,\mathrm{GeV}^{-1}\) for
\(M_{Z'}=100--500\) GeV. This is the only calculated non-resonant branch
with a broad resolved-mediator size window.

**Why it remains first:** it has the desired sterile sink topology, a resolved
propagator, and exceeds 5% away from the old contact-only extrapolation.

**Why it is not certified:** the derivative neutrino vertex is itself
higher-dimensional/loop-generated in weakly coupled renormalizable completions;
the 2019 PKR IceCube normalization is not a modern mass-dependent bound; and
monojet constraints are competitive or stronger over part of the 40--220 GeV
region. A tau-specific monojet plus IceCube recast and an explicit UV origin of
the derivative vertex are decisive.

**Verdict:** most promising non-resonant candidate; **SIZE PASS, CONDITIONAL**.

### Priority 2 — tau-aligned scalar leptoquark resonance

The cleanest renormalizable comparator is a strange-aligned
\(\widetilde R_2\sim(3,2,1/6)\) coupling
\(\bar s_R\widetilde R_2^T i\tau_2L_\tau\). It produces an on-shell
\(\nu_\tau\bar s\) resonance and has an active NC daughter. The repository
narrow-width screen gives, for \(M=1.5\) TeV and \(\lambda=2\),

\[
\Delta\sigma/\sigma_{\rm SM}=7.17\%\;(100\,\mathrm{PeV}),
\qquad 9.71\%\;(300\,\mathrm{PeV}).
\]

This is the most conventional gauge-invariant, renormalizable interaction in
the portfolio, but it is resonant rather than PKR-like. Finite width, the
charged \(s\tau\) partner, pair and single production, ditau tails,
jet+missing-energy searches, and loop/flavour constraints must be combined.
Existing R2/S1 studies find collider constraints dominant for the textures
they examined.

**Verdict:** best UV-clean resonant candidate; **SIZE PASS, CONDITIONAL**.
Use \(S_1\) only when a tau-producing CC channel is wanted and include its
CKM-linked charged-current constraints.

### Priority 3 — tau-neutrinogluon resonance

For a vector-like colour-octet doublet with the dimension-five transition
\(\bar L_\tau\sigma^{\mu\nu}G^A_{\mu\nu}L^A_{8R}\), the repository benchmark
\(M=3\) TeV, \(\Lambda=9\) TeV, \(a=4\) gives

\[
6.14\%\;(100\,\mathrm{PeV}),
\qquad 12.95\%\;(300\,\mathrm{PeV}).
\]

The gluon luminosity makes this the strongest stored resonant size benchmark.
Its vector-like mass sector is gauge invariant and anomaly free, but the
transition dipole is not a renormalizable UV completion and \(a=4\) is large
for a loop-generated dipole. Charged-octet partner searches, tau+jet,
missing-energy+jet, nonresonant ditau, validity and detector topology remain
open.

**Verdict:** strong size candidate but less UV-clean than the scalar
leptoquark; **SIZE PASS, CONDITIONAL/strong matching**.

### Priority 4 — strong/composite non-resonant benchmarks

Two sterile-sink operators pass only near their strong-coupling endpoint:

| branch | 100 PeV | 300 PeV | reason not promoted |
|---|---:|---:|---|
| derivative-tensor sterile | 5.52% | 7.32% | \(C_5\Lambda^3\simeq4\pi\); weakly coupled UV is loop suppressed |
| sterile gluon/Rayleigh | 3.11% | 6.42% | endpoint only; gauge-invariant field-strength parent is loop generated |

They remain useful composite/strong-UV targets and attenuation benchmarks, but
must not be advertised as generic perturbative EFT successes.

**Verdict:** **SIZE PASS only at STRONG/ENDPOINT**.

### Priority 5 — resolved charged scalar, retained as a flavour-model challenge

For \(M_H=90\) GeV and the tuned example \(g_\tau=1\),
\(g_L=g_R=0.4\), the stored screen gives about 7.6% and 7.4% at 100 and
300 PeV. The equality cancels the simple \(D_s\to\tau\nu\) pseudoscalar
combination, but a gauge-complete second-doublet texture produces CKM-related
\(\tau\to K\nu\), neutral-scalar interactions, and collider/electroweak
constraints. The cancellation is not a certificate.

**Verdict:** size works, but gauge/flavour correlations make it lower priority;
**CONDITIONAL, likely tuned**.

## Branches that failed the project gate

### Gauge-invariant sterile dimension-six scalar/tensor contacts

This is the decisive v0.32 correction to v0.30/v0.31. With current HighPT
90% bounds and the weak-coupling \(c=1\) domain, the exact CT18/NNPDF folds are

| operator | 100 PeV | 300 PeV |
|---|---:|---:|
| scalar \(O_{lNqd}^{S}[\tau23]\) | 0.0055--0.0057% | 0.0067--0.0069% |
| tensor \(O_{lNqd}^{T}[\tau23]\) | 0.540--0.568% | 0.678--0.706% |

Even the `4 pi` contact domain gives only about 1.0%/1.8% for the CT18 tensor
benchmark. The preliminary 8--11% tensor result was a factor-16 convention
error from mixing low-energy \(\epsilon_{S,T}\) and SMEFT normalizations; it
is withdrawn. There is no independent unsuppressed active-to-sterile vector
transition at dimension six.

**Verdict:** **FAILED**. Do not run transport. A colourless scalar-doublet UV
completion of the scalar operator cannot rescue the size at fixed bounds: its
resolved spacelike propagator suppresses the contact envelope, while its
charged partner adds constraints.

### Active quark twist-two \(d=8\), \(c_R/s_R\)

After exact interference, positivity, unitarity and dynamic-\(Q\) folding, the
largest total-cross-section shifts are below 0.25% through 300 PeV. Collider
bounds can only reduce them.

**Verdict:** **FAILED/CLOSED** for the 5% goal.

### Active gluon twist-two \(d=8\)

At \(C_G=-8\pi\), \(\Lambda=4\) TeV it reaches 7.32% only at 300 PeV with the
maximal \(\hat s<\Lambda^2\) domain. Tightening to
\(\hat s<0.5\Lambda^2\) gives 1.48%; 72--77% of the maximal signal is generated
above \(\hat s/\Lambda^2=0.5\). The collider recast remains incomplete, and
retaining 5% requires \(|C_G|\gtrsim20.8\), close to the unitarity ceiling.

**Verdict:** not robust; **ENDPOINT/COLLIDER-CONDITIONAL**.

### Sterile electroweak dipole portal

The photon/\(Z\) propagator drives the fold toward low \(Q^2\), preventing
sustained hard-DIS growth; transition-moment constraints further reduce the
coefficient.

**Verdict:** **FAILED** for an inclusive >5% hard-UHE effect; retain only for
low-recoil or distinctive decay signatures.

### Universal \(U(1)_{B-3L_\tau}\) vector mediator (additional suggestion)

This is a genuine anomaly-free, renormalizable non-resonant model, but the
oscillation bound implies approximately
\(g_X/M_X\lesssim4.7\times10^{-4}\,\mathrm{GeV}^{-1}\), corresponding to a
small neutral-current coupling shift. Its universal quark charge also causes
substantial cancellation between the SM up- and down-quark NC interference
terms. The analytic expectation is sub-percent in the inclusive CC+NC cross
section, not >5%.

**Verdict:** theoretically clean but **FAILED on size**. No expensive fold is
justified unless a non-universal anomaly-free charge assignment with a complete
spectator/collider analysis is supplied.

### Literal PKR/ATZ contact interpretations

The old Case-6 hard kernel demonstrates momentum growth, not a fully certified
theory. ATZ's Dirac derivative branch produces a right-handed/sterile daughter
and is a sink. PKR Eq. (4.7) read literally is off-diagonal/LNV and requires
flavour, antineutrino and secondary-species transport. Neither should be
silently interpreted as a diagonal active NC operator, and the historical
collider coefficient cannot be transferred without matching the external
states and normalization. The useful part of this branch is now represented
by the resolved derivative sterile model above.

## Final decision matrix

| Rank | Interaction | Mechanism | 5% size | Present status |
|---:|---|---|---|---|
| 1 | resolved derivative \(Z'\), \(\nu_\tau q\to Nq\) | non-resonant sterile sink | yes, 5--500 GeV window | best non-resonant; modern recast + UV vertex required |
| 2 | strange-aligned \(\widetilde R_2\) | renormalizable quark resonance | yes for strong benchmark | best UV-clean resonance; collider/width audit required |
| 3 | tau-neutrinogluon | gluon resonance | 6.1--13.0% | strong size; dipole UV and coloured-partner constraints open |
| 4 | derivative-tensor / sterile Rayleigh | contact sterile sink | yes only near \(4\pi\) | composite/strong benchmark only |
| 5 | charged scalar | resolved CC | yes in tuned example | gauge/flavour viability unresolved |
| -- | sterile \(d=6\) S/T | non-resonant contact | no: <=0.71% controlled | failed |
| -- | active \(c_R/s_R\) \(d=8\) | non-resonant contact | no: <0.25% total | closed |
| -- | active gluon \(d=8\) | non-resonant contact | only endpoint at 300 PeV | not robust |
| -- | sterile EW dipole | low-\(Q\) exchange | no | failed for this goal |
| -- | \(B-3L_\tau\) vector | renormalizable \(t\)-channel | analytic estimate <1% | failed on size |

## Next justified work

No further broad operator scan or Earth-transport run is justified now.

1. For the resolved derivative mediator, specify a UV completion of the
   derivative vertex and obtain a mass-dependent tau-specific monojet plus
   IceCube bound \(\kappa_{\max}(M)\). Compare it directly with the stored
   \(\kappa_{5\%}(M)\) curve.
2. In parallel, perform one current constraint audit for the strange-aligned
   \(\widetilde R_2\) benchmark, including finite width and both weak-doublet
   components. It is the cleanest fallback if resonance is acceptable.
3. Retain the tau-neutrinogluon model only if an explicit matching mechanism
   can justify the large transition coefficient after coloured-partner bounds.
4. Run Earth attenuation only after one benchmark satisfies
   \(\kappa_{\max}>\kappa_{5\%}\) or an equivalent resonance constraint gate.

The present honest conclusion is:

\[
\boxed{
\text{A plausible non-resonant window exists only for the resolved derivative
sterile mediator, but no all-constraints survivor is yet certified.}
}
\]
