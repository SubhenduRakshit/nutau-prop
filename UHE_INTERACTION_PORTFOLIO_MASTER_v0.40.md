# nutau-prop v0.40 — authoritative UHE interaction portfolio

Date: 2026-10-30
Status: authoritative portfolio checkpoint through final PKR Eq. (4.3)/(4.4) certification.

> **Audit policy:** `AUDIT_PORTFOLIO_v0.33_REVIEW.md` is **SUPERSEDED—DO NOT USE**. Preserve it and its failed reasoning; do not delete it.

## Decision standard
Target: >=5% modification of inclusive UHE neutrino–nucleon scattering around 100–300 PeV, with gauge/UV interpretation, partonic validity, collider/flavour/neutrino constraints, mediator width and transport topology all kept explicit.

## Key correction introduced in v0.40
For the two PKR dimension-5 resolved vertices,

kappa_i = (c_i/Lambda) gq,

so the actual vertex scale is

Lambda_vert = Lambda/c_i = gq/kappa_i,

not 1/kappa_i unless gq=1. The correct validity condition is therefore

sqrt(shat) < gq/kappa_i.

This correction closes much of the apparent medium/high-mass window once current direct dijet bounds force gq << 1.

## Current leading non-resonant candidate
### PKR Eq. (4.4), resolved derivative Z'
Topology: active-to-sterile/off-diagonal neutrino scattering through a resolved quark-coupled mediator.

Peak-capacity size+validity critical quark couplings at 100 PeV are approximately:

| M [GeV] | gq,crit for any 5% solution |
|---:|---:|
| 5 | 0.0051 |
| 40 | 0.0183 |
| 100 | 0.0527 |
| 220 | 0.1426 |
| 500 | 0.433 |
| 1000 | >1 |
| 2000 | >1 |

Modern direct dijet information closes the 500-GeV and heavier region. The 5–100-GeV region remains genuinely open; 220 GeV is narrow/conditional and requires an operator-specific modern monojet recast.

Downward-energy screen for the 100-PeV 5% benchmark: at M=40 GeV the NP/SM fractions are roughly 2.3%, 3.2%, 3.6%, 4.1% at 0.1,1,3,10 PeV; at M=100 GeV they are ~0.8%,2.0%,2.7%,3.4%. These do not constitute an obvious IceCube cross-section exclusion by themselves.

Verdict: **LEAD NON-RESONANT CANDIDATE / OPEN at ~5–100 GeV / CONDITIONAL near 220 GeV / CLOSED >=~500 GeV.**

Preferred benchmark class for the next real recast: M_Z' ~40–100 GeV with gq explicit, above the size+validity critical value and below direct-dijet bounds.

## Second non-resonant survivor
### PKR Eq. (4.3), resolved transition-dipole Z'
Off-diagonal only (i != j). In the validated massless DIS interpretation,

d2sigma_4.3 = 2 d2sigma_4.4

at equal M and kappa.

Peak-capacity critical quark couplings:

| M [GeV] | gq,crit for any 5% solution |
|---:|---:|
| 5 | 0.0036 |
| 40 | 0.0129 |
| 100 | 0.0373 |
| 220 | 0.1008 |
| 500 | 0.306 |
| 1000 | 0.819 |
| 2000 | >1 |

It is more efficient in raw UHE size than Eq. (4.4), but carries additional transition-moment/kinetic-mixing and width issues. Current direct dijet constraints close the 500-GeV and 1-TeV windows once split-coupling validity is imposed. The 5–100-GeV window remains open; ~220 GeV remains conditional.

Verdict: **STRONGER SIZE FALLBACK / OPEN at ~5–100 GeV / CONDITIONAL near 220 GeV / CLOSED >=~500 GeV.**

## Resonant candidates
### Strange-aligned R2-tilde leptoquark
Renormalizable, gauge-complete resonance. Stored benchmark M=1.5 TeV, lambda=2 gives about 7.2% at 100 PeV and 9.7% at 300 PeV. Finite-width, both weak-doublet components, pair/single production and flavour constraints remain binding.

Verdict: **best UV-clean resonant fallback / CONDITIONAL.**

### Tau-neutrinogluon / leptogluon resonance
Stored benchmark M=3 TeV, Lambda=9 TeV, a=4 gives about 6.1% at 100 PeV and 13% at 300 PeV. Main problems: large chromomagnetic matching coefficient and coloured charged partner limits.

Verdict: **SIZE PASS / strong matching / CONDITIONAL.**

## Strong/endpoint non-resonant benchmarks
- derivative-tensor sterile: >5% only near C Lambda^3 ~4pi;
- sterile gluon/Rayleigh: passes 5% only near strong/endpoint coupling at high energy;
- active gluon twist-2 d=8: >5% only under aggressive validity assumptions and collapses with tighter shat cuts.

Verdict: **retain only as strong/composite benchmarks.**

## Other conditional branch
### Charged scalar
Stored tuned example ~7.6% at 100 PeV and ~7.4% at 300 PeV. Gauge completion/flavour/collider requirements make it lower priority.

Verdict: **CONDITIONAL / tuned.**

## Failed / closed branches
1. Gauge-invariant sterile d=6 scalar: ~0.005–0.007% at current bounds.
2. Gauge-invariant sterile d=6 tensor: ~0.54–0.71%; below 5% even with aggressive contact domain.
3. Preliminary 8–11% d=6 tensor result: **WITHDRAWN** — factor-16 convention error.
4. Active quark twist-2 d=8: sub-percent after exact interference/positivity/unitarity/dynamic-Q treatment.
5. Sterile electroweak dipole portal: low-Q dominated; fails hard-UHE >5% target.
6. Universal U(1)_(B-3Ltau): clean but sub-percent in allowed region.
7. Literal PKR contact reinterpretations: hard-kernel diagnostics only, not final models.
8. Generic weakly-coupled derivative-tensor/Rayleigh UV completions: loop suppression makes endpoint >5% benchmarks unnatural without strong/composite dynamics.
9. PKR Eq. (4.4) and Eq. (4.3) with M_Z' >= ~500 GeV: **CLOSED in perturbative resolved form by direct-gq constraints plus correct split-coupling validity.**

## Modern experimental status
- Full Run-2 ATLAS monojet (139 fb^-1) is the correct collider dataset, but no published bound exists for these exact PKR dimension-5 operators. A modern kappa(M) curve requires actual signal generation, acceptance and likelihood reproduction. Luminosity rescaling is not certification.
- IceCube all-flavour nuN cross-section data cover 60 TeV–10 PeV and are SM-compatible. The v0.37 energy screen shows the 5-GeV points are most exposed; 40–100 GeV benchmarks give only few-percent effects in the measured PeV range. An event/zenith likelihood is still needed.
- Direct dijet/quark-coupling constraints are already sufficient, together with Lambda_vert=gq/kappa, to close the medium/high-mass resolved options.

## Final non-resonant frontier

The earlier near-no-go statement is sharpened to:

**Weakly coupled non-resonant models are effectively closed at medium/high mediator mass. The only remaining perturbative non-resonant loophole is a light/intermediate resolved PKR mediator, principally Eq. (4.4) or Eq. (4.3), with M_Z' roughly 5–100 GeV (and a narrow conditional extension toward ~220 GeV).**

Eq. (4.4) is the preferred lead because it is phenomenologically cleaner. Eq. (4.3) is the stronger raw-size fallback.

## Binding next gate
To finish experimental certification rather than merely model-space certification:
1. generate Eq. (4.4) and Eq. (4.3) signals for M~5,40,100,220 GeV with (gq,c_i/Lambda) separated;
2. reproduce ATLAS 139-fb^-1 monojet acceptance/bin likelihood;
3. reproduce the relevant IceCube event/zenith likelihood for the low-mass benchmarks;
4. apply specialized low-mass direct-Z' / quarkonium constraints and exact Gamma_Z'/M;
5. choose one surviving explicit benchmark.

No Earth propagation before that explicit benchmark passes.

Detailed two-option analysis: `PKR_FINAL_TWO_OPTION_CERTIFICATION_v0.40.md`.
