# nutau-prop v0.34 — authoritative UHE interaction portfolio

Date: 2026-10-30
Status: authoritative portfolio checkpoint, updated through v0.35 self-consistent derivative gate.

> **Audit policy:** `AUDIT_PORTFOLIO_v0.33_REVIEW.md` is **SUPERSEDED—DO NOT USE**. Preserve that file and its failed reasoning if/when it is present in the repository; do not delete it. At this commit the named file is not present on the default GitHub branch, so this notice is recorded here until the local-only copy is pushed and can be bannered directly.

## Decision standard

Target: an interaction that changes inclusive UHE neutrino–nucleon scattering by at least 5% around 100–300 PeV while remaining consistent with gauge invariance, partial-wave/EFT control, collider/flavour/neutrino bounds, and the correct transport daughter topology.

Labels: SIZE PASS, CONDITIONAL, STRONG/ENDPOINT, FAILED/CLOSED.

## Current leading candidate

### 1. Resolved derivative active-to-sterile Z' interaction — PRIORITY 1

PKR Eq. (4.4):

\[
\mathcal L\supset {c_2\over\Lambda}(\bar\nu^c i\overleftrightarrow{\partial_\mu}\nu)Z'^\mu+g_q\bar q\gamma_\mu q Z'^\mu,
\qquad \kappa\equiv {c_2g_q\over\Lambda}.
\]

For the sterile interpretation relevant here the daughter is a sterile state N, so the topology is \(\nu_\tau q\to Nq\), an active-flux sink.

Historical PKR anchors:

- at \(M_{Z'}=5\) GeV: \(\kappa_{\rm LHC}<3.3\times10^{-3}\,\mathrm{GeV}^{-1}\);
- at \(M_{Z'}=5\) GeV: \(\kappa_{\rm IceCube}<2.5\times10^{-4}\,\mathrm{GeV}^{-1}\);
- PKR Fig. 2(c): LHC dominates approximately \(40\! -\!220\) GeV;
- PKR Fig. 4: explicit 5-GeV cross-section comparison.

### v0.35 decisive size update

The required coupling-dependent validity condition is now imposed **inside every PDF fold**:

\[
\sqrt{\hat s}<1/\kappa.
\]

No fixed-cut quadratic rescaling is used. The completed CT18 decision scan finds the first 5% crossing:

| \(M_{Z'}\) [GeV] | \(\kappa_{5\%}\), 100 PeV [GeV\(^{-1}\)] | \(\kappa_{5\%}\), 300 PeV [GeV\(^{-1}\)] |
|---:|---:|---:|
| 5 | 3.0943e-5 | 3.1251e-5 |
| 40 | 5.4254e-5 | 5.1974e-5 |
| 100 | 7.9706e-5 | 7.4021e-5 |
| 220 | 1.2008e-4 | 1.0898e-4 |
| 500 | 2.1643e-4 | 1.9065e-4 |
| 1000 | no crossing in scanned controlled domain | no crossing |
| 2000 | no crossing | no crossing |

At 5 GeV the required \(\kappa_{5\%}\simeq3.1\times10^{-5}\) is about eight times below the historical PKR IceCube upper bound and more than two orders below the historical LHC upper bound. Thus the historical 5-GeV constraints do not close the 5% window.

The complete PKR Fig. 2(c) curve has **not** been numerically invented from textual anchors: the exact 5-GeV values and the stated 40–220 GeV crossover are verified; a point-by-point digitization remains a separate figure-reconstruction task.

Modern data identified for the next constraint gate:

- ATLAS full Run-2 monojet/jets+MET, 13 TeV, 139 fb^-1;
- IceCube modern all-flavour neutrino–nucleon cross-section measurements through the multi-PeV range.

Neither publishes a ready-made bound for this exact derivative active-to-sterile operator. Therefore modern monojet and IceCube limits remain **screening/recast pending**, not certified curves. Signal acceptance and likelihood reproduction are required.

Keep independently: \(c_2/\Lambda\), \(g_q\), \(\Lambda\), \(M_{Z'}\), \(\Gamma_{Z'}\), branching fractions and direct-Z' constraints. The one-parameter \(\kappa\) scan is a cross-section diagnostic, not a full UV model.

**Why it remains #1:** it is the only non-resonant branch with (i) a self-consistent >5% size window, (ii) a resolved propagator, (iii) sterile-sink attenuation topology, and (iv) a historical mass-dependent LHC/IceCube constraint framework that can be recast directly.

Verdict: **MOST PROMISING NON-RESONANT / SIZE PASS / MODERN-CONSTRAINT CONDITIONAL / UV CONDITIONAL**.

Detailed v0.35 record: `RESOLVED_DERIVATIVE_SELFCONSISTENT_v0.35.md`.

## Resonant fallbacks

### 2. Strange-aligned \(\widetilde R_2\sim(3,2,1/6)\)

Renormalizable, gauge-complete quark resonance. Stored benchmark \(M=1.5\) TeV, \(\lambda=2\): about 7.17% at 100 PeV and 9.71% at 300 PeV. Finite width and both weak-doublet components must be included; collider pair/single production and flavour constraints remain binding.

Verdict: **best UV-clean resonant fallback; SIZE PASS, CONDITIONAL**.

### 3. Tau-neutrinogluon / leptogluon resonance

Gauge-invariant colour-octet doublet with transition dipole. Stored benchmark \(M=3\) TeV, \(\Lambda=9\) TeV, \(a=4\): 6.14% at 100 PeV and 12.95% at 300 PeV. Main problems: large dipole matching coefficient and coloured charged-partner limits.

Verdict: **SIZE PASS, CONDITIONAL/strong matching**.

## Non-resonant strong/endpoint benchmarks

- derivative-tensor sterile: 5.52% / 7.32% at 100/300 PeV only near \(C\Lambda^3\simeq4\pi\);
- sterile gluon/Rayleigh: 3.11% / 6.42%, passing only near the endpoint at 300 PeV;
- active gluon twist-2 d=8: can exceed 5% only near the validity/strong-coupling edge and collapses under tighter \(\hat s\) cuts.

Verdict: **retain as strong/composite benchmarks, not perturbative successes**.

## Other conditional branch

### Charged scalar

Stored tuned example gives about 7.6% / 7.4% at 100/300 PeV, but gauge completion, CKM-linked flavour observables, neutral-scalar partners and collider constraints make it lower priority.

Verdict: **CONDITIONAL, likely tuned**.

## Failed / closed branches

1. Gauge-invariant sterile d=6 scalar: controlled current-bound result only ~0.005–0.007%.
2. Gauge-invariant sterile d=6 tensor: controlled current-bound result only ~0.54–0.71%; even the 4π contact domain remains below 5%.
3. Preliminary 8–11% d=6 tensor result: **withdrawn**; factor-16 normalization mistake from mixing low-energy epsilon and SMEFT Wilson conventions.
4. Active quark twist-2 d=8: <0.25% after exact interference/positivity/unitarity/dynamic-Q treatment.
5. Sterile electroweak dipole portal: low-Q dominated; fails the hard-UHE >5% objective.
6. Universal \(U(1)_{B-3L_\tau}\): theoretically clean but sub-percent because allowed coupling is small and up/down NC interference largely cancels.
7. Literal PKR/ATZ contact reinterpretations: useful as hard-kernel diagnostics only; external states/gauge completion and historical collider normalizations prevent direct use as final models.
8. Generic weakly-coupled derivative-tensor/Rayleigh UV completions: loop suppression makes the >5% endpoint benchmarks unnatural without strong/composite dynamics.

## Monojet status

The strong-contact monojet indication is retained only as a **screening estimate** until signal acceptance, cuts, detector response and likelihood are reproduced. Do not label it a recast/certification before that reproduction.

## Remaining justified work

1. Reconstruct/digitize PKR Fig. 2(c) point-by-point without inventing intermediate values.
2. Recast ATLAS 139-fb^-1 monojet/jets+MET for the resolved derivative signal and derive modern \(\kappa_{\max}^{\rm LHC}(M)\).
3. Recast the modern IceCube cross-section/event likelihood to derive \(\kappa_{\max}^{\rm IC}(M)\).
4. Compare those curves with the v0.35 self-consistent \(\kappa_{5\%}(M)\).
5. Only if \(\kappa_{\max}>\kappa_{5\%}\) somewhere, choose a benchmark and proceed to Earth attenuation.

## Current decision

\[
\boxed{\text{Most promising candidate: resolved derivative }Z'\text{ active-to-sterile mediator.}}
\]

This remains a prioritization decision, not an all-constraints certification. If the modern recasts close the 5–500 GeV window, the branch is closed and strange-aligned \(\widetilde R_2\) becomes the lead model.