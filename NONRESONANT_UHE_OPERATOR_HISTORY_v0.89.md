# v0.89 — canonical history of non-resonant UHE operator exploration

Date: 2026-10-31

## Purpose

This file is the canonical project record of the non-resonant / momentum-growing interaction classes explored for the `nutau-prop` project, why each was considered, what was actually calculated, which theoretical/experimental gates were applied, and why each branch was closed, demoted, or left conditional.

It supersedes older ranking statements when they conflict with later normalization or collider corrections, while preserving the historical files and failed approaches.

The working triage target was a modification of the inclusive UHE neutrino-nucleon cross section at the several-percent level around 100--300 PeV, nominally >=5%. This is a screening criterion, not an observability theorem: Earth attenuation can amplify or suppress the detector-level effect, and a final claim must be based on propagated flux and detector response.

Status labels used here:

- **CLOSED**: best controlled/current-bound realization is well below the project target; additional known constraints can only reduce it.
- **STRONG/ENDPOINT ONLY**: >=5% occurs only near strong coupling or near the EFT validity boundary; not a robust perturbative-EFT signal.
- **CONDITIONAL / NOT CERTIFIED**: the hard cross section can be large enough, but a required collider/UV/transport gate is incomplete.
- **WITHDRAWN RESULT**: an earlier numerical result was invalidated by a normalization/convention error.

## Supersession map

1. `STERILE_UHE_SCREEN_v0.30.md` and `UHE_INTERACTION_PORTFOLIO_v0.33.md` are historical screening documents.
2. `STERILE_D6_RESOLVED_CERTIFICATION_v0.32.md` supersedes the preliminary d=6 scalar/tensor numbers and withdraws the old 8--11% tensor result.
3. `D8_CONTACT_BRANCH_CLOSED_v0.27.md` is the detailed closure note for the active c_R/s_R two-derivative d=8 contacts.
4. `GLUON_D8_CERTIFICATION_v0.28.md` and `AUDIT_UHE_OPERATORS_v0.29_2026-10-30.md` define the current status of the active gluon d=8 branch.
5. `UHE_INTERACTION_PORTFOLIO_MASTER_v0.40.md` is useful for the broad portfolio, but its PKR Eq. 4.3/4.4 normalization and ranking are superseded by v0.86--v0.87.
6. `PKR_CONSTRAINT_DISPOSITION_v0.87.md` is the current collider-level disposition for corrected PKR Eq. 4.3/4.4 at M=50,100,220 GeV.

---

# 1. Gauge-invariant d=6 active-to-sterile scalar and tensor contacts

Representative nuSMEFT class: `O_lNqd` scalar/tensor transitions, e.g. tau23 flavour.

The corrected direct-chiral SMEFT kernels used in v0.32 are

\[
\frac{d\hat\sigma_S}{dy}=\frac{G_S^2\hat s\,y^2}{64\pi},
\qquad
\frac{d\hat\sigma_T}{dy}=\frac{G_T^2\hat s(2-y)^2}{2\pi}.
\]

Inputs:

- HighPT 90% bounds used in v0.32:
  - `|C^S_lNqd[tau23]|/Lambda^2 <= 0.41 TeV^-2`,
  - `|C^T_lNqd[tau23]|/Lambda^2 <= 0.25 TeV^-2`.
- weak-coupling contact domain with `sqrt(shat) < 1/sqrt(G)`;
- exact dynamic-Q CT18NLO and NNPDF31 folds.

Results:

| operator | PDF | 100 PeV | 300 PeV |
|---|---|---:|---:|
| scalar tau23 | CT18NLO | 0.00550% | 0.00667% |
| scalar tau23 | NNPDF31 | 0.00573% | 0.00689% |
| tensor tau23 | CT18NLO | 0.540% | 0.678% |
| tensor tau23 | NNPDF31 | 0.568% | 0.706% |

Even the much looser `4 pi` contact-domain exercise gives only about 1.0% at 100 PeV and 1.8% at 300 PeV for the tensor benchmark.

**Why it failed:** the amplitude grows only as a conventional d=6 contact interaction and present collider bounds force the coefficient too small. PDF choice changes the result by only a few percent relative, far too little to rescue it.

**Important withdrawn result:** an earlier 8--11% tensor estimate was wrong by a factor 16 because a low-energy GNI epsilon normalization containing `2 sqrt(2) G_F` was mixed with the SMEFT Wilson coefficient convention. That result is withdrawn permanently.

**Verdict:** **CLOSED** for the >=5% UHE objective.

Evidence: `STERILE_D6_RESOLVED_CERTIFICATION_v0.32.md`, results under `results/v0.32/`.

---

# 2. Active two-derivative quark d=8 contacts: c_R and s_R

Operators:

\[
Q_{l^2u^2D^2}^{3322}=(\bar L_3\gamma_\mu i\overleftrightarrow D_\nu L_3)
(\bar c_R\gamma^\mu i\overleftrightarrow D^\nu c_R),
\]

\[
Q_{l^2d^2D^2}^{3322}=(\bar L_3\gamma_\mu i\overleftrightarrow D_\nu L_3)
(\bar s_R\gamma^\mu i\overleftrightarrow D^\nu s_R).
\]

For an unpolarized right-handed quark PDF:

\[
\frac{d\hat\sigma_q}{dy}=\frac{C_q^2\hat s^3}{8\pi\Lambda^8}(2-y)^2(1-y)^2,
\]

\[
\frac{d\hat\sigma_{\bar q}}{dy}=\frac{C_q^2\hat s^3}{8\pi\Lambda^8}(2-y)^2.
\]

The antiquark J=0 partial wave binds:

\[
|C_q|\left(\frac{\hat s}{\Lambda^2}\right)^2\le \frac{8\pi}{3}.
\]

Requiring unitarity through `shat=Lambda^2` gives

\[
|C_q|\le 8\pi/3=8.37758.
\]

Forward positivity fixes `C_q <= 0`. The strict benchmark used was

\[
C_q=-8\pi/3,\qquad \Lambda=4\ \mathrm{TeV}.
\]

The exact CT18NLO dynamic-Q fold included the pure term, SM-Z interference and the partonic EFT cut.

For the most permissive `xi=1`:

| operator | E | Delta sigma_NC / sigma_NC | Delta sigma_tot / sigma_tot |
|---|---:|---:|---:|
| c_R | 100 PeV | 0.4269% | 0.1196% |
| c_R | 300 PeV | 0.8071% | 0.2277% |
| s_R | 100 PeV | 0.3580% | 0.1003% |
| s_R | 300 PeV | 0.6501% | 0.1834% |

Tighter validity cuts reduce these values further.

**Why it failed:** once the Wilson coefficient is fixed by exact partial-wave unitarity and positivity, the apparently rapid `shat^2/Lambda^4` amplitude growth does not translate into a large inclusive UHE rate after the PDF fold and validity cut. The total cross-section shift stays below 0.25% through 300 PeV.

**Verdict:** **CLOSED** for these active L3-c_R/s_R single-coefficient directions.

Evidence: `D8_CERTIFICATION_STATUS_v0.26.md`, `D8_CONTACT_BRANCH_CLOSED_v0.27.md`, `results/v0.25/d8_uhe_dynamicQ_v0_25.csv`.

---

# 3. Active gluon d=8 twist-2 / energy-momentum-tensor family

Operator:

\[
Q_G^{33}=i(\bar L_3\gamma^\mu\overleftrightarrow D^\nu L_3)
G^A_{\mu\rho}G^{A\rho}{}_{\nu}.
\]

Independent Jacob--Wick projection gives the binding helicity partial wave

\[
a_{1/2}=-\frac{C_G}{16\pi}\left(\frac{\hat s}{\Lambda^2}\right)^2,
\]

hence

\[
|C_G|\left(\frac{\hat s}{\Lambda^2}\right)^2\le 8\pi.
\]

The benchmark used `C_G=-8pi`, `Lambda=4 TeV`.

Exact CT18NLO fold:

| xi in shat < xi Lambda^2 | E | Delta sigma_tot / sigma_tot |
|---:|---:|---:|
| 1.0 | 100 PeV | 3.46% |
| 1.0 | 300 PeV | 7.32% |
| 0.5 | 100 PeV | 0.871% |
| 0.5 | 300 PeV | 1.48% |
| 0.25 | 300 PeV | 0.314% |

At `xi=1`, roughly 72--77% of the stored 100--300 PeV signal comes from `shat/Lambda^2 > 0.5`; 93--95% comes from `>0.25`.

**Why it did not become a robust solution:** the only point exceeding 5% is the most aggressive full-domain endpoint benchmark at 300 PeV. The signal collapses under modestly tighter EFT cuts. A collider recast was attempted but the v0.29 Contur workflow failed at environment initialization before the event/recast stage, so this branch is not experimentally excluded by that failure.

A threshold near `|C_G|~20.8` is required for a 5% total-CS effect at the stored 4-TeV/300-PeV/xi=1 benchmark, already close to the unitarity ceiling `8pi`.

**Verdict:** **STRONG/ENDPOINT ONLY / NOT FULLY CERTIFIED**, not a robust weakly-coupled nonresonant solution. It should not be described as universally excluded.

Evidence: `GLUON_D8_CERTIFICATION_v0.28.md`, `AUDIT_UHE_OPERATORS_v0.29_2026-10-30.md`, `results/v0.27/d8_gluon_uhe_v0_27.csv`.

---

# 4. Sterile derivative-tensor momentum-growing contact

Historical PKR Case-5 / ATZ-Q10-like transition structure, interpreted as an active-to-sterile or chirality-flipped daughter.

The repository endpoint screen gave approximately

- 100 PeV: 5.52%,
- 220 PeV: 6.78%,
- 300 PeV: 7.32%,

with mean y around 0.5.

However the benchmark reaches this only near the perturbative/contact endpoint, summarized in the portfolio as `C Lambda^3 ~ 4pi`.

**Why it was demoted:** a weakly coupled UV completion naturally generates the corresponding higher-dimensional tensor/derivative interaction with loop and/or mediator suppression; the >=5% point relies on a coefficient close to the strong-coupling ceiling and has unresolved gauge-completion/collider issues.

**Verdict:** **STRONG/ENDPOINT ONLY**. Retain as a composite/strong-UV benchmark, not a generic perturbative EFT success.

Evidence: `STERILE_UHE_SCREEN_v0.30.md`, `UHE_INTERACTION_PORTFOLIO_v0.33.md`, `UHE_INTERACTION_PORTFOLIO_MASTER_v0.40.md`.

---

# 5. Sterile gluon / Rayleigh transition

Representative gauge-invariant parent:

\[
\frac{C_G}{\Lambda^4}(\bar N\tilde H^\dagger L_\tau)
G^A_{\mu\nu}G^{A\mu\nu},
\]

or the CP-odd dual-field-strength analogue.

At the repository perturbative/contact endpoint the stored screen gave

- 100 PeV: 3.11%,
- 220 PeV: 5.38%,
- 300 PeV: 6.42%,

with a hard inelasticity distribution, mean y about 0.8.

A separate common-scale normalization example gave only about 0.71% at 100 PeV and 1.42% at 300 PeV.

**Why it was demoted:** the >5% result appears only near the strong/contact endpoint. Gauge-invariant field-strength parents are naturally higher-dimensional and often loop generated; a generic weakly coupled common-scale realization is much smaller.

**Verdict:** **STRONG/ENDPOINT ONLY**. Useful as a strong/composite attenuation benchmark.

Evidence: `STERILE_UHE_SCREEN_v0.30.md`, v0.29/v0.33 portfolio notes.

---

# 6. Sterile electroweak dipole portal

Gauge-invariant operators:

\[
O_{NB}=(\bar L_\tau\sigma^{\mu\nu}N)\tilde H B_{\mu\nu},
\qquad
O_{NW}=(\bar L_\tau\sigma^{\mu\nu}\tau^I N)\tilde H W^I_{\mu\nu}.
\]

These generate `nu_tau q -> N q` through photon/Z exchange.

**Why it failed:** although the dipole vertex carries momentum, the t-channel propagator drives the integral toward low Q^2. It therefore does not produce the sustained hard-DIS polynomial growth needed for a >=5% inclusive 100--300 PeV effect at realistic transition-moment coefficients. Existing dipole constraints further reduce the available coefficient.

**Verdict:** **CLOSED for the hard-UHE inclusive objective**; still potentially interesting for low-Q or decay signatures.

Evidence: `STERILE_UHE_SCREEN_v0.30.md`, portfolio summaries.

---

# 7. Universal U(1)_(B-3L_tau) vector mediator

This is a theoretically clean anomaly-free non-resonant vector model.

The allowed `g_X/M_X` is small, and universal quark charges produce substantial cancellation between SM up- and down-quark neutral-current interference terms.

Stored portfolio conclusion: inclusive UHE effect is sub-percent in the allowed region.

**Why it failed:** collider/oscillation constraints force too small a coupling and the interference structure further suppresses the total effect.

**Verdict:** **CLOSED on size** for the project objective.

Evidence: `UHE_INTERACTION_PORTFOLIO_v0.33.md`, `UHE_INTERACTION_PORTFOLIO_MASTER_v0.40.md`.

---

# 8. Resolved PKR derivative and transition-dipole Z' interactions

These are not ordinary local contact operators because the mediator propagator is kept explicitly, but they were the main nonresonant loophole after the contact scans.

The current corrected hard kernels under matched conventions are

\[
\overline{|M_{4.4}|^2}=4\frac{\kappa^2}{(M^2+Q^2)^2}\hat s^3y(1-y),
\]

and

\[
\overline{|M_{4.3}|^2}=16\frac{\kappa^2}{(M^2+Q^2)^2}\hat s^3y(1-y),
\]

so at equal kappa the corrected Eq. 4.3 rate is four times Eq. 4.4. Older repository statements using a factor 2 are superseded.

The v0.87 collider-compatible maxima at 100 PeV, using exact CMS dijet inputs plus the calibrated monojet implementation, are:

| M | Eq. 4.4 max | Eq. 4.3 max |
|---:|---:|---:|
| 50 GeV | 0.3395% | 0.7953% |
| 100 GeV | 0.0903% | 0.1084% |
| 220 GeV | 0.0456% | 0.0581% |

At the same 100-PeV maximizing couplings, the 220-PeV values are approximately

- M=50 GeV: 0.3373% / 0.7578%,
- M=100 GeV: 0.0951% / 0.1162%,
- M=220 GeV: 0.0497% / 0.0647%,

for Eq. 4.4 / Eq. 4.3 respectively.

Direct 5% tests:

- Eq. 4.4 cannot reach 5% at 50,100,220 GeV after exact dijet+width screens even before monojet: maxima 2.729%, 1.581%, 0.3615%.
- Eq. 4.3 can mathematically reach 5% at 50 and 100 GeV after dijet+width, but the minimum calibrated monojet r values are 5.832 and 54.71, so those target points fail the collider screen.
- At 220 GeV Eq. 4.3 reaches only 2.483% even without monojet.

**Why this branch presently fails at the tested masses:** corrected normalization plus direct dijet limits plus monojet acceptance reduce the available coupling far below what is needed for a several-percent UHE effect.

**Important limitation:** v0.87 does **not** close 5 GeV, 40 GeV, every other mediator mass, every flavour assignment, or every UV/cutoff prescription. The full IceCube detector reweight also lacks a public true-y response needed for an exact arbitrary-y signal fold. Therefore the current claim must remain mass-restricted.

**Verdict:** **CLOSED for the current 50/100/220-GeV implementations; low-mass 5/40-GeV region not fully recertified.**

Evidence: `PKR_CONSTRAINT_DISPOSITION_v0.87.md`; older v0.40 PKR ranking is superseded where inconsistent.

---

# 9. Literal contact reinterpretations of PKR/ATZ momentum-growing kernels

These were useful diagnostics for energy scaling but are not independent fully certified models.

Problems include external-state interpretation, active vs sterile daughter, charge-conjugated fields, gauge completion, validity scale, and transfer of historical collider bounds between conventionally different operators.

**Verdict:** diagnostic only. Do not use old contact numbers as an allowed physical benchmark.

---

# 10. Why simply going to d=10 or d=12 does not automatically help

For a generic local dimension-d operator with coefficient `C_d/Lambda^(d-4)`, the amplitude scales schematically as

\[
M_d\sim C_d\left(\frac{\sqrt{\hat s}}{\Lambda}\right)^{d-4}.
\]

Going from d to d+2 introduces another factor approximately

\[
\frac{\hat s}{\Lambda^2}.
\]

Inside a controlled EFT domain this factor is below one. Therefore higher dimension gives faster formal energy growth only as the partonic energy approaches the cutoff, exactly where the derivative expansion becomes least reliable. A new d=10/d=12 operator is not automatically a loophole unless it gains an independent enhancement: large gluon luminosity, special flavour alignment, symmetry-forbidden lower-dimensional operators, strong/composite matching, or a resolved mediator.

This is why the project should not claim a mathematical no-go theorem over all possible higher-dimensional operators. The evidence instead supports a practical statement: **the ordinary weakly coupled nonresonant structures tested so far do not give a robust several-percent UHE effect.**

---

# Compact decision table

| interaction class | 100 PeV | 300 PeV | current reason it does not provide a robust weakly-coupled solution | status |
|---|---:|---:|---|---|
| sterile d=6 scalar | ~0.0055% | ~0.0067% | collider-bound coefficient too small | CLOSED |
| sterile d=6 tensor | ~0.54% | ~0.68% | below target even at current bound; 4pi exercise still <2% at 300 PeV | CLOSED |
| active d=8 quark c_R | 0.1196% total | 0.2277% | unitarity+positivity+PDF fold suppress signal | CLOSED |
| active d=8 quark s_R | 0.1003% | 0.1834% | same | CLOSED |
| active d=8 gluon, xi=1 | 3.46% | 7.32% | only passes at UV endpoint; strong cutoff dependence; collider gate incomplete | STRONG/ENDPOINT / CONDITIONAL |
| active d=8 gluon, xi=1/2 | 0.871% | 1.48% | controlled-domain signal too small | below target |
| sterile derivative tensor | 5.52% | 7.32% | >=5% only near ~4pi/endpoint | STRONG/ENDPOINT ONLY |
| sterile Rayleigh/gluon | 3.11% | 6.42% | >=5% only near endpoint; weakly coupled common-scale example much smaller | STRONG/ENDPOINT ONLY |
| sterile EW dipole | <target | <target | low-Q dominated; transition constraints | CLOSED for hard-UHE goal |
| B-3L_tau vector | sub-percent | sub-percent | allowed coupling small + interference cancellation | CLOSED |
| corrected PKR Eq.4.4 | 0.34%,0.09%,0.046% at 50,100,220 GeV | similar sub-percent at same optima | dijet+monojet+normalization | CLOSED at tested masses |
| corrected PKR Eq.4.3 | 0.80%,0.11%,0.058% at 50,100,220 GeV | similar sub-percent at same optima | 5% points fail monojet or size | CLOSED at tested masses |

---

# What is not yet justified to claim

Do **not** write that all possible nonresonant BSM interactions are excluded.

The defensible current statement is:

> After explicit tests of gauge-invariant d=6 sterile scalar/tensor contacts, active d=8 quark and gluon momentum-growing operators, sterile derivative-tensor and Rayleigh/gluon transitions, electroweak dipoles, a clean B-3L_tau vector model, and corrected PKR resolved derivative/dipole interactions, no weakly coupled nonresonant benchmark studied so far gives a robust several-percent UHE effect while simultaneously remaining safely inside EFT/partial-wave control and the implemented collider constraints. The remaining >=5% contact examples are strong/endpoint dominated; the corrected PKR branch is sub-percent at the tested 50/100/220-GeV masses. Low-mass unresolved mediator windows and specially protected higher-dimensional UV structures are not ruled out by this statement.

This is a project-level near-no-go conclusion, not a theorem over all EFT operator bases or UV completions.

---

# Resonant models are tracked separately

The following are not counted as nonresonant failures:

- strange-aligned `R2-tilde` leptoquark resonance: stored benchmark about 7.2% / 9.7% at 100/300 PeV, constraints still to be fully certified;
- tau-neutrinogluon / leptogluon resonance: stored benchmark about 6.1% / 13% at 100/300 PeV, with strong matching and coloured-partner issues;
- tuned charged-scalar resolved branch: about 7.6% / 7.4%, gauge/flavour viability unresolved.

These should be treated in their own certification notes rather than used to weaken the nonresonant conclusion.
