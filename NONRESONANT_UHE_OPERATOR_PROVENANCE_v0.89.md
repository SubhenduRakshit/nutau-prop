# v0.89 — provenance map for the nonresonant UHE-neutrino operator programme

Date: 2026-10-31

## Purpose

`NONRESONANT_UHE_OPERATOR_HISTORY_v0.89.md` is the concise synthesis.  This file
is the **provenance map**: it identifies where the detailed derivations,
normalization checks, numerical folds, failed branches, superseded statements,
code and machine-readable outputs are preserved.

Nothing listed here is to be deleted merely because a later checkpoint
supersedes its conclusion.  Historical failures and normalization mistakes are
part of the project record.

The current decision target throughout the programme was a several-percent
modification of the inclusive UHE neutrino-nucleon interaction around
100--300 PeV, with the final judgement required to respect gauge/UV structure,
partial-wave/EFT validity, collider/flavour/neutrino constraints and the
correct Earth-transport daughter topology.

---

# 1. Global development history

## `PROGRESS.log`

This is the chronological execution log beginning with the original transport
validation.  It preserves:

- the Rakshit--Reya/Dutta--Reno--Sarcevic historical validation sequence;
- the original BSM charge-radius/derivative scans;
- detector observability stages;
- later operator scans;
- corrections to transport, normalization and interpretation;
- PKR v0.86 normalization-resolution entry.

Use it when the question is *why did the project move from one model to the
next?* rather than merely *what is the final result?*

---

# 2. Active d=8 quark derivative/twist-two branch

## Main detailed derivation

### `D8_CERTIFICATION_STATUS_v0.26.md`
Contains:

- exact operators
  `l2u2D22[3322]` and `l2d2D22[3322]`;
- exact q and qbar differential partonic cross sections;
- angular integrals `31/30` and `7/3`;
- Jacob--Wick/Legendre partial-wave decomposition;
- binding antiquark J=0 channel;
- strict unitarity bound `|C| <= 8 pi/3` for unitarity through
  `shat=Lambda^2`;
- positivity sign `C<=0`;
- discussion of Drell--Yan/HighPT and why the numerical UHE fold is the
  decisive gate;
- exact benchmark definition `C=-8pi/3, Lambda=4 TeV`.

### `D8_CONTACT_BRANCH_CLOSED_v0.27.md`
This is the detailed closure note.  It preserves the normalization
reconciliation including the otherwise easy-to-miss initial-spin factor, the
fixed-helicity amplitude reconstructed from the unpolarized cross section,
partial waves, positivity, dynamic-Q CT18 fold and the numerical closure:

- cR: total shift +0.1196% at 100 PeV, +0.2277% at 300 PeV;
- sR: total shift +0.1003% at 100 PeV, +0.1834% at 300 PeV.

It explicitly states that this closes those two active L3-cR/sR derivative
contact directions, **not all momentum-growing new physics**.

## Numerical code/results

- `scripts/d8_uhe_cert_v0_25.py`
- `scripts/run_d8_cert_v0_26.py`
- `scripts/highpt_d8_tautau_v0_26.wls`
- `results/v0.25/d8_uhe_dynamicQ_v0_25.csv`

The dynamic-Q fold includes exact SM-Z interference, the pure d=8 term,
`Q^2=xys`, and the partonic validity cut `shat<xi Lambda^2`.

---

# 3. Active d=8 gluon operator

## Main detailed derivations

### `GLUON_D8_SCREEN_v0.27.md`
Contains:

- the Murphy-basis operator
  `Q_l2G2D^33 = i(L3bar gamma^mu <->D^nu L3) G^A_{mu rho} G^{A rho}_nu`;
- the QCD-coupling convention and why there is no explicit gs in the
  two-external-gluon contact vertex in the canonical basis;
- exact unpolarized matrix element and `d sigmahat/dy`;
- explicit gluon-helicity amplitudes;
- partial waves J=1/2 and J=3/2;
- strict bound `|C_G| rho^2 <= 8pi`;
- positivity sign `C_G<=0`;
- CT18NLO dynamic-Q fold and endpoint-support diagnostics.

### `GLUON_D8_CERTIFICATION_v0.28.md`
Independent spin/color/unitarity certification.  Preserves:

- color-space eigenvalue check (no spurious factor 8 or sqrt(8));
- polarization normalization check;
- corrected choice of **total CC+NC** denominator for the observability proxy;
- endpoint sensitivity.

Canonical stored benchmark `C_G=-8pi, Lambda=4 TeV`:

| xi | 100 PeV total shift | 300 PeV total shift |
|---:|---:|---:|
| 1.0 | 3.46% | 7.32% |
| 0.5 | 0.871% | 1.48% |
| 0.25 | 0.211% | 0.314% |

At xi=1, roughly 72--77% of the signal at 100--300 PeV is generated from
`shat/Lambda^2 > 0.5`, so the only >5% point is explicitly UV-endpoint
dominated.

### `AUDIT_UHE_OPERATORS_v0.29_2026-10-30.md`
Preserves the independent audit of this branch.  Important points:

- the v0.27 PDF screen completed;
- the v0.29 collider workflow failed for an environment/CONTUR initialization
  reason, not because physics excluded the model;
- the collider result therefore remained incomplete;
- the 7.32% number is not to be described as a robust allowed prediction;
- boundary-aware quadrature, crossed coupled-channel unitarity and collider
  validation remained outstanding.

## Numerical code/results

- `scripts/d8_gluon_uhe_v0_27.py`
- `results/v0.27/d8_gluon_uhe_v0_27.csv`

Keep the result labelled **endpoint/conditional**, not universally excluded.

---

# 4. Sterile d=6 scalar/tensor semileptonic operators

## Main detailed source

### `STERILE_D6_RESOLVED_CERTIFICATION_v0.32.md`
This is the authoritative detailed correction for the gauge-invariant sterile
scalar/tensor contact branch.  It contains:

- the HighPT bounds used for `O_lNqd^S,T[tau23]`;
- exact direct-chiral SMEFT partonic kernels;
- the weak-coupling contact-domain definition;
- CT18NLO and NNPDF31 dynamic-Q folds;
- the flavor/sterile-mass caveat;
- PDF robustness check;
- the explicit factor-16 normalization mistake in the preliminary run.

Final controlled values:

| operator | CT18 100 PeV | CT18 300 PeV |
|---|---:|---:|
| scalar tau23 | 0.00550% | 0.00667% |
| tensor tau23 | 0.540% | 0.678% |

NNPDF31 changes these only slightly.  Even the looser `4pi` contact-domain
screen leaves the tensor below the 5% project gate.

### Withdrawn result preserved
The earlier 8--11% tensor result is explicitly preserved as **wrong**.  The
mistake was a factor-16 normalization error from inserting low-energy GNI
`epsilon_S,T` kernels directly with an SMEFT `C/Lambda^2` coefficient while
also retaining the low-energy `2 sqrt(2) G_F` normalization.

## Numerical code/results

- `scripts/sterile_d6_resolved_uhe_v0_32.py`
- `results/v0.32/sterile_d6_resolved_dynamicQ_v0_32.csv`
- `results/v0.32/sterile_d6_resolved_dynamicQ_v0_32.json`
- `results/v0.32/sterile_d6_resolved_summary_v0_32.md`

Corrected numerical-result commit recorded in that note:
`196f2996413cb51c0edb5701023161758c5c5144`.

---

# 5. Sterile derivative-tensor and sterile gluon/Rayleigh branches

## Screening record

### `STERILE_UHE_SCREEN_v0.30.md`
Preserves the operator-by-operator sterile-transition survey, including:

- derivative-vector transition;
- derivative-tensor transition;
- scalar/CP-odd gluon Rayleigh transition;
- d=6 sterile semileptonic classes;
- electroweak dipole portal;
- resolved derivative mediator;
- active d=8 comparators.

The useful stored strong/endpoint numbers are:

- derivative-tensor sterile: 5.52% at 100 PeV, 7.32% at 300 PeV at its
  perturbative endpoint;
- sterile gluon/Rayleigh: 3.11% at 100 PeV, 6.42% at 300 PeV near its
  endpoint;
- electroweak dipole portal: rejected for the hard-UHE objective because the
  photon/Z propagator drives the fold toward low Q rather than sustained hard
  momentum growth.

These should be described as **strong/composite or endpoint benchmarks**, not
as generic weakly coupled EFT successes.

The later consolidated interpretation is in `UHE_INTERACTION_PORTFOLIO_v0.33.md`
and the historical v0.40 master, with v0.40 superseded for PKR normalization by
v0.86/v0.87.

---

# 6. Universal B-3L_tau vector branch

The consolidated portfolio files preserve the reasoning:

- anomaly-free and renormalizable;
- allowed `g_X/M_X` small enough that the neutral-current correction is
  sub-percent;
- up/down SM NC interference partly cancels for universal quark charge.

Primary sources in this repository:

- `UHE_INTERACTION_PORTFOLIO_v0.33.md`
- `UHE_INTERACTION_PORTFOLIO_MASTER_v0.40.md`

This branch was screened out on size before an expensive full transport fold
was justified.

---

# 7. PKR Eq. (4.4) derivative and Eq. (4.3) transition-dipole branches

This branch has the most complicated supersession chain and must be read in
order.

## Historical development (preserve, but not current normalization)

- `PKR_EQ4P3_WIDTH_CERTIFICATION_v0.82.md`
- `PKR_5GEV_ICECUBE_RECONSTRUCTION_v0.85.md`
- historical result directories under `results/v0.83_*`

These are retained because they document how the project reached the later
normalization audit.  They are **not** current DIS normalizations.

## Current normalization authority

### `PKR_NORMALIZATION_RESOLUTION_v0.86.md`
This is the detailed analytic/numerical normalization resolution.  It records:

- exact decay and DIS traces for the explicitly specified massless,
  distinct-state chiral transitions;
- the corrected absolute coefficients;
- the correct Eq. 4.3 / Eq. 4.4 relative rate of four in the matched
  convention;
- legacy absolute DIS rates too large by factors 8 and 16 respectively;
- separation between a project normalization error and any claim about what
  the original published PKR implementation actually did.

Supporting files:

- `docs/PKR_NORMALIZATION_PLAN_v0.86.md`
- `STATUS_AUDIT_v0.86.md`
- `manuscript/v0.86/normalization.tex`
- supersession marker
  `results/v0.83_eq4p3_full_envelope_collider/STATUS_SUPERSEDED_v0.86.md`

## Current collider/UHE disposition

### `PKR_CONSTRAINT_DISPOSITION_v0.87.md`
This is the current experimental disposition for the tested masses.  It
contains:

- exact CMS dijet table entries for 50, 100, 220 GeV;
- invisible branching treatment;
- calibrated mass-specific monojet shapes;
- explicit width/validity scan;
- direct dynamic-Q CT18 evaluation;
- limitations of the public IceCube response;
- the important active<->sterile reconversion transport issue.

Nominal collider-compatible 100-PeV maxima:

| M (GeV) | Eq. 4.4 | Eq. 4.3 |
|---:|---:|---:|
| 50 | 0.3395% | 0.7953% |
| 100 | 0.0903% | 0.1084% |
| 220 | 0.0456% | 0.0581% |

This does **not** close 5 or 40 GeV, all flavors, or all UV assumptions.

Machine-readable experimental inputs/results are under
`results/v0.87_constraints/`, including the preserved CMS HEPData JSON/CSV.

---

# 8. Portfolio files and their status

## `UHE_INTERACTION_PORTFOLIO_v0.33.md`
Useful consolidated snapshot after the sterile d=6 correction.  Contains the
four labels SIZE PASS / CONDITIONAL / STRONG-ENDPOINT / FAILED and a detailed
ranking.

## `UHE_INTERACTION_PORTFOLIO_MASTER_v0.40.md`
Historically authoritative through the then-current PKR certification.  Keep it
for project history and the non-PKR branches.

**Important:** its PKR Eq. 4.3/4.4 statements are superseded by
`PKR_NORMALIZATION_RESOLUTION_v0.86.md` and
`PKR_CONSTRAINT_DISPOSITION_v0.87.md`.

## `AUDIT_UHE_OPERATORS_v0.29_2026-10-30.md`
Independent audit supplement.  It is especially useful for distinguishing
"not certified" from "excluded" and for recording failed workflow execution.

---

# 9. What is preserved versus what is not claimed

Preserved in the repository:

1. operator definitions and gauge/flavour assignments;
2. partonic amplitudes and differential kernels where derived;
3. spin/color averaging conventions;
4. partial-wave projections and unitarity ceilings;
5. positivity signs where used;
6. EFT/cutoff prescriptions and endpoint-support diagnostics;
7. PDF choice and dynamic-Q treatment;
8. numerical cross-section tables;
9. collider-screen assumptions and exact experimental inputs where obtained;
10. transport topology distinctions (active regeneration versus sterile sink,
    and later inverse-transition/reconversion caveat);
11. failed executions and withdrawn numerical claims;
12. scripts and machine-readable outputs for the major executed branches.

Not claimed:

- a mathematical no-go theorem covering every conceivable nonresonant model;
- a complete all-constraints exclusion of the endpoint gluon branch;
- an exact IceCube likelihood for arbitrary changed y distributions;
- closure of untested PKR masses such as 5/40 GeV;
- that a historical superseded normalization should be reused.

---

# 10. Reading order for future reconstruction

For a future audit or manuscript reconstruction, use:

1. `PROGRESS.log` — chronology and why the strategy changed;
2. `D8_CERTIFICATION_STATUS_v0.26.md` +
   `D8_CONTACT_BRANCH_CLOSED_v0.27.md` — active quark d=8;
3. `GLUON_D8_SCREEN_v0.27.md` + `GLUON_D8_CERTIFICATION_v0.28.md` +
   `AUDIT_UHE_OPERATORS_v0.29_2026-10-30.md` — active gluon d=8;
4. `STERILE_UHE_SCREEN_v0.30.md` — broad sterile-transition landscape;
5. `STERILE_D6_RESOLVED_CERTIFICATION_v0.32.md` — corrected d=6 sterile and
   resolved-derivative screen;
6. `UHE_INTERACTION_PORTFOLIO_v0.33.md` — consolidated pre-PKR ranking;
7. `PKR_NORMALIZATION_RESOLUTION_v0.86.md` — current PKR normalization;
8. `PKR_CONSTRAINT_DISPOSITION_v0.87.md` — current tested-mass collider/UHE
   disposition;
9. `NONRESONANT_UHE_OPERATOR_HISTORY_v0.89.md` — concise final synthesis;
10. this file — provenance map tying the synthesis back to the detailed record.

This reading order is intended to prevent any future summary from erasing the
physics reasoning or reviving a superseded normalization.
