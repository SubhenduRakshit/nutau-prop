# nutau-prop — canonical UHE BSM / QCD / propagation exploration record v0.90

Date: 2026-10-31
Status: CANONICAL RESEARCH-HISTORY INDEX AND DECISION RECORD
Purpose: preserve in one place what was explored, why, what worked, what failed, what remains conditional, and where the full derivations/results live.

This document is intentionally broader than an operator summary. It records:

1. nonresonant momentum-growing operators;
2. resolved nonresonant mediators;
3. resonant models;
4. charged-scalar / charged-Higgs mechanisms;
5. CGC/small-x QCD as a control and possible competing explanation;
6. full nu_tau–tau–nu_tau propagation, NC feed-down and tau-decay secondary flavours;
7. detector observability (IceCube, ARCA, GRAND, POEMMA, Trinity, Gen2-like concepts);
8. novelty/literature checks;
9. withdrawn/incorrect intermediate results and why they were wrong;
10. the current decision tree.

The short summaries in this file do NOT replace the detailed derivation files, code, CSV/JSON outputs or failed checkpoints. Source files are listed explicitly under each branch.

---

# 1. Physics target and decision rules

Primary physics target:

- identify BSM interactions capable of producing a measurable distortion of UHE neutrino propagation through Earth;
- emphasis on nu_tau because CC production of tau, tau energy loss/decay and regeneration generate nontrivial transport;
- distinguish changes in the inclusive neutrino-nucleon cross section from detector-level observability.

Historical triage gate:

- nominal >=5% modification of the inclusive UHE nu-N cross section around 100–300 PeV.

Updated observability criterion:

- >=5% cross-section change is only a screening criterion, not a discovery criterion;
- Earth attenuation can amplify a small cross-section difference when optical depth is appreciable;
- very opaque bins can show enormous relative flux ratios while carrying negligible absolute event rate;
- the physically relevant quantities are energy-angle-flavour-topology distributions folded with detector acceptance.

Status labels used here:

- CLOSED: controlled calculation shows the branch cannot meet the project goal under the tested assumptions, and additional constraints can only weaken it;
- STRONG/ENDPOINT: several-percent effect appears only near strong coupling or the edge of the EFT-validity region;
- CONDITIONAL: size is adequate but one or more binding model/constraint gates remain open;
- SIZE PASS: calculated benchmark exceeds the historical size gate, but this does not imply an allowed or observable model;
- NOT CERTIFIED: a calculation exists, but a decisive normalization, collider, flavour, detector or UV-completion check remains incomplete.

---

# 2. Transport physics that is already implemented

The propagation is not a simple exp(-tau) attenuation model.

For e and mu flavours the transport includes:

- CC removal;
- NC energy redistribution/feed-down;
- separate neutrino and antineutrino evolution.

For tau flavour the coupled transport includes:

- nu_tau N -> tau X (CC);
- tau energy loss;
- tau decay;
- nu_tau regeneration;
- antinu_tau regeneration;
- secondary nu_mu / antinu_mu from leptonic tau decay;
- in the all-flavour CGC runner, the corresponding e secondary is included identically in the massless daughter approximation with the same leptonic branching fraction;
- NC energy migration;
- PREM chords.

Important consequence:

Phi_alpha(out) = surviving primary + NC feed-down + tau-decay regeneration/secondaries.

Thus a change in the high-energy interaction rate can alter a lower-energy bin even when the local cross section at that lower energy is unchanged.

Primary implementation references:

- `src/nutau_prop/chirality_transport_expm.py`
- `src/nutau_prop/kernels.py`
- `src/nutau_prop/tau_evolution.py`
- `PROGRESS.log`

Historical validation lineage:

- Rakshit–Reya reproduction;
- Dutta–Reno–Sarcevic secondary-flavour checks;
- corrected inclusive tau-decay kernel;
- exact sparse-exponential interaction evolution.

---

# 3. Nonresonant higher-dimensional operator program

## 3.1 Gauge-invariant active-to-sterile d=6 scalar

Representative structure: nuSMEFT semileptonic scalar `O_lNqd^S`.

Current-bound exact dynamic-Q result:

- 100 PeV: about 0.0055–0.0057% of total SM;
- 300 PeV: about 0.0067–0.0069%.

Reason for failure:

- energy growth is too weak;
- current HighPT bounds force the coefficient far below what would be needed;
- resolving the scalar mediator cannot increase the contact-envelope rate in ordinary spacelike exchange and generally adds constraints.

Verdict: CLOSED for this project.

Detailed source:

- `STERILE_D6_RESOLVED_CERTIFICATION_v0.32.md`
- `results/v0.32/sterile_d6_resolved_dynamicQ_v0_32.csv`

## 3.2 Gauge-invariant active-to-sterile d=6 tensor

Corrected current-bound result:

- CT18: ~0.540% at 100 PeV, ~0.678% at 300 PeV;
- NNPDF31: ~0.568% at 100 PeV, ~0.706% at 300 PeV.

Even a much looser 4pi contact domain gives only about 1.0% / 1.8% in the CT18 benchmark.

Critical withdrawn result:

- an earlier 8–11% tensor effect was WRONG;
- cause: factor-16 normalization error from inserting low-energy epsilon-normalized kernels directly with SMEFT Wilson coefficients;
- the withdrawn result is preserved in history and must not be reused.

Verdict: CLOSED for the several-percent target.

Detailed source:

- `STERILE_D6_RESOLVED_CERTIFICATION_v0.32.md`

## 3.3 Active two-derivative d=8 quark operators

Operators:

- `l2u2D22[3322]` ~ (L3bar gamma_mu i<->D_nu L3)(cRbar gamma^mu i<->D^nu cR)
- `l2d2D22[3322]` ~ (L3bar gamma_mu i<->D_nu L3)(sRbar gamma^mu i<->D^nu sR)

Key analytic results:

- amplitude grows as shat^2 / Lambda^4;
- binding partial wave is antiquark J=0;
- strict unitarity through shat=Lambda^2 gives |C| <= 8pi/3;
- positivity fixes C <= 0.

Exact dynamic-Q CT18 result at C=-8pi/3, Lambda=4 TeV, most permissive xi=1:

- cR: +0.1196% total at 100 PeV; +0.2277% at 300 PeV;
- sR: +0.1003% total at 100 PeV; +0.1834% at 300 PeV.

Tighter EFT cuts reduce the effect further.

Reason for failure:

- strict unitarity ceiling is too restrictive;
- positivity fixes the interference sign;
- exact PDF fold leaves total effect far below several percent;
- collider limits cannot increase the allowed coefficient.

Verdict: CLOSED for the tested active L3-cR/sR derivative directions. This is NOT a theorem against every d=8 operator.

Detailed sources:

- `D8_CERTIFICATION_STATUS_v0.26.md`
- `D8_CONTACT_BRANCH_CLOSED_v0.27.md`
- `results/v0.25/d8_uhe_dynamicQ_v0_25.csv`

## 3.4 Active gluon d=8 twist-2 / energy-momentum-like operator

Operator:

Q_G^33 = i(L3bar gamma^mu <->D^nu L3) G^A_{mu rho} G^{A rho}_{ nu}.

Key analytic results:

- genuine SMEFT d=8 basis element;
- no explicit g_s in the canonical two-external-gluon contact vertex;
- helicity amplitudes give J=1/2 and J=3/2 partial waves;
- J=1/2 binds;
- strict unitarity gives |C_G| <= 8pi at shat=Lambda^2;
- positivity gives C_G <= 0.

Exact CT18 dynamic-Q results for C_G=-8pi, Lambda=4 TeV:

xi=1:

- 100 PeV: +3.46% total;
- 300 PeV: +7.32% total.

xi=1/2:

- 100 PeV: +0.871% total;
- 300 PeV: +1.48% total.

xi=1/4:

- 300 PeV: +0.314% total.

Endpoint diagnostic:

- at xi=1, ~72–77% of the stored 100–300 PeV BSM contribution comes from shat/Lambda^2 > 1/2;
- therefore the only >5% point is strongly UV-endpoint dominated.

Collider status:

- exact operator-specific collider certification was not completed;
- a Contur workflow failure was environmental and must not be interpreted as a physics exclusion.

Verdict: STRONG/ENDPOINT / NOT FULLY CERTIFIED; not a robust weakly-coupled several-percent solution.

Detailed sources:

- `GLUON_D8_SCREEN_v0.27.md`
- `GLUON_D8_CERTIFICATION_v0.28.md`
- `AUDIT_UHE_OPERATORS_v0.29_2026-10-30.md`
- `results/v0.27/d8_gluon_uhe_v0_27.csv`

## 3.5 Sterile derivative-tensor transition

Historical endpoint result:

- ~5.52% at 100 PeV;
- ~6.78% at 220 PeV;
- ~7.32% at 300 PeV.

But the effect appears only near C Lambda^3 ~ 4pi / the perturbative endpoint.

Weakly-coupled UV completions tend to suppress such operators, often loop-wise.

Verdict: STRONG/ENDPOINT benchmark only.

Detailed sources:

- `STERILE_UHE_SCREEN_v0.30.md`
- `STERILE_UHE_CERTIFICATION_v0.31.md`
- `UHE_INTERACTION_PORTFOLIO_v0.33.md`

## 3.6 Sterile gluon / Rayleigh transition

Representative parent:

(C_G/Lambda^4)(Nbar Htilde^dagger L_tau) G G

or CP-odd G Gtilde.

Endpoint result:

- ~3.11% at 100 PeV;
- ~5.38% at 220 PeV;
- ~6.42% at 300 PeV;
- hard y distribution, mean y around 0.8.

A common-scale 4-TeV benchmark gave much smaller ~0.71% / 1.42% at 100/300 PeV.

Reason for demotion:

- >5% requires endpoint/strong normalization;
- generic weakly-coupled UV matching is not naturally this large.

Verdict: STRONG/ENDPOINT benchmark only.

Detailed source:

- `STERILE_UHE_SCREEN_v0.30.md`

## 3.7 Sterile electroweak dipole portal

Gauge-invariant operators:

- O_NB = (Lbar sigma^{mu nu} N) Htilde B_mu nu
- O_NW = (Lbar sigma^{mu nu} tau^I N) Htilde W^I_mu nu.

Why it fails as a hard-UHE mechanism:

- although the vertex contains momentum, photon/Z t-channel exchange drives the integral toward low Q^2;
- there is no sustained hard-DIS polynomial growth comparable to the derivative d=8 classes;
- transition-moment constraints further suppress the coefficient.

Verdict: CLOSED for the inclusive >5% hard-UHE target; may remain interesting for low-recoil or decay signatures.

Detailed source:

- `STERILE_UHE_SCREEN_v0.30.md`

## 3.8 Universal B-3L_tau vector

Theoretically clean/anomaly-free as a control model, but:

- allowed coupling-to-mass ratio is small;
- universal quark charge leads to substantial cancellation between SM up/down NC interference contributions.

Expected inclusive effect: sub-percent in the allowed region.

Verdict: CLOSED on size for present goal.

Detailed sources:

- `UHE_INTERACTION_PORTFOLIO_v0.33.md`
- `UHE_INTERACTION_PORTFOLIO_MASTER_v0.40.md`

---

# 4. PKR Eq. (4.3)/(4.4) resolved derivative/dipole Z' program

This branch underwent a major normalization correction and must be read only through the latest checkpoints.

## 4.1 Exact normalization correction

For the explicit massless distinct-state chiral transition used in the project:

Eq. 4.4 derivative vertex:

|M|^2 = 4 kappa^2 shat^3 y(1-y) / D^2

dsigma/dy = kappa^2 shat^2 y(1-y)/(4pi D^2)

with D=M_Z'^2+Q^2.

The inherited repository coefficient was too large by factor 16 in cross section.

Eq. 4.3 transition dipole:

- matched scattering result is four times Eq. 4.4 at equal coefficient under the specified external-state convention;
- inherited Eq. 4.3 absolute DIS rate was too large by factor 8.

Do not translate this automatically into a claim that the original published PKR paper is wrong without reconstructing its exact field/conjugation/channel conventions.

Detailed source:

- `PKR_NORMALIZATION_RESOLUTION_v0.86.md`
- `STATUS_AUDIT_v0.86.md`
- `docs/PKR_NORMALIZATION_PLAN_v0.86.md`

## 4.2 Corrected collider-compatible effects

Latest tested exact-dijet + calibrated monojet implementation at 100 PeV:

Eq. 4.4 maxima:

- 50 GeV: 0.3395%
- 100 GeV: 0.0903%
- 220 GeV: 0.0456%

Eq. 4.3 maxima:

- 50 GeV: 0.7953%
- 100 GeV: 0.1084%
- 220 GeV: 0.0581%

Relaxing the monojet screen to r<=1.2 gives at 50 GeV only ~0.406% (4.4) and ~0.957% (4.3).

Direct 5% target test:

- Eq. 4.4 cannot reach 5% at 50/100/220 GeV after exact dijet + width alone;
- Eq. 4.3 can reach raw 5% at 50/100 GeV before monojet, but the calibrated monojet metric fails badly;
- 220 GeV Eq. 4.3 reaches only ~2.48% even without monojet.

Scope limitation:

- this closes the tested 50/100/220 GeV collider implementation;
- it does NOT close 5 GeV, 40 GeV, all other masses/flavour structures or all UV assumptions.

Detector/transport caveat:

- nu_tau -> N transitions are not necessarily irreversible sinks if the Hermitian interaction allows N q -> nu_tau q;
- exact IceCube reweighting requires y-dependent detector response not available in the public HESE MC schema used here.

Verdict: tested collider-compatible points are sub-percent; broader low-mass space is not universally excluded.

Detailed source:

- `PKR_CONSTRAINT_DISPOSITION_v0.87.md`
- `PKR_5GEV_ICECUBE_RECONSTRUCTION_v0.85.md`

---

# 5. Resolved charged scalar / charged-Higgs loophole — REOPENED

This is an important branch and must not be lost among the EFT failures.

## 5.1 Mechanism

A light charged scalar can mediate

nu_tau + s -> tau^- + c

through t-channel exchange.

This is qualitatively different from our momentum-growing EFT operators:

- it is resolved rather than contact-like;
- it can be renormalizable;
- mediator mass can be weak-scale (~90 GeV);
- it directly produces a tau, so it modifies the CC conversion probability relevant to Earth-skimming tau telescopes;
- s/c sea PDFs are suppressed at the LHC relative to light valence production, but become much less disadvantaged at UHE small x.

## 5.2 Our stored size results

Broken-phase parameter A = |g_tau|^2 (|g_L|^2+|g_R|^2).

For M_H=90 GeV, A=1:

- 100 PeV: ~1658 pb, ~24% of the historical SM reference;
- 300 PeV: ~2394 pb, ~23%.

For the tuned example g_tau=1, g_L=g_R=0.4, A=0.32:

- 100 PeV: ~531 pb, ~7.6%;
- 300 PeV: ~766 pb, ~7.4%.

Therefore our own analysis DID find a substantial nonresonant resolved effect.

## 5.3 Why it was demoted

The issue was never raw UHE size. The issue was model viability.

Potential constraints/correlations include:

- D_s -> tau nu;
- tau and K/meson decays;
- CKM-related charged currents;
- gauge-related neutral scalar partners;
- electroweak precision constraints;
- collider tau+MET / dijet / charged-Higgs searches;
- model-dependent neutrino-mass/flavour structure.

A simple cancellation such as g_L=g_R in one pseudoscalar observable is not by itself a gauge-complete certificate.

## 5.4 Literature connection: Huang et al.

Huang, Jana, Lindner & Rodejohann (JCAP 2022; arXiv:2112.09476) studied a flavour-structured charged-scalar scenario with cs and tau-nu_tau couplings and propagated nu_tau/tau through Earth for GRAND, POEMMA and Trinity.

Their basic reason for a large effect is the same loophole identified above:

- weak-scale resolved charged scalar;
- order-one-ish flavour-selected couplings;
- cs sea quarks become important at UHE;
- direct additional tau-producing CC channel;
- collider production can be less severe than for first-generation couplings.

This resolves the apparent conflict with our EFT results: the charged scalar is not a higher-dimensional momentum-growing contact interaction. It is a different mechanism.

## 5.5 Current action

REOPEN this branch for a 2026 certification.

Required gate:

1. choose an explicit gauge-complete scalar model/texture;
2. include charged and neutral scalar partners;
3. apply current flavour/meson/tau constraints;
4. apply current charged-Higgs/tau+MET/cs searches according to the actual production couplings;
5. identify the allowed envelope in (M_H, g_tau, g_cs);
6. only then propagate through Earth and fold with emerging-tau detector acceptance.

Verdict: SIZE PASS / MODEL VIABILITY CONDITIONAL / high-priority revisit.

Detailed sources:

- `STERILE_UHE_SCREEN_v0.30.md` (F3)
- `PKR_ATZ_UV_REASSESSMENT_v0.29_2026-10-30.md`
- `CHAT_UHE_CROSS_SECTIONS_v0.29_c1_2026-10-30.md`
- `UHE_INTERACTION_PORTFOLIO_v0.33.md`

---

# 6. Resonant models

## 6.1 Strange-aligned R2-tilde leptoquark

Stored benchmark:

- M ~1.5 TeV;
- lambda ~2;
- ~7.2% at 100 PeV;
- ~9.7% at 300 PeV.

Advantages:

- renormalizable, gauge-complete interaction;
- on-shell resonance can evade the contact-EFT endpoint problem;
- strange alignment reduces some first-generation constraints.

Open constraints:

- finite width;
- both weak-doublet components;
- pair/single production;
- flavour and tau-decay constraints;
- ditau and jet+MET tails;
- comparison with literature that often finds collider reach stronger than neutrino reach for related R2/S1 textures.

Verdict: SIZE PASS / CONDITIONAL; best conventional resonance candidate.

Detailed sources:

- `CHAT_LEPTOQUARK_RESONANCE_v0.29_c2_2026-10-30.md`
- `UHE_INTERACTION_PORTFOLIO_v0.33.md`

## 6.2 Tau-neutrinogluon / leptogluon resonance

Stored benchmark:

- M ~3 TeV;
- Lambda ~9 TeV;
- a ~4;
- ~6.1% at 100 PeV;
- ~13% at 300 PeV;
- Gamma/M ~0.035 in the stored benchmark.

Advantage:

- on-shell coloured resonance + rapidly growing gluon luminosity.

Problems:

- transition chromomagnetic coefficient is large for a loop-generated interaction;
- coloured charged partner limits can be severe;
- matching/UV origin remains unresolved.

Verdict: SIZE PASS / CONDITIONAL / strong matching.

Detailed sources:

- `UHE_PORTFOLIO_FOLLOWUP_v0.34.md`
- `UHE_INTERACTION_PORTFOLIO_v0.33.md`

---

# 7. CGC / small-x QCD control calculation

CGC is not BSM, but it is essential because it can mimic or exceed BSM propagation effects at UHE.

## 7.1 v0.88b surrogate

Implementation:

- literature-calibrated rcBK+DLA-2015 total-cross-section suppression surrogate;
- both SM CC and NC totals multiplied by R_CGC(E);
- conditional SM y distributions retained;
- full PREM propagation;
- all flavours;
- tau regeneration and secondary flavours included.

Anchor ratios:

- 1e8 GeV: 1.0;
- 1e9 GeV: 2/3;
- 1e12 GeV: 0.4;
- 1e14 GeV: 1/4.5;
- log-log interpolation.

This is deliberately a control model, not a precision implementation of a CGC differential y-kernel.

## 7.2 Completed v0.88b propagation result

For incident 1:1:1, nu=antinu, gamma=2.5, upgoing-angle average:

100 PeV:

- R_CGC=1;
- e/mu flux change ~ -0.0075%;
- tau ~ -0.018%.

300 PeV:

- R_CGC ~0.824;
- e/mu flux change ~ +15.75%;
- tau ~ +15.72%.

1 EeV:

- R_CGC=2/3;
- e/mu flux change ~ +34.19%;
- tau ~ +34.18%.

Interpretation:

- a reduced cross section makes Earth more transparent;
- energy migration means even a bin with local R_CGC=1 can change through feed-down from higher energy;
- enormous vertical ratios at 300 PeV–EeV (up to millions or much larger percent) occur because the SM denominator is essentially zero in extremely opaque bins and are not large event yields.

Therefore future plots must show absolute survival/event support and mask/flag opaque bins.

Result commit from user run:

- `92b1f5d` — results: v0.88b CGC all-flavour Earth propagation.

Detailed sources:

- `CGC_EARTH_PROPAGATION_v0.88.md`
- `src/nutau_prop/cgc.py`
- `scripts/cgc_earth_allflavor_v0_88.py`
- `scripts/cgc_earth_allflavor_v0_88b.py`
- `results/v0.88b_cgc_allflavor/summary.csv`
- `results/v0.88b_cgc_allflavor/flux_grid.csv`
- `results/v0.88b_cgc_allflavor/upgoing_angle_average.csv`

## 7.3 Novelty status

Not novel by itself:

- QCD-model dependence of Earth attenuation/angular distributions has prior literature;
- UHE nu_tau propagation with alternative/saturation-suppressed cross sections and tau regeneration has prior literature;
- detector event-rate effects of small-x QCD have prior literature.

Potentially more interesting directions:

- full modern nonlinear nuclear-BK/CGC with PREM-composition dependence R(E,A(r));
- compare mantle-only versus core-crossing trajectories;
- combine all-flavour regeneration with detector topology;
- profile CGC/small-x uncertainty against explicit BSM models in the same likelihood.

Do not claim novelty until an exhaustive literature check confirms the specific combination.

---

# 8. What GRAND / POEMMA / Trinity actually measure

For Earth-skimming UHE tau experiments the key observable is not simply outgoing nu_tau flux.

The relevant chain is:

nu_tau -> propagation through Earth -> final CC conversion to tau -> tau exits Earth -> tau decay air shower.

A useful schematic quantity is:

P_exit(E_nu,E_tau,beta) = probability that an incident nu_tau produces an emerging tau with given energy and angle.

BSM can change:

- neutrino penetration;
- NC feed-down;
- tau production cross section;
- tau regeneration;
- energy/angular distribution of the emerging tau.

Therefore:

- a smaller total cross section can increase penetration but reduce the probability of the final CC conversion;
- an NC-like BSM interaction may distort flux without directly increasing the emerging-tau yield;
- a charged scalar or leptoquark CC channel can directly increase tau production and is much more favourable for GRAND/POEMMA/Trinity.

This distinction explains why the tiny controlled nonresonant EFT effects are unlikely to be visible even though other BSM models can be observable in those experiments.

Literature benchmark:

- Huang et al. explicitly found future sensitivity in charged/neutral scalar and leptoquark scenarios after solving coupled tau-neutrino propagation and studying energy/angular distributions.

Our current controlled nonresonant effects:

- d=6 scalar/tensor, active d=8 quark, corrected PKR tested points: generally sub-percent -> almost certainly too small;
- endpoint d=8 gluon / derivative-tensor / Rayleigh: potentially amplified but not robust controlled benchmarks;
- charged scalar, LQ resonance, leptogluon resonance: large enough to justify dedicated emerging-tau calculations if current constraints leave an allowed point;
- CGC: propagation effect can be tens of percent and is relevant as a QCD systematic/control.

---

# 9. Novelty lessons from literature checks

Statements that should NOT be made:

- “first BSM study of nu_tau-tau-nu_tau propagation” — false;
- “first demonstration that BSM can be observable at GRAND/POEMMA/Trinity” — false;
- “first CGC Earth-propagation calculation” — false;
- “all nonresonant BSM is excluded” — too strong;
- “5% cross-section shift is automatically observable” — false.

More defensible possible novelty targets:

1. all-flavour consequences of tau regeneration under a specific BSM interaction, including secondary nu_mu and nu_e;
2. systematic discrimination between small-x QCD/CGC and BSM using joint (E, cos theta, flavour, topology) information;
3. PREM-composition-dependent nuclear saturation with full flavour transport;
4. a 2026 recertification of a flavour-structured cs–tau-nu charged scalar, if current data reopen/close a parameter region differently from the 2021/2022 literature;
5. a distinctive surviving resonance topology not already covered by existing UHE-tau detector studies.

Any future novelty claim must be tied to the exact observable/model combination, not a broad phrase.

---

# 10. Current research decision tree

## 10.1 Do NOT open another generic high-dimension operator scan first

Reason:

- extra derivatives increase energy growth but also drive the signal toward the EFT cutoff;
- d=8 quark and gluon examples already show the central tension;
- d=10/d=12 alone are not a generic rescue because each additional derivative pair brings powers of shat/Lambda^2, which remain <=1 in the controlled EFT domain.

A new higher-dimensional operator is worthwhile only if it introduces a genuinely new enhancement mechanism:

- symmetry forbids lower-dimensional operators;
- unusual flavour alignment;
- gluon enhancement with a controlled UV origin;
- distinct resolved threshold/resonance;
- protected coefficient/matching structure.

## 10.2 Highest-priority physics tasks

A. Reopen and fully certify the 90-GeV-ish flavour-structured charged-scalar branch.

B. Fully certify the strange-aligned R2-tilde leptoquark benchmark.

C. Certify the leptogluon resonance only if a plausible matching origin for the large dipole coefficient survives coloured-partner bounds.

D. Use CGC as a mandatory QCD control/systematic, not as a new-paper result by itself.

E. For any surviving BSM point, calculate the actual emerging-tau distribution and detector-level observability for GRAND/POEMMA/Trinity rather than inferring observability from Delta sigma/sigma.

F. Preserve all failures and withdrawn results.

---

# 11. Canonical detailed-source map

Nonresonant synthesis/provenance:

- `NONRESONANT_UHE_OPERATOR_HISTORY_v0.89.md`
- `NONRESONANT_UHE_OPERATOR_PROVENANCE_v0.89.md`

Active d=8 quark:

- `D8_CERTIFICATION_STATUS_v0.26.md`
- `D8_CONTACT_BRANCH_CLOSED_v0.27.md`

Active d=8 gluon:

- `GLUON_D8_SCREEN_v0.27.md`
- `GLUON_D8_CERTIFICATION_v0.28.md`
- `AUDIT_UHE_OPERATORS_v0.29_2026-10-30.md`

Sterile / operator portfolio:

- `STERILE_UHE_SCREEN_v0.30.md`
- `STERILE_UHE_CERTIFICATION_v0.31.md`
- `STERILE_D6_RESOLVED_CERTIFICATION_v0.32.md`
- `UHE_INTERACTION_PORTFOLIO_v0.33.md`
- `UHE_INTERACTION_PORTFOLIO_MASTER_v0.40.md` (historical portfolio; PKR conclusions superseded by v0.86/v0.87 where inconsistent)

Broader model-building assessment:

- `PKR_ATZ_UV_REASSESSMENT_v0.29_2026-10-30.md`
- `CHAT_UHE_CROSS_SECTIONS_v0.29_c1_2026-10-30.md`
- `CHAT_LEPTOQUARK_RESONANCE_v0.29_c2_2026-10-30.md`
- `UHE_PORTFOLIO_FOLLOWUP_v0.34.md`

PKR normalization/constraints:

- `PKR_NORMALIZATION_RESOLUTION_v0.86.md`
- `STATUS_AUDIT_v0.86.md`
- `PKR_CONSTRAINT_DISPOSITION_v0.87.md`
- `PKR_5GEV_ICECUBE_RECONSTRUCTION_v0.85.md` (historical/reconstruction only; read with v0.86 warning)

CGC:

- `CGC_EARTH_PROPAGATION_v0.88.md`
- `results/v0.88b_cgc_allflavor/*`

Transport history:

- `PROGRESS.log`
- `VALIDATION.md`
- `src/nutau_prop/chirality_transport_expm.py`
- `src/nutau_prop/kernels.py`

---

# 12. Present concise conclusion

The broad program has not found a robust several-percent weakly-coupled momentum-growing contact interaction after proper validity and constraint treatment.

However, that does NOT mean all nonresonant new physics is dead.

The key surviving lesson is that large UHE effects are more naturally obtained through resolved weak-scale or resonant physics than by simply increasing EFT dimension.

The most important nonresonant resolved loophole to revisit is the flavour-structured charged scalar producing nu_tau s -> tau c. Our own calculation already showed a 7–24% size effect for representative 90-GeV benchmarks; the unresolved question is whether a gauge-complete 2026 parameter point survives flavour and collider constraints.

Resonant leptoquark/leptogluon candidates remain conditional alternatives.

CGC/small-x QCD can itself generate O(10–30%) propagated-flux differences in the 300-PeV–EeV regime in the current surrogate and therefore must be profiled as a competing QCD effect/systematic in any future BSM claim.

Future observability claims must be made at the detector level, especially through the emerging-tau energy/angle distribution for GRAND/POEMMA/Trinity-like experiments.
