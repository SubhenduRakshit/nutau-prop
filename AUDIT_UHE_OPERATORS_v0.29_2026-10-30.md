# UHE operator audit — v0.29, 30 October 2026

Audit base: a75a77aabf09c3f9267b1d9e53847a1a1bfc0e98 (GitHub main).
Scope: analytical review, committed documentation/results, static inspection of the gluon fold and collider generator, and existing CI logs. No simulations, reruns, code edits or local-working-directory inspection. This is an audit supplement, not a new physics release. Uncommitted Mac Studio work is outside this evidence.

## Finding

No inspected benchmark is yet demonstrated to produce a several-percent UHE cross-section increase after all relevant constraints. This is an absence of certification, not a universal exclusion.

The active gluon operator is the leading *existing* candidate:
\[
Q_G^{33}=i(\bar L_3\gamma^\mu\overleftrightarrow D^\nu L_3)
G^A_{\mu\rho}G^{A\rho}{}_\nu .
\]
It is gauge invariant, dimension eight, and preserves the active daughter. Its charged-tau partner is compulsory. At tree level its positive elastic gluon contribution has no SM elastic neutrino-gluon amplitude with which to interfere. This statement does not remove higher-order matching/mixing corrections.

## Latest execution status

- v0.27 gluon workflow 40919163210 completed successfully and committed the PDF screen.
- v0.29 collider workflow [40925881750](https://github.com/SubhenduRakshit/nutau-prop/actions/runs/40925881750) failed in “Install Contur and PDF”; the event-generation/recast and output-commit steps were skipped.
- Job 123631117996 logs show Contur 3.1.4 installed successfully, but its startup used os.getenv("CONTUR_ROOT") as a path while that variable was unset, raising TypeError. This is an environment-initialization failure, not an exclusion or failed physics model.
- There are no results/v0.29 files in the inspected tree. FeynRules/UFO sources and a MadGraph card are implementations awaiting validation, not evidence of completed independent normalization checks.
- No workflow was started by this audit.

## Stored numerical evidence (not recomputed)

Values below are percentages of the SM total CC+NC cross section for the affected tau flavour, using the repository's equal neutrino/antineutrino average. They are NOT all-flavour detector event enhancements.

| Operator / benchmark | 100 PeV | 300 PeV | Status |
|---|---:|---:|---|
| active c_R, C=-8pi/3, Lambda=4 TeV, xi=1 | +0.1196% | +0.2277% | Small at this benchmark |
| active s_R, same benchmark | +0.1003% | +0.1834% | Small at this benchmark |
| active gluon, C=-8pi, Lambda=4 TeV, xi=1 | +3.4579% | +7.3220% | Collider and detector constraints incomplete |
| active gluon, same C/Lambda, xi=1/2 | +0.8707% | +1.4820% | Strong cutoff dependence |
| active gluon, same C/Lambda, xi=1/4 | +0.2114% | +0.3140% | Strong cutoff dependence |

Evidence: D8_CONTACT_BRANCH_CLOSED_v0.27.md; results/v0.27/d8_gluon_uhe_v0_27.csv; GLUON_D8_CERTIFICATION_v0.28.md.

For the last gluon energy at xi=1, approximately 77% of the stored signal comes from shat/Lambda^2>1/2. The quoted 7.32% is consequently an endpoint-sensitive EFT screen, not a robust allowed prediction. The coefficient threshold |C_G| approximately 20.8 for a 5% total-CS increase follows from C_G^2 scaling only at this fixed mass, energy, cutoff and normalization.

## Analytic checks and corrections to the strategy

### 1. The quark-branch closure is benchmark-specific

The stored 4-TeV, single-coefficient point justifies deprioritizing that point. It does not establish that all masses, coefficient combinations or mediator completions fail. In particular, “collider bounds can only lower the coefficient” is correct at fixed mass and fixed spectrum; it does not justify fixing Lambda=4 TeV as a universal optimum.

Write G=C/M^4 and keep the physical heavy scale M separate from |G|^(-1/4). At fixed C, lowering M both changes G and moves the partonic cut. A global closure needs an allowed envelope in (M,G), not just rescaling C at one M. Existing HighPT outputs should be inspected before any extra numerical scan; projected flavour-universal bounds are not substitutes for an actual tau-specific constraint.

### 2. Strong momentum growth weights the upper allowed x region

For the pure dimension-eight term,
\[
d^2\sigma\propto G^2(2m_NE)^3 x^3 f(x,Q^2)F(y)\,dx\,dy,\qquad
x_{\max}=\min(1,\xi M^2/(2m_NE)).
\]
The logarithmic-x weight is x^4 f. Under the illustrative approximation xf proportional to x^{-lambda}, with scale evolution neglected and x_max<1, the cut integral scales as E^{lambda}, not E^3. This explains why raw small-x luminosities and uncut energy growth cannot rank candidates reliably.

### 3. A 5% total-CS gate is a triage convention, not a flux theorem

Even a pure absorber has delta ln F=-X delta sigma. A smaller cross-section change can yield a larger fractional transmitted-flux change at appreciable optical depth, although the absolute event rate can be tiny. Active NC regeneration makes the relationship still less direct.

For the stated gluon kernel K(y)=(1-y)[1+(1-y)^2],
\[
\int_0^1K\,dy=3/4,\qquad \langle y\rangle=13/45.
\]
These exact massless partonic moments are inexpensive checks for future transport. Number conservation requires the active NC loss and daughter-gain terms to cancel after integration over all daughter energies. Cross-section size alone does not prove observable track enhancement.

### 4. Elastic unitarity is necessary, not a full coupled-channel certificate

The displayed gluon helicity amplitudes do reproduce the stated spin-averaged kernel and elastic bounds algebraically. The absence of a color enhancement in nu g -> nu g is correct for that block. It does not establish the eigenvalues of the crossed color-singlet gg <-> L3 L3bar block. Include both tau and neutrino channels, Bose normalization and any interactions required by a mediator model before labeling the ceiling fully certified.

A stricter EFT cut at fixed G is an uncertainty diagnostic. Raising C by xi^{-2} to saturate unitarity at the smaller cut changes the benchmark and cannot be called the same robustness check. Unitarity saturation itself does not guarantee a controlled perturbative or derivative expansion.

### 5. Static numerical-review limitations

scripts/d8_gluon_uhe_v0_27.py uses tree-level DIS coefficient functions with CT18NLO PDFs: this is a LO hard-process screen with NLO PDFs, not an NLO calculation. Its fixed logarithmic x grid crosses a sharp EFT boundary without inserting that boundary as an integration endpoint. The endpoint-support diagnostic uses a different midpoint integration. Before precision certification, check boundary-aware quadrature and convergence of both quantities.

Both inspected PDF wrappers silently convert exceptions to zero. Out-of-range requests and real PDF failures must be distinguished before trusting a certified fold. This audit has not established that an exception actually occurred in the stored run.

The crossed collider angular kernel and C_G^2 scaling are internally consistent with the supplied elastic normalization; independence still requires the proposed separate amplitude check. Tau polarization must be verified through decay, not inferred solely from helicity tags in the LHE file.

## Constraint status

1. **Ditau:** pending. [ATLAS 2503.19836](https://arxiv.org/abs/2503.19836) provides the relevant high-mass measurement. Its published limits for other operator choices cannot be copied onto Q_G. Validate the specific Rivet/Contur configuration, SM prediction and correlations. The current environment error is identified above.
2. **Monojet:** pending for this direction. A gg -> tau tau two-gluon projection is insufficient for a gauge-complete hard-jet calculation. Use the non-Abelian field-strength vertices and consistent radiation matching; assess EFT validity for the actual multi-leg process.
3. **IceCube:** pending as a tau-only NC constraint. [IceCube 2011.03560](https://arxiv.org/pdf/2011.03560) fixes CC/NC and flavour-independent cross-section rescalings. Its intervals cannot directly certify this operator. The stored 10-PeV gluon correction is only 0.0702% of total CS at xi=1, suggesting that this energy range may be less decisive than collider/UHE data, but not proving that it is nonbinding.
4. **EHE/ARCA:** no complete operator-specific detector likelihood demonstrated by the inspected artifacts. Preserve the previous effective-area rescalings as historical screens.
5. **Theory/laboratory:** retain electroweak partners, coefficient running/mixing, lower-dimensional terms or their justified suppression, tau energy-loss estimates, and the model's correlated laboratory constraints. Positivity sign conventions must match the exact operator definition; [Chala–Li, PRD 109, 065015](https://doi.org/10.1103/PhysRevD.109.065015) supplies relevant dimension-eight positivity analysis, not an experimental bound.

## Judicious operator shortlist

**Priority 1: Q_G^{33}.** Best existing numerical candidate, but not yet allowed. Finish its minimal analytic normalization/coupled-channel checks and existing collider gate before expanding the scan. Include monojet and tau-specific neutrino constraints before claiming survival.

**Priority 2: c_R/s_R active two-derivative directions, conditional only.** Do not spend transport time at the failed 4-TeV benchmark. Revisit only if the already available collider information permits a genuinely different physical-mass region that could pass a conservative analytic envelope. Preserve the failed benchmark; narrow the claimed exclusion rather than erase it.

**Contingency: a specified mediator realization.** A resolved spin-two propagator is not an automatic rescue: at matched low-energy coefficient, a simple spacelike exchange supplies M^2/(M^2+Q^2), suppressing the contact amplitude within its validity region. Any gain must come from newly accessible kinematics and a consistent model. Massive-spin-two interactions also require their own cutoff and correlated-channel checks; adding a propagator alone is not a UV completion.

**Do not open a broad speculative operator scan yet.** A neutrino-selecting Higgs-dressed gluon construction can remove the tree-level charged-tau term, but raises the dimension and brings v^2/M^2 suppression plus Higgs/electroweak channels. It is not a demonstrated loophole. ATZ sterile-daughter branches remain comparisons, with reconversion included if appreciable.

## Minimum next test; no run authorized or performed here

First fix the Contur initialization and validate a few analytic/crossed matrix-element points; only then repeat the existing small collider pilot. At the stored 4-TeV benchmark, an upper bound below |C_G| approximately 20.8 rules out the chosen 5% CS target at 300 PeV, xi=1. A weaker collider bound is only a provisional pass: full constraints and cutoff stability remain required.

Maintain separate labels: “CS increase”, “passes selected screens”, “all specified constraints checked”, and “observable detector effect”. None should be promoted to the next label merely because a runner or model file exists.
