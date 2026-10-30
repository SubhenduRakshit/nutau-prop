# nutau-prop v0.32 — exact d=6 sterile + resolved-derivative UHE certification

Date: 2026-10-30
Status: exact dynamic-Q cross-section screen completed; resolved-mediator phenomenological certification remains conditional on a modern tau-specific likelihood/recast.

## Executive result

The two calculations requested after v0.31 are now complete at the hard-cross-section level:

1. exact dynamic-Q CT18NLO and NNPDF31 folds for the least-constrained gauge-invariant d=6 active-to-sterile scalar/tensor semileptonic operators at current HighPT bounds;
2. exact dynamic-Q resolved derivative-mediator scan with the propagator kept explicitly and the historical PKR kappa bound used as an external normalization anchor rather than matching to the old Case-6 contact coefficient.

The main result is negative for ordinary d=6 sterile contact operators and positive-but-not-yet-certified for the resolved derivative mediator.

### d=6 sterile contact verdict

Using the 2024 HighPT 90% CL bounds

- |[C^S_lNqd]_{tau23}|/Lambda^2 <= 0.41 TeV^-2,
- |[C^T_lNqd]_{tau23}|/Lambda^2 <= 0.25 TeV^-2,

and imposing the weak-coupling c=1 domain sqrt(shat)<1/sqrt(G), the exact folds give:

| PDF | operator | 100 PeV | 300 PeV |
|---|---|---:|---:|
| CT18NLO | scalar tau23 | 0.00550% | 0.00667% |
| CT18NLO | tensor tau23 | 0.540% | 0.678% |
| NNPDF31 | scalar tau23 | 0.00573% | 0.00689% |
| NNPDF31 | tensor tau23 | 0.568% | 0.706% |

Thus no current-bound d=6 S/T point reaches the 5% project gate. Even the much looser 4pi contact domain gives only ~1.0% at 100 PeV and ~1.8% at 300 PeV for the tensor benchmark in CT18. The scalar channel is negligible.

**Conclusion: ordinary gauge-invariant d=6 active-to-sterile scalar/tensor contact interactions are demoted for this project.**

### Resolved derivative mediator verdict

For

L superset kappa Z'_mu (bar N i<->partial^mu P_L nu_tau) + g_q Z'_mu bar q gamma^mu q,

with kappa denoting the product coefficient in the repository convention, the exact kernel is

 d2sigma/dxdy = (4/pi) kappa^2 x s^2 y(1-y) Sum[x f(x,Q2)] /(M_Z'^2+Q2)^2.

The scan uses the historical PKR 5-GeV IceCube benchmark bound

 kappa <= 2.5e-4 GeV^-1

as the normalization anchor and imposes the c=1 derivative-vertex domain sqrt(shat)<Lambda_v with Lambda_v=1/kappa=4 TeV.

At kappa=2.5e-4 GeV^-1, CT18 gives:

| M_Z' | Delta sigma/SM, 100 PeV | Delta sigma/SM, 300 PeV | status of kappa=2.5e-4 |
|---:|---:|---:|---|
| 5 GeV | 320% | 309% | historical PKR IceCube-bound anchor |
| 40 GeV | 102% | 107% | shape point; mass-specific current bound not certified |
| 100 GeV | 45.7% | 50.6% | shape point; mass-specific current bound not certified |
| 220 GeV | 19.3% | 22.4% | shape point; mass-specific current bound not certified |
| 500 GeV | 6.26% | 7.60% | shape point; mass-specific current bound not certified |
| 1 TeV | 1.86% | 2.33% | fails size at fixed kappa |
| 2 TeV | 0.384% | 0.495% | fails size at fixed kappa |

The coupling needed to reach exactly 5% follows from sigma proportional to kappa^2. For the same CT18/c=1-vertex setup:

| M_Z' | kappa_5%(100 PeV) [GeV^-1] | kappa_5%(300 PeV) [GeV^-1] |
|---:|---:|---:|
| 5 GeV | 3.12e-5 | 3.18e-5 |
| 40 GeV | 5.55e-5 | 5.40e-5 |
| 100 GeV | 8.27e-5 | 7.86e-5 |
| 220 GeV | 1.27e-4 | 1.18e-4 |
| 500 GeV | 2.23e-4 | 2.03e-4 |
| 1 TeV | 4.10e-4 | 3.66e-4 |
| 2 TeV | 9.02e-4 | 7.94e-4 |

This is the useful output of the resolved scan: a broad 5--500 GeV mass window can in principle clear the 5% gate with kappa between roughly 3e-5 and 2e-4 GeV^-1. However only the 5-GeV point has been normalized here to an explicit published PKR IceCube bound. PKR also found that LHC monojet becomes stronger than IceCube over an intermediate mass interval (roughly 40--220 GeV), so the fixed-kappa rows at those masses must not be labeled allowed.

A search for a later paper directly superseding the exact derivative operator and providing a modern mass-dependent kappa likelihood did not identify one. The published IceCube cross-section fits after 2019 assume different flavor/CC-NC structures and are not a literal bound on this tau-specific sterile transition. A dedicated modern recast is therefore still required before any resolved-mediator point is called allowed.

## 1. Important normalization correction

The first v0.32 run incorrectly inserted the low-energy GNI epsilon_S,T kernels directly using the HighPT SMEFT coefficient G=C/Lambda^2. The low-energy convention contains the overall 2 sqrt(2) G_F normalization. Combining it directly with the SMEFT coefficient overestimates the S/T cross sections by a factor 16.

The corrected direct-chiral SMEFT kernels used in the final run are

 d sigmahat_S/dy = G_S^2 shat y^2 /(64 pi),

 d sigmahat_T/dy = G_T^2 shat (2-y)^2 /(2 pi).

The preliminary 8--11% tensor number from the first run is therefore withdrawn. The corrected result is 0.54--0.71% in the c=1 domain.

The repository script explicitly records this convention to prevent recurrence.

## 2. Operator basis point: no independent d=6 active-to-sterile vector transition

A gauge-invariant d=6 operator with a single active doublet L and sterile singlet N naturally gives scalar/tensor semileptonic transitions such as O_lNqd. The vector structures O_Nq, O_Nd involve sterile bilinears and require active-sterile mixing to turn an incoming active neutrino into N. They are therefore not an independent unsuppressed active->sterile neutral-current vector analogue of O_lNqd at d=6.

Accordingly the exact d=6 transition screen contains S and T. The momentum-growing derivative-vector candidate remains a higher-dimensional/resolved-mediator problem and is treated separately.

## 3. Flavor and sterile-mass caveat

O_lNqd[tau,2,3] contains both the neutral-current component relevant here (schematically nu_tau s -> N b and conjugates) and the charged-current crossed channel used by HighPT. The PDF fold therefore uses a charge-averaged strange/bottom luminosity.

For a very light sterile state, B->K(*)+invisible provides additional constraints, particularly on the tensor operator. A benchmark with m_N above the B-decay kinematic region (for example around 6 GeV) removes those B-decay production constraints while remaining effectively massless at 100--300 PeV. This does not remove the HighPT bound and does not by itself certify cosmology/lifetime/mixing.

Since the d=6 tensor fails the UHE size gate even at the collider bound, no further flavor/cosmology optimization of this branch is warranted for the present project.

## 4. PDF robustness

CT18NLO vs NNPDF31 changes the c=1 tau23 tensor result only from

- 0.540% to 0.568% at 100 PeV,
- 0.678% to 0.706% at 300 PeV.

This is far too small to affect the <5% decision.

## 5. Current project ranking after v0.32

1. **Resolved derivative active->sterile mediator:** survives the size test over a potentially useful mass/coupling window, but needs a current tau-specific IceCube + monojet recast and a UV completion of the derivative vertex.
2. **Derivative-tensor d=7 sterile/contact:** already gives 5.5--7.3% only near c~4pi; retain only as strong/composite benchmark.
3. **Sterile gluon/Rayleigh:** reaches >5% only near its strong-coupling endpoint; retain as strong/composite benchmark.
4. **d=6 sterile tensor:** FAIL, 0.54--0.71% at current HighPT bound with c=1 control.
5. **d=6 sterile scalar:** FAIL decisively, O(10^-2)% or less.
6. **EW dipole:** remains demoted for hard-UHE enhancement.

The non-resonant weakly-coupled question is therefore narrowed sharply: **the resolved derivative mediator is the only branch from this v0.32 pair that still merits a dedicated phenomenological recast.**

## 6. Reproducibility

Final corrected numerical results:

- `results/v0.32/sterile_d6_resolved_dynamicQ_v0_32.csv`
- `results/v0.32/sterile_d6_resolved_dynamicQ_v0_32.json`
- `results/v0.32/sterile_d6_resolved_summary_v0_32.md`

Code:

- `scripts/sterile_d6_resolved_uhe_v0_32.py`
- `.github/workflows/sterile-d6-resolved-v032.yml`

Corrected numerical result commit: `196f2996413cb51c0edb5701023161758c5c5144`.

Relevant literature anchors:

- Pandey, Karmakar, Rakshit, JHEP 11 (2019) 046, arXiv:1907.07700.
- EPJC 84 (2024) 795, DOI 10.1140/epjc/s10052-024-13104-x, for HighPT bounds on O_lNqd^S,T and B->K(*)+invisible discussion.

## Decision

Do not run Earth propagation for d=6 S/T.

For the resolved derivative mediator, the next and only justified gate is a **modern mass-dependent constraint recast**. If an allowed point satisfies kappa >= kappa_5%(M), then run sterile-sink Earth attenuation with the repaired transport solver. If all allowed kappa lie below that curve, close the non-resonant weakly-coupled branch and retain only strong/composite and resonance models.
