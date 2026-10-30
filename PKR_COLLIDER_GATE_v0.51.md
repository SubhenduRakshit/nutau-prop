# PKR collider gate v0.51

> **FINAL STATUS: FAILED NORMALIZATION / NO COLLIDER CONCLUSION.**
>
> The four committed v0.51 benchmarks reported the same MG5 cross section,
> `30.81509 pb`, despite `BR_invisible` spanning approximately `0.370--0.590`.
> For these extremely narrow mediators the NWA requires
> `sigma(pp->Z'j->inv) = sigma(pp->Z'j) * BR_inv`; therefore `sigma/BR_inv`
> should be constant (up to small integration/NWA corrections), which v0.51
> does not satisfy.  The resulting CheckMATE single-bin values `r ~ 14--27`
> have **no physical interpretation** and must not be quoted even qualitatively.
> Eq. (4.3) and Eq. (4.4) remain open pending a correctly normalized collider
> screen.  v0.51 is retained only as a failed-normalization audit artifact.

## Purpose

Test the conditional PKR benchmark

- `M_Zp = 50 GeV`
- `g_q = 0.03`
- Eq. (4.3) dipole and Eq. (4.4) derivative couplings fixed by the converged v0.42d 5%-at-100-PeV benchmark

against the current CMS jets+MET / monojet implementation in CheckMATE2.

CMS EXO-20-004 is arXiv:2107.13021. CheckMATE2 internally labels this implementation `cms_2107_13201` (the last digits are transposed relative to the arXiv identifier). Its `analyses.json` maps INSPIRE 1894408 to this internal name. v0.51 therefore uses `cms_2107_13201` exactly as CheckMATE expects.

## Intended width-matched logic

The v0.45 benchmarks have `Gamma_Zp/M_Zp < 9e-4`. In the narrow-width limit the visible recoil process factorizes,

`pp -> Z' + j -> invisible + j`,

and the measured jet/MET kinematics are controlled by the vector quark production current and `BR(Z'->invisible)`. The internal Lorentz structure of an entirely invisible on-shell decay is not measured.

v0.51 therefore attempted to use the standard `DMsimp_s_spin1` UFO with the same vector quark coupling and a nearly massless Dirac invisible state, choosing its vector coupling so that

`Gamma(Y1 -> Xd Xdbar) = Gamma_invisible^PKR`,

while `WY1` was explicitly intended to be fixed to the v0.45 PKR total width.

The generator-level normalization diagnostic shows that this intended matching was not demonstrated by the produced samples.  Consequently the downstream detector results are invalid as physics constraints.

## Event chain attempted

1. MG5_aMC LO: `p p > xd xd~ j`, 13 TeV, CT18NLO, `pT(j)>100 GeV`.
2. Pythia8 ISR/FSR/MPI/hadronization.
3. HepMC2 event record.
4. CheckMATE2 detector simulation and internal `cms_2107_13201` analysis (CMS EXO-20-004 / arXiv:2107.13021).
5. Multibin post-processing also failed on macOS/Spey, leaving only single-bin diagnostic values; these are invalid anyway because the upstream normalization failed.

## Replacement certification gate

v0.52 must pass all of the following before any new detector run:

1. read back each persisted `param_card.dat` and verify `M`, `g_q`, `g_chi`, and total width;
2. verify independent generated-event payload hashes (no sample reuse);
3. compare MG5's persisted banner cross section with an independent LHE `<init>` extraction;
4. require `sigma/BR_inv` to be constant within 2% across three fixed-`M`, fixed-`g_q` test points;
5. only after a PASS, run exactly one CheckMATE benchmark.

## Audit status

v0.47--v0.50 are retained as failed/deprecated direct-likelihood attempts. v0.51 is retained as the failed-normalization attempt described above. No PKR exclusion or survival conclusion is taken from any of these versions.
