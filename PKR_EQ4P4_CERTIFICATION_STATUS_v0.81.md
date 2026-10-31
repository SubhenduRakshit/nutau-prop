# PKR Eq. (4.4) certification status through v0.81

Date: 2026-10-31

## Canonical conclusion

PKR Eq. (4.4), in the project active-to-sterile resolved-Z' realization, is **not ruled out** and is **not yet certified collider-allowed**.

The currently least-visible screened mass point is M=50 GeV.  Its deterministic same-shape CMS EXO-20-004 surrogate recast gives

- gq = 0.038
- a = c/Lambda = 2.379776052e-3 GeV^-1
- BR_inv = 0.246335804
- gq^2 BR_inv = 3.557089009e-4
- normalization scale relative to v0.75 = 0.907532728
- multibin mu95 = 0.968238709
- multibin r = 1.032803162
- single-bin robscons = 0.404549024

With the predeclared project thresholds r<0.8 (robust survival), 0.8<r<1.2 (boundary unresolved), r>1.2 (robust exclusion), the canonical monojet status is

**MARGINAL_BOUNDARY_UNRESOLVED**.

## Audit trail

- v0.73 fixed the physical O+h.c. Eq. (4.4) invisible-width normalization.
- v0.74 found a narrow 5%-UHE M=50 screen window, but kept only the first a crossing.
- v0.75 generated one fresh 100k-event CT18/CheckMATE point and found r=1.1380; normalization passed at 0.051%, while selected-bin shape differences were sizable.  It was correctly classified boundary-unresolved.
- v0.76 introduced mass-dependent dijet screening but still used only the first a crossing and incorrectly ranked masses by gq^2 BR_inv alone.
- v0.78 fixed the full-a problem: all connected >=5%-UHE a intervals are retained.  M=50,100,220 GeV all contain some portion passing the current *screening* dijet prescription.
- v0.79 fixed cross-mass ranking by adding mass-dependent production x acceptance.  The 50-GeV representative is clearly the least collider-visible of 50,100,220 GeV; the ranking is separated well beyond pilot MC uncertainty.  No exclusion claim is taken from the 20k-event 100/220-GeV pilots.
- v0.80 changed the CheckMATE/Delphes random seed and therefore was not a pure normalization reweight.  Its r=0.7419 is diagnostic-only and must not be used as a survival certificate.
- v0.81 repeats the reweight with the identical v0.75 seed 7501 and explicitly verifies all 66 signal bins and ds values rescale by the common normalization factor.  The fixed-shape check passes at <5e-7 and r=1.032803 agrees with the linear diagnostic at 1.3e-6 relative accuracy.  v0.81 supersedes v0.80 for interpretation.

## Remaining Eq. (4.4) collider issue

The low-mass dijet input used in v0.76/v0.78 is still a digitized screening prescription.  CMS EXO-24-007 assumes BR(Z'->qq)=1 and publishes the mass-dependent universal-gq limit in Figure 2a, with tabulated data in HEPData DOI 10.17182/hepdata.165379.  Our NWA visible-rate mapping uses gq_equiv = gq sqrt(BR_qq).  For the v0.81 representative, gq_equiv is approximately 0.033, close enough to the published 50-GeV curve that the exact tabulated value should be used before a final dijet pass/fail statement.

Therefore the present Eq. (4.4) status is:

**UHE size: PASS at screened points. Width normalization: PASS. Generator NWA normalization: PASS. Monojet surrogate: MARGINAL. Dijet: exact 50-GeV tabulated comparison still required. Overall: OPEN.**

No general PKR or non-resonant no-go follows from Eq. (4.4).

## Next branch

Proceed independently to PKR Eq. (4.3) transition-dipole certification.  Do not transfer Eq. (4.4) monojet conclusions to Eq. (4.3).  First fix the exact dipole width/h.c. counting and coupling separation, then build the complete >=5%-UHE (M,gq,a) envelope before any new expensive detector simulation.
