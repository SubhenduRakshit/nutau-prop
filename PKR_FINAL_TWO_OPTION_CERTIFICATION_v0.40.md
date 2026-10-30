# PKR final two-option certification v0.40

Date: 2026-10-30

## Scope
Final targeted certification of the two remaining non-resonant PKR resolved-mediator candidates for a >=5% modification of the inclusive UHE neutrino-nucleon cross section near 100--300 PeV:

1. PKR Eq. (4.4), resolved derivative Z':
   kappa_2 = c2*gq/Lambda.
2. PKR Eq. (4.3), resolved transition-dipole Z' (off-diagonal i != j):
   kappa_1 = c1*gq/Lambda.

The common quark interaction is gq qbar gamma_mu q Z'^mu. The resolved propagator is kept exactly.

## Critical correction: split couplings and the validity scale
The previous self-consistent scans imposed sqrt(shat) < 1/kappa. That condition is correct only for the special choice gq=1. In general,

kappa_i = (c_i/Lambda) gq,

so the dimension-5 neutrino-vertex scale is

Lambda_vert = Lambda/c_i = gq/kappa_i.

Therefore the correct partonic validity condition is

sqrt(shat) < gq/kappa_i.

This split-coupling correction is binding once direct dijet searches force gq below unity.

## Exact kernel relation used
For the same off-diagonal/sterile external-state interpretation and massless DIS limit, the validated hard kernels obey

d2sigma_4.3/dxdy = 2 d2sigma_4.4/dxdy

at equal M_Z' and kappa. This relation is used only for that specified interpretation.

## UHE size benchmarks (gq=1 proxy)
The earlier v0.35/v0.36 scans found 5% thresholds at 100 PeV of approximately:

- Eq. (4.4): 3.09e-5, 5.43e-5, 7.97e-5, 1.20e-4, 2.16e-4 GeV^-1 at M=5,40,100,220,500 GeV.
- Eq. (4.3): 2.19e-5, 3.84e-5, 5.64e-5, 8.46e-5, 1.44e-4, 3.30e-4 GeV^-1 at M=5,40,100,220,500,1000 GeV.

These values are diagnostics, not full-model benchmarks, because gq must be separated.

## IceCube-range downward energy screen
Using each 100-PeV 5% benchmark and the same coupling at lower energy:

- M=5 GeV gives roughly 5--8% from 0.1 to 10 PeV for both kernels.
- M=40 GeV gives about 2.3%, 3.2%, 3.6%, 4.1% at 0.1,1,3,10 PeV.
- M=100 GeV gives about 0.8%,2.0%,2.7%,3.4%.
- M=220 GeV gives about 0.2%,1.1%,1.8%,2.8%.
- M>=500 GeV is still smaller in the directly measured IceCube range.

Thus the low-mass point is the one most exposed to a modern IceCube likelihood. The higher-mass points do not produce an obviously excluded O(1) PeV-scale distortion. A proper event/zenith likelihood remains required for a certified IceCube bound.

## Split-coupling validity result
At fixed old kappa_5% benchmark, reducing gq lowers Lambda_vert=gq/kappa and removes the high-shat part of the UHE integral. Examples at 100 PeV:

- M=500 GeV, gq=0.07:
  Eq. (4.4) -> 0.085%; Eq. (4.3) -> 0.235%.
- M=1000 GeV, gq=0.04:
  Eq. (4.3) -> 0.001% at the old kappa_5% benchmark.

Therefore one must optimize over kappa at each gq rather than merely rescale the old benchmark.

## Absolute peak-capacity test
Let k' = kappa/gq = c/Lambda = 1/Lambda_vert. At fixed gq,

R_4.4^phys(k') = gq^2 R_deriv(k'),
R_4.3^phys(k') = 2 gq^2 R_deriv(k').

Scanning k' therefore gives the absolute maximum UHE effect attainable at that gq while enforcing sqrt(shat)<1/k'. CT18, 100 PeV results:

| M [GeV] | Eq4.4 max @ gq=.04 | @ .07 | @ .10 | @ .30 | Eq4.3 max @ .04 | @ .07 | @ .10 | @ .30 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 5 | 306.5% | 938.7% | 1915.7% | 17241% | 613.0% | 1877% | 3831% | 34483% |
| 40 | 23.91% | 73.22% | 149.44% | 1345% | 47.82% | 146.45% | 298.88% | 2690% |
| 100 | 2.88% | 8.82% | 18.01% | 162.1% | 5.76% | 17.65% | 36.01% | 324.1% |
| 220 | 0.393% | 1.205% | 2.459% | 22.13% | 0.787% | 2.410% | 4.918% | 44.26% |
| 500 | 0.043% | 0.131% | 0.267% | 2.400% | 0.085% | 0.261% | 0.533% | 4.800% |
| 1000 | 0.006% | 0.018% | 0.037% | 0.335% | 0.012% | 0.037% | 0.075% | 0.671% |
| 2000 | 0.001% | 0.002% | 0.004% | 0.040% | 0.001% | 0.004% | 0.009% | 0.081% |

The approximate minimum gq required for any 5% solution follows from the gq^2 scaling:

| M [GeV] | gq,crit Eq4.4 | gq,crit Eq4.3 |
|---:|---:|---:|
| 5 | 0.0051 | 0.0036 |
| 40 | 0.0183 | 0.0129 |
| 100 | 0.0527 | 0.0373 |
| 220 | 0.1426 | 0.1008 |
| 500 | 0.433 | 0.306 |
| 1000 | >1 | 0.819 |
| 2000 | >1 | >1 |

These critical values are a size+validity criterion, not experimental limits.

## Modern direct-quark constraints and mass-region verdict
Current public dijet resonance results strongly constrain a flavor-universal leptophobic quark coupling. The exact gq limit is analysis- and width-dependent, so the following uses only robust published/current anchors and does not invent a continuous curve.

### 500 GeV
Current ATLAS low-mass dijet constraints reach approximately gq < 0.07 in the 250--650 GeV interval. Since the 5% requirement is gq >~0.433 (Eq4.4) or >~0.306 (Eq4.3), both branches are closed by orders of magnitude in coupling-space once the correct vertex-validity scale is imposed.

Verdict: CLOSED for both.

### 1 TeV
Current CMS dijet scouting constraints are at the few-times-10^-2 to O(0.1) level; a representative current 1-TeV limit is around gq~0.04. Eq4.4 cannot reach 5% even at gq=1 in the peak scan, and Eq4.3 requires gq~0.82.

Verdict: CLOSED for both.

### 220 GeV
A conservative direct-dijet anchor below 250 GeV is gq of order 0.2. The size+validity requirement is ~0.143 for Eq4.4 and ~0.101 for Eq4.3. Thus direct dijet information alone does not robustly close this mass; the modern monojet likelihood is binding.

Verdict: OPEN BUT NARROW / MODERN MONOJET RECAST REQUIRED. Eq4.3 has the larger margin.

### 100 GeV
The required couplings are only ~0.053 (Eq4.4) and ~0.037 (Eq4.3). Existing direct-dijet constraints at this mass do not by themselves eliminate such values.

Verdict: OPEN / MODERN MONOJET + DIRECT-Z' RECAST REQUIRED.

### 40 GeV
Required gq is only ~0.018 (Eq4.4) or ~0.013 (Eq4.3). This mass lies below the standard CMS low-mass dijet range used for the simple direct anchor; specialized low-mass/quarkonium/electroweak constraints must be used.

Verdict: OPEN / SPECIALIZED DIRECT-Z' + MONOJET RECAST REQUIRED.

### 5 GeV
Required gq is only ~0.005 (Eq4.4) or ~0.004 (Eq4.3), so size+validity is easy. Historical PKR IceCube limits are much weaker than the required kappa, but the benchmark generates a 5--8% effect throughout much of the modern IceCube cross-section range. Modern IceCube event/zenith likelihood and low-mass Z' constraints are therefore the binding tests.

Verdict: OPEN / ICECUBE + LOW-MASS DIRECT CONSTRAINTS REQUIRED.

## Collider monojet status
The full Run-2 ATLAS monojet analysis (139 fb^-1) is the correct current collider dataset. It publishes model-independent visible-cross-section limits but not the PKR dimension-5 operators. A modern kappa(M) bound cannot be obtained by luminosity rescaling. Signal generation, acceptance through the ATLAS MET bins, detector efficiency validation, and the likelihood must be reproduced. Existing modern Z' monojet recasts with renormalizable vector interactions demonstrate the workflow but cannot be copied numerically onto Eqs. (4.3)/(4.4).

Therefore no modern operator-specific monojet curve is claimed here.

## Width and perturbativity
The quark width is controlled separately by gq and is narrow for the small couplings relevant to the surviving low-mass window. The neutrino derivative/dipole width depends on c_i/Lambda = kappa/gq and therefore on the split couplings, not kappa alone. The validity scan already prevents use of partonic energies above Lambda_vert, but a final benchmark must additionally check Gamma_Z'/M explicitly. PKR already noted a severe high-mass width issue for Eq. (4.3); the high-mass branch is independently closed here by direct-gq plus validity.

## Final verdict
The two branches are NOT globally closed, but their viable territory has collapsed to a light/intermediate mediator loophole:

- M >= ~500 GeV: CLOSED for both in perturbative resolved-mediator form by direct-gq constraints plus the correct split-coupling validity condition.
- M ~220 GeV: conditional narrow window; Eq. (4.3) easier to keep above 5%.
- M ~5--100 GeV: genuine non-resonant loophole remains for both; modern monojet/IceCube/direct-Z' likelihoods are the remaining certification gate.

Between the two, Eq. (4.4) remains the cleaner model target because it avoids the extra transition-dipole/kinetic-mixing constraints; Eq. (4.3) wins only on raw cross-section efficiency by a factor of two in the hard kernel.

### Recommended lead benchmark class

Resolved derivative Eq. (4.4), M_Z' ~ 40--100 GeV, with gq kept explicit and chosen above the size+validity critical value but below direct-dijet bounds. This region has a modest effect in the currently measured IceCube PeV range while retaining a 5% UHE effect at 100--300 PeV.

Eq. (4.3) should be retained as the stronger fallback, especially near 100--220 GeV, but requires the additional transition-moment/kinetic-mixing and width audit.

## No Earth transport yet
No Earth-propagation run is justified until one explicit (M_Z', gq, c_i/Lambda) benchmark passes the modern operator-specific monojet likelihood, modern IceCube event/zenith likelihood, direct-Z' constraints, and explicit width test.
