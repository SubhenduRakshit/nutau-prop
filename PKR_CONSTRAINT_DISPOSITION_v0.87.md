# v0.87 — exact collider inputs and the actual IceCube completion boundary

31 October 2026. Continues the user-authorized Eq. 4.3/4.4 task after v0.86.

## Responsibility and scope

v0.86 finished the normalization correction but not the assigned full collider/
IceCube assessment. Describing that checkpoint as completion of the requested
job was premature. This follow-up finishes an accessible input that was left
approximate and documents the precise experimental/model information still
needed for a faithful IceCube fit. It does not label an approximate fit exact.

## Exact CMS dijet table obtained and applied

The browser failed to retrieve HEPData previously, but the Python HTTP client
successfully obtained its JSON endpoint. The table is preserved verbatim in
`results/v0.87_constraints/cms_dijet_figure2a.json` and parsed into CSV.
Source: [CMS Figure 2a, HEPData DOI 10.17182/hepdata.165379.v1/t3](https://www.hepdata.net/record/165379).
Observed universal-vector limits for BR(qq)=1 are:

| M (GeV) | observed gq limit |
|---:|---:|
| 50 | 0.035670 |
| 100 | 0.082596 |
| 220 | 0.114380 |

These are actual table entries; no graph digitization or mass interpolation is
used. Invisible branching is included through gq sqrt(BR_qq), in the narrow-
width approximation. We retain the calibrated mass-specific monojet shapes
from v0.81/v0.79 and the explicit v0.86 kernels/widths. The validity prescription
is sqrt(shat)<1/a, with 2e-5<=a<=0.999/M; the convention/UV qualification in
v0.86 remains. This is not a simultaneous experimental confidence construction.

Analytic monotonicity in z=gq^2 eliminates a brute-force coupling scan. All
one-dimensional a intervals are searched. Each maximizing point is directly
folded with dynamic-Q CT18 at orders72 and112 using14 worker processes, with
one numerical-library thread per process. Interpolation and quadrature errors
pass the1% gate; metadata gives their actual much smaller values.

### Nominal collider-compatible maxima at100 PeV

Here R=mean(sigma_NP/sigma_SM_nu, sigma_NP/sigma_SM_antinu), SM=CC+NC,
with nominal monojet r<=1 and the exact dijet bound.

| M (GeV) | Eq. 4.4 maximum R | Eq. 4.3 maximum R |
|---:|---:|---:|
| 50 | 0.3395% | 0.7953% |
| 100 | 0.0903% | 0.1084% |
| 220 | 0.0456% | 0.0581% |

At the same maximizing couplings, the220-PeV values are0.3373/0.7578%,
0.0951/0.1162%, and0.0497/0.0647% respectively (44/43). These are evaluations
at the100-PeV optima, not separately optimized220-PeV upper envelopes.

Even relaxing monojet to r<=1.2, the50-GeV maxima are0.4063% (44) and0.9568%
(43). The loose global dijet bound0.13 used in v0.86 is no longer needed at
these three masses.

### Direct5% target test

At fixed a the minimal coupling needed for the target is
z_required=0.05/F(a), where R=z F(a). Since both dijet and monojet yields
increase with z, testing this coupling is sufficient; larger z cannot rescue
a failed collider test.

- Eq. 4.4 cannot reach5% anywhere in the tested domain after just the exact
  dijet and width screens: maxima without monojet are2.729%,1.581%,0.3615%
  at50,100,220GeV respectively.
- Eq. 4.3 can reach5% after those two screens at50/100GeV, but its minimum
  calibrated monojet r is5.832/54.71 respectively. Direct CT18 evaluation of
  these target points gives4.99961%/4.99993%, validating the interpolation.
- At220GeV Eq. 4.3 reaches only2.483% even without monojet.

Thus the5% target fails the current collider implementation at these three
masses by appreciable margins. An IceCube fit is not necessary to make that
restricted screening decision. This does not close5GeV,40GeV, other masses,
different flavours, different UV/cutoff assumptions or the whole propagation
program. M100/M220 still use pilot monojet shapes. Finite-width, higher-order
and model-specific experimental validation are not supplied by this scan.

## What a complete IceCube calculation actually needs

The [official HESE7.5-year release](https://icecube.wisc.edu/data-releases/2021/12/hese-7-5-year-data/)
provides observed events, MC, detector-systematic splines and example likelihood
code. Those resources are already present locally. Therefore saying simply
"IceCube has not released a likelihood" would be incorrect for HESE.

We inspected the complete keys of all three local MC JSON files, each with
821764 raw entries, and saved their names, lengths and hashes. Interaction
type is present in the raw truth file, although the example loader used by
v0.43 drops that field. True inelasticity y, Bjorken x, and true hadronic shower
energy are absent. Reconstructed deposited energy is provided and is not
interchangeable with true yE. The public schema is also documented in the
[collaboration repository](https://github.com/icecube/HESE-7-year-data-release).

This matters physically: nu_tau q -> N q still deposits hadronic energy yE.
Calling the outgoing neutrino sterile does not make the detector interaction
invisible. v0.43 modified incident astrophysical fluxes only; it did not add
this BSM shower contribution with its changed y distribution. Profiling the
available detector splines cannot repair a missing signal response.

For a reconstructed bin i, the released SM NC response contains a marginal

    A_i(E) = integral dy R_i(E,y) f_SM(y|E).

The desired prediction instead uses f_NP. Knowing A_i does not determine that
integral unless f_NP/f_SM is constant or R_i is known. The script records an
exact two-bin counterexample: f_SM=(.5,.5), f_NP=(.75,.25), and admissible
responses R_A=(.2,.8), R_B=(.8,.2). Both give SM acceptance.5, but BSM acceptance
.35 and.65. This proves that an exact arbitrary-y reweight cannot be recovered
from the marginal response alone. It is not a claim that real IceCube has
these toy response values.

The [12-year HESE release](https://icecube.wisc.edu/data-releases/2023/07/icecube-hese-12-year-data-release/)
contains directional probability maps, not the missing hadronic-y response.
The [12.6-year EHE release](https://icecube.wisc.edu/data-releases/2025/02/a-search-for-extremely-high-energy-neutrinos-and-first-constraints-on-the-ultra-high-energy-cosmic-ray-proton-fraction-with-icecube/)
points to DOI10.7910/DVN/JHK49D. Today's metadata requests returned403 through
multiple public API forms. The local repo contains effective area and a
differential limit, not the reconstructed-bin response/background needed to
rebuild the full EHE likelihood. We do not claim to have reverified the entire
current remote EHE inventory while that endpoint is inaccessible.

## A necessary transport correction before any long fit

For the specified stable, massless N realization, the Hermitian transition
also allows N q -> nu_tau q. The active collision term is a sink, but the full
transport problem is not an irreversible sink unless reconversion is negligible
or an additional removal/decay mechanism is specified. The old dirac_sink
generator does not carry an N population.

This is apparent analytically in the fixed-energy limiting system

    d/dX [phi_nu,phi_N]^T = [[-S-T,T],[T,-T]] [phi_nu,phi_N]^T,

where S is an SM removal rate and T a bidirectional transition rate (the column
normalization absorbs N_A). The attenuation eigenvalues are

    lambda_slow,fast = (S+2T -/+ sqrt(S^2+4T^2))/2.

For incident (1,0), the active solution has a slow component with coefficient
[1-S/sqrt(S^2+4T^2)]/2 ~ (T/S)^2. At small T/S it can still dominate an opaque
chord whose direct SM flux is exponentially tiny. A pure-sink approximation
therefore needs a quantified domain; it cannot be declared exact globally.

The closed form was checked against a2x2 matrix exponential, including T=0,
S=0 number conservation, and thin/opaque cases. This is a diagnostic, not an
Earth prediction: the full problem needs inelastic redistribution, nu/antinu,
tau production/loss/decay and secondary flavours. The literal active-flavour
PKR realization needs its own off-diagonal active/antineutrino system, not this
sterile model. N-initiated detector showers also matter if an N flux reaches it.

## The attached ATZ/PKR reconstruction

The input is archived as an attachment hash in this checkpoint. Its relative
normalization argument supports the v0.86 trace result under matched external
states/conventions. The old5-GeV bound pair does not itself supply a modern
tau-specific detector likelihood or identify an absolute kernel uniquely.

The [full2021 ATZ erratum](https://doi.org/10.1007/JHEP11(2021)113) corrects
dimension-six DIS bounds in Eqs.(6.22),(6.23) and recalculates Rayleigh bounds.
Its MadGraph reference is explicitly for tree-level Rayleigh nu-e scattering
with photon emission. It does not state a replacement for the derivative DIS
Eqs.(77),(78). Thus it neither resolves that discrepancy nor establishes that
MadGraph uncovered those particular coefficients. A plausible lineage from
printed ATZ formulas to later code is not proof of PKR's original implementation.

## Precise remaining dependency, not a generic "unresolved" label

To produce the requested faithful IceCube exclusion one must obtain a
channel-resolved response in true E,y,direction (and any necessary x/hadronic
dependence), or MC truth sufficient for an exact hard-scattering reweight;
then add new hadronic showers and use a consistent bidirectional transport
solution. For EHE one also needs reconstructed-bin observations/backgrounds
and the corresponding response. `docs/ICECUBE_RESPONSE_REQUIREMENTS_v0.87.md`
specifies the deliverable required from a release or detector simulation.
No message was sent to collaborators and no missing response was fabricated.

The accessible exact-dijet calculation is now complete. A full IceCube result
is not claimed: replacing its missing response with a proxy would not satisfy
the user's requirement for a reliable calculation. The project remains short
of the original all-constraints propagation/publication goal, despite a much
stronger disposition for the three tested collider benchmarks.
