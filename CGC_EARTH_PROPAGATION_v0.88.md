# v0.88 — CGC/rcBK all-flavour Earth-propagation benchmark

Date: 2026-10-31

## Physics question

Quantify whether small-x saturation/non-linear QCD can produce an observable UHE
Earth-transmission distortion relative to the CT18/DGLAP-like SM baseline, for
all three active neutrino flavours.

This is a QCD control model, not BSM.

## CGC input

Primary reference:

- J. L. Albacete, J. I. Illana and A. Soto-Ontoso, *The neutrino-nucleon cross
  section at UHE and its astrophysical implications*, Phys. Rev. D 92, 014027
  (2015), arXiv:1505.06583.

That work finds DGLAP and running-coupling BK predictions similar through
approximately E_nu ~ 1e8 GeV, followed by systematically slower BK growth.  It
states that at E_nu=1e9 GeV the rcBK+DLA cross section is about 1/1.5 of the
DGLAP result, and that the difference reaches about a factor 4.5 at 1e14 GeV.
The reported 1e12--1e14 GeV flux-limit shifts correspond to cross-section
suppression factors of roughly 2.5--4.5 over that range.

v0.88 encodes those published statements as a log-log interpolation of
R_CGC(E)=sigma_CGC/sigma_DGLAP through

- 1e4 GeV: 1,
- 1e8 GeV: 1,
- 1e9 GeV: 1/1.5,
- 1e12 GeV: 1/2.5,
- 1e14 GeV: 1/4.5.

The 1e12 value is a representative anchor within the published 1e12--1e14
2.5--4.5 range, not a digitized exact curve point.  Therefore the object is
named a literature-calibrated **rcBK+DLA-2015 transport surrogate**.

A useful comparison is V. P. Goncalves, D. R. Gratieri and A. S. C. Quadros,
Eur. Phys. J. C 81, 496 (2021), which independently finds that CGC/IIMS modifies
the neutrino cross section and consequently the Earth absorption function.

## Differential-kernel qualification

The cited papers do not provide a machine-readable full CGC d sigma/dy grid.
v0.88 therefore multiplies both SM CC and NC totals by R_CGC(E) while retaining
the conditional SM y distributions.  This preserves the published total
cross-section effect and the validated energy-migration machinery without
inventing an unpublished y distribution.

This assumption must be varied or replaced by a direct dipole-model d sigma/dy
calculation before claiming precision CGC phenomenology.

## Baseline and transport

- CT18NLO isoscalar nucleon baseline through the existing `Case6DIS` SM kernel.
- PREM density profile.
- Separate neutrino and antineutrino propagation.
- nu_e, nu_mu: CC absorption + NC down-scattering.
- nu_tau: CC tau production, continuous tau loss, tau decay, nu_tau
  regeneration and NC down-scattering using the validated mixed-grid sparse
  exponential solver.
- tau leptonic decay secondaries are added to both e and mu flavours with
  branching fraction 0.18 per flavour and the existing massless daughter
  spectrum.
- UHE oscillations inside Earth are neglected.
- Glashow resonance is omitted because the principal report range is >=10 PeV.

The incident astrophysical mix is e:mu:tau=1:1:1 with equal neutrino and
antineutrino fluxes.  Both E^-2 and E^-2.5 spectra are propagated; E^-2.5 is the
primary printed benchmark.

## Outputs

`results/v0.88_cgc_allflavor/` contains

- `cross_section_ratio.csv` — CGC/DGLAP cross-section surrogate,
- `flux_grid.csv` — fixed-chord all-flavour fluxes and CGC/SM ratios,
- `upgoing_angle_average.csv` — isotropic upgoing angle average,
- `summary.csv` — target-energy/chord summary,
- `metadata.json` and `run.txt`.

Target energies are 10, 30, 100, 300 PeV, 1, 3 and 10 EeV.  Fixed report
chords include |cos(theta_z)|=0.05,0.10,0.20,0.30,0.50,0.70,1.00; a 12-point
Gauss-Legendre integration gives the upgoing-angle average.

## Interpretation rule

The primary observable is not Delta sigma/sigma by itself but

    Delta phi_alpha(E,theta)/phi_alpha^SM(E,theta)

for alpha=e,mu,tau, together with the angle-averaged upgoing ratio.  Very large
ratios in bins where the SM transmitted flux is exponentially tiny should not
be interpreted as high-statistics observability; the angle-averaged result and
horizon chords are the physically useful quantities.
