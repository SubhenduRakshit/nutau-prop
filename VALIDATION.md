# Validation protocol — v0.3

The BSM result is frozen until the historical transport gate passes.

## Why v0.2 was rejected

The first real run showed a 74.6% third-vs-second iterative change, so v0.2
was not accepted.  The failure led to an audit that found:
- incomplete tau->nu_tau decay physics;
- neutrino-only historical comparisons instead of nu+antinu;
- an ad hoc Picard update that was not the Rakshit-Reya Z iteration.

Those are corrected in v0.3.

## Rakshit-Reya gate

v0.3 uses equal incident nu_tau and antinu_tau fluxes, propagates them
separately, and combines them at the end.

Two independent solvers are compared:
1. direct stepwise coupled nu_tau/tau transport;
2. RR-style effective-Z iteration with a tau-decay-aware Z^(0).

Mandatory checks:
- max |F^(3)/F^(2)-1| < 5% in relevant bins;
- the second iteration must be visibly non-negligible for the hard vertical
  E^-1 spectrum;
- direct transport must agree with Z^(2) to 10%;
- the hard E^-1 spectrum must exhibit the 10^4-10^5 GeV regeneration bump;
- E^-2 must not exhibit a comparable bump.

Outputs:
- `results/historical/rr_reproduction.png`
- `results/historical/rr_reproduction.csv`
- `results/historical/rr_gate.json`

## Dutta-Reno-Sarcevic gate

The inclusive tau decay distribution is taken from PRD 62, 123001,
Appendix A/Table I.  The secondary nu_mu/antinu_mu channel uses the polarized
leptonic spectrum with branching fraction B=0.18.

Checks:
- vertical E^-1 flux near 1 PeV must reproduce the historical scale:
  nu_tau ~6% and nu_mu/e ~1% of incident flux, using broad ranges because
  the historical work used older CTEQ5-era cross sections;
- direct secondary flux must be compatible with their approximate Eq. (5)
  `F_secondary ~= B [F_nutau - F_numu]`;
- the secondary contribution must be much smaller for E^-2 than E^-1.

Outputs:
- `results/historical/drs_fig1_reproduction.png`
- `results/historical/drs_fig1_reproduction.csv`
- `results/historical/drs_gate.json`

## Combined gate

`results/historical_gate.json` must contain `"PASS": true`.

`run_all.sh` stops before the BSM scan if either historical branch fails.

## Remaining work after a historical PASS

A historical PASS validates the transport machinery. It does not yet establish
paper-level sensitivity. Remaining gates include:
- NLO/NNLO cross-section systematics;
- nuclear PDFs;
- event weighting / detector response;
- DRS/RR upward-muon event-rate reproduction;
- dedicated charge-radius antineutrino convention audit.


## v0.4 detector-level gate

The official 2015 IceCube contained-event total tau effective area is folded
with the published per-flavor `nu+nubar` power-law fit.

This stage is intentionally downstream of the historical transport gate.

The code reports both:
- **propagation-only**: public SM Aeff × `T_BSM/T_SM`;
- **propagation+detector proxy**: additionally multiply by the inclusive
  CC+NC BSM/SM cross-section ratio.

The main conservative statistic is an Asimov Poisson shape test after profiling
both an overall astrophysical normalization and spectral tilt.  This avoids
mistaking a nearly pure rate rescaling for a discoverable BSM signature.

A later paper-level refinement should replace the inclusive detector proxy with
the public per-channel CC/NC/topology response tables or the newer HESE Monte
Carlo release.


## v0.5 future publishability gate

### External-bound scan
Charge-radius benchmarks are extended through `|r_tau^2|=20e-32 cm^2`,
the 2025 global-fit 1σ interval.  The result is flagged as a phenomenological
high-Q benchmark because the modern charge-radius literature stresses
momentum dependence.

### Gen2 scaling
No unavailable response matrix is fabricated.  The already validated IceCube
binning is multiplied uniformly by 5 and 10 for a 10-year exposure.  Hence
the nuisance-profiled Asimov `q` scales linearly and `Z` as the square root of
the exposure/acceptance factor.

### Energy-localization gate
The post-profile likelihood contribution is reported in four true-energy
bands.  This is essential for deciding whether a larger coupling genuinely
moves sensitivity toward ~100 TeV or only increases a PeV effect with very
poor statistics.

### Generic future tau telescopes
The package reports the acceptance multiplier needed for 2σ and 3σ in ten
years.  This avoids assigning unsupported effective areas to proposed
Earth-skimming detectors whose energy response differs substantially from the
IceCube contained-event sample.


## v0.6 KM3NeT/ARCA gate

The ARCA calculation uses the full-detector performance paper (EPJC 84, 885)
and the associated official Zenodo record 10949703.

The selected shower response is decomposed into:
- electron-neutrino CC;
- tau-neutrino CC shower decays;
- all-flavour NC;
and into upgoing/down-going hemispheres.

The track control response uses the four published true-zenith ranges for
selected muon-neutrino CC events.

The digitised response shapes are separately normalized to the exact published
cosmic-neutrino rates of 220 tracks/year and 60 showers/year for the reference
E^-2 flux.  This is used only until the exact Zenodo numerical response is
machine-identified; `zenodo_inventory.json` is retained for that audit.

The likelihood profiles a common astrophysical normalization and spectral
index.  Results are given for showers-only and tracks+showers.

Because bin-by-bin atmospheric backgrounds are not included in the fallback
response, the significance is an optimistic upper bound.  A result below
2 sigma is therefore a robust negative gate.  A result above 2 sigma is only a
motivation for a second-stage exact-response/background analysis.
