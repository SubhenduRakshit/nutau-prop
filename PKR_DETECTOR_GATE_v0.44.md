# PKR detector gate v0.44

## Scope

This note closes the detector-sensitivity branch for the constrained PKR Eq. (4.3) and Eq. (4.4) benchmarks studied in v0.42d-v0.44. It does **not** establish a detector no-go theorem for every BSM interaction in the broader `nutau-prop` portfolio.

Benchmarks use the v0.42d coupled transport at `M=50 GeV`, `g_q=0.03`, normalized to a 5% NP/SM neutrino-nucleon cross-section contribution at 100 PeV, with the pointwise validity cut `sqrt(shat) < g_q/kappa`.

## Corrected observable

The useful propagation signal is the broad absorption-horizon strip

`1 <= tau_SM <= 3`,

not the nearly opaque vertical chord. v0.42d established that the coupled tau-regeneration calculation retains several-percent local `nu_tau` distortions in this band.

v0.44 corrects an important detector-fold issue: the published ARCA upgoing effective area already contains SM Earth absorption. Therefore the detector calculation uses the incident astrophysical flux times the published SM effective area and applies only the coupled `F_BSM/F_SM` transport modifier. It does not multiply an already propagated flux by an absorption-containing effective area.

The v0.44 fold includes:

- reconstructed-energy migration with `sigma_log10E = 0.3`;
- fine zenith slicing over the available v0.42d absorption-horizon region;
- showers and tracks;
- the change in total muon flavor, including tau-decay secondary `nu_mu`;
- direct PKR hadronic shower production through an NC-response proxy;
- profiling of the astrophysical normalization and spectral index;
- no bin-by-bin atmospheric background, making the ARCA significance an optimistic signal-only sensitivity.

## ARCA 10-year result

For the horizon-selected track+shower likelihood:

| model | topology | Z (10 yr ARCA, optimistic) |
|---|---|---:|
| Eq. 4.3 | active redistribution proxy | 0.00491 |
| Eq. 4.3 | sterile sink | 0.01978 |
| Eq. 4.4 | active redistribution proxy | 0.00625 |
| Eq. 4.4 | sterile sink | 0.02154 |

The strongest case is therefore Eq. (4.4) with the sterile-sink completion, but even that gives only `Z = 0.0215` before atmospheric backgrounds.

The reason is physical rather than a numerical convergence failure. The horizon-selected shower sample contains only about `1.15` expected SM events in ten years for the adopted astrophysical flux. The several-percent propagated tau-flavor distortion is then diluted in the observable shower sample by electron flavor and all-flavor NC contributions, while direct PKR hadronic showers partially compensate the propagation loss. Tracks have substantially more events (`~26.1` in the selected horizon bins), but the total muon-flavor modification from tau-decay secondaries is only at the `10^-4` level.

## Generic Gen2-like acceptance scaling

No Gen2 detector response is fabricated. Scaling the same public-response likelihood by acceptance/exposure only gives, for the strongest Eq. (4.4) sterile-sink case:

- `1 x` ARCA acceptance: `Z = 0.02154`;
- `5 x`: `Z = 0.04817`;
- `10 x`: `Z = 0.06812`.

Statistics alone would require approximately

- `8.62 x 10^3` times the ARCA exposure/acceptance for `2 sigma`;
- `1.94 x 10^4` times for `3 sigma`.

The active-redistribution cases require roughly `10^5`-scale multipliers for `2 sigma`.

## Final status

The correct status is

`Absorption-horizon propagation mechanism: PASSED.`

`Conditional sterile-sink coupled-transport screen: PASSED.`

`IceCube HESE sensitivity to these benchmarks: CLOSED / insensitive.`

`ARCA 10-year detector sensitivity to these benchmarks: CLOSED, already at the optimistic signal-only level.`

`Generic 5x-10x Gen2-like acceptance scaling: CLOSED.`

Thus the several-percent propagation distortion is real but does not translate into an observable event-level signal for the present constrained PKR benchmarks because the useful absorption-horizon tau sample is too sparse and the inclusive detector channels dilute the effect.

## Important limitation of the conclusion

This is a no-go result for the **specific constrained PKR Eq. (4.3)/(4.4) branch and its active-redistribution/sterile-sink topology choices**. It must not be quoted as a model-independent no-go theorem for all momentum-growing, resonant, gluonic, leptoquark, or other BSM neutrino interactions.

The remaining model-side caveats (exact mediator widths/branching ratios, dedicated monojet/dijet recasts, and an operator-specific IceCube constraint recast) cannot rescue detectability of these benchmark points because the detector sensitivity is already orders of magnitude below one sigma in the optimistic signal-only calculation. They remain relevant only if the PKR branch is retained for completeness or theoretical comparison.

## Reproducibility

Primary result commit: `852d007`.

Key files:

- `results/v0.44/arca_horizon_summary.csv`
- `results/v0.44/arca_horizon_bins.csv`
- `results/v0.44/detector_gate.csv`
- `results/v0.44/gen2_acceptance_scaling.csv`
- `results/v0.44/metadata.json`
- `scripts/pkr_arca_gen2_horizon_v0_44.py`

Earlier pure-sink vertical suppression results remain useful only as mechanism diagnostics and must not be presented as detector sensitivity.
