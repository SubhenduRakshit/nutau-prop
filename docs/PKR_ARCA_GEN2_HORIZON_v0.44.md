# PKR absorption-horizon detector gate — v0.44

## Status entering v0.44

- IceCube HESE public-data recast (v0.43b): insensitive at the present benchmark.
- Coupled transport (v0.42d): local reconstructed distortions survive at the few-percent level in the absorption-horizon region.
- The remaining question is whether an analysis deliberately concentrated on the moving Earth-absorption horizon has enough event statistics.

## Main correction relative to the v0.42d ARCA event-density proxy

The published ARCA upgoing effective area already contains SM Earth absorption. Therefore one must **not** multiply it by an already propagated SM flux. v0.44 uses

`incident astrophysical flux × published SM effective area × (BSM/SM arriving-flux ratio)`

for the BSM hypothesis. This avoids double-counting Earth attenuation.

## Analysis definition

- Models: PKR Eq. (4.3) and Eq. (4.4).
- Topologies: active-redistribution proxy and explicit sterile sink.
- Transport input: converged `results/v0.42d/coupled_true_energy_map.csv`.
- Zenith region: `|cos(theta)| = 0.25–0.55`, subdivided into `Delta |cos(theta)| = 0.05` bins.
- Reconstructed-energy bands: 0.1–0.3, 0.3–1, 1–3, 3–10, 10–30 and 30–100 PeV.
- Energy migration: Gaussian in `log10(E)` with width 0.30 dex.
- Keep only bins with SM-event-weighted `1 <= tau_SM,eff <= 3`.
- Channels: ARCA showers and tracks.
- Secondary `nu_mu`: included through the change in the **total** primary + tau-decay-secondary muon-flavour flux.
- Direct PKR detector interaction: included as an NC-like hadronic-shower response proxy using `sigma_NP/sigma_NC` from the v0.42d tables.
- Nuisance profiling: common astrophysical normalization and spectral index.

## Public-response limitation

The public ARCA Fig. 7 response is coarser than the desired horizon binning. The transport is evaluated/interpolated in fine zenith bins, while the shower effective area is hemisphere-averaged and the track effective area is held constant inside its published parent zenith bin. Therefore the ARCA result remains a **public-response forecast**, not an official KM3NeT likelihood.

No bin-by-bin atmospheric background model is publicly implemented here. Consequently the ARCA significance is deliberately interpreted as a **signal-only optimistic upper bound**:

- `Z_ARCA < 2`: ARCA route is robustly closed at this benchmark.
- `Z_ARCA >= 2`: promising only; exact backgrounds/response are mandatory before a positive sensitivity claim.

## Gen2 treatment

No IceCube-Gen2 detector response is fabricated. v0.44 reports generic `x5` and `x10` acceptance/exposure scalings of the optimized horizon analysis and the acceptance multipliers required for 2 sigma and 3 sigma. These are scaling forecasts only.

The public IceCube-Gen2 material states that the optical array is expected to increase instrumented volume to about 8 km^3 and yield roughly 5–10 times as many neutrinos as current IceCube above 30 TeV; this motivates the scaling bracket but does not turn it into a detector likelihood.

## Outputs

- `results/v0.44/arca_horizon_summary.csv`
- `results/v0.44/arca_horizon_bins.csv`
- `results/v0.44/gen2_acceptance_scaling.csv`
- `results/v0.44/detector_gate.csv`
- `results/v0.44/metadata.json`
- `results/v0.44/run.log`

## Runner

```bash
bash run_pkr_arca_gen2_horizon_v0_44.sh
```

The runner commits and pushes the result files and then returns to the shell prompt.
