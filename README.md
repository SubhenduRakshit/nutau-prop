# nutau-prop v0.6

Reproducible **tau-neutrino Earth transport** for the project

> forward-blind / charge-radius / derivative neutral-current interactions  
> → modified \(\nu_\tau N\) DIS  
> → PREM Earth propagation  
> → \(\nu_\tau \to \tau \to \nu_\tau\) regeneration.

This package is deliberately transparent. The earlier chat-level percentage tables are **not** embedded as reference answers. The code must generate the result.

## What is implemented

- LHAPDF-based LO \(\nu N\) CC and NC double-differential DIS.
- Isoscalar target constructed from proton PDFs by isospin.
- Exact SM/BSM interference at the **partonic LO amplitude/coupling** level.
- Diagonal \(\nu_\tau\) charge-radius benchmark:
  \[
  \delta g_V^q=-\frac{2}{3}Q_q M_W^2\sin^2\theta_W\langle r_{\nu_\tau}^2\rangle .
  \]
- Forward-blind derivative benchmark:
  \[
  \delta g_V^q(Q^2)=Q_q\,\epsilon_{\rm em}^{\infty}
       \frac{Q^2}{Q^2+m_X^2}.
  \]
- Full ten-layer PREM density parametrization.
- Coupled \(\nu_\tau\)-\(\tau\) transport:
  - CC tau production;
  - NC energy regeneration;
  - continuous tau energy loss;
  - tau decay back to \(\nu_\tau\).
- Rakshit-Reya-style SM benchmark output.
- convergence check and a machine-readable publishability gate.

## Important precision statement

`v0.1` is a **research-grade LO baseline**, not yet an NNLO neutrino DIS library.

The modern LHAPDF set is used directly, but coefficient functions are LO.
For a final paper:
1. validate total CC/NC cross sections against modern CT18/CSMS predictions;
2. quantify LO→NLO/NNLO effects;
3. if needed apply validated K-factors or replace the structure functions with a modern NLO implementation;
4. add nuclear PDFs / nuclear corrections.

The key BSM observable is primarily a **ratio** to the SM computed with the same PDF/perturbative setup, so much of the common normalization uncertainty cancels.

## Install on the Mac Studio

Copy/unpack this folder as

```bash
/Users/subhendurakshit/Research/nutau-prop
```

Then

```bash
cd /Users/subhendurakshit/Research/nutau-prop
chmod +x setup_macos.sh run_all.sh
./setup_macos.sh
conda activate nutau-prop
python -m pip install -e .
lhapdf install CT18NLO
```

Check:

```bash
python scripts/check_environment.py
```

## One-command run

```bash
cd /Users/subhendurakshit/Research/nutau-prop
conda activate nutau-prop
./run_all.sh
```

The first cross-section computation is the expensive step. Cached tables are written under `cache/`.

## Faster first run

Edit `config/default.yaml`:

```yaml
energy_grid:
  n: 51
integration:
  nx: 48
  y_log_bins: 40
  y_linear_bins: 40
```

For final numbers restore the defaults and then increase them for convergence.

## Outputs

- `results/sm_xsec_validation.txt`
- `results/rr_benchmark.csv`
- `results/rr_benchmark.png`
- `results/bsm_scan.csv`
- `results/gate.json`
- `results/gate.txt`

## Benchmark logic

The final gate is deliberately mechanical. In bins satisfying the configured energy and minimum-SM-transmission cuts, a benchmark passes if

\[
\left|\frac{\Phi_{\rm BSM}-\Phi_{\rm SM}}{\Phi_{\rm SM}}\right|\ge 0.05.
\]

Change this in `config/default.yaml` if desired.

## Physics conventions

### Neutral current

At LO,

\[
\frac{d^2\sigma_\nu^{NC}}{dxdy}
=\frac{2G_F^2m_NE_\nu}{\pi}
\frac{M_Z^4}{(Q^2+M_Z^2)^2}
\sum_q \left[
g_{qL}^2(xq+x\bar q(1-y)^2)
+g_{qR}^2(xq(1-y)^2+x\bar q)
\right].
\]

A charge radius shifts the quark vector coupling. Since

\[
g_V=g_L+g_R,\qquad g_A=g_L-g_R,
\]

the code implements

\[
g_L\to g_L+\frac{\delta g_V}{2},\qquad
g_R\to g_R+\frac{\delta g_V}{2}.
\]

### Antineutrino sign

The 2025 FPF charge-radius treatment notes the opposite relative charge-radius sign for antineutrinos. The setting is explicit in YAML (`antinu_charge_radius_sign`) rather than hidden in the code.

### Tau decay

The baseline uses an analytic polarized inclusive high-energy kernel. A final precision paper should compare it to TAUOLA/PROPOSAL or the exact decay treatment used in a modern tau-transport code.

## Parallelism

The v0.1 package caches one benchmark at a time and therefore avoids large RAM use. The expensive benchmark jobs can safely be launched independently.

On the Mac Studio, for example:

```bash
python scripts/precompute_xsec.py --out cache/cr_p6e-32.npz --kind charge_radius --r2 6e-32 &
python scripts/precompute_xsec.py --out cache/cr_m6e-32.npz --kind charge_radius --r2 -6e-32 &
wait
```

Do not oversubscribe LHAPDF jobs. Four simultaneous processes is a reasonable initial choice on a 36–64 GB Apple-Silicon machine.

## Primary references

- S. Rakshit and E. Reya, *On the transport equations of cosmic neutrinos passing through Earth and secondary nu_mu fluxes*, Phys. Rev. D 74, 103006 (2006), arXiv:hep-ph/0608054.
- K. Xie et al., *High-energy neutrino deep inelastic scattering cross sections*, Phys. Rev. D 109, 113001 (2024), arXiv:2303.13607.
- P. Weigel et al., *Cross sections and inelasticity distributions of high-energy neutrino deep inelastic scattering*, Phys. Rev. D 111, 043044 (2025).
- R. M. Abraham et al., *Neutrino electromagnetic properties and the weak mixing angle at the LHC Forward Physics Facility*, Phys. Rev. D 111, 015029 (2025).

See `VALIDATION.md` before using any number in a manuscript.


## v0.2: mandatory historical validation

Run the complete gated pipeline:

```bash
export NUTAU_NPROC=8
./run_all.sh
```

Run only the historical gate:

```bash
export NUTAU_NPROC=8
./run_historical_validation.sh
```

All terminal output is appended automatically to `PROGRESS.log`.
The BSM stage is blocked unless both the Rakshit-Reya and
Dutta-Reno-Sarcevic gates pass.

See `VALIDATION.md` before interpreting any BSM result.


## v0.3 correction after the first real RR run

v0.2 was rejected when the user's Mac run found non-convergent historical
iterations.  v0.3 fixes the historical physics rather than weakening the gate:
full inclusive tau decay, separate nu/antinu transport, and an RR-style
decay-aware Z iteration.

For the fastest next test, run only:

```bash
export NUTAU_NPROC=8
./run_historical_validation.sh
```

Only if that passes should `./run_all.sh` be used for the BSM stage.


## v0.4: IceCube event-weighted observability

v0.4 adds a detector-level fold using the official IceCube 2015 contained-event
tau-neutrino effective area in true energy and zenith.

The effective area already includes **SM Earth absorption**.  Therefore the BSM
Earth effect is applied only as `T_BSM/T_SM`; multiplying by the attenuated
flux again would double-count absorption.

Run only the new stage after the v0.3 historical gate has passed:

```bash
export NUTAU_NPROC=8
./run_detector_gate.sh
```

The script downloads the official 2.3 MB IceCube public data release if needed.

Main outputs:

- `results/bsm_scan_combined.csv`
- `results/event_weighted_counts.csv`
- `results/event_weighted_summary.csv`
- `results/event_weighted_summary.json`
- `results/event_weighted_best_model.png`

Use `Z_profile_norm_tilt` as the conservative first observability diagnostic.
The detector cross-section reweight is still an inclusive proxy because the
total public effective area does not expose the exact CC/NC acceptance mixture.


## v0.5: larger allowed couplings + future-detector gate

This release answers the final publishability question.

The externally bounded charge-radius scan is extended to

```text
±6, ±10, ±15, ±20 × 10^-32 cm^2
```

where `±20 × 10^-32 cm^2` is the 2025 global-fit 1σ tau-neutrino
charge-radius interval.

The derivative forward-blind scan is also extended, but those points are
**exploratory** and are not called globally allowed.

Future sensitivity is projected in three ways:

1. IceCube-Gen2, 10 years, response ×5;
2. optimistic Gen2-like yield, 10 years, response ×10;
3. generic future tau-neutrino telescope acceptance multipliers up to ×1000.

Run:

```bash
export NUTAU_NPROC=8
./run_future_gate.sh
```

Main outputs:

- `results/future/future_sensitivity_summary.csv`
- `results/future/future_sensitivity_energy_bands.csv`
- `results/future/generic_tau_telescope_scaling.csv`
- `results/future/future_publishability_gate.json`

The code does **not** invent a Gen2 effective-area table.  The ×5 and ×10
scenarios are explicit scalings of the already validated IceCube response,
based on public Gen2 design performance statements.


## v0.6: KM3NeT/ARCA track-control + shower test

v0.6 tests whether full KM3NeT/ARCA can break the astrophysical-flux nuisance
degeneracy that suppresses the IceCube and simple Gen2 projections.

The key new idea is to fit simultaneously:

- the tau-sensitive **shower** sample; and
- the largely BSM-insensitive, high-statistics **nu_mu CC track** sample as a
  common flux-normalisation/spectral-index control.

The official 2024 ARCA simulation predicts 220 cosmic-neutrino tracks/year and
60 showers/year for its reference E^-2 flux.  Its associated datasets and code
are public at Zenodo record 10949703.

Run:

```bash
export NUTAU_NPROC=8
./run_arca_gate.sh
```

The script automatically tries to download the official Zenodo record and saves
its inventory.  The current likelihood uses a coarse published Fig.7 response
digitisation calibrated exactly to the published track/shower rates; it is
therefore explicitly labelled a **signal-only optimistic upper bound**.

Main result:

- `results/arca/arca_publishability_gate.json`

If the track+shower profiled significance is below 2 sigma even in this
optimistic calculation, ARCA cannot rescue the benchmark.

## v0.8 additions
See `README_v0.8.md` and run `./run_v0_8.sh`. The original v0.6 validation/ARCA machinery is retained.
