# PKR Eq.4.4 v0.80 seed audit — v0.81

Date: 2026-10-31

## Finding

The v0.80 same-HepMC reweight was intended to hold the detector-level shape fixed while changing only the signal normalization from the v0.75 M=50 GeV sample to the v0.78 representative coupling factor.

However, the CheckMATE card in `scripts/pkr_eq4p4_M50_reweight_checkmate_v0_80.py` used

- v0.75: `RandomSeed: 7501`
- v0.80: `RandomSeed: 8001`

CheckMATE reruns Delphes from the HepMC input, so changing the random seed changes the detector stochastic realization and can change the selected-bin shape. Therefore v0.80 was not a pure same-shape normalization reweight.

The observed change

- v0.75: `r_MB = 1.1380325`
- v0.80: `r_MB = 0.7419261`

for only a normalization factor `P_rep/P_v075 = 0.907533` is consequently not interpretable as evidence for robust survival.

## Canonical status

`v0.80 r=0.7419` is preserved as an invalid/diagnostic result and **must not be used as a collider-survival certificate**.

The corrected test must reuse the v0.75 detector seed 7501 and verify that the detector-bin vector is identical up to the common normalization factor. Only then may the rescaled multibin result be interpreted.

## Expected scaling check

For a truly fixed detector-bin shape and fixed relative signal uncertainty, the signal-strength limit should approximately scale with the overall normalization. Using the v0.75 value and the v0.78 normalization factor gives the rough expectation

`r_new ~ 1.1380325 * 0.907533 ~ 1.033`,

which lies in the project's `0.8 < r < 1.2` boundary-unresolved band.

This expectation is a diagnostic only; the corrected CheckMATE/Spey rerun decides the recorded value.
