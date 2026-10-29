# v0.28 — Collider gate for the active gluon d=8 operator

## Goal

Constrain

\[
Q_{l^2G^2D}^{33}=i(\bar L_3\gamma^\mu\overleftrightarrow D^\nu L_3)
G^A_{\mu\rho}G^{A\rho}{}_{\nu}
\]

using the crossed charged-tau process

\[
gg\to\tau^+\tau^-.
\]

The quantity needed by the UHE analysis is the collider-allowed value of `|C_G|/Lambda^4`. Since the UHE pure term scales as `C_G^2`, any collider bound tighter than the strict unitarity ceiling `|C_G|=8pi` reduces the UHE effect quadratically.

## HighPT suitability

HighPT is designed for high-pT semileptonic Drell–Yan observables with quark-initiated semileptonic operators. The present operator is a two-gluon/two-lepton contact interaction and is not in that class. Therefore **do not assume HighPT can perform this recast** without a dedicated extension.

## Correct public ATLAS input

Use ATLAS arXiv:2503.19836 / JHEP 10 (2025) 054, based on 140 fb^-1 at 13 TeV.

The paper provides an unfolded differential measurement of the high-mass hadronic-tau-pair visible invariant-mass distribution. HEPData record: **155470**. The public record includes covariance information; Table 2 is the combined covariance matrix and Table 4 the systematic covariance matrix.

A validated Rivet implementation exists:

`ATLAS_2025_I2905252`.

Its particle-level selection is sufficiently explicit for a model-independent recast:

- hadronic tau decays;
- visible tau acceptance `|eta|<2.47`, excluding `1.37<|eta|<1.52`;
- exactly two accepted visible hadronic taus;
- leading visible tau `pT > 90 GeV`;
- subleading visible tau `pT > 60 GeV`;
- visible ditau mass `m_vis > 100 GeV`.

Because the measurement is unfolded to particle level, a full ATLAS detector simulation is not required for this cross-check. Generate the BSM signal, decay/shower it, run the validated Rivet analysis, and compare the resulting fiducial spectrum with the unfolded measurement and covariance matrix.

## Preferred recast route

1. Implement the `gg tau tau` contact interaction in a minimal UFO/FeynRules model using the same canonical Murphy normalization as the UHE fold.
2. Validate the UFO matrix element analytically against the crossed amplitude and its `C_G^2 shat^3/Lambda^8` scaling.
3. Generate `pp -> tau+ tau-` through the BSM `gg` channel at 13 TeV for a reference point, preferably `Lambda=4 TeV`, `C_G=-8pi`.
4. Shower and decay taus with Pythia/Herwig while retaining spin correlations as supported by the generator chain.
5. Run `ATLAS_2025_I2905252` in Rivet.
6. Build the binned likelihood from the HEPData unfolded spectrum plus covariance, or use Contur, which already includes this validated measurement and SM predictions.
7. Extract the 95% CL bound on `|C_G|/Lambda^4`.
8. Rescale/re-run the UHE CT18NLO fold at that allowed coefficient.
9. Run Earth propagation only if the collider-allowed UHE effect still reaches the several-percent target.

## Why the unfolded measurement is preferable to copying the published EFT fit

The published ATLAS new-physics interpretations use categories including b-jet multiplicity and target operator/model classes different from this gluon-initiated contact term. Reusing those published Wilson-coefficient bounds directly would therefore be unsafe.

The unfolded inclusive `m_vis` measurement is more appropriate for this new production mechanism because it is detector-corrected and accompanied by a validated Rivet implementation and covariance information.

## Contur option

Contur is explicitly designed to constrain new BSM theories with unfolded particle-level collider measurements and can ingest BSM Monte Carlo passed through Rivet. `ATLAS_2025_I2905252` is already in the Contur measurement set with SM theory predictions available. Once the custom interaction model produces events, this provides a cleaner route than hand-building an approximate detector-level likelihood.

## Gate decision before running

The exact UHE result at the strict unitarity ceiling gives only one robust-looking total-rate point above 5%: `xi=1`, 300 PeV, with `Delta sigma_tot/sigma_tot = 7.32%`. Therefore even a modest collider reduction of `|C_G|` can kill the contact branch. Since the UHE pure term scales as `C_G^2`, retaining a 5% total-rate effect at that point requires approximately

\[
|C_G| \gtrsim 8\pi\sqrt{5/7.32}\approx0.827\,(8\pi)\approx20.8.
\]

Thus the collider recast has a very sharp decision criterion:

\[
\boxed{|C_G|_{95} \gtrsim 20.8\quad\text{at}\quad\Lambda=4\,\mathrm{TeV}}
\]

is needed even to preserve the 5% total-rate proxy at 300 PeV for `xi=1`.

If the ATLAS/Contur bound is below this, close the contact gluon branch without Earth propagation and proceed to the resolved spin-2 mediator.
