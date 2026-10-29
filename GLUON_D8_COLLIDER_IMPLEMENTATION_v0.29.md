# v0.29 — Exact-operator collider implementation

## Literature status

A targeted search found earlier dimension-8 tau–gluon studies, but their operators are scalar/chirality-flipping structures such as `G^2 (bar L tau_R) H` (or CP-odd analogues). They are not the derivative twist-2 SMEFT direction

\[
Q_{l^2G^2D}^{33}=i(\bar L_3\gamma^\mu\overleftrightarrow D^\nu L_3)
G^A_{\mu\rho}G^{A\rho}{}_{\nu}.
\]

No public/published recast of this exact operator against the 2025 ATLAS high-mass ditau spectrum was identified. The exact direction is therefore implemented here rather than borrowing a scalar-operator limit.

## Two independent implementations

### A. FeynRules/UFO route

`models/gluon_d8/GLUON_D8.fr` implements the two-gluon projection of the canonical Murphy operator for tau and tau-neutrino fields. `models/gluon_d8/export_ufo_v0_29.wls` loads the standard FeynRules SM, checks hermiticity, requires a nonzero `gg tau tau` vertex, and exports `SM_GLUON_D8_UFO`.

For the collider process with exactly two external gluons, replacing the lepton covariant derivative by the ordinary derivative is sufficient: electroweak gauge pieces generate vertices with additional electroweak bosons and do not contribute to `gg -> tau+ tau-`. The file is therefore a targeted collider projection, not a complete all-process gauge-completion.

The UFO route is intended as an independent matrix-element validation and for a future MadGraph production sample.

### B. Analytic-LHE route

To avoid making the recast depend on proprietary Mathematica/FeynRules execution, `scripts/generate_gluon_d8_lhe_v0_29.py` generates the same hard process directly from the crossed certified matrix element.

For the single chiral final state `tau-_L tau+_R`, after averaging over both initial gluon helicities and colors,

\[
\overline{|\mathcal M|^2}
=\frac{C_G^2}{\Lambda^8}\frac{\hat s^4}{64}(1-z^4),
\qquad z=\cos\theta,
\]

and therefore

\[
\frac{d\hat\sigma}{dz}
=\frac{C_G^2\hat s^3}{2048\pi\Lambda^8}(1-z^4),
\qquad
\boxed{\hat\sigma=\frac{C_G^2\hat s^3}{1280\pi\Lambda^8}}.
\]

The `1-z^4` angular form is the standard gluon-fusion spin-2 dilepton kernel. The normalization follows by crossing the already certified `nu g -> nu g` amplitude and accounting for the second initial-gluon color/polarization average. It is also consistent with writing the operator as the chiral fermion/gluon stress-tensor product (up to the fixed normalization factor).

The generator convolves the partonic rate with `NNPDF31_lo_as_0130`, imposes the same EFT-validity cut `sqrt(shat)<sqrt(xi) Lambda`, generates unweighted LHE events, and marks the tau helicities (`tau-` left, `tau+` right).

## Shower and official ATLAS analysis

`scripts/shower_gluon_d8_pythia8_v0_29.cc` reads the LHE sample, performs Pythia8 ISR/FSR/MPI/hadronisation/tau decays, and writes HepMC3.

Rivet then runs the validated official analysis

`ATLAS_2025_I2905252`

on the HepMC sample. The analysis measures the unfolded visible ditau mass spectrum with the published particle-level hadronic-tau cuts.

## Contur

The current Contur database contains `ATLAS_2025_I2905252` in the `ATLAS_13_TAU` pool and provides an SM theory prediction. Consequently no ad-hoc covariance likelihood is needed: the BSM Rivet YODA result can be passed directly to Contur.

The BSM shape is coefficient-independent and the rate scales exactly as `C_G^2` at this order. `scripts/scan_gluon_d8_contur_v0_29.sh` therefore generates one high-statistics reference YODA file at `C_G=-8pi` and rescales it by `(C/C_ref)^2` for the coefficient scan before running Contur.

The scan includes the phenomenologically decisive neighborhood around

\[
|C_G|=20.8\quad(\Lambda=4\;\mathrm{TeV}),
\]

because the UHE total-rate proxy falls below 5% for smaller coefficients.

## One-shot execution

`run_gluon_d8_collider_v0_29.sh` performs

1. analytic LHE generation;
2. Pythia8 shower/decay;
3. official Rivet `ATLAS_2025_I2905252` analysis;
4. Contur at `C_G=-8pi`;
5. coefficient-rescaled Contur scan.

The GitHub Actions workflow `.github/workflows/gluon-d8-collider-v029.yml` runs an initial 20k-event pilot and commits only compact outputs; LHE/HepMC files are removed before committing.

## Decision rule

- if the 95% collider upper limit satisfies `|C_G|_95 < 20.8` at `Lambda=4 TeV`, close the contact-gluon branch without Earth propagation;
- if `|C_G|_95 >= 20.8`, increase MC statistics around the boundary, validate against the FeynRules/UFO matrix element, then reconsider Earth propagation;
- in either case the strong endpoint sensitivity remains a motivation for the resolved spin-2 mediator calculation.

Earth propagation is not part of v0.29.
