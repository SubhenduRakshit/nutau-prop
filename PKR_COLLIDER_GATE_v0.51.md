# PKR collider gate v0.51

## Purpose

Test the conditional PKR benchmark

- `M_Zp = 50 GeV`
- `g_q = 0.03`
- Eq. (4.3) dipole and Eq. (4.4) derivative couplings fixed by the converged v0.42d 5%-at-100-PeV benchmark

against the current CMS jets+MET / monojet implementation in CheckMATE2.

The analysis identifier is `cms_2107_13021`, corresponding to CMS EXO-20-004. Current CheckMATE2 provides a multibin likelihood for this analysis.

## Why a width-matched vector invisible state is sufficient for this gate

The v0.45 benchmarks have `Gamma_Zp/M_Zp < 9e-4`. In the narrow-width limit the visible recoil process factorizes,

`pp -> Z' + j -> invisible + j`,

and the measured jet/MET kinematics are controlled by the vector quark production current and `BR(Z'->invisible)`. The internal Lorentz structure of an entirely invisible on-shell decay is not measured.

Therefore v0.51 uses the standard `DMsimp_s_spin1` UFO with the same vector quark coupling and a nearly massless Dirac invisible state. Its vector coupling is chosen so that

`Gamma(Y1 -> Xd Xdbar) = Gamma_invisible^PKR`,

while `WY1` is explicitly fixed to the v0.45 PKR total width. This matches production times invisible branching ratio for each of the four width conventions.

This is **not** claimed to validate off-shell Eq. (4.3)/(4.4) amplitudes. If the CheckMATE result is close to the exclusion boundary, a dedicated charge-conjugate derivative/dipole UFO is required.

## Event chain

1. MG5_aMC LO: `p p > xd xd~ j`, 13 TeV, CT18NLO, `pT(j)>100 GeV`.
2. Pythia8 ISR/FSR/MPI/hadronization.
3. HepMC2 event record.
4. CheckMATE2 detector simulation and `cms_2107_13021` analysis.
5. Primary result: observed multibin signal-strength upper limit `mu95`; define `r = 1/mu95`.

The default is 100,000 generated events per benchmark and a 10% signal normalization uncertainty.

## Gate

- `r < 0.8`: robust CheckMATE survival for this collider screen.
- `r > 1.2`: robust CheckMATE exclusion.
- `0.8 <= r <= 1.2`: marginal; require a dedicated operator-level UFO / higher-fidelity collider recast.

The ordinary CheckMATE single-bin `robscons` is also recorded as a cross-check.

## Audit status of earlier attempts

v0.47--v0.50 are retained as failed/deprecated attempts to reconstruct or reproduce the CMS likelihood directly. No PKR exclusion or survival conclusion is taken from those versions.
