# Concrete input contract for an Eq. 4.3/4.4 IceCube recast

This is a saved technical specification, not a message sent to anyone.

## Signal definition

The current explicit branch is nu_tau q -> N q and its conjugate/inverse,
with a light stable N and the normalized derivative/tensor kernels of v0.86.
The final-state hadronic system deposits energy; N is invisible but its
subsequent transition can regenerate an active neutrino. An active-flavour
LNV realization must be treated as a separate flavour-transport hypothesis.

## Either of two sufficient response products

1. Selected MC with true incident energy/direction/species, CC/NC label,
   Bjorken x and y (or an adequately validated marginalized y response),
   nominal differential-cross-section/generator weights, true hadronic
   energy, reconstructed bin observables and their detector-systematic weights.
   Include definitions of whether primary truth is before or after Earth
   propagation and how SM tau-regenerated events are weighted.
2. A validated response R_i(E,y,cos(theta),channel;eta), with generator and
   exposure normalization, applicable to an otherwise SM-like hadronic
   shower with an invisible outgoing fermion. Supply any x-dependent effects
   or evidence that marginalizing them is adequate at the stated precision.

HESE observations, likelihood formula and detector splines are already available.
True y is not available in the inspected JSON; substituting reconstructed
deposited energy/E is not an exact recovery of it.

For the EHE likelihood additionally provide the reconstructed-bin observations,
background templates/correlations and signal response used in the12.6-year
fit, or an equivalent validated public likelihood interface.

## Acceptance tests once supplied

- Reproduce SM bin yields, selected flavour/topology totals, nuisance best fit
  and at least one published benchmark likelihood/limit.
- At zero coupling recover that result; at equal amplitude/common cutoff
  Eq.4.3 and4.4 produce equal signal templates.
- Include the new NC-like hadronic shower term, not only an absorption weight.
- Solve active/sterile forward and inverse transitions and tau regeneration;
  verify conservation with SM removal switched off and energy-grid convergence.
- Profile the documented nuisance parameters and calibrate the stated test
  statistic; do not relabel a fixed-template screening ratio as a new CL.
