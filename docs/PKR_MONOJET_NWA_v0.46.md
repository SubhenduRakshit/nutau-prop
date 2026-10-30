# PKR monojet narrow-width public-limit screen — v0.46

## Why this step

v0.45 shows the 50 GeV mediator is extremely narrow and has a substantial invisible branching fraction. Therefore, before building a full detector recast, test the benchmark against the public spin-1 monojet coupling limits in the narrow-width approximation (NWA).

For an on-shell narrow spin-1 mediator produced from quarks,

`N_monojet ∝ sigma(pp -> Z' + j) BR_inv ∝ g_q^2 BR_inv`,

up to the mass-dependent production/acceptance factor. The invisible derivative/dipole operator mainly fixes the decay width/BR; the hard recoil is controlled by production of the on-shell Z'.

## Public CMS anchor

CMS EXO-20-004 / JHEP 11 (2021) 153 uses 137 fb^-1 and derives coupling limits by rescaling the production cross section and mediator branching fractions. For a spin-1 mediator with `m_DM=m_med/3`, the published observed reach includes `g_q=0.018` at `m_med=100 GeV` with `g_chi=1`.

v0.46 computes the corresponding standard-vector reference `BR_inv` and maps each v0.45 PKR width convention onto the equivalent monojet yield.

The decisive output is not an invented 50-GeV limit. It is

`critical_gq95_at_50GeV = gq_PKR sqrt(BRinv_PKR / BRinv_CMSref)`.

If the actual observed CMS (or validated recast) 50-GeV limit lies below this number, the PKR benchmark is excluded in NWA.

## Why the width convention remains explicit

v0.45 retained both `one_Weyl_pair` and `operator_plus_hc` invisible-width conventions. Their monojet thresholds differ enough to matter near the public coupling sensitivity. v0.46 therefore carries both rather than silently choosing one.

## What follows

If the entire convention bracket is clearly above/below an exact low-mass public limit, the collider gate can be decided without expensive simulation. If it is marginal, the next justified step is the validated CMS EXO-20-004 MadAnalysis/simplified-likelihood recast using generated `pp -> Z' + j`, followed by the PKR invisible decay.
