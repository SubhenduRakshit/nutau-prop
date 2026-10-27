# Status and audit — v0.9 (2026-10-27)

## Verified inputs and decisions

- The official IceCube EHE area is in m² and summed across flavors; v0.8's
  first run omitted the factor 10⁴ and reported 0.001767 rather than about
  17.675 SM Ehlert events. Its later fix used a writable copy. v0.9 reads
  the first numeric row and confirms the flavor-column sum.
- The released KM3NeT bright-track response reproduces the published Ehlert
  expectation: 0.242754 versus 0.242767 in 335 days; central 90% 0.211081
  versus 0.211117.
- Case 6 is NC-like, so the `KM3_NC_equivalent_proxy` column is not a
  bright-track prediction. The direct new-track contribution is zero in this
  calculation. Indirect propagation effects remain uncomputed.
- The primary EFT benchmark uses √ŝ < Λ. At 3 TeV its coefficient exceeds
  the simple c₄/(4π) < 1 diagnostic; 2.44 TeV is near the boundary.

## Current result

The strongest released IceCube differential-limit window for the tested
Ehlert spectrum is centered at 316 PeV. For Λ = 1.5, 2.0, 2.44, 3.0 TeV,
the all-flavor **response-proxy** limiting common-flux scales are 0.1549,
0.1280, 0.1025, 0.0747. Applying them to the unchanged direct track fold
gives 0.0376, 0.0311, 0.0249, 0.0181 events. Tau-initial proxy choices
give somewhat weaker scales; see the CSV. These are not 90% confidence
limits and should not be presented as such.

## Blocking data for a publishable likelihood

IceCube's published EHE likelihood uses reconstructed energy and direction,
an astrophysical nuisance component, backgrounds, and pseudoexperiments.
Its public release contains only true-energy effective areas, an observed
differential-limit curve, a README, and a demo. The reconstructed-bin counts,
signal migration/response, backgrounds, and nuisance distributions needed
for a new-physics likelihood are not public in that release. A formal
coverage-correct limit and topology-aware ARCA/IceCube prediction require
those detector inputs or a collaboration-level simulation. The next
scientific checkpoint is to obtain those inputs and implement the
Case-6 flavor/charge-conjugation transport kernel consistently.

## Alternative mechanisms: analytic gates

- A muon-producing charged-current interaction could alter the track count;
  derive gauge embedding, collider limits, EFT validity, and both detectors'
  angle/energy response before scanning its parameters.
- A tau-specific NC interaction may affect secondary muons through Earth
  transport. Evaluate the coupled ντ–τ–νμ kernels with its differential
  inelasticity and neutral final-state assignment before claiming a track.
- A transient or directional flux can evade a stationary isotropic exposure
  ratio, but requires time- and direction-dependent IceCube/ARCA analyses.
- A resonance or narrow spectral feature needs an explicit UV-complete
  model; the contact EFT cannot be extrapolated through a resonance.

Failed and rejected approaches are retained in `PROGRESS.log` and Git history.
The JHEP draft still needs the author's current affiliation and correspondence
address before journal submission; neither was supplied in this task.
