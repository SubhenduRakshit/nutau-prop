# PKR collider gate v0.52 — MG5 normalization certification

## Scope

Generator-level certification only.  No detector or collider exclusion claim is allowed at this stage.

The certification keeps `M=50 GeV`, `g_q=0.03`, `m_X=1 MeV`, PDF `CT18NLO` (`lhaid=14400`) and `pT(j)>100 GeV` fixed, and generates three fresh `pp -> Xd Xd~ j` samples with target invisible branching ratios

- `BR_inv = 0.25`
- `BR_inv = 0.45`
- `BR_inv = 0.60`

The visible width is reconstructed from the common v0.45 `Gamma_total-Gamma_nu` value.  The invisible width and `gVXd` are then chosen to realize each target branching ratio while remaining in the narrow-width regime.

## Mandatory checks

A PASS requires all of the following simultaneously:

1. **Persisted parameter card:** the actual `Cards/param_card.dat` contains the intended `MY1`, `MXd`, all six diagonal vector quark couplings, `gVXd`, and `WY1`.
2. **No sample reuse:** the three parameter-card hashes and, more importantly, the event-payload hashes (the LHE content after `</init>`) are all distinct.
3. **Cross-section source consistency:** MG5's persisted banner `Integrated weight (pb)` agrees with an independent sum of the subprocess cross sections in the LHE `<init>` block to relative tolerance `1e-6`.
4. **NWA scaling:** `sigma/BR_inv` is constant to within 2% across the three points.
5. **Monotonicity:** the generated cross section increases with the target invisible branching ratio.

The script writes `cache/v0.52_mg5_norm/NORMALIZATION_CERTIFIED.json` only on PASS.  Both PASS and FAIL result tables/logs are preserved under `results/v0.52_mg5_norm/`.

## Next step

Only after this gate passes should exactly one CheckMATE benchmark be run.  A full mass scan is forbidden until that single detector benchmark is internally consistent.

v0.51 remains a preserved failed-normalization test; its `r ~ 14--27` values have no physical interpretation.
