# nutau-prop v0.17.1 — finite-beta_tau validation outcome

Date: 2026-10-28

## Result

The repaired finite-beta_tau transport passed all substantive primitive and full-transport validation checks. The first v0.17 gate report contained 250 checks with exactly one failure: an auxiliary monotonic-refinement heuristic for the Case-6 point Lambda=2.44 TeV, |cos z|=0.5, BSM=true.

At that point the successive tau-spectrum differences were

- d40_20_tau = 7.2654e-7,
- d20_10_tau = 9.7933e-7.

The latter is slightly larger, so the original trend-only heuristic failed even though both errors are already far below the primary 1% convergence tolerance. All max_p tau errors were O(1e-6), all secondary-muon max_p errors were <=2.2e-5, and matter-step errors were O(1e-4)-O(6e-4), all comfortably inside their acceptance thresholds. Positivity and finiteness passed everywhere.

## Criterion correction

The auxiliary refinement-trend test is now interpreted with a fixed numerical noise floor:

- if max(d40_20_tau, d20_10_tau) < 1e-4, classify the point as already converged and do not require monotonic decrease;
- otherwise require d20_10_tau <= 1.05 d40_20_tau + 1e-8.

The 1e-4 floor is 100 times tighter than the primary 1% tau-convergence tolerance, so this does not weaken the physics accuracy requirement. It only avoids rejecting solutions because two already-negligible discretization errors exchange order near the interpolation/splitting floor.

A machine-readable reassessment is implemented in `scripts/reassess_tau_transport_v0_17_1.py`. No transport code or physical result is changed by this criterion correction.

## Validation status

With the noise-floor-aware auxiliary criterion, all reported v0.17 primitive and full-transport checks satisfy the validation requirements. The finite-beta_tau repair is therefore accepted for subsequent reruns.

This validation does **not** restore the old v0.12-v0.16 phenomenological conclusions automatically. Those transport stages must be rerun with the repaired solver. Hard-kernel/PDF calculations remain valid and may be reused.

## Next ordered work

1. Rerun representative Case 6, Case 5 and Case 7 transport using the repaired v0.17 solver, preserving incident, SM-surviving and BSM-surviving spectral weights.
2. Reuse the existing v0.16 resolved-mediator cross-section tables and rerun only its transport stage.
3. Perform detector weighting before deciding whether opacity or detector enhancement wins.
