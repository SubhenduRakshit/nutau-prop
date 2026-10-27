# nutau-prop v0.13 — end-to-end Case-6 propagation + detector fold

This stage closes the remaining Case-6 consistency gap identified after v0.12.

## What is now combined

- the exact v0.11 Case-6 `d sigma/dy` detector vertex;
- the v0.12 PREM Earth propagation with SM CC, SM+Case6 NC energy degradation, tau regeneration, and SM-only secondary muon flavor;
- ARCA shower topology (`e CC`, `tau CC`, all-flavor `NC`) and zenith-separated track responses;
- IceCube EHE public flavor effective areas and released differential-limit curve;
- the KM3NeT/KM3-230213A public bright-track effective area.

## Secondary nu_mu treatment

Case 6 has no direct muon-track hard vertex.  The only new track contribution comes from the change in tau-decay secondary `nu_mu + antinu_mu` generated during Earth propagation.  The code adds only the **BSM-SM secondary difference** to published SM detector responses.  This avoids double counting a standard secondary contribution that may already be embedded in collaboration simulation.

For ARCA, the published zenith track curves allow a spectrum-specific correction for the SM Earth propagation already contained in the effective area.  For IceCube/KM3NeT the public UHE products are sky averaged, so the result is reported as a transparent public-response approximation/envelope, not an exact collaboration likelihood.

## Changed-y detector systematic

Published effective areas have already integrated over the SM inelasticity distribution.  Exact selection efficiency for the Case-6 `y` distribution therefore cannot be reconstructed from public data.  ARCA results are reported both as

1. `propagation_only`: modified Earth flux with the SM detector vertex;
2. `joint`: modified Earth flux plus the Case-6 NC detector cross section.

The difference is retained as a public-response systematic rather than hidden.

IceCube reports `central`, `low`, `high`, and `propagation_only` response folds.  The `central` result assumes equal up/down solid-angle weight for the sky-averaged flavor area; `low/high` span the available propagated zenith factors and CC/NC detector-vertex ambiguity.

## IceCube constraint

The v0.9 released sliding-window differential-limit screen is recomputed using the **joint propagated response**, not the obsolete `A_eff*(1+r)` detector-only rescaling.  The output therefore directly tests whether the earlier IceCube-vs-KM3NeT tension survives once Case-6 Earth opacity is included.

The script also reports the fraction of the SM Ehlert IceCube expectation lying above the v0.12 transport grid (20 EeV).  The QA gate requires this omitted fraction to be below 5%; otherwise the Case-6 transport grid must be extended before interpreting the public-limit result.

## Run after v0.12

The v0.13 package must be placed on top of a directory containing the generated v0.11/v0.12 results, or those stages must be rerun first.

```bash
python -m pip install -e .
bash run_case6_joint_v0_13.sh
```

Outputs:

- `results/v0.13/case6_arca_joint_v0_13.csv`
- `results/v0.13/case6_arca_secondary_track_bins_v0_13.csv`
- `results/v0.13/case6_icecube_joint_v0_13.csv`
- `results/v0.13/case6_icecube_public_windows_joint_v0_13.csv`
- `results/v0.13/case6_km3_joint_v0_13.csv`
- `results/v0.13/case6_joint_status_v0_13.json`
- `results/v0.13/case6_joint_gate_v0_13.json`
- `results/v0.13/run_case6_joint_v0_13_summary.txt`

A successful QA gate validates the numerical bookkeeping of this public-response theory fold.  It does **not** claim access to IceCube's unreleased reconstructed-bin likelihood or to collaboration detector MC for a changed inelasticity distribution.
