# v0.18.1 response to external audit of v0.18

Date: 2026-10-28

## Verdict accepted

The audit verdict is accepted: the analytic finite-beta_tau repair is sound, but v0.18 production was not sufficiently certified for quantitative phenomenology. v0.18 remains a provisional historical artifact. Detector conclusions remain withheld.

## Findings and corrective actions

### 1. Production solver differed from validated solver — fixed

A new validation script, `scripts/validate_chirality_transport_v0_18_1.py`, tests the exact production classes:

- `TauHelicityFlipTransport`;
- `TauMuonTransitionTransport`.

It uses the Ehlert spectrum, both external-state brackets, shallow/intermediate/vertical chords, Case 5, Case 6, Case 7 and representative resolved-mediator masses. Production resolution (`max_p=0.02`, `max_step_km=12`) is compared directly with a jointly refined reference (`max_p=0.01`, `max_step_km=6`). The gate includes pointwise errors at 1, 10, 100, 220 and 300 PeV, species-resolved outputs, integrated L1 errors, positivity and finiteness. The acceptance threshold is 0.5%.

This empirical convergence check is what supports production accuracy; no claim of globally second-order Strang accuracy is made.

### 2. Production gate was incomplete — fixed

`scripts/certify_transport_v0_18_1.py` requires all of the following before production:

1. Rakshit-Reya historical gate PASS;
2. Dutta-Reno-Sarcevic historical gate PASS;
3. every substantive v0.17 gate PASS except the explicitly superseded `refinement_trend_tau` heuristic;
4. v0.17.1 refinement reassessment PASS;
5. v0.18.1 chirality-solver validation PASS.

The production runner refuses to proceed unless `production_certified=true`.

### 3. nu and antinu information was lost — fixed

`scripts/run_certified_transport_v0_18_1.py` exports neutrino and antineutrino spectra separately. Depending on the model this includes:

- `SM_nu_tau`, `SM_antinu_tau`, `BSM_nu_tau`, `BSM_antinu_tau`;
- primary or secondary `nu_mu` and `antinu_mu` separately.

This preserves the information required for a later species-dependent detector response.

### 4. Transmission outside source support and interpolation gaps — fixed

The Ehlert tabulated support is carried explicitly as `in_incident_flux_support`.

- Incident flux is set to zero outside its tabulated support.
- Transmission is undefined (`NaN`) outside source support.
- Interpretation ratios are masked unless the target energy lies inside support and the matched-SM transmission exceeds the stated floor.
- Summary interpolation uses only the two local bracketing grid points and refuses interpolation if either required value/support flag is invalid. It does not delete invalid points and bridge across gaps.

Regenerated flux below source support remains stored as a physical propagated flux, but no `flux/incident` transmission is assigned there.

### 5. Symmetric splitting does not by itself prove second-order accuracy — accepted

The audit's equal-rate exchange example is valid. The neutrino interaction map is not an exact exponential of the full coupled exchange generator. Therefore the code no longer relies on the phrase "second-order Strang accuracy" as a certification argument.

The production accuracy claim is instead empirical and model-specific, from the v0.18.1 convergence tests of the actual exchange solvers at the production resolution and UHE target energies.

### 6. Table provenance and version drift — fixed for production

Production now requires each NPZ table to have a sidecar JSON and checks the physical metadata before loading it:

- Cases 5/7: v0.15, correct operator, `shat` cutoff, correct Lambda, CT10nlo;
- Case 6: v0.11, `shat` cutoff, correct Lambda, CT10nlo;
- resolved mediator: v0.16, correct mediator mass and low-Q coefficient, CT10nlo, and `vertex_cut_TeV = null` for the primary physical-propagator scan.

A SHA256 checksum and full sidecar metadata for every consumed table are written to `results/v0.18.1/input_provenance_manifest_v0_18_1.json`.

Package metadata and `nutau_prop.__version__` are synchronized to `0.18.1` in GitHub. Existing local checkouts must fetch those files to synchronize locally.

The manuscript is intentionally not bumped to v0.18.1 before certified numerical results exist.

### 7. Duplicate SM runs and core control — fixed

Each model/parameter/zenith task now computes the identical matched-SM propagation once and reuses it for the two external-state modes.

`run_certified_transport_v0_18_1.sh` caps internal threaded libraries to one thread per worker:

- `OMP_NUM_THREADS=1`;
- `OPENBLAS_NUM_THREADS=1`;
- `MKL_NUM_THREADS=1`;
- `VECLIB_MAXIMUM_THREADS=1`;
- `NUMEXPR_NUM_THREADS=1`.

The multiprocessing worker count remains controlled by `NUTAU_NPROC` (default 12).

## Physics caveats retained

- `dirac_sink` is an absorbing-state bracket, not a complete Dirac RH-state transport model if the BSM interaction can repopulate the active state.
- The resolved mediator softens the propagator but does not by itself establish UV validity of the dimension-five neutrino vertex.
- No laboratory/collider/astrophysical viability conclusion follows from the transport calculation alone.
- No general event-level no-go statement is justified.
- Detector weighting and the proper IceCube likelihood remain separate unfinished tasks.

## Ordered status

1. Run v0.18.1 production-solver validation.
2. Require composite certification.
3. Run certified Cases 5/6/7 and resolved-mediator propagation.
4. Inspect absolute and species-resolved transmitted spectra.
5. Only after that perform detector-response/likelihood work.
