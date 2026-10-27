# nutau-prop v0.12 progress checkpoint — 2026-10-27

Author/project owner: Subhendu Rakshit

## Ordered task requested in Chat

The required order is locked as:

1. **Task 3 first:** make the Case-6 detector interaction treatment channel-consistent.
2. **Task 1 second:** propagate the same differential Case-6 kernel through the Earth.
3. Reassess distance to submission only after 1–2 have numerical outputs.
4. Then test other momentum-growing operators.
5. Only after that optimize/UV-complete the physical turnover scale.

No earlier file or result has been deleted or replaced.

## Historical transport validation already established

The repository's v0.3 validation protocol documents and previously passed both historical branches:

- **Rakshit–Reya:** separate nu/antinu propagation; RR-style effective-Z iteration; max `|F3/F2-1| < 5%` in relevant bins; non-negligible second iteration for a hard `E^-1` spectrum; independent direct transport agreeing with `Z^(2)` to 10%; the expected `10^4–10^5 GeV` regeneration bump for `E^-1`; no comparable bump for `E^-2`.
- **Dutta–Reno–Sarcevic:** inclusive tau-decay distribution from PRD 62, 123001; vertical hard-spectrum historical scale near 1 PeV (`nu_tau ~6%`, secondary `nu_mu/e ~1%`, with broad tolerance for modern-vs-CTEQ5 inputs); consistency with the approximate `B(F_nutau-F_numu)` relation; smaller secondary contribution for `E^-2` than `E^-1`.

`PROGRESS.log` records both branches as **PASS**.  The generated JSON/CSV historical outputs were not committed on the current GitHub main branch, so `run_case6_v0_11_v0_12.sh` refuses to silently assume them on a clean machine: it reruns the historical gate whenever `results/historical_gate.json` is absent locally.

## v0.11 — detector NSI consistency

New code:

- `src/nutau_prop/case6.py`
  - computes the Case-6 `d sigma/dy` kernel with CT10nlo;
  - retains the v0.7/v0.8 coefficient convention;
  - primary domain prescription `sqrt(shat)<Lambda`;
  - keeps SM CC unchanged;
  - adds Case 6 only to NC;
  - returns the true combined NC redistribution distribution rather than only a total cross section.

- `scripts/build_case6_tables_v0_11.py`
  - parallelizes the PDF calculation over energy and nu/antinu using `NUTAU_NPROC` (default 14);
  - stores matched SM `p_nc_sm`, Case-6 `p_np`, and combined `p_nc` separately;
  - spans 10 TeV–20 EeV so the same tables cover propagation and the UHE stress test.

- `scripts/validate_case6_differential_v0_11.py`
  - checks probability normalization;
  - checks that CC is independent of the Case-6 cutoff;
  - reconstructs the nu+antinu averaged `sigma_NP/sigma_SM` and requires agreement with the frozen v0.8 cross-section result to 10%.

- `scripts/run_case6_detector_v0_11.py`
  - reweights the public ARCA shower components channel by channel;
  - `nu_e CC` unchanged;
  - `nu_tau CC` unchanged at the detector vertex;
  - only the tau one-third of all-flavor NC gets the Case-6 detector interaction;
  - direct muon-track increment is identically zero, as required by the NC-like topology;
  - reports the Case-6 mean inelasticity as a visible-energy diagnostic;
  - retains IceCube EHE tau-initial reweighting only as an explicitly labelled response proxy because the released IceCube area does not separate tau CC and NC selection.

Residual detector limitation: public effective areas have integrated over the SM inelasticity distribution.  Exact altered-y reconstruction efficiency requires detector MC and cannot be manufactured from the released one-dimensional effective areas.

## v0.12 — full Case-6 Earth transport

New code:

- `src/nutau_prop/transport_case6.py`
  - new flavor-aware transport path, leaving the legacy validated transport untouched;
  - primary tau flavor uses SM CC + SM/Case-6 NC as appropriate;
  - tau decay and regeneration use the established RR/DRS machinery;
  - secondary muon flavor propagates with a matched **SM-only** table, preventing a tau-specific Case-6 NSI from being incorrectly inherited by secondary `nu_mu`.

- `scripts/run_case6_transport_v0_12.py`
  - propagates nu_tau and antinu_tau separately through PREM;
  - includes NC energy degradation with the actual combined `d sigma/dy` kernel;
  - includes tau production, energy loss, polarized decay/regeneration and secondary muon flavor;
  - compares against a matched CT10 SM baseline;
  - scans `Lambda = 1.5, 2.0, 2.44, 3.0 TeV`, `|cos(theta)| = 0.2,0.5,0.8,1.0`, and `E^-2`, `E^-2.46`, and Ehlert incident shapes;
  - records the propagated ratios at 1, 10, 100, 220 and 300 PeV.

The numerical submission-distance judgment is deliberately not written before these stages are run.

## Ordered execution

Use:

```bash
export NUTAU_NPROC=14
bash run_case6_v0_11_v0_12.sh
```

The run writes the compact terminal record to `results/v0.12/run_case6_v0_11_v0_12_summary.txt` and versioned numerical products under `results/v0.11/` and `results/v0.12/`.
