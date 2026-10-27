# nutau-prop v0.8 — UHE momentum-growing Case-6 event forecast

This package is built directly on the supplied v0.6 project and preserves the earlier validation and ARCA machinery.

## New in v0.7
- Genuine dimension-7 momentum-growing neutrino-quark Case-6 operator.
- CT10nlo PDF integration over Bjorken x and inelasticity y.
- Two EFT validity prescriptions:
  - **primary/conservative:** `sqrt(shat) < Lambda_UV` (`r_shatcut`)
  - **aggressive envelope:** `Q < Lambda_UV` (`r_Qcut`)
- Benchmarks Lambda = 1.0, 1.5, 2.0, 2.44, 3.0, 4.0 TeV.
- The supplied `results/v0.7/case6_pdf_cut_v0_7_summary.txt` records the successful Mac Studio validation supplied by the user.

## New in v0.8
- Extends the Case-6 grid to 20 EeV.
- Uses a cosmogenic Ehlert spectrum rather than the old TeV-PeV E^-2.46 flux.
- Downloads and folds the official 12.6-year IceCube EHE effective-area product.
- Downloads and validates the official KM3NeT KM3-230213A bright-track response against the collaboration's Ehlert expectation (0.242767 events in 335 days; central-90% expectation 0.211117).
- Computes IceCube inclusive-response proxy counts and clearly-labelled KM3NeT NC-equivalent acceptance proxies.

## Important topology caveat
Case 6 is NC-like. It does not directly generate a new muon track. Therefore `KM3_brighttrack_direct_NP` remains equal to the SM bright-track expectation before a dedicated propagation/topology treatment. `KM3_NC_equiv_*` columns are diagnostic acceptance proxies, not claimed physical shower or track counts.

## Run on the Mac Studio
```bash
conda install -c conda-forge lhapdf   # only if LHAPDF is absent
lhapdf install CT10nlo                 # only if CT10nlo is absent
./run_v0_8.sh
```

Outputs are written under `results/v0.8/`.

## Scientific status
v0.8 is a detector-level screening calculation, not yet the final publication likelihood. The next refinement is topology-aware IceCube/KM3NeT detector response plus Earth transport/regeneration using the same differential Case-6 kernel.
