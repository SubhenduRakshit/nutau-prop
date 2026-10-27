# nutau-prop v0.12 — Case-6 detector + Earth-transport consistency gate

This stage implements the user-requested order: **detector NSI first (v0.11), then full Case-6 Earth transport (v0.12)**.

## Physics changes

- Case 6 is treated as a tau-specific NC-like interaction.
- The exact `d sigma/dy` kernel is computed with CT10nlo and the primary `sqrt(shat)<Lambda` EFT-domain cut.
- SM CC is unchanged.
- Detector-level ARCA showers are reweighted channel by channel; direct muon tracks receive no Case-6 vertex enhancement.
- Earth transport uses the same differential NC kernel for energy degradation.
- Tau regeneration and polarized secondary muon-neutrino production are retained.
- Secondary muon flavor propagates with SM-only interactions, not the tau-specific Case-6 table.

## Historical gate

The repository previously passed the Rakshit-Reya and Dutta-Reno-Sarcevic validation gates.  Because their generated JSON/CSV outputs are not tracked on current `main`, the ordered runner reruns historical validation if `results/historical_gate.json` is absent locally.

## Run

```bash
cd ~/Research/nutau-prop
git pull --ff-only
python -m pip install -e .
export NUTAU_NPROC=14
bash run_case6_v0_11_v0_12.sh
```

Key outputs:

- `results/v0.11/case6_differential_gate_v0_11.json`
- `results/v0.11/case6_detector_channel_fold_v0_11.csv`
- `results/v0.11/case6_icecube_tau_detector_proxy_v0_11.csv`
- `results/v0.12/case6_earth_transport_v0_12.csv`
- `results/v0.12/case6_earth_transport_summary_v0_12.csv`
- `results/v0.12/run_case6_v0_11_v0_12_summary.txt`

No submission claim should be made until these outputs are inspected.  In particular, the public ARCA/IceCube effective areas have already integrated over the SM inelasticity distribution; exact altered-y selection efficiency still requires collaboration MC or an equivalent response release.
