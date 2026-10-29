# nutau-prop v0.26 — active d=8 certification status

Date: 2026-10-29
Status: analytic gates closed; exact dynamic-Q + Run-II di-tau execution packaged; Earth propagation blocked pending numerical gate.

## Executive decision

The preferred active operators remain

- `l2u2D22[3322]`: `(L3bar gamma_mu i<->D_nu L3)(cRbar gamma^mu i<->D^nu cR)`;
- `l2d2D22[3322]`: `(L3bar gamma_mu i<->D_nu L3)(sRbar gamma^mu i<->D^nu sR)`.

They are gauge invariant, conserve lepton number, keep `nu_tau -> nu_tau` active, and generate amplitudes growing as `shat^2/Lambda^4`.

The certification sequence has produced a stronger result than v0.24: **the binding elastic partial wave is the antiquark J=0 channel.** With the strict convention that the EFT is unitary through `shat=Lambda^2`,

`|C| <= 8 pi / 3 = 8.37758`.

Forward-limit positivity for both single-coefficient directions is

`-4 C >= 0`, hence **C <= 0**.

Therefore the canonical strict benchmark is

`C = -8 pi/3`, `Lambda = 4 TeV`.

At this point the pure `C^2` contribution is small; the survival of the contact-operator branch is decided by the SM x d=8 interference in the exact dynamic-Q UHE fold.

## Gates closed

### 1. Exact partonic normalization — PASS

For fixed LH neutrino and an unpolarized q_R PDF,

`d sigma_q^NP/dy = C^2 shat^3/(8 pi Lambda^8) (2-y)^2 (1-y)^2`,

`d sigma_qbar^NP/dy = C^2 shat^3/(8 pi Lambda^8) (2-y)^2`.

The angular moments are

`Iq = 31/30`, `Iqbar = 7/3`, `Iqbar/Iq = 70/31 = 2.2581`.

Thus sea antiquarks are intrinsically important.

### 2. Partial-wave unitarity — PASS

Quark channel:

`a1(q) = 7 C rho^2/(96 pi)`, `a2(q)=C rho^2/(160 pi)`.

Antiquark channel:

`a0(qbar)=3 C rho^2/(16 pi)`, `a1(qbar)=C rho^2/(48 pi)`.

The antiquark J=0 channel binds:

`|C| rho^2 <= 8 pi/3`.

For `rho=1`: `|C|<=8.37758` under `|Re a_J|<=1/2`.

If the trusted EFT domain ends instead at `Mcut^2=xi Lambda^2`, the invariant form is

`|G8| Mcut^4 <= 8 pi/3`, `G8=C/Lambda^4`.

Do not mix this `unitary-to-cut` convention with the stricter `unitary-to-Lambda` default.

### 3. Positivity — PASS

Li et al., JHEP 10 (2022) 107, arXiv:2204.13121, Table 1 gives

- `-4 C_8,lu partial2 >=0`;
- `-4 C_8,ld partial2 >=0`.

Hence the single-coefficient branch is `C<=0`.

This fixes the interference sign rather than allowing it to be selected for phenomenological convenience.

### 4. Published Drell-Yan comparison — NON-BINDING at strict unitarity point

The same positivity/Drell-Yan study quotes projected HL-LHC individual bounds at `Lambda=2 TeV` for the flavor-universal directions. With the full angular information the approximate negative-side limits are

- `C_8,lu partial2 > -1.2`;
- `C_8,ld partial2 > -1.1`.

For a fixed dimensionful coefficient, rescaling to `Lambda=4 TeV` multiplies dimensionless C by `(4/2)^4=16`, giving rough projected limits

- `C_lu > -19.2`;
- `C_ld > -17.6`.

The strict unitarity value `C=-8.3776` is inside both. These projected constraints are based on a far more collider-visible flavor-universal/light-lepton setup; the exact `tau tau c c` and `tau tau s s` Run-II bounds should therefore not be stronger in ordinary circumstances. This is a conservative literature argument that Drell-Yan is unlikely to bind before unitarity. The exact HighPT flavor recast is nevertheless packaged and should be retained as validation.

## Pure-term size at strict unitarity

Rescaling the v0.23 CT18 prescreen from `C=26` to `C=8.37758` gives, for `xi=1`, before interference:

| E | c_R pure / SM NC | s_R pure / SM NC | c_R pure / total SM | s_R pure / total SM |
|---|---:|---:|---:|---:|
|100 PeV|0.34%|0.46%|0.11%|0.14%|
|300 PeV|0.70%|0.87%|0.22%|0.28%|

Therefore any phenomenologically interesting contact signal must now come mainly from interference.

For the positivity sign `C<0`:

- charm has `g_R^u<0`, so interference is constructive in the project convention;
- strange has `g_R^d>0`, so interference is destructive.

This can reverse the old raw-PDF ranking and make charm the preferred candidate.

## Gates still requiring numerical execution

### A. Exact dynamic-Q CT18 UHE fold — DECISIVE

Committed code: `scripts/d8_uhe_cert_v0_25.py`.

v0.26 execution uses

- CT18NLO central member;
- `C=-8 pi/3`, `Lambda=4 TeV`;
- `xi=1,1/2,1/4`;
- `E=1,10,100,300 PeV`;
- exact `Q^2=x y s` PDFs;
- SM CC+NC reference;
- exact SM x d8 interference plus pure d8 term;
- endpoint-support diagnostics.

This gate decides whether either contact operator can still produce a few-percent NC/propagation effect.

### B. Exact Run-II tau-tau HighPT recast — validation

Committed script: `scripts/highpt_d8_tautau_v0_26.wls`.

It uses the exact HighPT WC labels and the ATLAS full Run-II `di-tau-ATLAS` likelihood. It bootstraps the official HighPT package if needed. The fit uses `EFTorder->4`, i.e. the consistent SM x d8 interference truncation.

### C. Monojet/dineutrino — not certified

No published exact dimension-8 monojet limit for the two `D22[3322]` directions was identified. Dimension-6 missing-energy/jet SMEFT literature is not transferable as an exact bound. Because the strict-unitarity UHE fold may already kill the contact branch, a dedicated monojet implementation should only be undertaken if Gate A leaves a useful signal.

### D. IceCube 60 TeV–10 PeV — published fit is not directly transferable

IceCube PRD 104, 022001 fits four cross-section scaling factors in the bins 60–100, 100–200, 200–500 TeV and 0.5–10 PeV, but assumes

- common scaling for all three flavors;
- simultaneous scaling of CC and NC with fixed CC/NC ratio.

Our operator is tau-only and NC-like. Therefore the published `x_i` likelihood cannot be imposed as an exact constraint. The proper recast requires propagation plus detector/PID response. A simple all-flavor scaling proxy may be reported only as a screen, not as a likelihood.

The published 68.3% frequentist intervals are very broad, especially in the highest bin, so a sub-percent to few-percent tau-only correction below 10 PeV is not expected to bind; the exact dynamic-Q output should verify the size before any dedicated reweight is attempted.

## Execution package

Run from the repo root:

`bash run_d8_cert_v0_26.sh`

The runner performs the analytic regression, exact dynamic-Q CT18 fold, exact HighPT Run-II tau-tau scan if Wolfram is available, writes all outputs under `results/v0.26/`, then commits and pushes those outputs automatically.

A GitHub Actions workflow was also added (`.github/workflows/d8-dynamicq-v025.yml`) to provide an LHAPDF runtime, but no workflow run was registered after its creation; it should not be treated as executed.

## Decision rule

1. If the strict `C=-8 pi/3` dynamic-Q result stays below a few percent in the NC channel at 100–300 PeV, **stop this contact branch** and promote the resolved spin-2/gluon/mediator branches.
2. If charm or strange remains at several percent, check its exact HighPT result and only then invest in monojet and a dedicated IceCube tau-only NC recast.
3. Earth propagation remains blocked until these gates are closed.
