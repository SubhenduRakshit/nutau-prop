# nutau-prop v0.34 — follow-up to UHE interaction portfolio

Date: 2026-10-30
Status: targeted follow-up to `UHE_INTERACTION_PORTFOLIO_v0.33.md`; no Earth transport run.

## Mandate

Follow only the four next-justified-work items in v0.33:
1. resolved derivative active-to-sterile mediator: UV origin + current monojet/IceCube gate;
2. strange-aligned `\widetilde R_2` resonance: current constraint/width audit;
3. tau-neutrinogluon: retain only if a credible origin for the large transition dipole remains after coloured-partner constraints;
4. Earth attenuation remains blocked until one benchmark passes its applicable constraint gate.

## 1. Resolved derivative active-to-sterile mediator

Low-energy interaction:

\[
\mathcal L \supset \kappa Z'_\mu\,\bar N i\overleftrightarrow{\partial^\mu}P_L\nu_\tau
+ g_q Z'_\mu\bar q\gamma^\mu q .
\]

The v0.32 fold established the size curve. The required coupling for a 5% inclusive UHE effect is approximately

- 5 GeV: `kappa_5 ~ 3.1e-5 GeV^-1`;
- 40 GeV: `~5.5e-5`;
- 100 GeV: `~8e-5`;
- 220 GeV: `~1.2e-4`;
- 500 GeV: `~2.0-2.2e-4`;
- 1 TeV: `~3.7-4.1e-4`.

### UV-completion screen

Electroweak gauge invariance requires an active-neutrino/Higgs insertion above EWSB. A representative gauge-completed scaling is

\[
\kappa \sim c\,{v\over \Lambda_{UV}^2}.
\]

Thus

\[
c\simeq \kappa\Lambda_{UV}^2/v.
\]

For `kappa = 1e-4 GeV^-1`:

- `Lambda_UV = 0.5 TeV` -> `c ~ 0.10`;
- `Lambda_UV = 1 TeV` -> `c ~ 0.41`;
- `Lambda_UV = 2 TeV` -> `c ~ 1.63`.

If the derivative vertex is loop generated in a weakly-coupled completion, `c ~ g_*^n/(16 pi^2)`; reproducing `c~0.4` requires an underlying coupling product of order `16 pi^2 c ~ 60`, and even `c~0.1` requires order `16`. Therefore a conventional perturbative loop origin is not credible for the 5%-size window. A tree-level completion, unusual symmetry structure, or strong/composite sector would be required.

The complete sterile-neutrino EFT literature confirms that derivative higher-dimensional structures exist and distinguishes tree- from loop-generated classes, but no simple renormalizable single-mediator model was found that reproduces the exact low-energy derivative vertex with the required normalization while preserving the same clean `nu_tau q -> N q` topology.

### Current monojet/IceCube evidence audit

No post-2019 publication was found that supplies a directly transferable mass-dependent `tau`-specific bound on this exact derivative sterile operator. Current IceCube sterile-neutrino searches mainly constrain active-sterile oscillation mixing at GeV--TeV energies, not this UHE scattering vertex. The 2025 EHE public release is suitable for a future dedicated likelihood/reweight but is not itself a published `kappa(M)` bound for this operator.

Therefore the historical PKR `kappa` bound must not be promoted to a 2026 certification. A dedicated recast is still necessary.

### Decision

- UHE size: PASS.
- Weakly-coupled loop UV completion: FAIL / highly implausible.
- Exact current collider+IceCube bound: unresolved because no directly reusable public mass-dependent result exists.
- Overall: **CONDITIONAL, now explicitly STRONG/UV-CHALLENGED** rather than a generic perturbative non-resonant success.

## 2. Strange-aligned `\widetilde R_2` scalar leptoquark

Take

\[
\widetilde R_2\sim(3,2,1/6),\qquad
\lambda\,\bar s_R\widetilde R_2^T i\tau_2 L_\tau.
\]

The two weak-doublet components necessarily couple as

\[
R^{2/3}\to s\tau,\qquad R^{-1/3}\to s\nu_\tau.
\]

The v0.33 UHE benchmark is `M=1.5 TeV`, `lambda=2`, with about 7.17% / 9.71% at 100/300 PeV.

### Width

For a single massless fermion channel,

\[
\Gamma\simeq {|\lambda|^2 M\over16\pi}.
\]

At `lambda=2`, `M=1.5 TeV`, this gives about `119 GeV`, i.e. `Gamma/M ~ 7.96%` for each component's dominant Yukawa decay. The narrow-width approximation is therefore only moderately accurate; the benchmark should be recomputed with finite width before publication, but it is not automatically invalid.

### Collider status

QCD pair production is coupling independent. Run-2 ATLAS searches exclude scalar LQs decaying 100% to `b tau` up to about 1.46 TeV. That is not directly the same as a strange-jet final state, so it does not by itself exclude the 1.5 TeV strange-aligned point. Conversely, the `s nu_tau` component gives jets+MET and must be confronted with squark/LQ-like searches; no exact 2026 recast for the two-component `s tau` / `s nu` doublet was found in the literature search.

Recent 2025 ATLAS resonant lepton-jet searches significantly strengthen electron/muon LQ limits, but do not provide a tau+strange result that can be transferred directly. Hence the 1.5 TeV, lambda=2 point remains conditional, not certified or excluded.

### Decision

- UHE size: PASS.
- Renormalizable/gauge-complete interaction: PASS.
- Width: moderate (`~8%`), finite-width correction required.
- Direct exact tau+strange + jets+MET recast: unresolved.
- Overall: **still the cleanest UV-complete candidate in the portfolio if resonances are allowed**.

## 3. Tau-neutrinogluon / leptogluon resonance

Benchmark retained from v0.33:

\[
M=3~\mathrm{TeV},\quad \Lambda=9~\mathrm{TeV},\quad a=4,
\]

with 6.14% / 12.95% at 100/300 PeV.

The vector-like colour-octet doublet itself is gauge consistent, but the transition operator

\[
{a g_s\over2\Lambda}\bar L_\tau\sigma^{\mu\nu}G^A_{\mu\nu}L^A_{8R}
\]

is dimension five. If generated perturbatively at one loop, a coefficient `a=4` is far above a natural loop-sized value; schematically `a~g_*^2/(16 pi^2)` would require `g_*~sqrt(4*16 pi^2)~25`, well outside perturbativity. Thus a strong/composite matching mechanism is required.

A targeted literature search did not find a modern LHC tau-leptogluon analysis that directly certifies the 3 TeV point. QCD pair production of the charged colour-octet partner and its `tau+jet` decay remains the decisive collider channel.

### Decision

- UHE size: PASS.
- Perturbative dipole matching: FAIL.
- Strong/composite origin: possible but unproven.
- Coloured-partner direct bound: unresolved for exact model.
- Overall: **retain only as a strong/composite resonance benchmark, below `\widetilde R_2` in credibility**.

## 4. Portfolio update / Earth-attenuation gate

No benchmark yet satisfies an all-constraints gate. Therefore the v0.33 instruction remains active:

**Do not run Earth attenuation yet.**

Updated ranking after the targeted follow-up:

1. strange-aligned `\widetilde R_2`: best UV-clean candidate, but resonant and collider/finite-width recast pending;
2. resolved derivative sterile mediator: only broad non-resonant size window, but weakly-coupled UV origin is strongly disfavoured and an exact modern recast is absent;
3. tau-neutrinogluon resonance: strong UHE size, but requires strong/composite dipole matching and coloured-partner certification;
4. derivative-tensor / sterile Rayleigh endpoint benchmarks: strong/composite only.

The key change relative to v0.33 is that Priority 1 is no longer regarded as naturally perturbative. The absence of a public exact 2026 `kappa_max(M)` result means the correct next technical task is a dedicated monojet + IceCube likelihood recast if one wants to keep this non-resonant branch alive.

## Evidence notes

- IceCube 2024--2025 sterile-neutrino searches constrain oscillation mixing and are not the same observable as the derivative UHE scattering vertex.
- IceCube's 2025 EHE data release provides public material for a dedicated UHE likelihood/reweight.
- ATLAS Run-2 scalar LQ pair production to `b tau` excludes masses below about 1.46 TeV for 100% branching, useful only as a nearby comparator for the strange-aligned `R2` benchmark.
- ATLAS 2025 resonant lepton-jet searches improve electron/muon LQ constraints but do not directly close the tau+strange benchmark.

