# PKR Earth propagation v0.41 — constrained Eq. (4.4) and Eq. (4.3)

**Status:** propagation physics gate passed at screening level.  This note supersedes propagation conclusions based on unconstrained/fixed-cut Case-6 benchmarks, but preserves those older records as audit history.

## 1. Question

Return to the original aim of the project: after imposing the updated constraints and the corrected EFT-validity condition, can either of the final non-resonant PKR interactions produce a significant Earth-propagation signal?

We study

- PKR Eq. (4.4), resolved derivative `Z'` interaction;
- PKR Eq. (4.3), resolved transition-dipole `Z'` interaction.

For both models

\[
\kappa_i \equiv {c_i g_q\over\Lambda},\qquad
\Lambda_{\rm vert}={\Lambda\over c_i}={g_q\over\kappa_i},
\]

and the validity cut is imposed inside every PDF integral,

\[
\boxed{\sqrt{\hat s}<g_q/\kappa_i}.
\]

The resolved propagator is retained exactly.  For equal `(M,kappa)`, the Eq. (4.3) hard kernel is twice Eq. (4.4).

## 2. Current constraint choice

CMS EXO-24-007 (23 March 2026) searches for 50--300 GeV narrow spin-1 dijet resonances with 138 fb^-1 and quotes observed universal-quark coupling limits spanning approximately

\[
g_q<0.03\text{--}0.13,
\]

assuming the resonance decays only to quarks.  Because the PKR mediator also has an invisible neutrino channel, using the strongest end, `gq=0.03`, uniformly is deliberately conservative for this first propagation gate; a mass-by-mass branching-ratio recast can only be less trivial than this screen.

Current IceCube HESE cross-section data cover 60 TeV--10 PeV and remain compatible with the SM, but the published energy-bin uncertainties are much larger than a 5% extra contribution.  There is no published modern operator-specific `kappa(M)` likelihood for either PKR dimension-5 vertex.  Therefore this note does **not** invent a modern monojet/IceCube `kappa` exclusion; it uses the direct CMS quark-coupling envelope plus the already documented PKR historical product bounds as screening constraints.

## 3. Why only 50 GeV is propagated in this conservative pass

The v0.39 peak-capacity calculation already proved that at `gq=0.03` the 100 and 220 GeV points cannot reach a 5% UHE cross-section modification for either interaction after the true vertex-validity cut is imposed.

Using the exact v0.39 scaling with fixed `gq`, the maximum 100-PeV fractions are approximately

| M | Eq. (4.4) max | Eq. (4.3) max |
|---:|---:|---:|
| 100 GeV | 1.62% | 3.24% |
| 220 GeV | 0.22% | 0.44% |

so these masses are not useful for a conservative 5% propagation benchmark.  The 50 GeV point remains open and lies inside the new CMS search range.

## 4. Self-consistent 50 GeV benchmarks

CT18NLO dynamic-Q folds with `gq=0.03` give

| model | kappa_5% [GeV^-1] | Lambda_vert=gq/kappa |
|---|---:|---:|
| Eq. (4.4) derivative | 1.05797e-4 | 283.6 GeV |
| Eq. (4.3) dipole | 5.88269e-5 | 510.0 GeV |

Both are genuine first crossings of `sigma_NP/sigma_SM=5%` at 100 PeV with the moving validity cut evaluated inside the fold.

The NP/SM fractions remain close to 5% over the relevant energy range:

| E | Eq. (4.4) | Eq. (4.3) |
|---:|---:|---:|
| 1 PeV | 5.43% | 4.80% |
| 3 PeV | 5.35% | 4.87% |
| 10 PeV | 5.24% | 4.94% |
| 30 PeV | 5.13% | 4.98% |
| 100 PeV | 5.00% | 5.00% |
| 300 PeV | 4.88% | 5.00% |

This flatness is a useful physical feature of the light resolved mediator: the propagator tames the momentum growth that made the contact operators pathological.

## 5. Earth sterile-sink propagation envelope

For the process `nu_tau N -> N_s X` with a sufficiently sterile/long-lived daughter, the new channel is a genuine active-flux sink.  Independently of SM regeneration, its exact disappearance optical depth is

\[
\tau_N(E,\theta)=N_A X(\theta)\,\sigma_{\rm NP}(E),
\qquad
P_{\rm sink}=e^{-\tau_N}.
\]

The PREM vertical column obtained by the repository code is

\[
X_{\rm vert}=1.09479\times10^{10}\ {m g\,cm^{-2}}.
\]

### Vertical chord

| E | Eq. (4.4) survival | Eq. (4.3) survival |
|---:|---:|---:|
| 1 PeV | 0.704 | 0.733 |
| 3 PeV | 0.562 | 0.592 |
| 10 PeV | 0.386 | 0.407 |
| 30 PeV | 0.231 | 0.241 |
| 100 PeV | 0.102 | 0.102 |
| 300 PeV | 0.0355 | 0.0328 |

Thus even a 5% interaction-level effect corresponds to a pure-sink suppression of roughly

- 27--30% at 1 PeV;
- 59--61% at 10 PeV;
- 90% at 100 PeV;
- 96--97% at 300 PeV

for a vertical Earth chord.

### Zenith dependence at 100 PeV

For both models (they are matched to 5% there), the extra sterile-sink suppression is approximately

| |cos theta| | pure-sink suppression |
|---:|---:|
| 0.1 | 7.5% |
| 0.3 | 23.6% |
| 0.5 | 41.2% |
| 0.7 | 56.9% |
| 1.0 | 89.8% |

This strong angular dependence is the central observable.

## 6. What this does and does not prove

### Proven by v0.41d

1. Both final PKR models have a conservative current-constraint benchmark at `M=50 GeV, gq=0.03` that gives a self-consistent 5% UHE cross-section contribution.
2. The corrected `sqrt(shat)<gq/kappa` validity cut is applied point-by-point in the dynamic-Q PDF fold.
3. The new interaction is a sterile sink, and its PREM optical depth is large: the Earth can amplify a 5% cross-section modification into an O(10--90%) zenith-dependent flux effect.
4. Eq. (4.3) needs the smaller product coupling and has the higher vertex-validity scale; Eq. (4.4) has fewer transition-moment-specific model complications.

### Not yet claimed as exact detector prediction

The table above isolates the BSM sink.  The final arriving `nu_tau` spectrum also contains SM NC downscattering and `nu_tau -> tau -> nu_tau` regeneration.  These refill lower-energy bins and must be solved together with the sink for a precision spectrum.  The repository's repaired mixed transport also tracks the secondary `nu_mu` from tau decay; because the BSM sink removes parent `nu_tau`, the secondary `nu_mu` component must change as well and cannot be set identically to the SM value.

A full mixed-grid PREM calculation was launched during this certification pass but the matrix-exponential solve remained numerically stiff on the UHE grid.  No result from that unfinished solve is used here.  The v0.41d optical-depth envelope is therefore the certified propagation **screen**, not a detector likelihood.

## 7. Physics interpretation

This substantially changes the paper outlook.  The surviving light resolved mediators do not need an enormous raw cross-section enhancement.  Because Earth attenuation is exponential, an experimentally allowed few-percent sterile sink can generate an order-one change in the upward `nu_tau` spectrum while leaving downward events essentially unaffected.

That provides a clean target observable:

\[
\boxed{\text{energy-dependent up/down and zenith distortion of the }\nu_\tau\text{ flux}}
\]

with a correlated reduction/modification of the tau-decay secondary `nu_mu` component.

The two operators are almost degenerate at the propagation level once matched to the same 5% cross-section fraction; their discrimination comes mainly from collider/width/transition-moment constraints and their different required `kappa`, not from radically different Earth attenuation shapes.

## 8. Decision

\[
\boxed{\text{PROPAGATION GATE: PASS}}
\]

The original paper aim is therefore still viable.  The correct next physics task is no longer an operator survey.  It is to stabilize/converge the coupled PREM `nu_tau-tau` transport for these 50 GeV constrained benchmarks and then fold the species-resolved spectra into IceCube/ARCA/Gen2 observables, including secondary `nu_mu` tracks.

## 9. Reproducibility

- `scripts/pkr_sink_envelope_v0_41d.py`
- `results/v0.41d/benchmark_gate.csv`
- `results/v0.41d/sink_envelope.csv`
- `results/v0.41d/metadata.json`
- PDF: CT18NLO
- benchmark: `M=50 GeV, gq=0.03`
- validity: `sqrt(shat)<gq/kappa` inside each PDF fold
- Earth: PREM chord profile

Earlier failed/contact/strong-coupling candidates remain preserved in the master portfolio and are not deleted.
