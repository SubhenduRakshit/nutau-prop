# PKR Earth propagation v0.41 — conditional sterile-sink / absorption-horizon screen

**Status:** **SUPERSEDED AS A PROPAGATION CERTIFICATION.**  v0.41 is retained as a conditional sterile-sink screen only.  The allowed-model and detector-level propagation gate remains open.

## Corrected decision

The three binding statuses are

\[
\boxed{\text{Absorption-horizon propagation mechanism: PROMISING}}
\]

\[
\boxed{\text{Conditional active-to-sterile sink screen: PASSED}}
\]

\[
\boxed{\text{Allowed-model + detector-level propagation gate: OPEN}}
\]

The earlier statement `PROPAGATION GATE: PASS` is **retracted**.

## 1. Scope

This screen asks whether the two surviving light resolved PKR interactions can generate a potentially observable Earth-propagation distortion after the updated split-coupling validity condition is imposed:

- PKR Eq. (4.4), resolved derivative `Z'` interaction;
- PKR Eq. (4.3), resolved transition-dipole `Z'` interaction.

For both models

\[
\kappa_i={c_i g_q\over\Lambda},\qquad
\Lambda_{\rm vert}={\Lambda\over c_i}={g_q\over\kappa_i},
\]

with

\[
\boxed{\sqrt{\hat s}<g_q/\kappa_i}
\]

imposed inside every PDF integral.

## 2. Conditional 50 GeV benchmark

The screen used

\[
M_{Z'}=50\ {\rm GeV},\qquad g_q=0.03,
\]

and found the self-consistent 100-PeV 5% interaction benchmarks

| model | `kappa_5%` [GeV^-1] | `Lambda_vert=gq/kappa` |
|---|---:|---:|
| Eq. (4.4) derivative | 1.05797e-4 | 283.6 GeV |
| Eq. (4.3) dipole | 5.88269e-5 | 510.0 GeV |

These are **conditional benchmarks, not established allowed points**.  A complete model must still satisfy the mediator-width/branching-ratio, monojet/missing-energy, dijet and IceCube likelihood constraints with the same parameter split.

## 3. Topology caveat: literal PKR versus sterile-N completion

The pure-sink calculation is **not** the generic literal PKR prediction.  Read literally, Eqs. (4.3)/(4.4) are off-diagonal/LNV active-flavour transitions.  Treating the daughter as a disappearing sterile state requires an explicit active-to-sterile completion.

Therefore the v0.41 sink result should be interpreted as the **maximal disappearance envelope** for such a sterile-N completion.  For a literal active/antineutrino daughter, part of the interaction is redistribution rather than loss and the `nu_tau` deficit must be computed by the coupled transport system.

## 4. Why the vertical headline was misleading

For a sterile sink with fractional extra interaction strength `r`, the additional suppression relative to SM transmission is approximately

\[
\Delta=1-e^{-r\tau_{\rm SM}}.
\]

Large fractional suppression at very large `tau_SM` is not automatically useful because the SM event density is already exponentially depleted there.  A Poisson-only, monochromatic, single-bin toy figure of merit behaves for `r << 1` as

\[
\Delta\sqrt{N}\propto r\tau_{\rm SM}e^{-\tau_{\rm SM}/2},
\]

which peaks at

\[
\tau_{\rm SM}=2.
\]

The maximum is broad; the relevant analysis region is approximately

\[
\boxed{1\lesssim\tau_{\rm SM}\lesssim3},
\]

not the vertical chord.  The exact detector optimum can move within this band because of effective area, solid angle, backgrounds, energy migration and regeneration.

## 5. Certified result from v0.41d

The exact dynamic-Q sterile-sink optical-depth calculation demonstrates that a few-percent allowed interaction contribution can generate an order-10% distortion near the SM absorption horizon.  It also demonstrates that the enormous vertical deficits previously highlighted occur mostly where the direct SM flux is already strongly attenuated.

Hence the robust physics statement from v0.41 is:

\[
\boxed{\text{A moving absorption-horizon distortion is potentially observable.}}
\]

It is **not** yet a detector-level sensitivity result.

## 6. Regeneration has no predetermined sign

No claim is made that the pure-sink result is a conservative floor.  In the full system:

1. `nu_tau -> tau -> nu_tau` regeneration refills lower-energy bins;
2. this can reduce the local fractional deficit at fixed `(E,theta)`;
3. the sterile sink removes parents that would otherwise generate the full descendant chain, which can increase the deficit elsewhere;
4. NC migration and reconstructed-energy migration further mix these effects.

The net sign is therefore energy- and zenith-dependent and must be computed.

## 7. Analysis strategy after correction

The next calculation must center the detector analysis on the broad absorption-horizon strip rather than on vertical events:

- use `(E_rec, cos(theta))` bins;
- bin more finely in zenith than in reconstructed energy;
- overlay `tau_SM = 1,2,3` contours;
- display the SM event density on the same axes as the BSM fractional distortion;
- include the coupled `nu_tau <-> tau` system, NC degradation, tau-decay secondary `nu_mu`, detector effective area, migration and backgrounds;
- compare the sterile-N disappearance completion against the literal active/LNV redistribution topology.

The approximate scaling `sigma_SM ~ E^0.36` is useful only for intuition.  Numerical work should use the local slope

\[
\alpha(E)={d\ln\sigma_{\rm SM}\over d\ln E}
\]

from the actual cross-section tables.

## 8. Reproducibility

The original v0.41d numerical screen remains preserved:

- `scripts/pkr_sink_envelope_v0_41d.py`
- `results/v0.41d/benchmark_gate.csv`
- `results/v0.41d/sink_envelope.csv`
- `results/v0.41d/metadata.json`
- CT18NLO PDFs
- PREM Earth model

These results are retained as the conditional sterile-sink envelope and are not deleted.

## 9. Binding next gate

\[
\boxed{\text{v0.42: coupled transport + detector-weighted absorption-horizon map}}
\]

Only after that calculation, together with the explicit allowed-model constraint check, may the project claim a propagation gate pass.
