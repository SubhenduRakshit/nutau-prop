# PKR absorption-horizon propagation v0.42

## Status

\[
\boxed{\text{Coupled absorption-horizon mechanism: NUMERICALLY PASSED}}
\]

\[
\boxed{\text{Allowed-model + official detector-likelihood gate: OPEN}}
\]

This document supersedes any interpretation of v0.41 as a full propagation certification.  v0.41 remains the conditional pure-sink screen; v0.42 adds coupled \(\nu_\tau\leftrightarrow\tau\) transport, NC degradation, tau-decay secondary \(\nu_\mu\), and a public-response detector fold.

## 1. Conditional benchmark used

The calculation reuses the v0.41d self-consistent 100-PeV 5% interaction benchmarks at

\[
M_{Z'}=50\ {\rm GeV},\qquad g_q=0.03,
\]

with the validity condition imposed inside the PDF fold,

\[
\boxed{\sqrt{\hat s}<g_q/\kappa}.
\]

| interaction | \(\kappa\) [GeV\(^{-1}\)] | \(g_q/\kappa\) |
|---|---:|---:|
| PKR Eq. (4.4) derivative | \(1.0579696\times10^{-4}\) | 283.6 GeV |
| PKR Eq. (4.3) dipole | \(5.8826895\times10^{-5}\) | 510.0 GeV |

These are **conditional transport benchmarks, not established allowed parameter points**.

## 2. Two transport topologies

The same hard interaction is propagated in two limiting realizations:

1. **sterile sink** — the daughter is an invisible sterile state; this is the maximal-disappearance envelope;
2. **active redistribution proxy** — implemented by the validated `majorana_flip` transport mode to quantify active \(\nu_\tau\leftrightarrow\bar\nu_\tau\) redistribution.

The second is a transport proxy and is **not** asserted to be the unique literal-PKR flavour completion.

## 3. UHE coupled transport result

The useful region is the absorption-horizon band

\[
1\lesssim\tau_{\rm SM}\lesssim3,
\]

not the vertical chord.  A UHE-only discovery grid, \(1\)–\(300\) PeV, was run at

\[
|\cos\theta|=0.12,\ 0.25,\ 0.45,\ 0.70.
\]

All eight operator/zenith transport tasks completed successfully.  They include

- SM CC and NC transport;
- \(\nu_\tau\to\tau\to\nu_\tau\) regeneration;
- tau-decay secondary \(\nu_\mu\);
- PREM Earth density;
- the PKR extra interaction with the coupling-dependent \(\hat s\) cut.

A public detector projection then applies the calibrated ARCA Fig. 7 upgoing tau-CC energy response and a controlled Gaussian migration of width

\[
\sigma_{\log_{10}E}=0.30.
\]

This is not an official ARCA likelihood.

## 4. Focused numerical convergence

The dominant discovery-grid chord, \(|\cos\theta|=0.70\), was rerun at three resolutions:

| setting | \(N_E\) | tau refinement | max PREM step |
|---|---:|---:|---:|
| coarse | 15 | 1 | 192 km |
| medium | 25 | 2 | 96 km |
| fine | 31 | 4 | 48 km |

The medium-to-fine relative changes are:

| model/topology | reconstructed cell | event-weighted horizon | true 1 PeV | true 1.5 PeV |
|---|---:|---:|---:|---:|
| Eq. (4.3) active redistribution | 2.89% | 3.05% | 2.76% | 1.12% |
| Eq. (4.3) sterile sink | 0.65% | 0.71% | 0.004% | 0.043% |
| Eq. (4.4) active redistribution | 2.03% | 2.15% | 2.31% | 0.96% |
| Eq. (4.4) sterile sink | 0.60% | 0.66% | 0.007% | 0.039% |

Thus every monitored observable satisfies the prechosen preferred convergence criterion of less than 5% medium-to-fine drift.

## 5. Fine-grid result at the discovery optimum

At the discovery-grid point

\[
|\cos\theta|=0.70,\qquad E_{\rm rec}\simeq1.259\ {\rm PeV},\qquad
\tau_{\rm SM}^{\rm eff}\simeq2.707,
\]

the fine-grid reconstructed distortions are

| model | sterile sink | active redistribution proxy |
|---|---:|---:|
| Eq. (4.4) derivative | \(-14.08\%\) | \(-8.86\%\) |
| Eq. (4.3) dipole | \(-12.69\%\) | \(-7.33\%\) |

The SM-event-weighted absolute distortion over the local \(1\le\tau_{\rm SM}^{\rm eff}\le3\) reconstructed horizon band is

| model | sterile sink | active redistribution proxy |
|---|---:|---:|
| Eq. (4.4) derivative | 14.05% | 8.91% |
| Eq. (4.3) dipole | 12.66% | 7.40% |

The important physical conclusion is that active redistribution weakens the maximal sterile-sink signal by roughly one third, but **does not erase the absorption-horizon feature**.

## 6. Secondary muon-neutrino signal

The coupled calculation also confirms that tau decay transfers the BSM change into the secondary \(\nu_\mu\) flux.  On the fine grid at \(|\cos\theta|=0.70\):

- Eq. (4.4) derivative sterile sink: secondary \(\nu_\mu\) deficit about 8.0% at 1 PeV and 8.9% at 1.5 PeV;
- Eq. (4.4) active redistribution proxy: about 4.3% and 4.9%;
- Eq. (4.3) dipole sterile sink: about 7.3% and 8.2%;
- Eq. (4.3) active redistribution proxy: about 3.5% and 3.9%.

Therefore the track channel is not merely an unaffected control sample; it must be included in the eventual detector likelihood.

## 7. Important cutoff caveat

The strongest reconstructed cell above lies close to the transport lower boundary \(E_{\min}=1\) PeV.  Because the illustrative migration width is 0.30 dex, events reconstructed near 1.26 PeV can receive contributions from true energies below the present grid.

Therefore the \(1.26\)-PeV detector-fold magnitude is **numerically converged with respect to the internal grid but not yet certified against the lower-energy boundary**.

This caveat does not invalidate the coupled-transport mechanism itself.  The next certification should use an interior absorption-horizon point, approximately

\[
E_{\rm rec}\sim3\ {\rm PeV},\qquad |\cos\theta|\sim0.45,
\]

for which the 1-PeV lower boundary is well separated from the detector migration kernel.

## 8. Binding interpretation

The result that is now justified is

\[
\boxed{\text{Regeneration and active redistribution do not kill the PKR absorption-horizon signal.}}
\]

The result that is **not** yet justified is

\[
\boxed{\text{an experimentally allowed PKR model is detectable by ARCA/IceCube.}}
\]

That stronger statement still requires both:

1. an interior-grid detector-fold certification with a proper reconstructed-energy/zenith response;
2. the modern allowed-model constraint gate: mediator width/branching fractions, monojet or missing-energy likelihood, dijet/direct-\(Z'\) constraints, and IceCube likelihood, with \(c_i\), \(g_q\), \(\Lambda\), and \(\Gamma_{Z'}\) kept separate.

## 9. Reproducibility / audit trail

The focused convergence calculation completed successfully in GitHub Actions run `41030476091`.  Its physics step produced the full coarse/medium/fine tables and the convergence numbers quoted above.  The workflow's final `git push` was rejected only because the remote `main` branch advanced concurrently; this was a repository race, not a calculation failure.  The failed push is retained as part of the audit trail.

The v0.42 code and earlier failed/slow workflow variants are intentionally preserved rather than deleted.
