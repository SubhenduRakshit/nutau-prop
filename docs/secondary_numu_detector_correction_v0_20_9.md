# Secondary muon-neutrino regeneration: required detector-level correction

**Status:** open physics correction identified after v0.20.9 transport validation work  
**Date:** 2026-10-29  
**Scope:** Case 6/7/resolved-mediator tau-sector BSM propagation, ARCA and IceCube event folding

## Finding

The Earth-transport calculation already includes secondary muon neutrinos and antineutrinos produced in polarized tau decay,

\[
\nu_\tau N \to \tau^- X,\qquad
\tau^-\to \nu_\tau+\nu_\mu+\cdots,
\]

with the charge-conjugate channel for \(\tau^+\). The implemented leptonic secondary channel uses branching fraction \(B\simeq0.18\).

The mixed combined-generator state is

\[
(\nu_\tau,\bar\nu_\tau,\nu_\mu^{\rm sec},\bar\nu_\mu^{\rm sec},\tau^-,\tau^+),
\]

and the tau-decay blocks inject \(\nu_\mu^{\rm sec}\) and \(\bar\nu_\mu^{\rm sec}\), which subsequently propagate with their Standard-Model CC/NC interaction blocks.

Therefore a tau-sector BSM interaction changes not only the transmitted/regenerated tau-neutrino flux but also the secondary muon-neutrino flux:

\[
\text{BSM}\to \Phi_{\nu_\tau,\bar\nu_\tau}
\to \tau^\mp\text{ production}
\to \tau^\mp\text{ decay}
\to \Phi_{\nu_\mu^{\rm sec},\bar\nu_\mu^{\rm sec}}.
\]

This effect has already been exposed at transport level through quantities such as

\[
R_\mu^{\rm sec}(E,\theta)=
\frac{\Phi_{\nu_\mu}^{\rm sec,BSM}+\Phi_{\bar\nu_\mu}^{\rm sec,BSM}}
{\Phi_{\nu_\mu}^{\rm sec,SM}+\Phi_{\bar\nu_\mu}^{\rm sec,SM}}.
\]

## Detector-level omission identified

The current ARCA projection code treats the selected \(\nu_\mu\) CC track sample as unchanged by tau-only BSM and uses it only as an astrophysical normalization/spectral-index control. That approximation is incomplete once tau regeneration is included.

The detector-level muon-neutrino flux must instead be decomposed schematically as

\[
\Phi_\mu^{\rm det}(E,\theta)
=
\Phi_\mu^{\rm primary}(E,\theta)
+
\Phi_\mu^{\tau\text{-decay}}(E,\theta),
\]

where the second term differs between the SM and BSM hypotheses. Consequently:

1. **ARCA tracks:** the secondary component contributes through \(\nu_\mu N\to\mu X\) and can change the track expectation.
2. **ARCA showers:** the same secondary component also contributes to NC showers and must be included consistently.
3. **IceCube:** any event-level likelihood using muon-flavor acceptance must include the same BSM-modified secondary component rather than treating the muon channel as exactly unchanged.
4. **Nuisance control:** the ARCA track sample cannot be assumed to be a perfectly BSM-blind control sample until the size of this regenerated contribution is explicitly calculated.

## Required implementation before a physics conclusion

After the current transport gate closes, the first detector-level correction is:

- propagate and save species-resolved
  \(\nu_\tau,\bar\nu_\tau,\nu_\mu^{\rm sec},\bar\nu_\mu^{\rm sec}\)
  spectra for SM and BSM benchmarks over the zenith/energy grid used by the detector fold;
- form the total muon-flavor flux entering the detector response, including the primary astrophysical \(\nu_\mu+\bar\nu_\mu\) component and the tau-decay secondary component;
- fold the BSM-modified secondary flux through the ARCA \(\nu_\mu\) CC track response and the appropriate NC shower response;
- apply the corresponding correction in the IceCube comparison;
- rerun the profiled ARCA/IceCube likelihoods before making any viability or no-go statement.

## Interpretation rule

A large value of \(R_\mu^{\rm sec}\) does **not** by itself imply a large track-rate shift, because the secondary component may still be small compared with the primary astrophysical muon-neutrino flux. The physically relevant quantity is the change in the **total detector-weighted muon-flavor event expectation**.

Therefore no statement that "tau-only BSM leaves ARCA tracks unchanged" should appear in the final analysis unless the secondary contribution has first been shown quantitatively negligible after detector folding.

## Publication status

This is a required correction to the detector-level analysis, not a failure of the current transport solver. The transport already contains the relevant secondary channel. The missing step is carrying that channel consistently into ARCA and IceCube event predictions.
