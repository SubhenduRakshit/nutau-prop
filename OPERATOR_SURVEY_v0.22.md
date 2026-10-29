# nutau-prop v0.22 — audit reconciliation and active d=8 UHE branch

Date: 2026-10-29

## Why this note exists

Audit(4) is internally consistent with the corrected partonic EFT-validity logic, but its work plan still belongs to the earlier ATZ/nuLEFT branch. Since then the project identified a separate and potentially better high-energy branch: purely active, gauge-invariant dimension-8 semileptonic twist-2 / spin-2 operators, plus an analogous gluon twist-2 possibility. This note prevents the two branches from being mixed.

## What Audit(4) correctly preserves

1. The old no-go based on E_nu < Lambda^2/(2 m_N) is withdrawn. The correct DIS validity condition is partonic, shat = 2 x m_N E_nu < Lambda^2.
2. The PKR/ATZ splice remains real: the draft wrote the charge-conjugated PKR operator but used an ATZ cross section belonging to a different operator.
3. The ATZ factor-64 / y(y-1) normalization-sign problem remains unresolved and must not be inherited silently.
4. The earlier all-flavour effective-area rescaling is not a detector prediction.
5. The v0.20.9 transport/grid-convergence work remains useful.
6. A proper energy-resolved recast of IceCube's 60 TeV-10 PeV cross-section information is a high-value constraint and should replace the crude single-energy proxy.

## ATZ branch: keep, but no longer the unique main line

Audit(4) proposes the clean ATZ families Q3,4; Q8,9; Q10,11 and reports a first-pass saturated screen in which all three families clear the few-percent UHE gate, with partial-wave unitarity binding Q8-Q11 and IceCube binding Q3/Q4.

This is worth preserving as a separate branch, but the quoted screen must be independently reproduced in code before being treated as a project result. For the Dirac ATZ derivative operators the outgoing neutrino is right-handed/sterile, so the BSM channel is a pure sink. That is specific to this branch.

## New preferred branch: active d=8 two-derivative semileptonic operators

The project should now prioritize gauge-invariant SMEFT operators of the schematic form

  (Lbar_3 gamma_mu i<->D_nu L_3)(qbar gamma^mu i<->D^nu q) / Lambda^4

with q chosen first as c_R and then s_R.

Primary benchmark:

  [O_{8,lu partial2}]_{3322}

with L_3 = (nu_tau, tau)_L and u_2 = c_R.

Physics advantages:

- fully SU(3)c x SU(2)L x U(1)Y invariant;
- keeps nu_tau -> nu_tau active;
- lepton number conserved;
- no sterile daughter and no PKR flavour-antisymmetry issue;
- amplitude grows as shat^2/Lambda^4, stronger than the old ATZ/PKR derivative case;
- heavy-sea-quark choice exploits small-x UHE PDFs while suppressing high-mass LHC parton luminosity;
- can arise naturally from heavy spin-2 exchange without requiring a dominant dimension-6 interaction.

The strange-right-handed analogue should be screened in parallel.

## Gluon branch

Two gluon-related mechanisms must be kept distinct.

1. Old scalar/gluonic neutrino operators such as ATZ/PKR-type (nu bilinear) G^2 interactions. These are established in the literature and already contributed directly to nu-N DIS in previous studies.
2. A new active spin-2/twist-2 gluonic structure, schematically

  (Lbar_3 gamma_{ {mu} } i<->D_{nu} L_3)
  [G^{a mu rho} G^a_{ nu rho} - trace]

or its canonical SMEFT representative.

The second option is especially interesting at UHE because the small-x gluon PDF is much larger than the charm or strange sea, but the same gluon luminosity can strengthen LHC constraints. It must therefore be compared on a bound-saturated basis rather than by raw cross section.

## Candidate hierarchy to test

A. active charm twist-2: best chance of weak collider bounds;
B. active strange twist-2: larger PDF but stronger LHC reach;
C. active gluon twist-2: largest raw UHE luminosity, probably strongest collider pressure;
D. ATZ Q8-Q11 sterile-daughter branch: retain as comparison;
E. ATZ/PKR scalar-gluon branch: retain as historical/robustness comparison;
F. Z' branch: secondary resolved-mediator benchmark, not presently the main momentum-growth mechanism.

## Constraints that must be applied to every active d=8 candidate

1. exact gauge-invariant operator normalization;
2. SM interference and sign/positivity constraints;
3. partial-wave unitarity;
4. EFT validity with partonic cuts shat < xi Lambda^2 for xi = 1, 1/2, 1/4;
5. high-pT tau+tau- Drell-Yan recast, preferably with HighPT where supported;
6. monojet + MET recast;
7. IceCube energy-resolved cross-section / attenuation constraints from 60 TeV-10 PeV;
8. IceCube EHE and KM3NeT UHE constraints;
9. direct nu_tau bounds, though derivative suppression is expected to make GeV constraints weak;
10. QCD evolution / PDF uncertainties, especially charm/strange/gluon at small x.

## Immediate calculation order

1. Derive exact d^2 sigma/dx dy for nu_tau c -> nu_tau c for [O_{8,lu partial2}]_{3322}, including SM interference.
2. Recast the LHC tau+tau- tail for the same flavor assignment.
3. Screen monojet + MET.
4. Compute the bound-saturated UHE PDF fold at 1, 10, 100, 300 PeV with xi = 1, 1/2, 1/4.
5. Repeat for the strange analogue.
6. Identify the canonical active gluon twist-2 SMEFT operator and repeat the same constraint/fold exercise.
7. Only operators with robust r_tau above a few percent enter the validated Earth-propagation solver.

## Important branch distinction for transport

ATZ derivative branch:
  nu_tau,L -> nu_R, so BSM = active-flux sink with no BSM regeneration.

Active d=8 twist-2 branch:
  nu_tau,L -> nu_tau,L, so BSM modifies a neutral-current-like active channel and must include an active daughter redistribution kernel.

Do not use the pure-sink transport conclusions from Audit(4) for the active twist-2 branch.

## Current project decision

The main UHE search is now:

  active d=8 charm / strange / gluon twist-2 operators
  -> all phenomenological constraints
  -> bound-saturated UHE cross-section screen
  -> Earth propagation only for survivors.

ATZ Q3,4 and Q8-Q11 remain useful comparison branches, but they are no longer the sole starting point.
