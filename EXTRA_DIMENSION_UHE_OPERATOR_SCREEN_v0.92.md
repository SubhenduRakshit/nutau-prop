# nutau-prop v0.92 — extra-dimensional / low-scale-gravity UHE operator screen

Date: 2026-11-01
Status: ANALYTIC TRIAGE; no exact LHAPDF fold yet
Purpose: identify extra-dimensional operators and related low-scale-gravity channels worth carrying into the UHE nu_tau Earth-propagation program.

## 1. Scope

Separate controlled effective operators from genuinely trans-Planckian processes. The relevant classes are:

1. ADD virtual KK-graviton exchange -> universal dimension-8 stress-tensor operator;
2. graviton-loop dimension-6 universal axial-current operator;
3. resolved Randall-Sundrum KK spin-2 exchange;
4. real KK-graviton emission;
5. bulk-sterile-neutrino KK-tower up-scattering;
6. trans-Planckian microscopic/quantum black holes and string balls (not local operators);
7. ordinary LED active-sterile oscillations (control only for the present UHE project);
8. torsion-induced four-fermion operators in extra-dimensional gravity (secondary branch).

## 2. ADD tree-level virtual graviton exchange: d=8

PDG convention:

L8 = +/- 4/M_TT^4 [ T_{mu nu} T^{mu nu} - 1/(delta+2) T^mu_mu T^nu_nu ].

For effectively massless neutrino external states, the trace contribution is absent at the neutrino vertex. The hard amplitude scales schematically as

M ~ c_T shat^2,

and therefore

sigmahat ~ c_T^2 shat^3,

up to spin/angular factors and interference with the SM neutral current.

Current PDG 2025 summary quotes the strongest virtual-graviton limit around

M_TT > 9.02 TeV (95% CL),

from dijet angular distributions.

At M_TT=9.02 TeV,

4/M_TT^4 = 6.04e-4 TeV^-4.

For comparison, our already-tested strict d=8 coefficients were

(C_q/Lambda^4) = (8pi/3)/(4 TeV)^4 = 3.27e-2 TeV^-4,

(C_g/Lambda^4) = (8pi)/(4 TeV)^4 = 9.82e-2 TeV^-4.

Thus the allowed ADD contact coefficient is smaller in amplitude by factors about 54 and 162, respectively. Pure d=8 contributions therefore carry an additional coefficient-squared suppression of about 2.9e3 to 2.6e4 relative to those benchmarks before angular/PDF differences.

Interpretation: the controlled ADD contact operator is unlikely to rescue the several-percent UHE goal. However the exact SM-interference term with light quarks is linear in the coefficient and must be folded before declaring closure.

Status: QUICK EXACT FOLD REQUIRED; likely sub-percent in its controlled domain.

Critical validity condition: do not extrapolate the contact operator into shat >= M_TT^2. The interesting high-energy region is precisely where the EFT ceases to be predictive.

## 3. Graviton-loop d=6 universal axial operator

Extra-dimensional graviton loops can generate

L_Upsilon = +/- (4 pi / Lambda_gamma^2) Upsilon,

Upsilon = 1/2 (sum_f fbar gamma_mu gamma5 f)^2.

This contains universal nu-q axial contact terms and can interfere with the SM. Historical LEP limits were Lambda_gamma ~15-17 TeV; a modern current contact-interaction recast is required before numerical use.

Unlike the d=8 tree operator, the partonic rate only grows linearly with shat, but the light-quark universality and SM interference can compensate somewhat. A naive low-Q comparison gives an amplitude well below the weak interaction for Lambda_gamma ~15 TeV, so an O(10%) inclusive UHE effect is not expected, but this must be checked with the UHE Q^2 distribution.

Status: MEDIUM-LOW priority; exact dynamic-Q fold plus current compositeness bound is cheap and should close or retain it quickly.

Important caveat: its coefficient is UV/cutoff dependent; it is not a model-independent ADD prediction.

## 4. Resolved Randall-Sundrum KK graviton

A KK spin-2 state couples schematically as

L ~ - G^{mu nu} T_{mu nu}/Lambda_pi.

For nu q -> nu q the exchange is t-channel, so

M ~ shat^2 / [Lambda_pi^2 (Q^2 + m_G^2)].

This is the resolved counterpart of the contact stress-tensor operator. Current CMS Run-3 dilepton limits exclude, for the standard RS benchmark,

m_G < 2.5, 4.4, 5.0 TeV

for k/Mbar_Pl = 0.01, 0.05, 0.1 respectively.

For the strongest-coupled standard benchmark, the low-Q effective coefficient is already of order a few 1e-4 TeV^-4, comparable to or below the allowed ADD contact coefficient; increasing Q^2 further suppresses t-channel exchange. Multiple KK modes and bulk-localization variants must be treated explicitly if revisited.

Status: LOW-MEDIUM priority for the standard universal RS benchmark; not yet closed for non-universal/bulk-localized variants.

## 5. Real KK-graviton emission

Processes with an escaping KK graviton are qualitatively different from virtual exchange. Summing the KK tower converts the Planck-suppressed single-mode rate into powers of the fundamental scale M_D. Collider monojet/monophoton searches directly constrain this branch.

Current published limits depend strongly on the number of extra dimensions. PDG summaries quote monojet limits roughly M_D > 5.9-11.2 TeV (ATLAS) and 5.5-10.7 TeV (CMS) for typical delta ranges; the 2025 CMS monophoton analysis reaches M_D up to about 3.2 TeV for n=3-6 but is not the globally strongest channel.

For UHE neutrinos, real-graviton emission can remove a significant fraction of energy while leaving a neutrino or charged lepton in the final state, so its y-distribution and active/sterile bookkeeping differ from ordinary NC scattering.

Status: MEDIUM-HIGH priority. Need exact partonic channels and current M_D(delta) bounds before Earth propagation.

## 6. Bulk sterile-neutrino KK tower: potentially important loophole

A gauge-singlet neutrino propagating in an extra dimension appears in 4D as a tower N_n. Existing LED work has shown that chirality-flipping dipole or scalar up-scattering can populate many heavy KK modes; published applications have focused mainly on CEvNS rather than UHE DIS.

For our project this is especially interesting because several single-sterile operators already failed on size. The tower sum can change the conclusion:

nu_tau + q -> N_n + q,

with all modes satisfying m_n < sqrt(shat) included. Schematically,

sigma_total = sum_{n accessible} sigma(nu q -> N_n q).

The multiplicity grows with the density of KK states, but so do perturbativity, cosmology, oscillation, invisible-width and unitarity constraints. One must not multiply a single-state cross section by an unrestricted number of KK modes.

Novelty target: no UHE nu-N / Earth-propagation calculation of this tower-enhanced up-scattering was identified in the first literature pass. This requires a dedicated search before any novelty claim.

Status: HIGH priority among genuinely new operator mechanisms.

Required first calculation:

1. choose one gauge-invariant parent (dipole and scalar separately);
2. derive mode-dependent couplings and density of states;
3. impose current LED/bulk-neutrino limits;
4. sum only kinematically accessible modes;
5. test perturbative unitarity;
6. exact PDF fold at 100 PeV, 300 PeV, 1 EeV;
7. if >5-10%, implement bidirectional/decay transport as required.

## 7. Trans-Planckian black holes / string balls: not operators, but high-priority gravity channel

Once sqrt(shat) approaches the fundamental gravity scale, local operators are not trustworthy. Microscopic black-hole/string-ball production is a distinct regime.

Current 2026 CMS limits exclude semiclassical black holes roughly below 8.4-11.4 TeV and string balls below 9.0-10.7 TeV depending on the model. ATLAS Run-3 quantum-black-hole searches reach about 9.4 TeV in lepton+jet channels.

Kinematic relevance for UHE neutrinos can be seen from

x_min = M_thr^2 / (2 m_N E_nu).

Representative thresholds:

| E_nu | x_min(M=5.5 TeV) | x_min(9 TeV) | x_min(11 TeV) |
|---|---:|---:|---:|
| 100 PeV | 0.161 | 0.431 | 0.644 |
| 300 PeV | 0.0537 | 0.144 | 0.215 |
| 1 EeV | 0.0161 | 0.0431 | 0.0644 |

Therefore modern collider-allowed thresholds are strongly PDF-suppressed at 100 PeV but become kinematically important at 300 PeV-EeV. Such events are high-inelasticity/absorptive and can strongly modify Earth opacity and emerging-tau probabilities even if they contribute only in the high-energy tail.

Status: HIGH priority, especially 300 PeV-EeV. Must use current collider-allowed threshold parameters and model-dependent formation factors; old 1-TeV-gravity benchmarks are obsolete.

## 8. Ordinary LED oscillation / KK mixing

The KK sterile tower can also modify oscillations/matter effects without changing the hard DIS cross section. IceCube's 2026 high-energy atmospheric-neutrino analysis finds no evidence and constrains the largest extra dimension to below about 0.17 micrometers in the benchmark studied.

At 100 PeV-EeV this is not the primary mechanism for our Earth-attenuation project; the relevant oscillation phases/matter resonances belong mainly to much lower energies.

Status: CONTROL / LOW priority for current UHE transport objective.

## 9. Extra-dimensional torsion

Gravity with torsion in extra dimensions can induce four-fermion contact operators at the fundamental gravity scale. After eliminating non-propagating torsion, vector/axial structures involving neutrinos and quarks can result.

This is phenomenologically closer to a dimension-6 NSI than to KK-graviton exchange. Precision electroweak/contact-interaction constraints are likely binding, but a tau-specific or non-minimal torsion realization may differ from universal benchmarks.

Status: SECONDARY; first derive the exact nu-q operator and map modern constraints before any UHE fold.

## 10. Priority decision

### High priority

1. bulk-sterile KK-tower up-scattering (dipole/scalar parents);
2. black-hole / string-ball / quantum-black-hole production at 300 PeV-EeV;
3. real KK-graviton emission with current M_D(delta) limits.

### Quick closure calculations

4. ADD virtual-graviton d=8 stress-tensor operator with exact SM interference;
5. loop-induced d=6 universal axial operator;
6. standard RS KK-graviton t-channel exchange.

### Low priority / controls

7. ordinary LED oscillation mixing for the 100 PeV-EeV project;
8. extra-dimensional torsion unless a weakly constrained tau-specific realization is found.

## 11. Immediate numerical gates

Before full Earth propagation, require for each branch:

- exact current constraint mapping;
- exact d sigma/dy or d2 sigma/dxdy;
- >=5% inclusive effect or a clearly stronger high-y/sterile-sink attenuation effect at 100 PeV, 300 PeV or 1 EeV;
- controlled validity / unitarity;
- no use of an EFT beyond its cutoff.

Only survivors proceed to the existing nu_tau-tau-nu_tau PREM solver and GRAND/POEMMA/Trinity emerging-tau fold.
