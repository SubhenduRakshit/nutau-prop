# nutau-prop — gravity / LIV / PBH / cosmic-propagation gap register v0.91

Date: 2026-10-31
Status: LIVING SUPPLEMENT TO `UHE_BSM_EXPLORATION_MASTER_v0.90.md` AND `UHE_BSM_TODO_MASTER_v0.90.md`

Purpose: record BSM/gravity-related directions not yet covered by the contact/resonance portfolio, and classify them by whether they affect the source, cosmic propagation, or Earth `nu_tau-tau-nu_tau` transport.

## A. High-priority Earth-interaction / Earth-transport directions

### A1. Low-scale gravity / extra dimensions / KK graviton exchange

Mechanisms:
- ADD/warped extra dimensions;
- virtual KK-graviton exchange in `nu q` / `nu g` scattering;
- low-scale string/quantum-gravity completions;
- possible microscopic/quantum black-hole or string-ball production above threshold.

Why relevant:
- unlike weakly coupled higher-dimensional EFT operators already closed, these can produce large high-energy cross sections through a resolved tower/strong-gravity threshold;
- they directly modify `nu N` attenuation, inelasticity and therefore `nu_tau -> tau -> nu_tau` transport.

Modern constraint caveat:
- old TeV-gravity benchmarks with ~1 TeV gravity scales are obsolete as publication benchmarks;
- current collider searches push semiclassical black-hole/string-ball thresholds to roughly the multi-TeV to ~10 TeV regime depending on model assumptions, and RS/KK-graviton constraints are also multi-TeV and coupling dependent.

TODO:
- [ ] perform a 2026 literature/constraint audit for ADD, RS/warped and low-scale-gravity neutrino scattering;
- [ ] derive/implement current-allowed `d^2 sigma/dx dy` for KK-graviton/eikonal exchange;
- [ ] separately screen microscopic/quantum-black-hole and string-ball production, with threshold and inelasticity distributions;
- [ ] compute the maximum allowed `Delta sigma/sigma_SM` at 100 PeV, 300 PeV and 1 EeV before any Earth transport;
- [ ] if >=10% survives, run full `nu_tau-tau-nu_tau` PREM propagation and emerging-tau detector folds;
- [ ] compare against current CMS/ATLAS black-hole/string-ball/KK-graviton limits and astrophysical bounds;
- [ ] novelty gate: determine whether modern collider-allowed low-scale gravity has already been propagated through full tau regeneration for GRAND/POEMMA/Trinity.

Completion criterion: modern allowed parameter point with detector-relevant effect, or documented closure.

### A2. Spacetime torsion / geometry-induced matter interactions

Recent work shows that non-dynamical torsion can induce an effective parity-violating four-fermion matter interaction, producing a PREM-density-dependent oscillation potential.

Relevance:
- directly uses Earth matter profile;
- qualitatively different from cross-section enhancement;
- likely affects flavour conversion more than absorption at UHE.

TODO:
- [ ] derive the UHE limit of the torsion-induced Hamiltonian and determine whether standard oscillation averaging removes the observable effect at 100 PeV-EeV;
- [ ] translate existing atmospheric-neutrino bounds into the parameter range relevant to this project;
- [ ] only if a non-negligible UHE flavour effect survives, couple it to the attenuation/regeneration transport equation.

Priority: LOW/MEDIUM pending cheap analytic screen.

## B. High-priority cosmic-propagation directions

### B1. Lorentz-invariance violation (LIV) / modified dispersion

Possible signatures:
- energy-dependent velocities;
- neutrino splitting;
- vacuum `e+e-` emission for superluminal neutrinos;
- modified interaction/decay thresholds;
- flavour- and direction-dependent SME propagation.

Current status from literature:
- the ~220 PeV KM3NeT event already gives very strong constraints on superluminal LIV and on quadratic energy-dependent LIV scales;
- 2026 flavour-anisotropy analyses constrain many SME coefficients of dimensions 2-8.

Interpretation for this project:
- not promising as a large unconstrained Earth cross-section modification;
- still valuable as a distinct source-to-Earth spectral/flavour distortion and as a quantum-gravity control.

TODO:
- [ ] screen the surviving LIV coefficient ranges after KM3NeT/IceCube 2025-2026 limits;
- [ ] test whether any allowed coefficient can measurably change the incident 100 PeV-EeV spectrum or flavour composition;
- [ ] if yes, propagate the modified incident flux through the standard Earth solver rather than modifying the `nu N` kernel unless the chosen LIV model explicitly changes scattering;
- [ ] assess whether tau decay length/threshold effects are relevant in models where charged-lepton LIV is correlated.

Priority: MEDIUM as a cosmic-propagation/quantum-gravity branch, LOW as an Earth-cross-section branch.

### B2. Quantum-gravity-induced decoherence / spacetime foam

Mechanisms:
- Lindblad-like loss of coherence;
- stochastic metric fluctuations;
- modified-dispersion-induced wave-packet decoherence.

Current status:
- recent IceCube analyses strongly constrain Planck-suppressed decoherence models, including energy-growing rates.

TODO:
- [ ] use latest IceCube constraints to identify any surviving high-energy scaling;
- [ ] quantify residual effects on flavour composition for astrophysical/cosmogenic neutrinos;
- [ ] compare with all-flavour Earth regeneration observables.

Priority: MEDIUM/LOW unless a surviving parameter region gives a large flavour effect.

### B3. Neutrino-dark-matter scattering

Mechanism:
- attenuation/cascade during propagation through source-host, intergalactic and Milky-Way dark-matter columns.

Current status:
- the KM3-230213A event gives strong limits near ~220 PeV, approximately `sigma_DM-nu/m_DM <~ 1e-22 cm^2 GeV^-1` from the Milky-Way halo alone under the assumptions of the recent analysis;
- simple models with rapidly rising cross sections are often squeezed by unitarity; richer dark sectors can remain interesting.

TODO:
- [ ] screen tau-philic and energy-growing simplified models against the 2026 KM3NeT-based bounds;
- [ ] determine whether any surviving model yields a spectral dip/cascade large enough to matter at Earth;
- [ ] feed the resulting incident spectrum/flavour composition into the Earth transport code;
- [ ] distinguish this source-to-Earth attenuation from CGC/BSM Earth attenuation by its lack of zenith dependence inside Earth.

Priority: HIGH as a complementary propagation discriminator, though not an Earth-interaction model.

### B4. Neutrino self-interactions on the cosmic neutrino background

Mechanisms:
- resonant absorption on CnuB;
- cascades/regeneration;
- secret neutrino interactions.

TODO:
- [ ] audit current 2026 bounds and viable tau-specific mediator windows;
- [ ] identify whether a PeV-EeV dip/cascade survives;
- [ ] use resulting spectrum as incident input to Earth propagation and test degeneracy with source composition and Earth-BSM effects.

Priority: HIGH for cosmogenic/UHE spectral physics, but separate from the Earth-cross-section paper unless used as a nuisance/background model.

### B5. Weak-equivalence-principle violation / flavour-dependent gravity

Mechanisms:
- flavour-dependent Shapiro delay;
- gravitationally induced phase differences;
- possible modified flavour conversion.

Current status:
- high-energy neutrino/blazar timing already constrains WEP violation strongly (recent analyses reach ~1e-6 to 1e-8 depending on event/source assumptions).

TODO:
- [ ] cheap analytic screen of whether any allowed flavour-dependent gravitational potential can alter UHE flavour ratios after oscillation averaging;
- [ ] if negligible, close for the present project.

Priority: LOW.

## C. Source-physics directions

### C1. Primordial black holes (PBHs)

PBHs are primarily a source model, not an Earth-propagation modification.

Relevant mechanisms:
- direct Hawking neutrino emission from PBHs evaporating today;
- early-Universe PBH evaporation followed by redshift/wash-out;
- PBH-emitted dark particles/HNLs that decay later to UHE neutrinos;
- modified evaporation from quantum-gravity effects (memory burden, remnants, etc.);
- higher-dimensional PBHs, whose temperature/lifetime and visible branching fractions differ from 4D PBHs.

Recent literature motivation:
- 2025 work showed direct/indirect PBH mechanisms can populate PeV-EeV neutrino energies, with some indirect scenarios potentially visible at IceCube/GRAND;
- a 2025 PRL proposed late-stage PBH evaporation as an explanation of the highest-energy neutrino events;
- 2025 work on 5D PBHs in extra-dimensional scenarios connects the ~100 PeV scale to a higher-dimensional gravity scale;
- 2026 neutrino data are being used to constrain PBH dark-matter abundance and evaporation scenarios.

Interpretation for this project:
- do NOT mix PBH source physics with Earth `nu_tau` propagation as if it were a new `nu N` interaction;
- PBHs can instead supply an incident UHE spectrum to be propagated through Earth;
- gravitational lensing/capture by PBHs does not constitute the kind of `nu_tau-tau` regeneration effect presently targeted.

TODO:
- [ ] decide whether PBH should be a separate source paper/benchmark rather than part of `nutau-prop`;
- [ ] if retained as a benchmark, import a published PBH neutrino spectrum and propagate it through the existing Earth solver;
- [ ] test whether the harder/peaked PBH spectrum changes the observability of BSM/CGC Earth effects;
- [ ] keep higher-dimensional PBHs linked to the extra-dimension branch but separate source vs interaction physics.

Priority: LOW for the current Earth-propagation paper; HIGH only if the project expands to UHE source phenomenology.

## D. Other new-physics directions worth a cheap screen

- neutrino decay / invisible decay / pseudo-Dirac neutrinos;
- long-range fifth forces or ultralight mediator matter potentials;
- ultralight-dark-matter refractive effects on flavour propagation;
- electroweak sphaleron/instanton production in UHE neutrino-parton collisions;
- low-scale string excitations/string balls distinct from semiclassical black holes;
- heavy sterile-neutrino thresholds/decays beyond the already tested stable-sterile contact models;
- nonstandard tau energy-loss mechanisms if a concrete model alters electromagnetic energy loss or decay.

For each: first perform an analytic/current-bound size gate. Do not run Earth transport unless the surviving effect is plausibly >= few-percent in an experimentally populated region.

## Recommended new priority ordering

1. charged scalar `cs-tau nu` recertification (existing P0);
2. modern low-scale gravity / extra dimensions / microscopic-BH screen;
3. resonant LQ/leptogluon branches;
4. precision CGC / nuclear saturation;
5. neutrino-DM and secret-neutrino cosmic-propagation controls;
6. LIV / quantum-gravity decoherence;
7. torsion / WEP / other gravity-induced flavour effects;
8. PBH source physics as a separate benchmark/project.

## Update rule

When the canonical master/ToDo advances beyond v0.90, merge this register into the next master ToDo rather than leaving these branches orphaned.