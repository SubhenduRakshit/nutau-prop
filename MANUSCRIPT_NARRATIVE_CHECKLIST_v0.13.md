# Manuscript narrative checklist — v0.13

This file is a preservation note. Do not delete the earlier failed/less-promising interaction cases from the manuscript history. The final paper must explain the physics progression, not only the successful/latest calculation.

## Canonical project goal
Find a physically viable BSM nu_tau interaction that produces an observable Earth-propagation/regeneration signature while remaining consistent with UHE neutrino constraints. This is not a no-go project and is separate from the cosmogenic-neutrino/UHECR proton-fraction project.

## Narrative that must be preserved in the manuscript

### 1. Historical transport validation
- Coupled nu_tau-tau transport validated against Rakshit-Reya and Dutta-Reno-Sarcevic benchmarks.
- Preserve the role of regeneration, NC energy degradation, tau decay, and secondary nu_mu.

### 2. Earlier smooth interactions
The manuscript must give these cases a proper section, not only a passing sentence. For each case include:
- operator / interaction structure;
- energy and Q^2 dependence of the cross section;
- the physical effect expected before running the calculation;
- propagation result;
- detector-fold result;
- nuisance profiling result;
- the specific reason it failed.

The key conclusion to preserve is:

    smooth interaction -> modest/smooth propagation distortion -> too weak and/or degenerate with flux normalization/spectral tilt, with some large fractional effects occurring only in bins with negligible absolute flux.

These failures motivated the search for genuinely momentum-growing interactions. Do not rewrite them as a universal no-go theorem.

### 3. Case 6: derivative vector dimension-7 contact interaction
Preserve the distinction from the earlier failures.
- |M|^2 grows as C4^2 shat^3 y(1-y): the interaction is strong enough to generate a large, distinctive energy/zenith effect.
- C4 is held fixed at the benchmark value; Lambda_UV in the calculation is only an EFT-validity boundary sqrt(shat)<Lambda_UV.
- Increasing Lambda_UV admits more high-shat phase space and strengthens the opacity.
- The hard EFT cut is NOT a physical turnover.
- c4 = C4 Lambda_UV^3 is a dimensionless-coupling/perturbativity diagnostic, not the reason the cross section turns over.
- Full differential transport gives a PeV regeneration enhancement followed by strong tens-to-hundreds-of-PeV opacity.
- The joint propagation+detector v0.13 fold shows Earth opacity wins over detector-vertex enhancement at UHE.
- Even the favorable public-response envelope requires the Ehlert benchmark flux to be reduced to about 0.31, and the corresponding KM3NeT bright-track expectation remains only about 0.07 events.
- Therefore contact Case 6 succeeds at producing strong propagation physics but fails as the desired viable positive-signal model because the momentum growth becomes self-defeating at UHE.

### 4. Why Case 5 is next
- Case 5 is momentum-growing but flavor off-diagonal (i != j).
- A tau-initiated interaction can produce another neutrino flavor, so the daughter need not remain subject to repeated tau-specific opacity.
- A tau->mu channel can indirectly feed a muon-track signal after propagation.
- This is a qualitatively different loophole from Case 6 and must be computed explicitly rather than inferred by rescaling Case 6.

### 5. Why a resolved mediator is next after the contact-operator survey
- A resolved derivative Z' or related mediator can keep momentum enhancement in the useful propagation window while physically softening the high-Q amplitude through the propagator.
- This is the correct way to study a physical turnover; do not use the hard EFT cutoff as a stand-in for UV completion.

## Physics progression to state explicitly

    Earlier smooth cases:
        too weak / too nuisance-degenerate
            -> need stronger energy dependence

    Case 6 contact:
        strong enough and distinctive
        but Earth energy degradation/opacity wins at UHE
            -> need a mechanism that avoids repeated opacity and/or softens high-Q growth

    Case 5:
        test flavor-conversion escape from repeated tau opacity

    Resolved mediator:
        test physical high-Q turnover

## Language discipline
- Do not call the earlier smooth-case result a general no-go theorem.
- Do not call the Case-6 hard sqrt(shat) cut a physical cutoff/turnover mechanism.
- Prefer 'energy degradation/opacity' over simple 'absorption' for NC-like Case 6, since the neutrino is usually redistributed to lower energy rather than destroyed.
- Distinguish detector interaction enhancement from Earth-propagation suppression: the observable event rate is a competition between arriving flux and detector interaction probability.
- Do not promote public effective-area folds to collaboration likelihoods; retain the public-y-response limitation explicitly.

## Required manuscript revision before submission
Promote the earlier smooth interactions to a dedicated section alongside Case 6, with enough detail that the reader sees why each failed case motivated the next step. The final manuscript should read as a physics progression, not as a list of benchmarks.
