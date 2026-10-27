# nutau-prop v0.15 — contact-operator transport findings

> **AUDIT STATUS — PROVISIONAL / SUPERSEDED FOR INTERPRETATION.**  The first v0.15 transport run used different numerical substep schedules for the matched SM and BSM solutions, and the helicity solver accumulated tau-decay secondary muon neutrinos without propagating their subsequent SM interactions.  Several saved Dirac-sink ratios above unity are therefore a hard numerical warning.  Do not use the transport conclusions below in the manuscript or as a physics decision gate until the corrected v0.15.1 rerun passes.  See `AUDIT_CORRECTION_v0.15.1.md`.  The v0.14 exact-PDF hard-kernel results remain valid.

Date: 2026-10-27

This file records the quantitative findings from the completed Cases 5, 6, and 7 momentum-enhanced contact-operator study. It supplements `OPERATOR_SURVEY_v0.15.md`, which contains the operator definitions, exact-PDF prescreen, external-state/chirality interpretation, and decision rule.

## Common physics question

The project asks whether a momentum-enhanced BSM neutrino interaction can increase detector sensitivity without making the Earth so opaque that the arriving UHE neutrino flux is destroyed or degraded out of the signal energy bin.

The three contact operators tested here probe distinct possibilities:

- Case 5: derivative tensor quark operator, flavor-off-diagonal in the Majorana interpretation;
- Case 6: derivative vector quark operator, the strongest contact benchmark studied;
- Case 7: neutrino-gluon Rayleigh operator with a hard y distribution.

The EFT condition `sqrt(shat) < Lambda_UV` remains only a validity boundary. It is not a physical high-energy turnover.

## Case 5 — derivative tensor quark operator

### Exact-PDF prescreen

At the published benchmark coefficient and the perturbative endpoint `c/(4 pi)=1`,

- `Lambda_UV = 4.1176 TeV`;
- `r_NP/SM = 0.0552, 0.0678, 0.0732` at `100, 220, 300 PeV`;
- `<y>_NP ~ 0.50`.

At equal cutoff and benchmark coefficients, the Case-5 NP cross section is exactly `0.0219037` times the Case-6 NP cross section, because Case 5 has the same `y(1-y)` shape and half the matrix-element normalization at equal coefficient.

### External-state brackets

Two physical interpretations were transported:

1. Dirac sink: the chirality-flipped right-handed daughter is inactive to ordinary weak detection to leading order.
2. Majorana tau-mu transition: `nu_tau <-> antinu_mu` and `antinu_tau <-> nu_mu`, including the reverse transitions required by the off-diagonal Hermitian interaction.

The Majorana case is therefore not treated as an artificial one-way `tau -> mu` escape.

### Transport findings

**The numbers in this subsection are retained only as the historical first v0.15 output and are not validated physics results.**

At `Lambda_UV = 3 TeV`, vertical propagation (`|cos z|=1`) gives approximately:

- Dirac sink tau ratios at `1, 10, 100, 220, 300 PeV`: `0.983, 0.808, 0.317, 0.164, 0.116`;
- Majorana tau-mu tau ratios: `1.015, 0.861, 0.376, 0.195, 0.138`;
- Majorana active-pair ratios at the same energies: `1.017, 0.891, 0.380, 0.196, 0.139`.

At the perturbative endpoint `Lambda_UV = 4.1176 TeV`, vertical propagation gives:

- Dirac sink tau ratios: `0.987, 0.748, 0.0649, 0.0120, 0.00492`;
- Majorana tau-mu tau ratios: `1.057, 0.855, 0.0965, 0.0182, 0.00743`;
- Majorana active-pair ratios: `1.056, 0.901, 0.0985, 0.0183, 0.00744`.

The original interpretation that active flavor conversion helps modestly but does not prevent severe UHE attenuation is **withdrawn pending the corrected v0.15.1 rerun and detector weighting**.

## Case 6 — derivative vector quark operator

Case 6 remains the strongest contact benchmark and already passed the separate v0.11-v0.13 differential transport and detector-response program.

The first v0.15 external-state audit tested:

1. Dirac active-to-right-handed sink;
2. Majorana `nu_tau <-> antinu_tau` helicity flip.

**The v0.15 chirality-audit ratios below are provisional because of the matched-step failure.**

Representative vertical ratios at `220 PeV` from that first run were:

- `Lambda=1.5 TeV`: Dirac `2.93e-5`, Majorana `6.84e-5`;
- `Lambda=2.0 TeV`: Dirac `1.68e-11`, Majorana `9.03e-11`;
- `Lambda=2.44 TeV`: Dirac `5.15e-20`, Majorana `6.23e-19`;
- `Lambda=3.0 TeV`: Dirac `1.31e-35`, Majorana `4.79e-34`.

For `|cos z|=0.5`, the corresponding first-run Majorana 220-PeV ratios were `0.152, 8.97e-3, 1.45e-4, 5.43e-8` for `Lambda=1.5, 2.0, 2.44, 3.0 TeV`.

The older v0.12/v0.13 Case-6 transport/detector calculation uses separate machinery and is not automatically invalidated by this v0.15 solver bug.  However, the claim that the new chirality audit strengthens the Case-6 conclusion must be revalidated with matched steps.

## Case 7 — neutrino-gluon Rayleigh operator

### Exact-PDF prescreen

At the perturbative endpoint `Lambda_UV = 4.2825 TeV`,

- `r_NP/SM = 0.0311, 0.0538, 0.0642` at `100, 220, 300 PeV`;
- `<y>_NP ~ 0.80`.

The gluon PDF helps the total cross section, but the `y^3` kernel makes each interaction much more destructive to the parent neutrino energy than Cases 5 or 6.

### External-state brackets

Two transport interpretations were tested:

1. Dirac sink;
2. Majorana `nu_tau <-> antinu_tau` flip.

### Transport findings

**The numbers in this subsection are retained only as historical first-run output.  The secondary-muon ratios are specifically invalid because the daughters were not propagated after production.**

At `Lambda_UV = 3 TeV`, vertical tau ratios at `1, 10, 100, 220, 300 PeV` were approximately:

- Dirac sink: `1.016, 0.947, 0.433, 0.196, 0.126`;
- Majorana flip: `1.041, 0.980, 0.451, 0.201, 0.129`.

At the perturbative endpoint `Lambda_UV = 4.2825 TeV`:

- Dirac sink: `1.062, 0.944, 0.196, 0.0274, 0.00872`;
- Majorana flip: `1.115, 1.010, 0.214, 0.0291, 0.00916`.

The saved secondary-muon ratios from this first run must not be used.

## Combined contact-operator conclusion

The previous combined conclusion that Cases 5/6/7 establish a robust negative result for momentum-enhanced contact operators is **withdrawn pending corrected transport and detector weighting**.

What is still established is narrower:

- the v0.14 Case-5/Case-7 hard kernels and PDF prescreen are internally consistent;
- the operators can generate percent-level to large cross-section modifications in the target energy range;
- the physical competition between detector interaction probability and Earth propagation loss must be evaluated with validated transport and detector weights;
- ratios in already-opaque vertical bins are not by themselves a detector-level exclusion.

## Decision

Pause phenomenological interpretation of v0.15 and any v0.16 transport run.  Complete the v0.15.1 correction gate in `AUDIT_CORRECTION_v0.15.1.md` first.  Resolved-mediator v0.16 cross-section tables may be retained as hard-kernel work, but no positive/negative model conclusion should be drawn from them until the corrected transport baseline is established.
