# nutau-prop v0.15 — contact-operator transport findings

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

At `Lambda_UV = 3 TeV`, vertical propagation (`|cos z|=1`) gives approximately:

- Dirac sink tau ratios at `1, 10, 100, 220, 300 PeV`: `0.983, 0.808, 0.317, 0.164, 0.116`;
- Majorana tau-mu tau ratios: `1.015, 0.861, 0.376, 0.195, 0.138`;
- Majorana active-pair ratios at the same energies: `1.017, 0.891, 0.380, 0.196, 0.139`.

At the perturbative endpoint `Lambda_UV = 4.1176 TeV`, vertical propagation gives:

- Dirac sink tau ratios: `0.987, 0.748, 0.0649, 0.0120, 0.00492`;
- Majorana tau-mu tau ratios: `1.057, 0.855, 0.0965, 0.0182, 0.00743`;
- Majorana active-pair ratios: `1.056, 0.901, 0.0985, 0.0183, 0.00744`.

Thus active flavor conversion helps modestly relative to the Dirac sink, especially at intermediate energy, but it does not prevent severe UHE attenuation. At 220 PeV for a vertical chord, even the track-favorable Majorana active-pair bracket retains only about `1.8%` of the matched SM active tau-plus-muon flux.

### Case-5 conclusion

Case 5 does not rescue the contact-operator route. Its special flavor-conversion structure softens the suppression somewhat, but repeated propagation plus energy degradation still wins over any detector-side benefit at UHE.

## Case 6 — derivative vector quark operator

Case 6 remains the strongest contact benchmark and already passed the full v0.11-v0.13 differential transport and detector-response program.

The v0.15 external-state audit tested:

1. Dirac active-to-right-handed sink;
2. Majorana `nu_tau <-> antinu_tau` helicity flip.

The Majorana-flip results reproduce the previous same-active-species transport to very high accuracy, while the Dirac-sink interpretation is even more suppressive.

Representative vertical ratios at `220 PeV` are:

- `Lambda=1.5 TeV`: Dirac `2.93e-5`, Majorana `6.84e-5`;
- `Lambda=2.0 TeV`: Dirac `1.68e-11`, Majorana `9.03e-11`;
- `Lambda=2.44 TeV`: Dirac `5.15e-20`, Majorana `6.23e-19`;
- `Lambda=3.0 TeV`: Dirac `1.31e-35`, Majorana `4.79e-34`.

For `|cos z|=0.5`, the corresponding Majorana 220-PeV ratios are `0.152, 8.97e-3, 1.45e-4, 5.43e-8` for `Lambda=1.5, 2.0, 2.44, 3.0 TeV`.

Therefore the Case-6 conclusion is robust to the external-state interpretation: the interaction becomes strong enough to be distinctive, but Earth opacity and energy degradation become overwhelming at UHE. The Dirac interpretation strengthens, rather than weakens, this conclusion.

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

At `Lambda_UV = 3 TeV`, vertical tau ratios at `1, 10, 100, 220, 300 PeV` are approximately:

- Dirac sink: `1.016, 0.947, 0.433, 0.196, 0.126`;
- Majorana flip: `1.041, 0.980, 0.451, 0.201, 0.129`.

At the perturbative endpoint `Lambda_UV = 4.2825 TeV`:

- Dirac sink: `1.062, 0.944, 0.196, 0.0274, 0.00872`;
- Majorana flip: `1.115, 1.010, 0.214, 0.0291, 0.00916`.

The secondary-muon ratios remain close to unity rather than developing a compensating enhancement. For example, at the endpoint and a vertical chord they are approximately `1.006, 0.973, 0.968, 0.981, 0.987` in the Dirac bracket and `1.053, 1.008, 0.979, 0.987, 0.993` in the Majorana bracket across `1, 10, 100, 220, 300 PeV`.

### Case-7 conclusion

Case 7 also fails as a positive contact-operator signal. Its high inelasticity makes an NP interaction especially effective at removing neutrinos from the UHE signal bin, while no large compensating secondary-track enhancement appears.

## Combined contact-operator conclusion

The three contact operators establish a common structural result:

- smooth interactions studied earlier were too weak or too nuisance-degenerate;
- momentum-enhanced contact interactions can become large enough to be observable in principle;
- once they become large enough, Earth propagation becomes self-defeating through absorption/energy degradation;
- flavor conversion in Case 5 helps, but not nearly enough;
- Case 7's hard inelasticity makes its propagation loss especially severe;
- the Case-6 conclusion survives the Dirac/Majorana external-state audit.

This is not claimed as a theorem covering all BSM neutrino interactions. It is a robust negative result for the contact-operator classes explicitly tested.

## Decision

Stop broadening the contact-operator survey. The next physics target is a resolved mediator that preserves momentum enhancement in the interesting energy range while introducing a real propagator suppression at high Q. This is the motivation for v0.16.
