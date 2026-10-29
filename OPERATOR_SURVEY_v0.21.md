# nutau-prop v0.21 — UHE operator strategy after PKR/ATZ audit

Date: 2026-10-29

Status: canonical research-history note. This file does **not** delete or overwrite the earlier v0.13–v0.15 surveys. It records why the project direction changed after the v0.20.9 convergence gate and the subsequent operator/EFT audits.

## 1. Numerical status carried forward

The v0.20.9 coupled transport convergence gate completed successfully:

- `ALL_PASS=True`
- `NEEDS_FURTHER_REFINEMENT=False`

The PREM + nu_tau/tau transport machinery, tau loss/decay treatment, and spatial/tau-grid convergence work are therefore retained. The current problem is theoretical input, not numerical convergence.

## 2. The PKR/ATZ splice that invalidated the old Case-6 interpretation

The draft had effectively combined two different interactions.

### Literal PKR Eq. (4.7)

\[
O_B=(\bar\nu_i^c\,i\overleftrightarrow{\partial_\mu}\nu_j)(\bar q\gamma^\mu q).
\]

The charge-conjugated derivative bilinear is antisymmetric in flavour,

\[
B^\mu_{ij}=-B^\mu_{ji},\qquad B^\mu_{ii}=0.
\]

Therefore a diagonal tau-tau interpretation of the literal PKR operator is not valid.

### ATZ derivative operator

Altmannshofer–Tammaro–Zupan (ATZ, arXiv:1812.02778) instead use the non-charge-conjugated derivative interaction

\[
O_A=(\bar\nu_\beta i\overleftrightarrow{\partial_\mu}P_L\nu_\alpha)(\bar q\gamma^\mu q).
\]

Diagonal flavour is allowed. For Dirac neutrinos the dimension-7 operator converts an incoming active LH neutrino into an RH neutrino. Thus its Earth-transport interpretation is an active-flux sink, not an NC-like active daughter.

The previous manuscript combined the written PKR operator with an ATZ-like cross section. Those results are therefore preserved only as historical numerical diagnostics.

## 3. ATZ normalization issue remains open

ATZ Eq. (77) prints a coefficient proportional to

\[
64|\widehat C_8^{(7)}|^2\hat s^3 y(y-1),
\]

with a sign that is negative over the physical 0<y<1 range. Independent trace checks reproduced the \(\hat s^3 y(1-y)\) shape but not the published overall normalization under the natural chiral/spin-average convention.

The earlier suggestion that the hat on \(\widehat C\) might absorb this discrepancy is withdrawn. ATZ define

\[
\widehat C_a^{(d)}=C_a^{(d)}/\Lambda^{d-4},
\]

so the hat is simply the dimensionful Wilson coefficient.

Conclusion: any operator retained for publication must have its amplitude and normalization rederived independently.

## 4. Retraction of the false UHE EFT no-go

An intermediate audit incorrectly imposed

\[
E_\nu<\Lambda^2/(2m_N)
\]

as the contact-EFT validity condition.

For DIS the correct invariant is partonic,

\[
\hat s=2xm_NE_\nu,
\]

so the validity restriction is

\[
\hat s<\Lambda^2
\quad\Longleftrightarrow\quad
x<x_{\max}=\Lambda^2/(2m_NE_\nu).
\]

At UHE, substantial small-x phase space survives. Therefore the previous claim that the contact programme is automatically closed above the PeV scale is retracted.

For future work the EFT-domain uncertainty should be quantified with several cuts, e.g.

\[
\hat s<\xi\Lambda^2,\qquad \xi=1,\;1/2,\;1/4.
\]

A signal that disappears when \(\xi\) is tightened is not robust.

## 5. Why a Z' mediator was explored earlier

The resolved-Z' line was pursued because the contact-operator calculation appeared to be dominated near the EFT cutoff. A resolved mediator provides a physical propagator,

\[
1/(Q^2+M_{Z'}^2),
\]

which automatically turns over at high momentum transfer and therefore avoids extrapolating a contact interaction beyond its domain.

That work was useful, but a conventional vector Z' is not now the preferred route for the original momentum-growth goal. At high Q the propagator removes the apparent energy growth, whereas the new dimension-8 derivative/twist-2 candidates have genuine amplitudes that grow with powers of \(\hat s\).

The Z' programme is therefore retained as a resolved-mediator comparison/robustness branch, not as the primary operator benchmark.

## 6. Correct high-energy framework: dimension-8 nuSMEFT

The preferred high-energy basis is now gauge-invariant dimension-8 nuSMEFT / SMEFT, rather than a dimension-7 LEFT operator extrapolated to TeV partonic energies.

Key basis reference:

- H.-L. Li, Z. Ren, M.-L. Xiao, J.-H. Yu, Y.-H. Zheng, arXiv:2105.09329, complete sterile-neutrino EFT bases through d<=9.

ATZ remains useful below the EW scale for matching and low-energy bounds. For UHE amplitudes, use the gauge-invariant parent theory first and match downward near the EW scale.

A nonlinear-EW robustness extension would use nuHEFT/HEFT, e.g. Sun–Xiao–Yu, arXiv:2206.07722.

## 7. Preferred new active-nu_tau candidate: fermionic twist-2 / spin-2 operator

A cleaner candidate than the sterile-N derivative operator is a purely active dimension-8 semileptonic two-derivative operator. The current priority benchmark is schematically

\[
O_{\tau c}^{(8)}=
(\bar L_3\gamma_\mu i\overleftrightarrow D_\nu L_3)
(\bar c_R\gamma^\mu i\overleftrightarrow D^\nu c_R),
\]

corresponding to the flavour choice commonly denoted by a structure such as

\[
[O_{8,lu\partial2}]_{3322}.
\]

Important properties:

- fully electroweak gauge invariant;
- no light sterile state is required;
- contains \(\nu_{\tau L}c\to\nu_{\tau L}c\), so the BSM daughter remains active;
- lepton number conserving and flavour diagonal;
- amplitude grows schematically as
  \[
  M_8\sim C\hat s^2/\Lambda^4;
  \]
- particularly natural in spin-2/twist-2 UV completions, which can generate the d=8 structure without a dominant d=6 contact term.

The physical attraction is that the interaction can be negligible at ordinary IceCube energies but become important around 100–300 PeV.

### Why charm first

Charm is chosen deliberately to exploit the mismatch between:

- small-x UHE DIS, where high-scale charm is generated radiatively and can be appreciable;
- large-x multi-TeV LHC tails, where the charm luminosity is much smaller than valence u,d.

This may permit a larger UHE effect for a collider-allowed coefficient.

The strange-right-handed analogue should be screened simultaneously as the second fermionic benchmark: larger UHE sea PDF, but likely stronger LHC parton luminosity.

## 8. Phenomenological constraints required for the charm/strange twist-2 candidates

The following must all be imposed before transport:

1. **High-pT tau-pair tails**
   - gauge invariance gives charged-lepton partners such as c cbar -> tau+ tau-;
   - use a dedicated Run-II recast, preferably HighPT where the d=8 semileptonic structures are implemented;
   - exact flavour choice [3322] must be evaluated rather than borrowing valence-quark bounds.

2. **Monojet + MET**
   - same operator gives c cbar -> nu_tau nubar_tau plus ISR;
   - dedicated screen/recast required.

3. **IceCube cross-section / Earth-attenuation data**
   - especially 60 TeV–10 PeV, where the operator should still be small if the desired UHE turn-on is realized.

4. **IceCube EHE and KM3NeT UHE constraints**
   - use event-level/response information only after proper propagation and flavour treatment.

5. **Low-energy nu_tau constraints**
   - expected to be weak because of derivative suppression, but should still be quoted.

6. **Partial-wave unitarity and EFT-domain cuts**
   - apply parton by parton.

7. **Positivity constraints**
   - the sign/correlations of the d=8 coefficients cannot be chosen solely to maximize the UHE signal.

8. **QCD running/matching**
   - required consistently for the canonical operator basis.

No current constraint identified in the literature search automatically excludes a few-TeV charm-twist-2 benchmark, but the dedicated tau-tau + monojet recast is load-bearing and must be completed before claiming an allowed 5% signal.

## 9. Gluon contribution: historical precedent and new opportunity

Gluons can contribute to nu-N scattering in two distinct ways.

### 9.1 Existing direct neutrino-gluon operators

PKR already considered a scalar gluon contact operator (their Eq. 4.8 / repository Case 7), schematically

\[
(\bar\nu^cP_L\nu)G_{\mu\nu}^aG^{a\mu\nu}.
\]

ATZ also includes low-energy neutrino-gluon operators and derives their coherent/DIS contributions.

Resolved-model nu-g scattering also exists in the literature, including leptoquark studies where nu_tau g -> LQ + q can be important, especially for heavy-quark couplings.

Therefore direct gluonic contributions to nu-N cross sections are established precedent.

### 9.2 New candidate: active spin-2/twist-2 gluon operator

The more interesting UHE direction is instead a purely active gluonic twist-2 structure, schematically

\[
O_{\tau G}^{(8)}\sim
(\bar L_3\gamma_{\{\mu}i\overleftrightarrow D_{\nu\}}L_3)
\left[G^{a\mu\rho}G^a{}^\nu{}_{\rho}-\text{trace}\right].
\]

This would generate

\[
\nu_\tau g\to\nu_\tau g
\]

with an active daughter and strong momentum growth.

Its attraction is the very large small-x gluon PDF at UHE. Its disadvantage is equally clear: the LHC gluon luminosity is also large, so collider bounds may be much stronger than for charm.

The project should therefore compare, at **bound-saturated coefficients**:

- charm twist-2;
- strange twist-2;
- gluon twist-2.

The key quantity is not the raw matrix element but

\[
r_\tau(E)=\sigma_{\rm BSM}/\sigma_{\rm SM}
\]

after all collider, low-energy, UHE-neutrino, unitarity, positivity and EFT-domain constraints.

A literature search did not identify a systematic PeV–EeV Earth-propagation study of this active neutrino–gluon twist-2 operator. Treat that as a provisional negative search result, not a proof of absence.

## 10. Current ranking of interaction families

For the original goal — a robust, significant UHE nu_tau-N cross-section modification with realistic Earth propagation — current priority is:

1. **Active d=8 fermionic twist-2 / spin-2 operator, charm initial state**
2. **Active d=8 strange twist-2 analogue**
3. **Active d=8 gluon twist-2 operator**
4. **Resolved spin-2 mediator UV completion of the twist-2 class**
5. **Resolved Z' comparison branch**
6. **Sterile-N derivative nuSMEFT / ATZ completion**
7. **Literal PKR charge-conjugated operator** only as historical/off-diagonal theory, not the main benchmark

This ranking may change after the collider and PDF screens.

## 11. Immediate work plan

Do **not** run another expensive Earth-propagation calculation until this screen is closed.

### A. Exact amplitudes

For charm, strange and gluon twist-2 candidates:

- derive the exact tree-level nu_tau-parton amplitudes;
- include SM interference where allowed;
- fix all spin/color averages;
- derive d2sigma/dxdy;
- map to a canonical d=8 SMEFT basis.

### B. Constraint saturation

For each operator:

- HighPT tau-tau recast;
- monojet + MET screen;
- low-energy neutrino constraints;
- IceCube/KM3NeT constraints;
- positivity;
- partial-wave unitarity;
- EFT validity with xi = 1, 1/2, 1/4.

### C. UHE PDF fold

Using modern PDFs/structure functions, calculate at least

- 1 PeV;
- 10 PeV;
- 100 PeV;
- 220 PeV;
- 300 PeV;
- extend toward EeV if the cosmogenic fold is sensitive there.

Report both total r(E) and the distribution of cross-section support in x, y and hat{s}/Lambda^2.

### D. Publication gate

Only operators with a **constraint-saturated, EFT-controlled** modification of at least about 5% in a region with appreciable flux proceed to the validated Earth solver.

### E. Earth propagation

For surviving active-daughter operators, retain full SM tau regeneration and modify the NC-like loss/regeneration kernel consistently. Include secondary nu_mu from tau decay and propagate species separately before detector folding.

## 12. Historical decision log

Preserve the following sequence explicitly:

- v0.12–v0.15 showed that several momentum-growing kernels could produce nominal O(5–10%) or larger raw UHE effects under old assumptions.
- v0.20.9 established numerical convergence of the coupled transport machinery.
- operator audit then showed that the literal PKR Case-6 diagonal interpretation is invalid and that the manuscript had spliced PKR and ATZ structures.
- ATZ normalization/sign remains unresolved and must be independently rederived.
- an intermediate EW-matching no-go was itself retracted because it cut on E_nu rather than the partonic invariant.
- the correct high-energy framework is a gauge-invariant d=8 basis; ATZ/nuLEFT is primarily the low-energy matching/bounds endpoint.
- Z' work remains useful as resolved-mediator comparison, but it is no longer required as the primary route.
- the current leading idea is an active d=8 spin-2/twist-2 semileptonic operator, with charm, strange and gluon initial states compared under all constraints.

## 13. Standing scientific rule

Do not make a broad no-go or positive-observability claim from a single operator family.

The project claim should be based on the **best allowed momentum-growing operators after a systematic bound-saturated screen**, followed by the already validated Earth-propagation calculation.
