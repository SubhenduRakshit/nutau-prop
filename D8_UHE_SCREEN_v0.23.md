# nutau-prop v0.23 — active d=8 twist-2 constraint survey and UHE cross-section prescreen

Date: 2026-10-29
Status: first quantitative prescreen; not yet a certified collider recast

## Executive result

The best active-UHE candidates identified so far are the gauge-invariant two-derivative SMEFT operators

\[
O_{8,lu\partial2}^{33\,22}=(\bar L_3\gamma_\mu i\overleftrightarrow D_\nu L_3)(\bar c_R\gamma^\mu i\overleftrightarrow D^\nu c_R)
\]

and

\[
O_{8,ld\partial2}^{33\,22}=(\bar L_3\gamma_\mu i\overleftrightarrow D_\nu L_3)(\bar s_R\gamma^\mu i\overleftrightarrow D^\nu s_R),
\]

with Lagrangian coefficient C/\Lambda^4. They keep nu_tau active, are lepton-number conserving, are fully SM-gauge invariant, and their amplitudes scale as shat^2/Lambda^4. The same operator class is established in the dimension-8 Drell-Yan literature (Alioli et al., arXiv:2003.11615; Boughezal et al., arXiv:2106.05337) and is subject to forward-limit positivity constraints (Li et al., arXiv:2204.13121).

A first CT18-based UHE prescreen shows that C=O(20-30) at Lambda=4 TeV can produce an O(5%) modification of the SM NC cross section at 100-300 PeV if the full shat<Lambda^2 domain is used. Requiring shat<Lambda^2/4 moves the required coefficient to O(100), demonstrating strong endpoint sensitivity. Therefore the operator remains promising but cannot yet be called phenomenologically allowed until exact tau-tau/monojet, positivity-sign and partial-wave unitarity recasts are done.

## 1. Literature constraints that genuinely apply

### 1.1 Drell-Yan d=8 literature

Alioli, Boughezal, Mereghetti and Petriello, Phys. Lett. B 809 (2020) 135703, arXiv:2003.11615, explicitly define O8,lu partial2 and O8,ld partial2. These two-derivative operators populate l=2 partial waves at amplitude level and generate s^2/Lambda^4 growth. Massive spin-2 exchange can generate this class without a simultaneous leading dimension-6 interaction.

Boughezal, Mereghetti and Petriello, Phys. Rev. D 104 (2021) 095022, arXiv:2106.05337, show that momentum-dependent d=8 four-fermion operators give the leading O(s^2/Lambda^4) corrections and that current Drell-Yan data probe such operators at several-TeV scales. Their published numerical fits are not the exact tau-charm/tau-strange flavor assignment used here.

HighPT implements general-flavor high-pT semileptonic likelihoods, including energy-enhanced SMEFT operators through dimension 8 and full Run-II tau-tau channels. This is the correct tool for the exact [L3 L3 c_R c_R] and [L3 L3 s_R s_R] collider recast.

No published table was found giving a direct 95% CL bound on the exact O8,lu partial2 [3322] or O8,ld partial2 [3322] coefficients. Therefore a numerical tau-tau recast is still mandatory.

### 1.2 Positivity

Li, Mimasu, Yamashita, Yang, Zhang and Zhou, JHEP 10 (2022) 107, arXiv:2204.13121, show that the same q q l l two-derivative operators are constrained by forward-limit positivity. The crossed elastic q l -> q l amplitude is precisely the amplitude relevant to UHE neutrino DIS. The single-coefficient sign therefore cannot be chosen arbitrarily to maximize interference; the final analysis must use the positivity-allowed sign/basis mapping. For the prescreen below the pure C^2 contribution is quoted separately because it is sign independent.

### 1.3 IceCube

IceCube PRD 104 (2021) 022001 measured the all-flavor nu-N DIS cross section from 60 TeV to 10 PeV in four energy bins (60-100, 100-200, 200-500 TeV, 0.5-10 PeV), including flavor/PID information and tau decay. It assumes a common all-flavor scaling and fixed CC/NC ratio, so its published scale factors cannot be imposed directly on a tau-specific NC-like operator. A dedicated energy-resolved reweight/recast is required.

### 1.4 UHE KM3NeT bound

Bertolez-Martinez and Hooper, accepted PRD 22 Sep 2026, constrain the total neutrino-nucleon cross section using KM3-230213A to sigma_nuN < 26 sigma_SM around E_nu~220 PeV / E_CM~20 TeV. This is far too weak to exclude a 5-50% modification but is a useful high-energy ceiling.

### 1.5 PKR trend translation: proxy only

PKR's old derivative coefficient C7=8.6e-10 GeV^-3 belongs to a different operator and cannot be transferred as a bound. As a dimensional trend check only, matching the size of a d=7 amplitude C7 shat^(3/2) to a d=8 spin-2 amplitude G8 shat^2 at a characteristic LHC partonic energy gives

G8 ~ C7/sqrt(shat_*).

At sqrt(shat_*)=1 TeV: G8~8.6e-13 GeV^-4, equivalent to Lambda_eff~1.04 TeV for unit coefficient.
At 2 TeV: G8~4.3e-13 GeV^-4, Lambda_eff~1.23 TeV.

At the bookkeeping reference Lambda=4 TeV these correspond to dimensionless C~220 and C~110 respectively. These are NOT certified bounds; they only show that the old PKR trend would permit coefficients large enough to generate an observable UHE effect.

## 2. Partonic cross section for active spin-2 operator

For

L_NP=(C_f/Lambda^4)(bar nu_tau gamma_mu i<->D_nu nu_tau)(bar f_R gamma^mu i<->D^nu f_R),  f=c,s,

the derivative contraction gives (p+p') dot (k+k') = shat(2-y). Using an unpolarized initial quark and a fixed left-handed neutrino,

\[
{d\hat\sigma_q^{NP}\over dy}={C_f^2\hat s^3\over 8\pi\Lambda^8}(2-y)^2(1-y)^2,
\]

\[
{d\hat\sigma_{\bar q}^{NP}\over dy}={C_f^2\hat s^3\over 8\pi\Lambda^8}(2-y)^2.
\]

The SM interference is proportional to C_f and must use the positivity-allowed sign in the final fit. The sign-independent pure term is sufficient for a conservative UHE size test.

For a nearly symmetric sea, the y-integrals are

Iq=31/30, Iqbar=7/3, average=101/60.

The nucleon-level pure contribution scales approximately as

\[
\sigma_f^{NP}(E)\simeq {C_f^2(2m_NE)^3\over 8\pi\Lambda^8}{101\over60}
\int_0^{x_{max}}dx\,x^3[f(x,Q^2)+\bar f(x,Q^2)],
\]

with

x_max = xi Lambda^2/(2 m_N E), xi=1,1/2,1/4.

The production code must evaluate Q^2=xys dynamically; the fixed-Q moments used below are a prescreen only.

## 3. CT18 prescreen

Reference: Lambda=4 TeV, CT18 central PDFs. Representative fixed Q~1.1 TeV moments were used only to estimate the size before the exact dynamic-Q run.

At 100 PeV, xi=1 gives x_max=0.0853. The x^3 moments are approximately

- charm: 3.38e-5
- strange: 4.60e-5
- gluon: 3.15e-4

At 300 PeV, x_max=0.0284:

- charm: 3.92e-6
- strange: 4.89e-6
- gluon: 3.69e-5

Thus the raw gluon twist-2 moment is roughly 7-9 times the c/s moment before operator normalization and alpha_s factors.

SM reference cross sections from Xie et al. (CT18 NNLO, PRD109 113001):

- 100 PeV: sigma_NC~1840 pb, sigma_CC~4080 pb, total~5920 pb.
- 300 PeV: log interpolation of their 200/500 PeV table gives sigma_NC~2785 pb, sigma_CC~5990 pb, total~8775 pb.

### 3.1 Sign-independent pure-NP result for C=26, Lambda=4 TeV

Full partonic domain xi=1:

| E | flavor | sigma_NP [pb] | sigma_NP/sigma_NC | sigma_NP/sigma_total |
|---|---|---:|---:|---:|
|100 PeV|c_R|60.1|3.26%|1.01%|
|100 PeV|s_R|81.8|4.44%|1.38%|
|300 PeV|c_R|188|6.75%|2.14%|
|300 PeV|s_R|234|8.41%|2.67%|

The interference can shift these values, but the pure term alone shows that the strange operator clears the 5% NC target by 300 PeV and charm is close.

With xi=1/2, the same C=26 pure-NP contribution falls to roughly 0.9-1.8% of NC.
With xi=1/4 it falls to roughly 0.2-0.34% of NC.

### 3.2 Coefficient required for a 5% NC modification from the pure term

At Lambda=4 TeV:

| xi | C_c needed at 100 PeV | C_s needed at 100 PeV | C_c needed at 300 PeV | C_s needed at 300 PeV |
|---:|---:|---:|---:|---:|
|1|32.2|27.6|22.4|20.0|
|1/2|61.9|54.8|48.4|43.9|
|1/4|129|117|109|99.7|

This is the key physics result of the prescreen: an observable signal exists naturally in the full EFT domain, but a signal robust to the conservative xi=1/4 truncation requires a much larger Wilson coefficient. Therefore partial-wave unitarity and the exact collider recast are decisive.

## 4. PKR-trend numerical proxy

Using the purely dimensional PKR-to-d8 proxy C~110 at Lambda=4 TeV (matching at sqrt(shat*)~2 TeV), the sign-independent pure contribution gives:

xi=1/4:

- 100 PeV: c_R ~66.6 pb =3.62% of NC; s_R ~82.0 pb =4.46% of NC.
- 300 PeV: c_R ~142 pb =5.10% of NC; s_R ~170 pb =6.09% of NC.

Thus a coefficient of the size suggested by the old PKR trend would still reach the few-percent gate even after the stringent shat<Lambda^2/4 cut. Again, this is a trend proxy, not a transferred experimental limit.

For the still looser 1-TeV matching proxy C~220, even xi=1/4 gives approximately 14-18% of NC at 100 PeV and 20-24% at 300 PeV, but such a large coefficient is much more likely to be excluded by the exact collider/unitarity analysis.

## 5. Current constraint verdict

We cannot yet truthfully label any C~20-100 point as satisfying ALL phenomenological constraints. The following certification steps remain:

1. exact Run-II tau+tau- HighPT recast for O8,lu partial2 [3322] and O8,ld partial2 [3322];
2. monojet/dineutrino constraint from the same L3 operator;
3. exact positivity sign/basis mapping for the chosen convention;
4. exact J=2 partial-wave unitarity limit as a function of xi;
5. energy-resolved IceCube 60 TeV-10 PeV recast for a tau-only NC-like contribution;
6. dynamic-Q CT18/NNPDF UHE integration and PDF uncertainty.

The literature does establish that this operator class is constrained at the several-TeV level for favorable light-flavor Drell-Yan channels, but no exact published [tau tau charm] or [tau tau strange] number was found. Heavy-sea initial states are expected to weaken the LHC reach, but the size must be computed, not assumed.

## 6. Gluon active operator

Cai, Schmidt and Valencia, JHEP05(2018)143, arXiv:1802.09822, considered dimension-8 lepton-gluon derivative operators and found LHC limits on LFV operator scales around 0.5-1.6 TeV in their flavor channels. Their operator class includes derivative lepton bilinears coupled directly to the gluon tensor and has |M|^2 proportional to alpha_s^2 s^4, hence sigma proportional to s^3. This is strong precedent that a direct active lepton-gluon momentum-growing channel is legitimate.

However, their published bounds are LFV and do not directly give the diagonal tau-tau coefficient required here. The exact diagonal active nu_tau/gluon normalization and tau-tau/monojet recast remain to be done. The CT18 x^3 gluon moment is 7-9 times the c/s moment, but the explicit alpha_s normalization can remove much of that raw PDF advantage. No gluon cross-section number is quoted here until its exact amplitude normalization is derived.

## 7. Current ranking

1. O8,ld partial2 [tau tau s_R s_R]: slightly larger UHE rate; collider constraint may be stronger.
2. O8,lu partial2 [tau tau c_R c_R]: slightly smaller UHE rate; probably more collider-friendly.
3. active d=8 lepton-gluon derivative/twist-2 operator: potentially very strong because of gluon PDF, but exact diagonal normalization/bounds not yet certified.
4. ATZ Q8-Q11 sterile-daughter branch: comparison only.

The c_R versus s_R winner should be decided only after the exact HighPT recast. The gluon branch should be promoted if its diagonal tau coefficient survives the collider screen.

## 8. Immediate next calculation

The next gate is not Earth propagation. It is a compact certification calculation:

A. exact dynamic-Q CT18/NNPDF cross section for c_R and s_R;
B. HighPT tau-tau likelihood for the exact flavor coefficients;
C. monojet screen;
D. J=2 partial-wave bound;
E. IceCube binned reweight screen.

Only the bound-saturated survivor should be inserted into v0.20.9 transport.
