# nutau-prop v0.25 — active d=8 certification pass

Date: 2026-10-29
Status: analytic gates closed; numerical collider/UHE recast scripts committed; Earth propagation deliberately not run

## Executive result

The certification sequence requested after v0.24 was

1. partial-wave unitarity;
2. exact normalization;
3. tau+tau- recast;
4. monojet/positivity;
5. IceCube energy-resolved constraint;
6. exact dynamic-Q UHE fold.

The first two gates and positivity can now be closed analytically. They produce a stronger correction than v0.24:

**the binding elastic partial wave is the J=0 antiquark channel, not the J=1 quark channel.**

For the strict criterion |Re a_J| <= 1/2 and demanding unitarity through shat=Lambda^2,

\[
|C| \le {8\pi\over3}=8.37758.
\]

The positivity-allowed single-coefficient branch is C<=0 in the Alioli/Li convention.

This greatly reduces the sign-independent d=8-squared UHE contribution. Using the v0.23 CT18 prescreen only as a rescaling reference, at Lambda=4 TeV, xi=1 and |C|=8.3776 the pure term is approximately

| E | c_R: pure / SM NC | s_R: pure / SM NC | c_R: pure / total SM | s_R: pure / total SM |
|---|---:|---:|---:|---:|
|100 PeV|0.34%|0.46%|0.11%|0.14%|
|300 PeV|0.70%|0.87%|0.22%|0.28%|

Therefore the viability of this contact branch no longer comes from the C^2 term. The decisive quantity is now the **SM x d=8 interference**, which is formally the leading O(Lambda^-4) correction and whose sign is fixed by positivity.

For C<0:

- charm/u-type has g_R^u<0, hence C*g_R^u>0: constructive interference in the convention used here;
- strange/d-type has g_R^d>0, hence C*g_R^d<0: destructive interference.

This can reverse the earlier raw-PDF ranking: charm may be the better phenomenological candidate even though strange has the larger sea PDF.

---

## 1. Exact operator convention

The two active operators remain

\[
O_{8,lu\partial2}^{3322}=
(\bar L_3\gamma_\mu i\overleftrightarrow D_\nu L_3)
(\bar c_R\gamma^\mu i\overleftrightarrow D^\nu c_R),
\]

\[
O_{8,ld\partial2}^{3322}=
(\bar L_3\gamma_\mu i\overleftrightarrow D_\nu L_3)
(\bar s_R\gamma^\mu i\overleftrightarrow D^\nu s_R),
\]

with

\[
\mathcal L_{NP}={C_f\over\Lambda^4}O_f.
\]

HighPT's exact labels are

- charm: `WC["l2u2D22",{3,3,2,2}]`;
- strange: `WC["l2d2D22",{3,3,2,2}]`.

This is the active spin-2/two-derivative branch, not the ATZ sterile-daughter operator.

---

## 2. Exact partonic normalization — PASS

For massless kinematics,

\[
(p+p')\cdot(k+k')=\hat s-\hat u=\hat s(2-y).
\]

For a fixed LH neutrino and an unpolarized quark PDF, the exact pure-NP contributions are

\[
{d\hat\sigma_q^{NP}\over dy}
={C_f^2\hat s^3\over8\pi\Lambda^8}
(2-y)^2(1-y)^2,
\]

\[
{d\hat\sigma_{\bar q}^{NP}\over dy}
={C_f^2\hat s^3\over8\pi\Lambda^8}
(2-y)^2.
\]

The integrated angular weights are

\[
I_q={31\over30},\qquad I_{\bar q}={7\over3},
\qquad {I_{\bar q}\over I_q}={70\over31}=2.2581.
\]

The fixed-chirality current traces are

\[
|J_q|^2=4\hat u^2,
\qquad
|J_{\bar q}|^2=4\hat s^2.
\]

The factor 1/2 between these fixed-helicity traces and the published DIS formula is the average over the initial unpolarized quark spin. This distinction is essential for partial-wave unitarity.

Analytic regression script: `scripts/check_d8_analytic_v0_25.py`.

---

## 3. Partial-wave unitarity — PASS, with a new correction

### 3.1 Quark channel

Undoing the initial spin average, the fixed-helicity q_R amplitude is

\[
\mathcal M_q=
2C\left({\hat s\over\Lambda^2}\right)^2
{(3+c)(1+c)\over4}.
\]

Using

\[
{(3+c)(1+c)\over4}
={7\over4}d^1_{11}+{1\over4}d^2_{11},
\]

and

\[
\mathcal M=16\pi\sum_J(2J+1)a_Jd^J,
\]

one obtains

\[
a_1^{(q)}={7C\over96\pi}\rho^2,
\qquad
a_2^{(q)}={C\over160\pi}\rho^2,
\qquad \rho\equiv{\hat s\over\Lambda^2}.
\]

At rho=1, |Re a_1|<=1/2 gives

\[
|C|\le {48\pi\over7}=21.54.
\]

### 3.2 Antiquark channel — binding

For the qbar state created by q_R the helicity difference is lambda=0. The fixed-helicity amplitude is

\[
\mathcal M_{\bar q}
=C\left({\hat s\over\Lambda^2}\right)^2(3+c)
=3C\rho^2P_0(c)+C\rho^2P_1(c).
\]

Therefore

\[
a_0^{(\bar q)}={3C\over16\pi}\rho^2,
\qquad
a_1^{(\bar q)}={C\over48\pi}\rho^2.
\]

The J=0 condition is the strongest:

\[
|\mathrm{Re}\,a_0|\le{1\over2}
\quad\Longrightarrow\quad
|C|\rho^2\le {8\pi\over3}.
\]

Thus at shat=Lambda^2,

\[
\boxed{|C|\le8.37758.}
\]

This supersedes the v0.24 C<30.46 statement, which inferred a fixed-helicity amplitude directly from an already spin-averaged DIS cross section and did not test the qbar J=0 channel.

### 3.3 Cutoff convention

If Lambda is the physical UV scale and the EFT is required to remain unitary through shat=Lambda^2, use |C|<=8.3776 independent of the numerical xi cut.

If instead the trusted contact domain ends at M_cut^2=xi Lambda^2, the invariant statement is

\[
|G_8|M_{cut}^4\le {8\pi\over3},
\qquad G_8={C\over\Lambda^4},
\]

or C<=8.3776/xi^2. Both conventions are retained in the numerical script; the default is the stricter `unitary-to-Lambda` interpretation.

---

## 4. Positivity — PASS

For the single-coefficient directions used here, the forward-limit positivity inequalities of Li et al. give

\[
-4C_{8,lu\partial2}\ge0,
\qquad
-4C_{8,ld\partial2}\ge0,
\]

hence

\[
\boxed{C\le0.}
\]

The exact UHE fold therefore no longer chooses the interference sign to maximize a signal; it uses this positivity-fixed sign.

---

## 5. SM interference — exact formula fixed

Define

\[
G_8={C\over\Lambda^4},
\qquad
P_Z(Q^2)={M_Z^2\over M_Z^2+Q^2}.
\]

The right-handed SM neutral-current coefficient is

\[
K_{SM}=2\sqrt2G_F g_R^f P_Z,
\]

while

\[
K_{NP}=G_8\hat s(2-y).
\]

After the unpolarized initial-quark average, the interference terms are

\[
{d\hat\sigma_q^{int}\over dy}
={\sqrt2G_Fg_R^fP_ZG_8\hat s^2\over2\pi}
(2-y)(1-y)^2,
\]

\[
{d\hat\sigma_{\bar q}^{int}\over dy}
={\sqrt2G_Fg_R^fP_ZG_8\hat s^2\over2\pi}
(2-y).
\]

For antineutrinos the q/qbar angular weights interchange.

This is implemented in `scripts/d8_uhe_cert_v0_25.py`.

A conceptual consequence is important: the interference is the controlled leading d=8 correction at O(Lambda^-4). The pure d=8 square is O(Lambda^-8), where unknown higher-dimensional interference terms also enter. Any claim that relies mainly on C^2 rather than interference has a much larger EFT-truncation burden.

---

## 6. Run-II tau+tau- recast — mapping complete; numerical likelihood not executed in this environment

The public HighPT package contains:

- the exact d=8 WC labels `l2u2D22` and `l2d2D22`;
- the full Run-II ATLAS di-tau likelihood (`di-tau-ATLAS`);
- dimension-8 SMEFT support;
- event-yield and chi-square machinery.

The exact recast script has been committed as

`scripts/highpt_d8_tautau_v0_25.wls`.

It scans Lambda=3,4,5 TeV on the positivity branch C<=0 and returns the one-parameter Delta chi2=3.84 interval.

This execution environment has no Wolfram kernel, so the likelihood could not be numerically run here. The repository script is ready for the user's Mac installation.

Important truncation note: `EFTorder -> 4` is the consistent SMEFT O(Lambda^-4) calculation and therefore constrains the SM x d=8 interference. It deliberately does not include the d=8 square, which is O(Lambda^-8).

Under the strict unitarity-to-Lambda bound |C|<=8.38, the external collider proxies previously quoted (which were substantially weaker) would not be binding, but they are not promoted to certified limits until the HighPT run is reproduced.

---

## 7. Monojet / dineutrino — NOT YET CERTIFIED

The same L3 operator produces q qbar -> nu_tau nubar_tau and therefore monojet+MET after ISR.

No published numerical bound was found for the exact flavor-diagonal d=8 `l2u2D22[3322]` or `l2d2D22[3322]` direction. Existing combined charged-lepton / MET semileptonic SMEFT studies are primarily dimension-6 and do not provide a directly transferable number for this spin-2 d=8 coefficient.

Therefore this gate remains `UNRESOLVED`, not silently passed.

This is less urgent in the strict unitarity-to-Lambda branch because |C|<=8.38 is already much smaller than the earlier collider-scale proxies, but a dedicated monojet recast is still required for publication and for the alternative `unitary-to-Mcut` interpretation.

---

## 8. IceCube 60 TeV–10 PeV — published measurement cannot be imposed directly

IceCube PRD 104, 022001 measures the all-flavour DIS cross section in four energy bins, but its fit:

- scales all flavours together;
- scales CC and NC together with a fixed CC/NC ratio;
- recalculates Earth propagation for that common rescaling;
- uses detector PID information.

Our operator is instead tau-specific and NC-like, and it changes only the tau NC redistribution kernel. The published x_i scale factors are therefore not a literal likelihood for this model.

A detector-level reweight or a faithful public-likelihood reconstruction is still required for a certified IceCube bound.

The measurement uncertainties are broad, especially in the 0.5–10 PeV bin, so the strict-unitarity d=8 model is not obviously excluded by the published numbers. This statement is only a screen, not a fit.

---

## 9. Exact dynamic-Q UHE fold — code complete; execution requires the installed LHAPDF environment

The production script

`scripts/d8_uhe_cert_v0_25.py`

implements:

- CT18NLO dynamic Q^2=xys evaluation;
- SM CC+NC and NC references;
- exact pure d=8 term;
- positivity-fixed SM interference;
- nu/nubar averaging;
- xi=1,1/2,1/4 partonic cuts;
- `unitary-to-Lambda` and `unitary-to-Mcut` conventions;
- endpoint support diagnostics.

This assistant runtime does not have the LHAPDF Python bindings, so the exact numerical table could not be generated here. The user's project Mac already has LHAPDF/CT18NLO from the validated repository environment.

The old fixed-Q v0.23 table is retained only as a cross-check. Rescaling its pure term from C=26 to the corrected strict bound C=8.3776 gives the table in the Executive Result. Those values are not the final result because the interference has not been folded dynamically yet.

---

## 10. Gate status

| Gate | Status | Result |
|---|---|---|
| exact operator normalization | PASS | cross sections and q/qbar weights verified |
| partial-wave unitarity | PASS | **qbar J=0 binds: |C|<=8pi/3 at shat=Lambda^2** |
| positivity | PASS | C<=0 |
| tau+tau- HighPT mapping | PASS | exact WC labels and likelihood script prepared |
| tau+tau- numerical recast | PENDING EXECUTION | requires Wolfram/HighPT runtime |
| monojet/dineutrino | UNRESOLVED | no exact published d=8 [3322] recast located |
| IceCube published-fit applicability | FAIL as direct bound | needs tau-specific NC reweight |
| exact dynamic-Q UHE code | PASS | production script committed |
| exact dynamic-Q numerical result | PENDING EXECUTION | requires LHAPDF runtime |
| Earth propagation | NOT RUN | intentionally waits for certification |

---

## 11. Current physics verdict

The active c_R/s_R dimension-8 contact branch is **much less promising than v0.23/v0.24 suggested** if Lambda is a genuine UV scale and perturbative unitarity is required through Lambda.

The corrected strict unitarity bound drives the pure d=8 square below 1% of the NC cross section at 100–300 PeV in the old prescreen. The only remaining way to obtain a significant controlled contact-EFT effect is through the positivity-fixed SM interference.

That makes the exact dynamic-Q interference fold the decisive numerical test.

If the interference also stays below a few percent of the relevant NC kernel, this active contact branch should be demoted and the project should move to either:

1. the active gluon derivative/twist-2 operator, which has much larger small-x luminosity but different collider bounds; or
2. a resolved spin-2 mediator, which is especially well motivated because the contact signal is endpoint dominated.

No Earth propagation should be run before the v0.25 numerical outputs close the remaining gates.

---

## 12. Reproducible files added in v0.25

- `scripts/check_d8_analytic_v0_25.py` — exact y integrals and partial-wave regression.
- `scripts/highpt_d8_tautau_v0_25.wls` — exact HighPT Run-II di-tau scan.
- `scripts/d8_uhe_cert_v0_25.py` — dynamic-Q CT18 UHE fold with interference, positivity and corrected unitarity.

The historical v0.23 and v0.24 notes are preserved unchanged so the correction chain remains auditable.
