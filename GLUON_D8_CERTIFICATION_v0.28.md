# v0.28 — Independent certification of the active gluon d=8 branch

## Operator

\[
Q_{l^2G^2D}^{33}=i(\bar L_3\gamma^\mu\overleftrightarrow D^\nu L_3)
G^A_{\mu\rho}G^{A\rho}{}_{\nu},\qquad
\mathcal L_{\rm NP}=\frac{C_G}{\Lambda^4}Q_{l^2G^2D}^{33}.
\]

The canonical Murphy-basis normalization contains no explicit `g_s` in the two-external-gluon contact vertex: both field strengths contribute their linear derivative pieces. Statements about an `alpha_s` cancellation under unitarity saturation correspond only to a rescaled operator convention.

## Independent Jacob–Wick partial-wave check

For incoming `nu_tau` helicity `-1/2`, the two diagonal gluon-helicity amplitudes are

\[
\mathcal M_-=-2C_G\rho^2\cos\frac\theta2,
\qquad
\mathcal M_+=-2C_G\rho^2\cos^3\frac\theta2,
\qquad \rho=\hat s/\Lambda^2.
\]

For gluon helicity `-1`, the incoming helicity difference has magnitude `1/2`, and

\[
\cos\frac\theta2=d^{1/2}_{1/2,1/2}(\theta).
\]

Thus this channel is pure `J=1/2`. With

\[
\mathcal M=16\pi\sum_J(2J+1)a_J d^J,
\]

one obtains

\[
\boxed{a_{1/2}=-\frac{C_G}{16\pi}\rho^2}.
\]

For gluon helicity `+1`, the helicity difference has magnitude `3/2`, and

\[
\cos^3\frac\theta2=d^{3/2}_{-3/2,-3/2}(\theta),
\]

so this channel is pure `J=3/2`, with

\[
\boxed{a_{3/2}=-\frac{C_G}{32\pi}\rho^2}.
\]

The `J=1/2` channel binds. Requiring `|Re a_J|<=1/2` gives

\[
\boxed{|C_G|\rho^2\le 8\pi}.
\]

Therefore, for unitarity required through `\hat s=\Lambda^2`,

\[
\boxed{|C_G|\le 8\pi=25.132741\ldots}.
\]

## Color check

The two-gluon vertex is proportional to `delta^{AB}`. In a fixed-helicity partial-wave block the color-space matrix is therefore

\[
(a_J)_{AB}=a_J\,\delta_{AB}.
\]

Its eight eigenvalues are all equal to `a_J`; there is no enhancement by `8` or `sqrt(8)` in the unitarity eigenvalue.

The same statement appears in the unpolarized cross section:

\[
\frac1{8}\sum_{A,B}|\delta^{AB}\mathcal M|^2=|\mathcal M|^2.
\]

Thus the color average cancels the final-color sum exactly.

## Polarization normalization check

Averaging the two fixed-gluon-helicity amplitudes gives

\[
\frac12\left(|\mathcal M_-|^2+|\mathcal M_+|^2\right)
=2C_G^2\rho^4(1-y)\left[1+(1-y)^2\right],
\]

which reproduces the matrix element used in the CT18NLO fold. Hence the `8pi` ceiling is not an artifact of a missing spin or color average.

## Flux/attenuation gate — corrected criterion

The project target is a several-percent propagated flux distortion. Therefore the first conservative proxy should be the BSM shift relative to the **total** `CC+NC` interaction cross section, not the NC denominator alone.

The exact v0.27 fold already stores both quantities. For `C_G=-8pi`, `Lambda=4 TeV`:

| xi | E_nu | Delta sigma_NC / sigma_NC | Delta sigma_tot / sigma_tot |
|---:|---:|---:|---:|
| 1.0 | 100 PeV | 12.34% | **3.46%** |
| 1.0 | 300 PeV | 25.96% | **7.32%** |
| 0.5 | 100 PeV | 3.11% | **0.871%** |
| 0.5 | 300 PeV | 5.25% | **1.48%** |
| 0.25 | 300 PeV | 1.11% | **0.314%** |

Thus the contact gluon branch clears the 5% total-rate proxy **only for `xi=1` at 300 PeV**. The previously quoted NC-based gate was too permissive.

The SM NC fraction in the exact fold is about 28% at 100–300 PeV, not a fixed 29.6%; the stored `r_total` values above should be used directly rather than converting with a constant NC fraction.

## Endpoint sensitivity

Tightening the partonic cut from `xi=1` to `xi=0.5` reduces the NC signal by factors

- about `3.97` at 100 PeV;
- about `4.94` at 300 PeV.

At `xi=1`, 93–95% of the BSM contribution at 100–300 PeV comes from `rho>0.25`, and 72–77% from `rho>0.5`.

Therefore the only point that clears the 5% total-rate proxy is also the most UV-endpoint dominated. This materially strengthens the motivation for a resolved mediator if the collider recast removes the `xi=1` contact window.

## Certification status

- exact spin/color normalization: **PASS**;
- independent Jacob–Wick projection: **PASS**;
- strict unitarity ceiling `|C_G|<=8pi`: **PASS**;
- positivity `C_G<=0`: **PASS**;
- exact CT18NLO fold: **PASS numerically**;
- 5% total-rate proxy: **only xi=1, 300 PeV survives**;
- collider recast: **binding next gate**;
- Earth propagation: **do not run yet**.
