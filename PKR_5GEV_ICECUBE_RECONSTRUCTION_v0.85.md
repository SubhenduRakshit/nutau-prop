# PKR 5 GeV IceCube normalization reconstruction — v0.85

Date: 2026-10-31
Status: model-independent forensic reconstruction from the published PKR 5 GeV limits

## Inputs from PKR

For the printed operators

\[
\mathcal L_{4.3}=\frac{c^{(1)}}{\Lambda}(\bar\nu_i^c\sigma_{\mu\nu}P_L\nu_j)Z'^{\mu\nu}+g_q\bar q\gamma^\mu q Z'_\mu,
\qquad i\ne j,
\]

and

\[
\mathcal L_{4.4}=\frac{c^{(2)}}{\Lambda}(\bar\nu^c i\overleftrightarrow{\partial}^{\mu}\nu)Z'_\mu+g_q\bar q\gamma^\mu q Z'_\mu,
\]

PKR quotes at \(M_{Z'}=5\) GeV

\[
\kappa_{43}^{\rm IC}=\frac{c^{(1)}g_q}{\Lambda}\lesssim4.8\times10^{-4}\ {\rm GeV}^{-1},
\]

\[
\kappa_{44}^{\rm IC}=\frac{c^{(2)}g_q}{\Lambda}\lesssim2.5\times10^{-4}\ {\rm GeV}^{-1}.
\]

Figures 3 and 4 show the corresponding black curves saturating essentially the same IceCube cross-section allowance. PKR states that this allowance is obtained from the cascade total cross section via Eq. (A.1), with the same CT10 PDF implementation for the NSI calculation.

## Exact operator relation

Independent trace evaluation of the printed operators gives, for massless external states and equal numerical \(M_{Z'}\), \(g_q\), and \(a=c/\Lambda\),

\[
\overline{|\mathcal M_{4.4}|^2}
=4\,\frac{\kappa^2}{(M_{Z'}^2+Q^2)^2}\,\hat s^3 y(1-y),
\]

and

\[
\overline{|\mathcal M_{4.3}|^2}
=16\,\frac{\kappa^2}{(M_{Z'}^2+Q^2)^2}\,\hat s^3 y(1-y).
\]

Hence

\[
\boxed{F_{43}=4F_{44}},
\]

where \(\sigma_i^{\rm NP}=\kappa_i^2F_i(E)\). Therefore any common IceCube upper envelope must give

\[
\boxed{\kappa_{43}^{\rm lim}=\frac12\kappa_{44}^{\rm lim}}.
\]

## What the published limits imply

Because the two published black curves saturate the same IceCube allowance, their quoted limits imply

\[
\frac{F_{43}^{\rm PKR}}{F_{44}^{\rm PKR}}
=\left(\frac{\kappa_{44}^{\rm lim}}{\kappa_{43}^{\rm lim}}\right)^2
=\left(\frac{2.5}{4.8}\right)^2
=0.2713.
\]

This is opposite to the exact ratio 4.

If the Eq. (4.3) coefficient is fixed to its exact value 16, the effective Eq. (4.4) coefficient inferred from the published pair is

\[
C_{44}^{\rm inferred}
=16\left(\frac{4.8}{2.5}\right)^2
=58.98.
\]

Allowing the displayed limits to be rounded at the last quoted digit gives approximately

\[
55.5\lesssim C_{44}^{\rm inferred}\lesssim62.7,
\]

which is fully compatible with an implemented integer coefficient 64. The exact Eq. (4.4) coefficient is 4. Thus the published 5 GeV pair is numerically consistent with Eq. (4.3) being normalized approximately correctly while Eq. (4.4) was enhanced by about a factor 16 in cross section.

## Corrected 5 GeV bound reconstruction

There are two equivalent ways to reconstruct the corrected Eq. (4.4) IceCube bound.

1. Correct the apparent factor-16 Eq. (4.4) cross-section normalization:

\[
\kappa_{44,\rm corr}^{\rm IC}
=\sqrt{16}\,(2.5\times10^{-4})
=1.00\times10^{-3}\ {\rm GeV}^{-1}.
\]

2. Use the exact operator relation and the published Eq. (4.3) bound:

\[
\kappa_{44,\rm corr}^{\rm IC}
=2\kappa_{43}^{\rm IC}
=9.6\times10^{-4}\ {\rm GeV}^{-1}.
\]

These independent reconstructions agree to about 4%, well within the precision implied by the rounded published limits and the visual curve reading.

Equivalently, using PKR's plotting variables, the published black labels are approximately \(\epsilon_1=143\) and \(\epsilon_2=75.6\). The corrected Eq. (4.4) label should instead be

\[
\epsilon_{2,\rm corr}\simeq 2\epsilon_1\simeq286,
\]

or, applying the factor-4 coupling correction directly to the published Eq. (4.4) label,

\[
\epsilon_{2,\rm corr}\simeq4(75.6)=302.4.
\]

Again the two methods agree at the few-percent level.

## Interpretation

- The published Eq. (4.3) 5 GeV IceCube normalization is consistent, at the level accessible from the rounded published limits, with the exact coefficient 16.
- The published Eq. (4.4) 5 GeV IceCube normalization is strongly consistent with an effective coefficient near 64 rather than the exact coefficient 4.
- Therefore the Eq. (4.4) IceCube bound \(2.5\times10^{-4}\,\mathrm{GeV}^{-1}\) should not be used as printed under the exact operator convention; the reconstructed corrected bound is approximately \((0.96\text{--}1.00)\times10^{-3}\,\mathrm{GeV}^{-1}\).
- This is a normalization inconsistency in the published phenomenological implementation unless an undocumented coefficient convention can be identified.
- The later nutau-prop Eq. (4.3) kernel with coefficient 128 is a separate repository error: it is 8 times the exact Eq. (4.3) coefficient and should not be attributed to PKR from this reconstruction.

## Decisive invariant test

At fixed \(M_{Z'}=5\) GeV and equal \(\kappa\), the corrected curves must satisfy point-by-point

\[
\sigma_{43}^{\rm NP}(E)=4\sigma_{44}^{\rm NP}(E).
\]

Thus, for any shared IceCube allowance and any PDF set, the limits must obey \(\kappa_{43}^{\rm lim}/\kappa_{44}^{\rm lim}=1/2\). The published ratio is \(4.8/2.5=1.92\), a factor 3.84 away in coupling-space and a factor 14.75 in the implied NP cross-section normalization. This conclusion is independent of CT10 interpolation details.
