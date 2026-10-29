# v0.27 — Active quark derivative-contact branch CLOSED

## Scope

This note closes the momentum-growing dimension-8 active-neutrino contact operators

\[
Q_{l^2u^2D^2}^{3322}=(\bar L_3\gamma_\mu i\overleftrightarrow D_\nu L_3)
(\bar c_R\gamma^\mu i\overleftrightarrow D^\nu c_R),
\]

\[
Q_{l^2d^2D^2}^{3322}=(\bar L_3\gamma_\mu i\overleftrightarrow D_\nu L_3)
(\bar s_R\gamma^\mu i\overleftrightarrow D^\nu s_R),
\]

with

\[
\mathcal L_{\rm NP}=\frac{C_q}{\Lambda^4}Q_q,
\qquad L_3=(\nu_{\tau L},\tau_L)^T.
\]

This conclusion applies to these active `L3-cR/sR` derivative contact directions. It is **not** a no-go theorem for all momentum-growing new physics.

## Exact partonic normalization

For a fixed left-handed neutrino and an unpolarized right-handed quark PDF,

\[
\frac{d\hat\sigma_q}{dy}
=\frac{C_q^2\hat s^3}{8\pi\Lambda^8}(2-y)^2(1-y)^2,
\]

\[
\frac{d\hat\sigma_{\bar q}}{dy}
=\frac{C_q^2\hat s^3}{8\pi\Lambda^8}(2-y)^2.
\]

The angular integrals are

\[
I_q=\frac{31}{30},\qquad I_{\bar q}=\frac73,
\]

so the antiquark channel is intrinsically important.

### Normalization reconciliation: why there is no residual sqrt(2)

The quoted partonic cross section is for an **unpolarized** antiquark PDF. Only one massless antiquark helicity participates in the chiral operator, so the partonic spin average contributes a factor `1/2`:

\[
\frac{d\hat\sigma_{\bar q}}{dy}
=\frac{1}{16\pi\hat s}\,\frac12\,|\mathcal M_h|^2.
\]

Combining this with the certified cross section gives

\[
|\mathcal M_h|^2
=4 C_q^2\left(\frac{\hat s^2}{\Lambda^4}\right)^2(2-y)^2,
\]

hence, choosing the real tree amplitude,

\[
\mathcal M_h
=2C_q\rho^2(2-y)
=C_q\rho^2(3+\cos\theta),
\qquad \rho=\frac{\hat s}{\Lambda^2}.
\]

Equivalently, if one writes

\[
\mathcal M_h=N\frac{3+\cos\theta}{2},
\]

then

\[
\boxed{N=2C_q\rho^2},
\]

not `sqrt(2) C_q rho^2`. The latter follows only if the unpolarized initial-spin average is inadvertently omitted when reconstructing a fixed-helicity amplitude from the cross section.

With

\[
\mathcal M=16\pi\sum_J(2J+1)a_J P_J(\cos\theta),
\]

this immediately gives

\[
a_0^{(\bar q)}=\frac{3C_q}{16\pi}\rho^2,
\qquad
a_1^{(\bar q)}=\frac{C_q}{48\pi}\rho^2.
\]

Thus the strict coefficient bound below is unchanged.

## Partial-wave unitarity

The crossed antiquark channel is the binding one:

\[
a_0^{(\bar q)}=\frac{3C_q}{16\pi}
\left(\frac{\hat s}{\Lambda^2}\right)^2.
\]

With the strict condition \(|\mathrm{Re}\,a_J|\le1/2\),

\[
|C_q|\left(\frac{\hat s}{\Lambda^2}\right)^2\le\frac{8\pi}{3}.
\]

Requiring unitarity through \(\hat s=\Lambda^2\) therefore gives

\[
|C_q|\le\frac{8\pi}{3}=8.37758.
\]

## Positivity

In the Murphy/Alioli convention used by the project, forward-limit positivity requires

\[
C_q\le0.
\]

The strict benchmark is therefore

\[
\boxed{C_q=-\frac{8\pi}{3}},\qquad \Lambda=4~\mathrm{TeV}.
\]

## Exact CT18NLO dynamic-Q UHE result

The executed `results/v0.25/d8_uhe_dynamicQ_v0_25.csv` fold includes the exact pure term, SM-Z interference, dynamic \(Q^2=xys\), and the partonic EFT cut \(\hat s<\xi\Lambda^2\).

For the most permissive \(\xi=1\):

| operator | Eν | Δσ_NC/σ_NC | Δσ_tot/σ_tot |
|---|---:|---:|---:|
| cR | 100 PeV | +0.4269% | +0.1196% |
| cR | 300 PeV | +0.8071% | +0.2277% |
| sR | 100 PeV | +0.3580% | +0.1003% |
| sR | 300 PeV | +0.6501% | +0.1834% |

Tighter EFT cuts reduce the effect further. The signal is also strongly endpoint dominated near the EFT cutoff.

A useful coherence check is the strange-channel sign flip. Positivity fixes `C_q<0`; the SM interference scales as `C_q shat^2`, while the positive BSM square scales as `C_q^2 shat^4`. The exact fold therefore changes from a tiny negative correction at 10 PeV to a positive correction at 100 PeV, exactly as expected.

## Decision

The project target is an observable Earth-propagation distortion at the several-percent level (nominal gate >=5%). Even when the positivity-allowed coefficient saturates strict partial-wave unitarity, these operators remain below 1% in the NC cross section through 300 PeV and below 0.25% in the total cross section.

Therefore:

**CLOSED / phenomenologically dead for the >=5% target.**

An exact ditau or monojet recast can only further reduce the allowed coefficient and cannot rescue this branch. Earth propagation must not be run for this operator class.

## Next momentum-growing direction

Move to the genuine SMEFT class-14 gluonic operator

\[
Q_{l^2G^2D}^{33}=i(\bar L_3\gamma^\mu\overleftrightarrow D^\nu L_3)
G^A_{\mu\rho}G^{A\rho}{}_{\nu},
\]

because it probes the much larger gluon PDF and has a weaker strict partial-wave coefficient ceiling. This is treated separately in `GLUON_D8_SCREEN_v0.27.md`.
