# v0.86: analytic prediction and numerical acceptance criteria

2026-10-31. User-authorized resolution of Eq. (4.3)/(4.4) on Mac Studio.
Historical results are preserved. This is a consistency audit of the repo's
specified massless, distinct-state chiral transition, not a declaration that
the original PKR analysis is incorrect.

## Definitions before computation

Metric (+---), sigma^{mu nu}=i[gamma^mu,gamma^nu]/2, PL=(1-gamma5)/2.
One incident active helicity; average the quark spin (1/2), not the neutrino
spin. Singlet-current quark colour sum/average cancels. One ordered transition
and its Hermitian conjugate, with two distinct conjugate decay channels.
The explicit derivative vertex is a(k+k')^mu PL and the tensor vertex is
2a q_nu sigma^{nu mu} PL. The latter factor 2 comes from Z'_{mu nu}.
These definitions remove field-sum and tensor-convention ambiguities.

Let k,p -> k',p', t=(k-k')^2, u=-s-t, kappa=gq*a. The predicted derivative
traces are S=-t and K_mu K_nu H^{mu nu}=-4su, hence

    |M44|^2 = 4 kappa^2 s t u / (M^2-t)^2
    d sigmahat44/dy = kappa^2 s^2 y(1-y)/(4 pi (M^2+sy)^2).

For massless external fermions the tensor current is, up to a phase, twice
the derivative current: |M43|^2=4|M44|^2. The predicted total conjugate widths
are Gamma44=a^2 M^3/(24pi), Gamma43=a^2 M^3/(6pi). Crossing must preserve
this factor. The unpolarized vector decay includes a polarization average 1/3.

The nucleon fold uses sum(x f_q+x f_qbar), and therefore the coefficient of
a^2 gq^2 x s_N^2 y(1-y)/(M^2+Q^2)^2 is 1/(4pi) (44), 1/pi (43).
The inherited 4/pi coefficient would overstate this specified 44 channel by
16; the inherited factor-two 43 prescription would overstate 43 by 8.
This is a prediction to test, not an assumed replacement of the ATZ formula.

## Minimum calculations and gates

1. Independent explicit 4x4 Dirac-matrix traces for scattering and decay,
   exact rational kinematics and Ward identities: symbolic equality required.
2. Implement explicit normalized vertices' closed-form kernel separately
   from the historical Case6 module. Check endpoints, positivity, contact
   integral and equal-amplitude/equal-cutoff equivalence of the two vertices.
3. Dynamic-Q CT18NLO folds at affected benchmarks and E=100,220,300 PeV,
   with tighter quadrature: target <1% change; otherwise refine and report.
   Single-thread numerical libraries, up to 14 worker processes.
4. Reuse prior full a envelopes only after direct recomputation checks.
   Analytically solve monotone coupling constraints rather than brute-force
   gq scans. Propagate the exact normalization changes; do not extrapolate
   collider contact bounds across mediator masses.
5. Report collider output as calibrated narrow-width *screens*, not an
   exact model likelihood or all-constraints certificate. No detector rerun
   is justified before this normalization gate.

## Domain and unresolved issues

Massless external fermions, LO hard process with CT18NLO PDFs, isoscalar
target, five quark flavours, Q>=1.3 GeV, existing small-x extrapolation policy.
Mass corrections are O(m^2/shat) away from thresholds; no massive sterile
prediction is made. PDF, Qmin, EW completion and collider systematics are
not controlled by a quadrature check.

The cutoff sqrt(shat)<1/a assumes a particular Wilson coefficient. In
general a=c/Lambda does not fix Lambda. Repeat fixed-coupling cuts at
0.5/a and 1/a, and compare the two vertices at a COMMON physical cutoff.
Since equal amplitudes require a44=2a43, separately imposing 1/a is not an
operator-basis invariant validity prescription.

The literal PKR charge-conjugated flavour operator can involve different
state multiplicities. Its mapping to active-to-sterile transport, and the
original UFO normalization, must not be inferred solely from this audit.
