# PKR-like Eq. (4.3)/(4.4): normalization resolution v0.86

Date: 2026-10-31. Author: Subhendu Rakshit. Baseline: commit
`5eb5c44fb01b95f34a4d202c44c578a31fa00046`.

Publication checkpoint: v0.86. During completion, remote v0.84/v0.85 work
arrived and was integrated without overwriting it; this independent audit was
renumbered from its local draft v0.84 to v0.86. No new physics run was needed
for that version-label change.

**Verdict:** the repo's specified distinct-state chiral transition had mutually
inconsistent DIS and decay normalizations. With its existing physical widths,
the inherited Eq. 4.4 DIS kernel is too large by 16 and the inherited Eq. 4.3
kernel by 8. The fixed-coupling benchmarks previously called 5% are about
0.31% and 0.62%, respectively. This is a demonstrated error in our implementation;
it is **not evidence sufficient to declare the published PKR analysis wrong**.

## 1. What theory is actually calculated

We specify one ordered, distinct-state, massless chiral transition, a universal
vector quark current, and its Hermitian conjugate. There is one incident active
helicity and no average over an absent incident sterile helicity. Quark spins
are averaged; the singlet-current colour sum cancels the colour average.

With metric (+---), sigma^{mu nu}=i[gamma^mu,gamma^nu]/2, PL=(1-gamma5)/2,
the explicit neutrino vertices are

    V44^mu = a (k+k')^mu PL,
    V43^mu = 2 a q_nu sigma^{nu mu} PL,   q=k-k',   kappa=gq*a.

The factor two in V43 is from the field strength. Any use of a different
coefficient convention must change scattering, widths, bounds and the assumed
cutoff together. The conjugate operator creates a second distinct decay
channel; it does not double the amplitude for a fixed incoming active state.

This is the project’s explicit active-to-sterile realization. The literal
charge-conjugated active-flavour operators in PKR must not automatically be
identified with a sterile sink. Flavour sums, identical-particle factors,
Majorana/Dirac assignments and the original UFO are not recovered by this audit.

There is also an exact flavour selection rule. With anticommuting fields,
`(C PL)^T=-C PL` and `(C sigma^{mu nu} PL)^T=+C sigma^{mu nu} PL` imply that
both the charge-conjugated derivative bilinear with an overleftrightarrow
derivative and the tensor bilinear are antisymmetric under interchange of
the two flavour labels. Both vanish for a literally identical diagonal
field. A nonzero Eq. 4.4 realization therefore also needs distinct fields or
an explicit off-diagonal flavour assignment. Suppressed indices in the
printed expression do not establish what the original numerical model used.
These charge-conjugation matrix identities are checked in the exact trace script.

## 2. External states, amplitude, traces and PDF fold

For nu(k) q(p) -> N(k') q(p'), set s=shat, t=(k-k')^2, u=-s-t and K=k+k'.
The derivative amplitude, up to a phase, is

    M44 = kappa/(t-M^2) [ubar(k') PL u(k)] K_mu [ubar(p') gamma^mu u(p)].

The longitudinal propagator term vanishes: q.K=0 and the vector quark current
is conserved. The independent exact Dirac traces give

    S = Tr[kslash' PL kslash PR] = -t,
    H^{mu nu} = (1/2)Tr[pslash' gamma^mu pslash gamma^nu],
    K_mu K_nu H^{mu nu} = -4 s u,
    |M44|^2 = 4 kappa^2 s t u/(M^2-t)^2.

With t=-s y, u=-s(1-y), the physical cross section is nonnegative:

    d sigmahat44/dy = kappa^2 s^2 y(1-y)/[4 pi (M^2+s y)^2].

For s_N=2m_N E, shat=x s_N and F=sum_q[x f_q+x f_qbar],

    d^2 sigma44/dxdy = kappa^2 x s_N^2 y(1-y) F/[4 pi (M^2+Q^2)^2].

The massless external-fermion equations imply the tensor current is twice
the derivative current, up to a phase. Therefore **sigma43=4 sigma44** at
equal numerical a, gq, M and integration domain. This is valid for the virtual
t-channel mediator too; it does not require q^2=M^2. The contact limit gives
sigmahat44=kappa^2 shat^2/(24 pi M^4)+O(shat/M^2) relative corrections.

For a vector at rest, average over its three polarizations. The one-channel
decay amplitudes squared are a^2 M^4/3 and 4a^2 M^4/3. The two conjugate
massless channels give

    Gamma44 = a^2 M^3/(24 pi),    Gamma43 = a^2 M^3/(6 pi).

Consequently the width-normalized DIS prediction is identical for the two:

    (d^2 sigma/dxdy)/Gamma_inv
      = (6 gq^2/M^3) x s_N^2 y(1-y) F/(M^2+Q^2)^2.

This last identity is the useful convention-independent consistency test,
for the specified two-conjugate-channel multiplicity. A relabelling of a cannot
repair a disagreement in sigma/Gamma. With extra invisible decay species the
width-to-one-incident-species mapping must be changed explicitly.

The legacy resolved kernel uses 4/pi, rather than 1/(4pi), for Eq. 4.4.
v0.82/v0.83 then used twice that kernel for Eq. 4.3 despite a fourfold width.
Those are the origins of the factors 16 and 8. The v0.82 width calculation
survives; its cross-section paragraph does not.

## 3. Fresh Mac Studio calculations

Host recorded in metadata: Subhendus-Mac-Studio-2.local, arm64 macOS.
14 worker processes; OMP/OpenBLAS/MKL/Accelerate/NumExpr threads set to one.
CT18NLO member 0 through the installed LHAPDF 6.5.6 binary, LO hard scattering,
isoscalar nucleon, five flavours, Q>=1.3 GeV and x>=1e-9. The y integration
starts at its physical Qmin-dependent boundary, replacing the historical
fixed 1e-5 floor. A ratio is the mean of NP/SM_nu and NP/SM_antinu, with SM=CC+NC;
it is not a ratio to SM NC alone, nor an arbitrary flux-weighted event ratio.

All numbers below hold fixed the listed couplings and cut sqrt(shat)<1/a.

| Benchmark | M (GeV) | gq | a (GeV^-1) | NP/SM at 100 PeV | at 220 PeV | at 300 PeV |
|---|---:|---:|---:|---:|---:|---:|
| v0.73 Eq. 4.4 | 50 | .030 | .00356306608 | 0.3106% | 0.3046% | 0.3022% |
| early Eq. 4.3 threshold | 50 | .030 | .00197258217 | 0.6217% | 0.6204% | 0.6196% |
| v0.81 Eq. 4.4 | 50 | .038 | .002379776 | 0.3179% | 0.3155% | 0.3144% |
| v0.83 Eq. 4.3 | 50 | .017 | .00674370848 | 0.6212% | 0.5972% | 0.5883% |

The corresponding 100-PeV NP cross sections for the first two rows are
2.17770e-35 and 4.35948e-35 cm^2. Values differ slightly from exactly 5%/16
and 5%/8 because we recomputed the integrals with improved quadrature/domain.
Complete numerical outputs, including absolute cross sections, are in
`results/v0.86_normalization/benchmarks.csv`.

### Cutoff dependence is physical input, not just a numerical error

Replacing 1/a by 0.5/a at fixed couplings gives 100-PeV ratios of 0.1232%
and 0.3416% for the first two rows. At 220 PeV they are 0.1183% and 0.3339%.
There is substantial cutoff sensitivity. In general a=c/Lambda does not
determine Lambda without c. The convention 1/a is retained only for comparison
with the historical program, not certified as a universal EFT boundary.

At equal amplitudes a44=2a43, widths and cross sections are equal when the
same physical cutoff is used. Applying 1/a separately gives different cutoffs
and makes an apparent operator hierarchy partly a validity-prescription effect.
The identity does not establish equivalent electroweak/UV completions or all
processes with extra gauge bosons, masses, or different flavour multiplicities.

## 4. Collider implications: what is and is not updated

At fixed mediator mass in the narrow-width approximation, define z=gq^2,
I=Gamma_inv and B=Gamma_qq/gq^2. The existing fully invisible recoil screen is

    P=z I/(Bz+I),    r=K_M P,    D=z-P=gq^2 BR_qq.

We solve the monotone constraints on z analytically at each a. For P<=Pmax,
z<=Pmax I/(I-B Pmax) only when I>B Pmax; otherwise there is no monojet bound
from this inequality. For D<=d, z<=[d+sqrt(d^2+4dI/B)]/2. We also impose
Gamma/M<=0.1, gq^2<=4pi and 1/a>M. All a intervals were optimized, not merely
the first crossing; no new detector events were necessary for this normalization
correction. Direct 72- and 112-point folds check every refined maximum.
The scanned coefficient range is 2e-5<=a<=0.999/M in GeV^-1; the quoted
maxima occur in its interior. These are numerical envelope maxima in that
domain, not rigorous bounds on an unspecified Wilson coefficient or UV model.

CMS publishes universal-vector dijet limits spanning gq=0.03--0.13 for
BR_qq=1. We use 0.03 only as a sufficient screening condition; failure does
not exclude. Using 0.13 at all three masses deliberately relaxes the real
mass-dependent dijet constraint. It is an approximate upper envelope subject
to the stated narrow-width reinterpretation, not an exact mass-specific recast.
See the [official CMS result](https://cms-results.web.cern.ch/cms-results/public-results/publications/EXO-24-007/).

The following maxima at 100 PeV use nominal r<=1 and that deliberately weak
gq sqrt(BR_qq)<=0.13 condition:

| Mass (GeV) | Eq. 4.4 maximum NP/SM | Eq. 4.3 maximum NP/SM |
|---:|---:|---:|
| 50 | 3.142% | 4.514% |
| 100 | 0.2197% | 0.2635% |
| 220 | 0.0583% | 0.0743% |

At M=50 the old approximate mass-specific 0.038 input instead gives 0.3695%
and 0.7953%. These more restrictive values are **sensitivity diagnostics**,
since we have not obtained the exact HEPData table. HEPData access returned
403/unavailable; the official CMS summary and its BR assumption were verified.

**Do not turn the nominal table into a robust exclusion.** With the globally
weak dijet envelope and r<=1.2, Eq. 4.3 at 50 GeV reaches 5.258%, illustrating
the sensitivity of the relaxed bound. The M=100/220 monojet shapes remain
20k-event pilot calibrations. M=50 uses the v0.81 100k same-seed normalization.
Exact operator-specific off-shell production, finite-width distortions,
mass-specific dijet limits and a tau-specific IceCube fit remain outstanding.
No new 95% CL model limit is asserted. Masses outside 50,100,220 GeV are not
closed by these calculations. Omitting dijets allows much larger 50-GeV
effects and is not a viable all-constraints result.

Correcting DIS does not change the existing collider r for the same a,gq,M:
the widths already had the audited normalization. It changes which points
meet the UHE target. Thus v0.73's r=1.239 remains its collider-screen value,
but describing that point as a 5% UHE benchmark is withdrawn.

## 5. Propagation and the original aim

For the specified sterile channel the extra collision term is a sink, with
no interference with the SM active final state. The unregenerated optical
depth is N_A X sigma_NP. Its legacy fixed-point value must be divided by 16
(44) or 8 (43) if the integration domain is unchanged. The previously reported
large vertical disappearance therefore cannot be retained unchanged.

The final active flux is not obtained by multiplying an SM regenerated
spectrum by this simple exponential: NC redistribution and tau regeneration
couple energies. No new coupled transport or detector likelihood is claimed
here. At 220 PeV the corrected fixed-point cross-section effects remain
subpercent; that establishes neither cosmogenic event enhancement nor an
interaction-independent no-go theorem. The goal remains a valid interaction
with a measurable tau-neutrino energy/zenith distortion, after constraints
and detector exposure. The former 5% claims do not establish that goal.

## 6. Reliability, provenance and history

- Exact symbolic scattering traces: three rational angles and a general
  rational parametrization of the scattering angle; Ward identities pass.
- Exact arbitrary-M decay traces and crossing factor pass.
- 20 tests cover contact normalization, endpoint positivity, equal-amplitude
  equivalence and feasibility/saturation of the analytic coupling bounds.
- 288 fresh base folds; higher-order benchmark and maximum checks. Largest
  refined-maximum relative quadrature change: 1.42e-7. Largest interpolation
  error: 1.93e-4 (0.0193%). These are numerical errors only, not PDF/systematic
  uncertainties. No evidence of a need for a long brute-force run.
- The first run imported the older installed package rather than this checkout;
  its traceback is preserved. Explicit PYTHONPATH=src fixed that. A pre-existing
  untracked py.py shadowed pytest's dependency; tests were run from /private/tmp.
  No user file was deleted or modified to fix either environment issue.
- Original v0.73--v0.83 artifacts are retained; the v0.82 cross-section statement
  and v0.83 UHE interpretation are marked superseded. New calculations live in
  a separate module and versioned results directory. Historical ATZ/Case6 code
  is not silently rewritten or reinterpreted.

The [PKR paper](https://arxiv.org/abs/1907.07700) distinguishes its resolved
mediator operators and mass-dependent LHC/IceCube bounds. Reproducing its
original model and conventions is needed to attribute a discrepancy to that
publication. Modern constraints and our implementation mistakes are separate
issues. This audit resolves the latter for the explicit vertices above; it
does not certify a full UV theory or an all-constraints allowed point.

### Audit of the concurrently added v0.85 reconstruction

The published 5-GeV IceCube coupling bounds, 4.8e-4 (43) and 2.5e-4 (44)
GeV^-1, were verified in PKR's text. **If** they refer to the same cross-section
allowance, same external-state/flavour sum and our common coefficient convention,
their ratio contradicts sigma43=4 sigma44. They imply F43/F44=(2.5/4.8)^2
=0.271267, instead of four, for sigma_i=kappa_i^2 F_i. At their respective
printed bounds the audited kernels give sigma43/sigma44=14.7456, not one.
This is a substantive conditional discrepancy in the printed normalization chain.

However, the pair alone cannot identify which implementation or convention is
wrong. Trusting the 43 anchor implies a 44 bound of9.6e-4; trusting the 44
anchor instead implies a 43 bound of1.25e-4 GeV^-1. Neither is a newly fitted
IceCube limit. Reproducing the common allowance or recovering the original UFO
is required to choose an absolute normalization for the published analysis.

The v0.85 assertion that inferred coefficient58.98 with rounding range55.5--62.7
is "fully compatible" with64 is not supported by rounding alone:64 lies outside
that interval. An additional allowance/curve-reading or convention uncertainty
would need to be established. The factor16 suggested for the published 44
implementation is therefore a hypothesis; the repo's own factor16/8 errors
are demonstrated directly from its kernel and widths. v0.85's claim of an
absolute, model-independent reconstruction is too strong. The concurrent
v0.84 threshold runner also retains the old Eq43 normalization and should not
be used to establish a corrected5% threshold.

Reproduce with `bash run_pkr_normalization_v0_86.sh`. The runner does not commit
or push automatically. See the matching JHEP-format technical note in
`manuscript/v0.86/normalization.pdf` and its LaTeX source.
