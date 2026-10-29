# Chat-only analytic/operator search and cross-section checkpoint — v0.29-c1

Date: 2026-10-30. Status: research checkpoint, **not an all-constraints-certified discovery or publication result**.

## Scope and provenance

User instruction: remain in the chat environment until explicitly authorized to return to work mode; identify interactions producing significant UHE neutrino–nucleon deviations despite other constraints; calculate cross sections and record progress in the repository.

This checkpoint used literature retrieval, the GitHub connector, and bounded remote Wolfram Language evaluations. No Mac Studio shell, local simulation, workflow, or production-code modification was used. This document preserves unsuccessful benchmarks. The preceding repository head was 10edafb41758dc3ea14bd29000613b7f28e5267f; the audited physics-code head was a75a77aabf09c3f9267b1d9e53847a1a1bfc0e98 (v0.29). This is a documentation/calculation supplement, not a new production-code version.

Prior context: v0.7/v0.8 momentum growth at the cross-section level remains a conditional result. The IceCube m²-to-cm² effective-area correction increased its benchmark from approximately 0.001767 to 17.674595 events. The KM3NeT check was 0.242754 versus 0.242767. The subsequent common-flux event screen strongly suppressed the proposed benefit. That screen is not the published IceCube likelihood. Neither the original NC topology assumption nor a cross-section enhancement proves a bright-track explanation.

Working size criterion here: an extra contribution of at least 5% of the repository's SM total at 100–300 PeV, followed by independent constraint and detector tests. This threshold is a research selection criterion, not experimental significance.

## 1. Analytic structure before numerical integration

Use s = 2 m_N E, shat = x s, t = -Q² = -x s y, with m_N = 0.939 GeV. Initial neutrinos have one physical helicity; average initial quark or gluon spin and colour, not neutrino spin. This matters for factors of two. Tree-level kernels are folded with central CT18NLO PDFs; using NLO PDFs does not make the hard cross sections NLO.

All DIS numbers below are per nucleon and use the equal neutrino/antineutrino average. Strange/charm/gluon densities are unchanged by proton-neutron isospin averaging. No flavour-population average is applied: a tau-only interaction must subsequently be weighted by the incident tau-neutrino flux.

Numerical question: evaluate the irreducible PDF convolution after obtaining analytic kernels and limiting forms. Pass criterion: less than 1% change under quadrature refinement. No propagation scan is justified by this checkpoint.

## 2. Resolved charged scalar: large rate, unresolved flavour viability

Define broken-phase vertices
\[
 {\cal L}\supset H^+\{\bar c(g_LP_L+g_RP_R)s+
 g_\tau\bar\nu_{\tau L}\tau_R\}+\mathrm{h.c.}
\]
and A = |g_tau|² (|g_L|²+|g_R|²).

For nu_tau s -> tau^- c, and the crossed antiquark process, the massless spin traces give
\[
 \overline{|{\cal M}|^2}=
 \frac{A\,t^2}{2(M_H^2-t)^2},\qquad
 \frac{d\hat\sigma}{dy}=\frac{A\hat s\,y^2}
 {32\pi(M_H^2+\hat s y)^2}.
\]
The two chiral quark contributions do not interfere in this limit. The scalar produces the opposite charged-lepton chirality to the SM CC interaction, so leading massless SM interference vanishes.

With r = shat/M_H²,
\[
 \hat\sigma=\frac{A}{32\pi\hat s}
 \left[1+\frac1{1+r}-\frac{2\log(1+r)}r\right].
\]
Limits: A shat/(96 pi M_H^4) for r << 1; A/(32 pi shat) for r >> 1, with logarithmic subleading corrections. A full scalar propagator therefore removes indefinite partonic contact growth. Hadronic growth can still arise through PDFs.

For neutrinos the luminosity is s(x)+cbar(x); the charge-averaged result uses [s+sbar+c+cbar]/2. External masses are neglected; their effects require a massive check for precision work, especially at low Q. We impose Q >= 1.3 GeV.

For M_H = 90 GeV and A = 1, the calculated extra cross sections are:

| E (PeV) | Extra sigma (pb) |
|---:|---:|
| 1 | 263.113 |
| 10 | 716.485 |
| 100 | 1658.394 |
| 300 | 2393.814 |

These are **normalization benchmarks, not allowed parameter points**. The older charged-scalar motivation is [Huang et al., arXiv:2112.09476](https://arxiv.org/abs/2112.09476), whose unprojected quark-vertex convention must not be silently equated with A = 1.

Constraint check: Ds -> tau nu measures the pseudoscalar combination g_R-g_L, while DIS measures the sum of squares. In a common low-energy normalization,
\[
 |\epsilon_P|=\frac{|g_\tau(g_R-g_L)|}
 {2\sqrt2 G_F |V_{cs}|M_H^2}
\]
before QCD evolution (phases/signs depend on the vertex convention). A unit single-chiral product at 90 GeV gives |epsilon_P| approximately 3.85. It cannot be defended by ignoring leptonic charm decays.

[Fuentes-Martin et al., arXiv:2003.12421](https://arxiv.org/abs/2003.12421), Table 2, quotes the real, single-coefficient SM-connected cs/tau interval epsilon_P in [-0.008,0.04] at 2 GeV. QCD evolution matters. This interval is not a global bound with arbitrary complex coefficients, cancellations, or disconnected amplitude-sign solutions; no global maximum is claimed here. Its high-pT contact bounds also cannot be applied unchanged to a 90-GeV mediator.

Setting g_L=g_R cancels the Ds pseudoscalar matrix element at matching. For example g_tau=1 and g_L=g_R=0.4 gives A=0.32 and 530.686/766.020 pb at 100/300 PeV. **This algebraic cancellation is not a validated loophole.**

Our gauge-completion check: in a Higgs-basis second-doublet model,
\[
-{\cal L}_Y=\bar Q_L\rho_d H_2d_R+\bar Q_L\rho_u\widetilde H_2u_R+
\bar L_L\rho_e H_2e_R+\mathrm{h.c.},
\]
the charged couplings contain V rho_d and -rho_u^\dagger V. A diagonal charm/strange texture cancelling the cs pseudoscalar leaves a us coupling proportional to V_us rho_ss, hence tau -> K nu; it also creates neutral-scalar interactions. CKM suppression does not by itself save a decay normalized to its CKM-suppressed SM amplitude. Additional flavour textures require their own FCNC audit. Tau/kaon decay, charged/neutral scalar collider limits, electroweak precision, vacuum stability and running remain mandatory. Reject the claim that the simple Ds cancellation certifies this benchmark; do not reject the entire multi-scalar family from this observation.

## 3. Scalar and dual-gluon operators: distinguish active LNV from sterile Dirac

A scalar left-Weyl bilinear is flavour symmetric and permits a tau-tau component. Its LNV active daughter is an antineutrino. This differs from an antisymmetric derivative LNV bilinear and from a Dirac chirality-flipping sterile-daughter operator.

A gauge-invariant LNV parent is dimension nine:
\[
 \frac{c}{M^5}(\overline{\chi_i^c}\chi_j)G^A_{\mu\nu}G^{A\mu\nu},
 \qquad \chi=\widetilde H^\dagger L.
\]
The dual field-strength replacement gives the CP-odd counterpart. A Dirac sterile completion has different field content and matching. Do not transfer collider limits between these theories.

To remove identical-fermion convention ambiguity, define the physical amplitude coefficient kappa by
\[
 {\cal M}=\kappa\,\delta^{ab}B_\nu
 [(k\cdot k')(\epsilon\cdot\epsilon'^*)-
 (k\cdot\epsilon'^*)(k'\cdot\epsilon)] ,
 \qquad \sum|B_\nu|^2=-t.
\]
Then the initial-gluon polarization average gives
\[
 \overline{|{\cal M}|^2}=\frac{|\kappa|^2(-t)^3}{4},\quad
 \frac{d\hat\sigma}{dy}=\frac{|\kappa|^2\hat s^2y^3}{64\pi},\quad
 \hat\sigma=\frac{|\kappa|^2\hat s^2}{256\pi}.
\]
The dual contraction gives the same unpolarized pure-BSM kernel. This equality does not equate their CP or neutrino-mass constraints. There is no tree-level SM elastic nu-g amplitude for interference.

Illustrate kappa = 6.4e-10 GeV^-3; this is **not a newly validated bound**. In a convention with the physical scalar-neutrino vertex normalized as g7 GG, kappa=4g7. Apply shat < xi M².

| E (PeV) | M=4 TeV, xi=1: extra pb | M=1 TeV, xi=1: extra pb |
|---:|---:|---:|
| 1 | 0.0136651 | — |
| 10 | 1.23332 | — |
| 100 | 49.4493 | 1.91711 |
| 300 | 147.480 | 3.68538 |

The 4-TeV column corresponds to only 0.71%/1.42% of the old SM total at 100/300 PeV. If g7=c v²/(2 M^5), this benchmark needs c approximately 5.4e3 at M=4 TeV, versus approximately 5.3 at 1 TeV. It is not a weakly-coupled common-scale UV benchmark at 4 TeV. The lower-scale cut removes most of the rate. Additional neutrino-mass and collider bounds have not been certified.

For a resolved scalar with a gluonic vertex generated by coloured particles, the schematic replacement is kappa(Q²)=y_nu c_SGG/(M_S²+Q²); the gluonic form factor must also be resolved at its own particle threshold. A propagator alone does not complete the whole theory. A vanishing CP-odd vacuum condensate does not prove absence of all neutrino-mass contributions.

The operator classification is compared against [ATZ, arXiv:1812.02778](https://arxiv.org/abs/1812.02778) and [PKR, arXiv:1907.07700](https://arxiv.org/abs/1907.07700); their coefficient bounds are not reused without matching the external states and vertex normalization.

## 4. Active gluon twist-2 control

Use the repository kernel
\[
 \frac{d\hat\sigma}{dy}=\frac{G_8^2\hat s^3}{8\pi}
 (1-y)[1+(1-y)^2],\qquad G_8=C/M^4.
\]
The integrated angular factor is 3/4 and mean y is 13/45. The active neutrino daughter requires redistribution, not a sterile sink.

For C=-8 pi, M=4 TeV:

| E (PeV) | xi=1 extra pb | xi=0.5 extra pb |
|---:|---:|---:|
| 1 | 0.00232191 | — |
| 10 | 1.95824 | — |
| 100 | 229.750 | 62.8280 |
| 300 | 733.984 | 164.157 |

Here xi means shat/M², not sqrt(shat)/M. The xi=1 column gives 3.30%/7.05% at 100/300 PeV; xi=0.5 gives 0.90%/1.58%. Coupled-channel unitarity and collider certification remain outstanding.

The independent xi=1 results differ by about 4–5% from stored v0.27 values at UHE. Our quadrature uses the exact cut endpoint, logarithmic Q² integration and linear interpolation of the PDF grid; the production script uses a fixed x trapezoid, finite y quadrature and LHAPDF interpolation. The discrepancy is unresolved between these implementation choices. It must not be silently assigned to only one cause or promoted to a production correction. The refinement check below does not measure interpolation error.

## 5. New priority: a tau-neutrinogluon resonance

The newer [Bai, Xie and Zhou, arXiv:2510.13948](https://arxiv.org/abs/2510.13948), including its updated appendix, supplies a useful resonant gluon-initiated direction and discusses charged partners. It motivates checking a different mechanism from contact growth. Its forecast is not an all-constraints certificate for our tau benchmark.

Our explicitly specified simplified interaction is
\[
 L_8=(N_8,E_8)^T\sim(8,2)_{-1/2},\quad
 {\cal L}\supset-M\bar L_8L_8+
 \frac{a g_s}{2\Lambda}\bar L_{\tau}\sigma^{\mu\nu}
 G^A_{\mu\nu}L_{8R}^{A}+\mathrm{h.c.}
\]
with vector-like L8, conserved lepton number and a tau-aligned coefficient. It is gauge invariant and anomaly free. Take approximately degenerate doublet masses with no open on-shell E8 W decay. The neutral state then dominantly returns to active nu_tau+g; its charged partner decays to tau+g. This is a resolved particle plus a dimension-five transition operator, **not a fully specified renormalizable UV completion**.

Independent detailed-balance check of the resonant normalization: a fixed-colour spin-averaged two-body decay gives sum|M_decay|²=16 pi M Gamma_in. The production gluon-spin average replaces the decay fermion-spin average, and the colour sum cancels the gluon-colour average. With one neutrino helicity, the one-body phase space gives
\[
 \hat\sigma_{\rm NWA}(\nu g\to N_8\to X)=
 \frac{16\pi^2\Gamma_{\rm in}}M\,\delta(\hat s-M^2)\,B_X,
 \quad
 \Gamma_{\rm in}=\frac{\alpha_s a^2 M^3}{4\Lambda^2}.
\]
Consequently
\[
 \sigma_{\rm res}(\nu N\to N_8\to X)=
 \frac{4\pi^2\alpha_s(M)a^2}{\Lambda^2}\,
 [xg(x,M^2)]_{x=M^2/(2m_NE)}\,B_X .
\]
This is a resonant contribution, not a prediction for the complete off-shell inclusive cross section. For the stated minimal decay spectrum B(nu g) is approximately one. The on-shell contribution vanishes below E=M²/(2m_N), while virtual effects need not vanish.

Targeted benchmark: M=3 TeV, Lambda=9 TeV, a=4. Here M/Lambda=1/3, Gamma_in/M=0.03446, and the resonance threshold is 4.79 PeV. CT18's alpha_s(3 TeV)=0.0775274.

| E (PeV) | Resonant extra sigma, 3 TeV (pb) | 5 TeV comparison (pb) |
|---:|---:|---:|
| 1 | 0 | 0 |
| 10 | 1.33717 | 0 |
| 100 | 427.718 | 30.2395 |
| 300 | 1347.508 | 157.514 |

The comparison uses M=5 TeV, Lambda=15 TeV, a=4. Both have M/Lambda=1/3. The weaker a=3 benchmarks give 240.591/757.973 pb for 3 TeV and 17.0097/88.6015 pb for 5 TeV at 100/300 PeV.

Relative to the repository SM reference, the 3-TeV a=4 point gives **6.14% at 100 PeV and 12.95% at 300 PeV**. It passes the size gate. Its extra cross sections are 4.277e-34 and 1.348e-33 cm², respectively. This is the strongest new candidate in this checkpoint.

Constraint status — provisional, not “allowed”:

- A vector-like degenerate doublet avoids gauge anomalies and a large arbitrary doublet mass splitting; tau alignment avoids tree-level e/mu flavour-changing vertices. This is structural, not a complete loop-level electroweak/flavour test.
- The charged partner and neutral coloured state require pair-production, tau+jet, missing-energy+jet, nonresonant ditau and single-production constraints, with the actual branching fractions. Electron/muon limits cannot be relabelled as a tau recast.
- [Almeida et al., arXiv:2212.06178](https://arxiv.org/abs/2212.06178) distinguishes available-data exclusions from prospective high-luminosity reach. Its approximately 3.5-TeV reach is not a universal present mass exclusion; conversely being outside an older electron benchmark does not certify this tau model.
- The coefficient a=4 is deliberately sizable. A weakly coupled renormalizable theory generally generates dipoles with loop suppression; an explicit strong/composite matching or another justified mechanism is needed. Calling the renormalizable vector-like mass sector a complete UV model would be misleading.
- M/Lambda=1/3 offers a separation of scales, but unknown higher operators are not quantitatively bounded. Gamma/M approximately 3.4% supports an NWA screen away from threshold, not a guaranteed 3.4% total error. PDF slopes, interference, nonresonant diagrams and width effects need assessment.
- Applying NWA does not restore high-energy unitarity of the dimension-five theory; it isolates its pole contribution. The complete theory still needs a validity/partial-wave analysis.
- IceCube and KM3NeT require the tau-specific response and changed inelasticity. The neutral branch is NC-like, not a new direct muon track. A charged-partner cascade could change this, but requires a different mass spectrum and branching calculation.
- No demonstrated compatibility with all current collider, laboratory and neutrino constraints is claimed. No transport run is authorized by finding this provisional point.

## 6. Reproducibility, precision and failure log

Source grids:
- https://lhapdfsets.web.cern.ch/current/CT18NLO/CT18NLO_0000.dat
- https://lhapdfsets.web.cern.ch/current/CT18NLO/CT18NLO.info

Remote Wolfram session: xtltUsYl. Grid: 161 x values, 37 Q values, 11 parton species. Q range approximately 1.295–100000 GeV. Linear interpolation of xf in log(x), log(Q), clipped to nonnegative values. This is not the default LHAPDF interpolator. Alpha_s is linearly interpolated in log(Q) from the set's own metadata.

The Wolfram expression hash Hash[pdfText,"SHA256","HexString"] was
488b284d3c6cbf1706a2956fa57195cf37ee7fd98942e35a40a3742378689a3e.
This is recorded as a Wolfram expression hash, **not a raw-file SHA256**.

Integration coordinates are log(x), log(Q²). Because dx dy = x y dlog(x) dlog(Q²), folding a density f=(xf)/x leaves weight y*(xf). Bounds: Qmin²/s <= x <= min(1,xi M²/s), Qmin² <= Q² <= xs; no EFT cut for the resolved charged scalar. Values outside the PDF grid are zero. Conversion: 1 GeV^-2 = 0.389379338e9 pb.

A 64x40 versus 128x64 Gauss-Legendre comparison over all three DIS kernels and four energies changed results by at most **0.0323%**, passing the 1% quadrature criterion. Additional 256x96 charged-scalar checks gave 1658.392255 and 2393.807468 pb at 100/300 PeV. This is numerical convergence only. PDF uncertainties, interpolation systematics, scale variation, heavy-particle mass corrections and higher hard orders remain unquantified.

SM denominator reused from results/v0.27, in pb for 1,10,100,300 PeV:
969.127310147587, 2778.2553968750053, 6961.783036135988, 10407.910381689187.
Percentages should therefore be read as screening comparisons, not precision ratios computed with a common newly validated integration implementation.

Failure log: one Wolfram result-formatting expression had mismatched brackets and returned a syntax error; no output from that call was used. The corrected separate evaluation succeeded. No failures or rejected models were deleted. The v0.29 collider workflow remained uncertified; this checkpoint did not repair or launch it.

### Minimal Wolfram reproduction

The following is documentation of the bounded remote calculation, not a local-job instruction.

~~~wolfram
pdfText=Import["https://lhapdfsets.web.cern.ch/current/CT18NLO/CT18NLO_0000.dat","Text"];
rows=Select[StringSplit[StringTrim[StringSplit[pdfText,"---"][[2]]],"\n"],
  StringLength[StringTrim[#]]>0&];
parse[s_]:=ToExpression[StringReplace[s,
  RegularExpression["([0-9.])E([+-]?[0-9]+)"]->"$1*^$2"]];
{xnodes,qnodes,pids}=(parse /@ StringSplit[#])& /@ Take[rows,3];
dat=Partition[(parse /@ StringSplit[#])& /@ Drop[rows,3],Length[qnodes]];
pidx=AssociationThread[pids,Range[Length[pids]]];
Do[ip[pid]=ListInterpolation[N[dat[[All,All,pidx[pid]]]],
 {Log[xnodes],Log[qnodes]},InterpolationOrder->1],{pid,pids}];
xf[pid_,x_?NumericQ,q_?NumericQ]:=
 If[x<First[xnodes]||x>1||q<First[qnodes]||q>Last[qnodes],
  0.,Max[0.,ip[pid][Log[x],Log[q]]]];

gq[n_]:=gq[n]=Module[{a=N[Table[k/Sqrt[4 k^2-1],{k,n-1}]],v,w,z},
 {v,w}=Eigensystem[Normal[SparseArray[
  {Band[{1,2}]->a,Band[{2,1}]->a},{n,n}]]];
 z=Ordering[v];Transpose[{v[[z]],2 w[[z,1]]^2}]];
pb=0.389379338*10^9; mn=.939; qmin=1.3;
fold[e_,kind_,mass_,cut_,nx_,ny_]:=Module[
 {ss=2 mn e,lo,hi,gx=gq[nx],gy=gq[ny],lx,x,ql,qu,lq,q2,y,lum,kern},
 lo=Log[qmin^2/ss];
 hi=Log[If[cut>0,Min[.999999,cut mass^2/ss],.999999]];
 If[hi<=lo,Return[0.]];
 pb (hi-lo)/2 Total[Table[
  lx=(hi+lo)/2+(hi-lo)/2 xx[[1]];x=Exp[lx];
  ql=Log[qmin^2];qu=Log[x ss];
  xx[[2]] (qu-ql)/2 Total[Table[
   lq=(qu+ql)/2+(qu-ql)/2 yy[[1]];q2=Exp[lq];y=q2/(x ss);
   lum=If[kind=="charged",
    (xf[3,x,Sqrt[q2]]+xf[-3,x,Sqrt[q2]]+
     xf[4,x,Sqrt[q2]]+xf[-4,x,Sqrt[q2]])/2,xf[21,x,Sqrt[q2]]];
   kern=Switch[kind,
    "charged",x ss y^2/(32 Pi (mass^2+q2)^2),
    "scalarG",(x ss)^2 y^3/(64 Pi),
    "twistG",(x ss)^3(1-y)(1+(1-y)^2)/(8 Pi)];
   yy[[2]] y lum kern,{yy,gy}]],{xx,gx}]]];
Table[{ep,fold[ep*10^6,"charged",90.,0,128,64],
 fold[ep*10^6,"scalarG",4000.,1.,128,64]*(6.4*10^-10)^2,
 fold[ep*10^6,"twistG",4000.,1.,128,64]*(8 Pi/4000.^4)^2},
 {ep,{1,10,100,300}}]

info=StringSplit[Import[
 "https://lhapdfsets.web.cern.ch/current/CT18NLO/CT18NLO.info","Text"],"\n"];
av=parse[StringReplace[StringTrim[StringSplit[
 SelectFirst[info,StringStartsQ[#,"AlphaS_Vals:"]&],":"][[2]]],
 {"["->"{","]"->"}"}]];
as=ListInterpolation[av,{Log[qnodes]},InterpolationOrder->1];
lg[ep_,m_,lam_,a_]:=Module[{x=m^2/(2 mn ep*10^6)},
 If[x>=1,0.,pb*4 Pi^2 as[Log[m]] a^2/lam^2*xf[21,x,m]]];
Table[{ep,lg[ep,3000.,9000.,4.],lg[ep,5000.,15000.,4.]},
 {ep,{1,10,100,300}}]
~~~

## Decision at this checkpoint

Continue chat-only. Prioritize the tau-neutrinogluon interaction for analytic matching and a specific constraint audit. It has a calculated >5% resonant UHE contribution and a concrete gauge-invariant field content, but it is not yet the requested all-constraints survivor. Keep the charged scalar as a flavour-sensitive comparison and the scalar/dual-gluon branches as distinct, presently unsuccessful benchmarks. The active twist-2 gluon branch remains cutoff-sensitive and collider-uncertified.

The original scientific goal is still a demonstrably viable UHE deviation, with correct daughter topology and observable propagation/detector consequences. The size gate has been passed conditionally; the viability and event-level gates have not.
