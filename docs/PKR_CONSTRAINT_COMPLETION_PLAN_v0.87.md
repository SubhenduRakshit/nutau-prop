# v0.87: finish accessible constraints and identify exact response requirements

2026-10-31. The user correctly objected that v0.86 completed normalization,
not the requested collider/IceCube assessment. No new permission is needed.

Analytic prediction before numerics:
1. Use the explicit v0.86 vertices and matched widths. At equal physical
   amplitude/common cutoff the two scattering channels agree. Retain each
   historical 1/a prescription only with that convention dependence explicit.
2. At fixed a,M, z=gq^2, I=Gamma_inv, B=Gamma_qq/gq^2:
   P=zI/(Bz+I), D=z-P. Both increase monotonically with z. Thus exact tabulated
   dijet limits and the existing monojet likelihood calibration impose analytic
   bounds on z; only the one-dimensional a optimization is numerical.
3. Replace approximate dijet entries by CMS HEPData Figure2a observed limits.
   Quote nominal and relaxed monojet caps, not a combined experimental CL.
   Fold the maximizing points directly with CT18; require <1% interpolation
   and quadrature changes. Use up to14 single-thread worker processes.
4. Audit the actual HESE JSON fields and EHE release inventory before calling
   another fit complete. A sterile outgoing neutrino does not remove the
   hadronic shower yE in the detector. A flux-only reweight is incomplete.
5. For an SM NC template, A_i(E)=integral dy R_i(E,y) f_SM(y|E). The desired
   BSM template uses f_NP. Without y truth or R_i(E,y), A_i alone does not
   identify the BSM rate if f_NP/f_SM is nonconstant. Document this with an
   exact counterexample rather than treating computational time as the obstacle.
6. Read the full ATZ erratum. Distinguish its actual corrections from speculation
   about why it cited MadGraph. Do not infer the original UFO from bound ratios.

No long detector run is justified until the signal response is specified.
An inaccessible exact response will be reported as a precisely evidenced
input limitation, while the public-data work above is completed and versioned.
