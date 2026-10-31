# Scope correction from v0.87

The saved v0.43/v0.43b likelihood numbers are historical flux-modifier diagnostics,
not a complete Eq. 4.3/4.4 detector recast. Their event-weight code modifies the
astrophysical flux but does not add the new interaction's hadronic-shower event
term. The inspected public MC lacks true y needed for an exact changed-inelasticity
reweight. A stable sterile state's inverse transition also requires treatment
or a quantified approximation. The detector-systematic splines being fixed is
therefore not the only limitation.

The original outputs remain unchanged. See `PKR_CONSTRAINT_DISPOSITION_v0.87.md`
and `docs/ICECUBE_RESPONSE_REQUIREMENTS_v0.87.md` for the concrete completion
requirements. An SM flux-fit check alone does not validate a BSM template.
