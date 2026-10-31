# v0.77 status after coupling-envelope audit

v0.77 is **not authorized for execution** at this stage.

Reason: its M=100 GeV benchmark was selected from v0.76 using only the first a-crossing and the cross-mass proxy `gq^2 BRinv`. The subsequent audit correctly notes that (i) the complete `sigma_NP/sigma_SM >= 5%` a-interval, including any second crossing, must be examined pointwise, and (ii) cross-mass monojet ranking requires `sigma_prod(M) x A_b(M)` in addition to `gq^2 BRinv`.

The v0.77 files are preserved for audit history but should not be run until the corrected full-a envelope (v0.78) and a mass-dependent production x acceptance comparison identify a justified benchmark.

This does not change the v0.75 conclusion: the 50-GeV point is boundary-unresolved, not a robust exclusion. It also does not establish a no-go for Eq. 4.4 or PKR.
