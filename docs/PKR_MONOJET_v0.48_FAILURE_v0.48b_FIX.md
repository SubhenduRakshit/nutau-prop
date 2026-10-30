# PKR monojet v0.48 numerical failure and v0.48b repair

## v0.48 failure

The first Mac Studio execution of v0.48 did not produce a physics result. The legacy public CMS simplified-likelihood implementation emitted invalid `slogdet` warnings and subsequently returned a value that made the fixed `g_q` root bracket fail with

`ValueError: f(a) and f(b) must have different signs`.

This is retained as a numerical/statistical failure, not interpreted as exclusion or survival of the PKR benchmark.

## v0.48b repair

v0.48b leaves all physics inputs and detector information unchanged and modifies only the numerical layer:

1. require finite covariance entries;
2. symmetrise the public covariance matrix;
3. diagonalise it and lift non-positive numerical eigenvalues to `1e-12` of the largest eigenvalue;
4. record the relative Frobenius-norm change and abort if it exceeds `1e-5`;
5. require a finite positive simplified-likelihood `mu95`;
6. replace the fixed `g_q in [1e-5,0.25]` root bracket by an adaptive physical-branch scan, recording the available maximum NWA scaling.

The v0.48 physics output is trusted only if v0.48b passes this numerical certification and the 100/250 GeV public-limit validation remains acceptable.
