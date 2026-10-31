# Canonical computation location

Effective 2026-11-01, all local Mac Studio computations for this project must use only

`/Users/subhendurakshit/Research/nutau-prop`

which is the local clone of

`https://github.com/SubhenduRakshit/nutau-prop`.

Do not start new computations from historical working directories such as

- `/Users/subhendurakshit/Research/nutau-prop-cert`
- `/Users/subhendurakshit/Research/nutau-prop-v0.*`

Historical result, log, provenance and audit files that record those paths are retained unchanged because they document where earlier computations actually ran. They must not be interpreted as instructions for new runs.

New run scripts should use paths relative to the repository root and, where practical, guard or `cd` to `${HOME}/Research/nutau-prop`. Mac Studio parallel runs should default to `NUTAU_NPROC=14` unless a numerical library or external program requires a different concurrency setting.
