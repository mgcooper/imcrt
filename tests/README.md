# Tests

Fast tier for `src/mcrt.m`. Run from the repo root in MATLAB:

    runtests('tests')

The suite must pass in under 30 seconds. It needs no toolbox.

## Seeding

Monte Carlo output is stochastic. The kernel `mcrt` takes no seed argument
and must not gain one. Every test that calls the kernel seeds the generator
itself with `rng(seed, 'twister')` immediately before the call. The shared
parameter sets and their seeds live in `mcrtcases.m`:

| case            | seed |
|-----------------|------|
| vdh_reflectance | 42   |
| absorbing       | 43   |

`kernelfixture.m` puts `src` and `src/derivative` on the path for one test
file through `matlab.unittest.fixtures.PathFixture`. `runtests` therefore
works without `Setup`. The fixture restores the caller's global random
stream when the file finishes. `testSetup.m` covers `Setup.m` itself.

## Files

- `testSetup.m`: `Setup.m` adds exactly the repo dirs and drops stale
  example dirs without a warning.
- `testSmoke.m`: N=1e3 on each case; output fields, grid-tied shapes, and
  finite non-negative values.
- `testDeterminism.m`: the same seed reproduces the output exactly; a new
  seed changes it.
- `testEnergyConservation.m`: reflectance + transmittance + absorption is
  within 1e-2 of 1.
