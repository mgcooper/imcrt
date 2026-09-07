# Tests

Fast tier for `src/mcrt.m`. Run from the repo root in MATLAB:

    runtests('tests')

The suite must pass in under 30 seconds. It needs no toolbox.

## Seeding

Monte Carlo output is stochastic. The kernel `mcrt` takes no seed argument
and must not gain one. Every test that calls the kernel seeds the generator
itself with `rng(seed, 'twister')` immediately before the call. The shared
parameter sets and their seeds live in `mcrtcases.m`:

| case            | albedo | N   | seed |
|-----------------|--------|-----|------|
| vdh_reflectance | 0.9    | 1e4 | 42   |
| absorbing       | 0.1    | 1e4 | 43   |
| mini_fluence    | 0.99   | 2e3 | 44   |

The albedo-0.999 fluence configuration from `mcrt_verify.m` costs about
11,500 steps per packet and never runs in this suite.

`kernelfixture.m` puts `src` and `src/derivative` on the path for one test
file through `matlab.unittest.fixtures.PathFixture`. `runtests` therefore
works without `Setup`. The fixture restores the caller's global random
stream when the file finishes. `testSetup.m` covers `Setup.m` itself.

## Files

- `testSetup.m`: `Setup.m` adds exactly the repo dirs and drops stale
  example dirs without a warning.
- `testSmoke.m`: each case at its own N; output fields, grid-tied shapes,
  and finite non-negative values.
- `testDeterminism.m`: the same seed reproduces the output exactly; a new
  seed changes it.
- `testEnergyConservation.m`: reflectance + transmittance + absorption is
  within 1e-2 of 1.
- `testGolden.m`: the digest of every output for every case matches
  `golden/mcrt_golden.txt` line for line.

## Golden digest

`mcrtgolden.m` runs every case and digests every output field. Scalars are
stored exactly. Each array is stored as its size, sum, sum of squares, and
first and last element. A polynomial hash of its bytes modulo a 25-bit prime
completes the record. Every number is printed with `%.17g`. The hash changes
with any bit or any element order.

The baseline was captured on the pre-fix kernel on purpose. A
behavior-neutral commit (style, dead code) must leave the file
bit-identical, which is the proof that it is neutral. Only a commit that
carries a physics fix may re-baseline, and it must quantify the delta in
its commit message.

Re-baseline with one command from the repo root:

    Setup; addpath('tests'); mcrtgolden('write')

`mcrtgolden('write', file)` writes to another file; the write test uses
it with a temporary folder.

The digest assumes one machine and one MATLAB version: floating-point sums
and the Mersenne Twister stream are reproducible there, and nowhere else is
promised. The current baseline was captured on macOS with R2025b Update 3.
