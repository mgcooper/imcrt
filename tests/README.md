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
11,500 steps per packet. The suite runs it only at 500 packets per run,
in `testVerify.m`.

`kernelfixture.m` puts `src`, `src/derivative`, `tests/oracle`, and
`tests/verify` on the path for one test file through
`matlab.unittest.fixtures.PathFixture`.
`runtests` therefore works without `Setup`. The fixture restores the
caller's global random stream when the file finishes. `testSetup.m` covers
`Setup.m` itself.

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
- `oracle/RotationMatrix.m`: the independent minimal rotation that checks
  the direction-cosine update. The production path excludes it.
- `testRotationMatrix.m`: the oracle is the minimal proper rotation that
  maps a onto b for generic, parallel, antiparallel, and near-antiparallel
  pairs.
- `testChgdirOracle.m`: the direction update keeps the unit norm, deflects
  by the sampled cosine, matches the oracle's azimuth spacing, and the
  inline copy in `mcrt.m` agrees with `chgdir`.
- `testDirectBeam.m`: direct transmittance is Beer-Lambert, direct
  reflectance is zero, and on-axis, grazing, and crafted exits land in
  valid bins.
- `testFluenceBalance.m`: ka times the volume integral of the fluence
  returns the absorbed weight, the direct beam's absorption sits in radial
  bin 1, and a purely absorbing slab is Beer-Lambert.
- `testBinMeasures.m`: the solid angles and annulus areas in `RT.grid`
  are the exact bin measures, resolved tallies sum back to the
  hemispherical fractions, and the diffuse angular tables match van de
  Hulst Table 35 at every nonzero mu: 20 runs of N=5e4 (seeds 1 to 20),
  run-spread standard error, |t| <= 3.5 with 19 degrees of freedom.
- `verify/vdhtable35.m`: van de Hulst Table 35 references and case
  parameters, the one source for the tests and `mcrt_verify.m`.
- `testVerify.m`: the verification layer: case table, van de Hulst and
  fluence verdict tables, and `mcrt_verify.m` running standalone for each
  case at its interactive size.
- `verify/verifycases.m`: the cases `mcrt_verify.m` runs, with kernel
  inputs, run count, packets, seeds, and the t cutoff.
- `verify/vdhverdict.m`: the 16-row PASS/FAIL verdict for the reflect
  case from M runs (hemispherical Rd, Tt, Tdr, Rdr and twelve angular
  rows). It is a struct of plain arrays, so it prints on Octave too.
- `verify/fluenceverdict.m`: the 12-row self-consistency verdict for the
  fluence case (energy, ten direct-beam depth bins, surface fluence).
- `verify/printverdict.m`: prints a verdict struct with `fprintf` and
  returns the number of passing rows.
- `verify/vdhovernight.m`: the multi-run driver (MATLAB only). It runs
  both cases at the full sizes in `verifycases.m` and checkpoints every
  run to a `.mat` file with the kernel hash, version, inputs, and run
  time. The fluence report includes the phi_z depth profile. It
  resumes from matching checkpoints, recomputes stale ones, and writes a
  dated PASS/FAIL report. Hemispherical rows get a 5e-4 absolute allowance
  and angular rows a 2% relative allowance for binning systematics. The
  report names any row that passed only by an allowance.
- `verify/variancecheck.m`: run-to-run variance of Rd and Tt against the
  binomial estimate; a ratio above 2 fails.
- `verify/reportname.m`: a dated report path that takes a numeric suffix
  when the name is taken, so no report is ever overwritten.
- `testOvernight.m`: the driver at scale 1e-3. It checks that the driver:
  - creates the output folder and writes 36 checkpoints and a report;
  - resumes from the checkpoints and recomputes stale or unreadable ones;
  - keeps every report and reaches OVERALL FAIL with tmax = 0;
  - errors on an unwritable folder and closes the report after an error.
- `verify/vdhangular.m`: interpolates each run's angular tallies (pchip
  between bin centers) to the table's mu values and returns z-scores from
  the spread over runs.
- `testRoulette.m`: roulette is terminate-or-boost and unbiased, and a
  deep absorbing slab conserves weight to 1e-6.
- `kernellines.m`: reads a block of `src/mcrt.m` between two marker
  comments so tests can evaluate the kernel's inline code.
- `runblock.m`: evaluates such a block with a struct as its workspace and
  returns every variable afterward.
- `testBuildgrid.m`: grid lengths, orientation, optimized centers, and
  widths for integral inputs; current floor-based lengths for non-integral
  inputs (characterization until hardening, defect K).

## Golden digest

`mcrtgolden.m` runs every case and digests every output field. Scalars are
stored exactly. Each array is stored as its size, sum, sum of squares, and
first and last element. A polynomial hash of its bytes modulo a 25-bit prime
completes the record. Every number is printed with `%.17g`. The hash changes
with any bit or any element order.

The baseline matches the kernel at HEAD. A behavior-neutral commit
(style, dead code) must leave the file bit-identical, which is the proof
that it is neutral. Only a commit that carries a physics fix may
re-baseline, and it must quantify the delta in its commit message. The
file's history therefore holds one version per physics fix; its first
version, from commit 3c99b8d, digests the pre-fix kernel.

Re-baseline with one command from the repo root:

    Setup; addpath('tests'); mcrtgolden('write')

`mcrtgolden('write', file)` writes to another file; the write test uses
it with a temporary folder.

The digest assumes one machine and one MATLAB version: floating-point sums
and the Mersenne Twister stream are reproducible there, and nowhere else is
promised. The current baseline was captured on macOS with R2025b Update 3.
