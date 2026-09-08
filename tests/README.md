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
  by the sampled cosine, and matches the oracle's azimuth spacing.
- `testHgcos.m`: the Henyey-Greenstein sampler's mean is g, every draw is
  a cosine, and g = 0 takes the isotropic branch.
- `testDirectBeam.m`: direct transmittance is Beer-Lambert, direct
  reflectance is zero, on-axis and grazing exits land in valid bins, and
  `binindex` clamps crafted coordinates into range.
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
  fluence verdict tables, the per-run errors below four runs, and
  `mcrt_verify` running for each
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
- `verify/vdhverify.m`: the multi-run driver (MATLAB only). It runs
  both cases at the full sizes in `verifycases.m` and checkpoints every
  run to a `.mat` file with the kernel hash, version, inputs, and run
  time. The fluence report includes the phi_z depth profile. It
  resumes from matching checkpoints, recomputes stale ones, and writes a
  dated PASS/FAIL report. Hemispherical rows get a 5e-4 absolute allowance
  and angular rows a 2% relative allowance for binning systematics. The
  report names any row that passed only by an allowance.
- `verify/impactreport.m`: writes `docs/impact-report.md`: every shared
  case through the pre-fix kernel and each fix tag (extracted from git),
  with attribution per quantity, plus the retrospective against the
  frozen 2021 verification script. MATLAB only.
- `verify/kernelfiles.m`: the one list of every source file that shapes an
  mcrt result, read by the driver's fingerprint, the source-clean check,
  and the impact report's kernel extraction.
- `verify/extractfiles.m`: extracts a list of repository files at a git ref
  into a scratch folder; it refuses an unknown ref or a ref without the
  kernel and skips a helper the ref predates.
- `verify/impactquantities.m`, `verify/impactattribution.m`,
  `verify/relchange.m`, `verify/pairedratio.m`, `verify/unpairedratio.m`,
  `verify/srcmatches.m`, `verify/impactargs.m`, `verify/srctext.m`,
  `verify/agreetext.m`: the quantity summary of one run, the rule that
  names the fixes at which a quantity changed, the change and ratio
  texts, the check that the worktree kernel matches a tag, the report
  defaults, and the two report sentences whose branches depend on the
  source state and the data. The source check is tested in a
  scratch git repository, not the live worktree.
- `testImpact.m`: the attribution rule, the quantity summary, and the
  report generator at a tiny size.
- `verify/verifyplot.m`: the figure of one run against its reference,
  shown by `mcrt_verify` and saved as PNG by `vdhverify`.
- `verify/variancecheck.m`: the run-to-run spread of Rd and Tt against
  the per-run standard errors, a ratio that passes between 0.5 and 2.
- `verify/reportname.m`: a dated report path that takes a numeric suffix
  when the name is taken, so no report is ever overwritten.
- `testVdhverify.m`: the driver at scale 1e-3. It checks that the driver:
  - creates the output folder and writes 36 checkpoints, a report, and one
    PNG figure per case;
  - resumes from the checkpoints and recomputes stale or unreadable ones;
  - keeps every report and reaches OVERALL FAIL with tmax = 0;
  - errors on an unwritable folder and closes the report after an error.
- `verify/vdhangular.m`: interpolates each run's angular tallies (pchip
  between bin centers) to the table's mu values and returns z-scores from
  the spread over runs.
- `testRoulette.m`: `roulette` is terminate-or-boost and unbiased, and a
  deep absorbing slab conserves weight to 1e-6.
- `perf/perfcases.m`, `perf/perfbench.m`, `perf/perfoverhead.m`,
  `perf/baseline.txt`: the timing cases, the timeit harness with its
  baseline reader and writer, the function-call overhead microbenchmark,
  and the baseline measured on one machine and release (see Performance).
- `perf/testPerf.m`: the on-demand perf suite, `runtests('tests/perf')`:
  the kernel's normalized time stays within a factor of two of the
  baseline on the baseline's host and release (filtered elsewhere), a
  written baseline reads back, and the overhead numbers are finite and
  positive. It is not in the fast suite because timing on a loaded laptop
  is not repeatable.
- `testCompute.m`: `computeReflectance`, `computeTransmittance`, and
  `computeAbsorption` normalize synthetic raw tallies by measure and N,
  add direct to diffuse, close the fluence balance on a small grid, and
  return standard errors through `mcstderr`.
- `testMcstderr.m`: `mcstderr` on known contributions, and the kernel's
  per-run `RT.se` against the spread of one hundred seeded runs.
- `testBuildgrid.m`: grid lengths, orientation, optimized centers, nominal
  widths and edge-based measures for integral inputs; a fractional bin
  count is refused by name, whole counts with floating-point error give
  round(n) bins (K), and the plot option draws four axes.
- `testInputValidation.m`: every bad mcrt argument raises `mcrt:input`
  before the loop, a fractional Z/dz raises `buildgrid:nonintegral`, the
  edge values g = -1, g = 1, ks = 0, and a one-bin slab run, and the
  options wmin, wrr, R, dr, and da take effect or are refused by name.

## Performance

`perfbench()` times one seeded run of each perf case with `timeit` three
times, times a fixed reference loop right after each, and keeps the least
case timing with its paired reference; `perfbench('write')` records them
with the architecture, host, and MATLAB release in `perf/baseline.txt`.
`testPerf` compares each case divided by its reference with the same ratio
from that file when the host and release match, with a tolerance of a
factor of two. The reference
cancels most load and thermal drift, and the tolerance is wide because a
laptop still swings by tens of percent, so the test catches a tripled loop,
not a small change. Run the perf suite and re-capture the baseline after a
deliberate speed change or a MATLAB upgrade:

    runtests('tests/perf')

    Setup; addpath('tests', 'tests/perf'); perfbench('write')

`perfoverhead()` times a bare function call and, written out and as a
call, the direction update, the radial index clamp, the Henyey-Greenstein
draw, and roulette, plus one tally event with and without a sum-of-squares
accumulate. It is the evidence behind the kernel's calls to `chgdir`,
`hgcos`, and `roulette`, the exit-only `binindex` call, and the absence of
per-bin variance tallies. The baseline records the architecture, the host
name, and the MATLAB release; `testPerf` compares only on that host.

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
version, from commit 3c99b8d, digests the pre-fix kernel; the last pre-fix
commit, b9d51d5, carries the tag `cooper2021-as-published`. Each physics-fix
commit carries a tag:
`fix-B-aliasing`, `fix-A-direct-tally`, `fix-C-fluence`, `fix-N-roulette`,
and `fix-S-bin-measures`. Any output can be compared across fixes by
checking out a tag and running `mcrtgolden`.

Re-baseline with one command from the repo root:

    Setup; addpath('tests'); mcrtgolden('write')

`mcrtgolden('write', file)` writes to another file; the write test uses
it with a temporary folder.

The digest assumes one machine and one MATLAB version: floating-point sums
and the Mersenne Twister stream are reproducible there, and nowhere else is
promised. The current baseline was captured on macOS with R2025b Update 3.
A MATLAB upgrade is the one exception to the physics-fix rule above: after
an upgrade, run the suite first, and if only `testGolden` fails, re-baseline
with the command above in a commit that changes nothing else and names the
new release. Re-capture the perf baseline the same way.
