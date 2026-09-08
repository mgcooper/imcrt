# Dispositions ledger

Items the cleanup campaign (DesignSpec of 2026-08-30, epic `imcrt-98z`)
considered and did not implement. Each entry records the decision, the
reason, and what would reopen it. Numbers cite `docs/impact-report.md`
and the perf measurements of 2026-09-08 (R2025b).

## Backlog

- **Angular-resolution input.** `mcrt` fixes 30 angular bins over the
  hemisphere (`da = A/30`); an earlier fork used 100 for its comparisons.
  A bin-count input is out of the campaign's scope. Reopen when a study
  needs angular densities finer than 3 degrees; the verdict functions
  interpolate between bin centers and would need no change.
- **Roulette threshold input.** `wmin = 1e-4` and `wrr = 10` are fixed
  in `mcrt`. Roulette is unbiased, so the threshold changes only the
  variance and the run time. The impact report compared 1e-4 with 1e-5
  for the van de Hulst case (albedo 0.9), where no packet reached either
  threshold; the fluence verification case (albedo 0.999, about 11,500
  steps per packet) was not compared. Reopen with the angular input, as
  one options argument.
- **Near-axis direction conditioning.** `chgdir` divides by sin(theta)
  and takes the on-axis branch only below 1e-12, so a direction at a
  polar angle theta between 1e-12 and 1e-6 loses about eps/theta^2 in its
  unit norm: 2e-4 at theta = 1e-6 and 1e-2 at theta = 5e-8. The chance of
  such an angle per scatter depends on the phase function: about 1e-12
  for a near-isotropic direction, about 5e-7 per scatter at g = 0.999.
  The verification and paper cases use g up to 0.9. Raising the branch
  threshold to 1e-6 is a physics change that moves the golden digest; do
  it in its own tagged commit before running strongly forward-peaked
  cases.
- **Per-run errors for absorption and fluence.** `RT.se` covers the
  reflectance and transmittance outputs, which each packet touches once.
  A packet deposits absorption in one bin many times, so the per-step
  squares would understate that variance; a per-packet accumulator per
  depth bin would cost about a fifth of a short packet's time. Reopen if
  a fluence verdict needs a single-run error; the run spread serves now.
- **Per-step index clamps as functions.** Two `binindex` calls per step
  measured 10 to 13 percent of the loop, above the 10 percent bar, so the
  radial and depth clamps stay written out in `mcrt` and `binindex` serves
  the exit angle. Reopen if the loop is restructured or a MATLAB release
  makes a call cheaper than the 6 to 11 ns measured.

## Declined

- **`log(rand)` for the path length.** MATLAB's `rand` returns values in
  the open interval (0, 1), so the logarithm is always finite. Changing
  the expression to `log(1 - rand)` would map every draw to a different
  path length and move the trajectory-dependent lines of the golden
  digest for no gain. The 2021 production kernel used `log(1 - rand)`
  and the mirrored Henyey-Greenstein coefficients; the impact report's
  copy with both mappings differs from the paired copy by under one
  standard error on the five hemispherical means.
- **Octave test matrix.** MATLAB release compatibility is preferred, not
  required (owner, 2026-09-07). The code keeps to core features where
  that costs nothing; the suite runs on R2025b, and Octave is not tested.
- **Overflow-bin normalization (defect L).** The radial bin beyond R
  pools every packet at r > R and is normalized as one more ring of width
  `dr`, so its per-area density is meaningless when much light leaves
  past R. The depth bin beyond Z stays empty: a packet that crosses Z is
  scored as transmitted before any absorption is tallied. The fluence
  verification case does populate the radial overflow bin (about 0.9
  percent of the absorbed weight), but every verdict row reads a
  hemispherical sum, a radially integrated angular total, or a depth
  total, which carry that weight correctly; only the last ring's density
  is affected. A known limitation, noted at
  the grid comment in `mcrt`, not a fix.
- **Precomputed Henyey-Greenstein coefficients.** `hgcos` recomputes its
  five coefficients from `g` on every call, about 7 ns of a 170 to 220 ns
  step. The simpler signature was kept.
- **Perf test in the fast suite.** Timing on the development laptop
  swings by a factor of two under load, so the perf test runs on demand
  (`runtests('tests/perf')`) with a paired reference loop and a factor-two
  tolerance rather than failing the fast suite at random.
