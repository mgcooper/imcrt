# Project-specific code style — imcrt

This file extends canonical `STYLE.md`.

## Naming

- Tally arrays encode quantity, scattering class, and resolved dimensions:
  `Xdc_dims` where `X` ∈ {`R`,`T`,`A`} (reflectance, transmittance,
  absorption), `c` ∈ {`f`,`r`} (diffuse, direct), and `dims` ⊆ {`r`,`a`,`z`}
  (radial, angular, vertical). Examples: `Rdf_ra`, `Tdr`, `Adf_rz`, `phi_z`
  (fluence). Keep this scheme for any new tally.
- Physical scalars are terse and match the paper notation: `w`
  (single-scattering albedo), `a` (co-albedo), `g` (asymmetry), `c` (mean free
  path), `wt` (packet weight), `ns` (scatter count). Direction cosines are
  `ux, uy, uz`. Grid/config parameters are capitalized: `R`, `A`, `Z`, `N`.
- Function files are lowercase with no separators: `buildgrid`, `chgdir`,
  `rodscat`, `rodintersect`. Loop counter is `n`.

## Formatting

- Indent function bodies by one level (3 spaces).
- Every `function` closes with `end`.
- Trailing comments carry units in square brackets:
  `c = 1/(ka+ks); % extinction path length [cm]`. Keep unit annotations
  when editing a line that has one — and fix the unit if it is wrong.

## Idioms and patterns

- Core MATLAB only — no toolbox functions in `src/`.
- Function calls inside the loop are decided by the perf harness under
  `tests/perf`: a block whose call overhead measures below 10% of
  the loop is a candidate for a function, and its comments move with
  it if it's moved to one. Measured on 2026-09-08 (R2025b): a bare call
  costs about 11 ns against a photon step of 170 to 220 ns, so one call per
  step is 5 to 7% and two calls per step are 10 to 13%. `chgdir`, `hgcos`,
  and `roulette` are calls; `binindex` is called at exits only, and the
  per-step radial and depth clamps stay inline.
- Keep citations attached to the lines that implement them.
- Randomness: bare `rand` in the kernel; seeding happens in callers and
  test fixtures via `rng(seed,'twister')`.
- In-place assignment of function outputs to input names is allowed ONLY when
  no right-hand side reads a variable the same statement block already
  overwrote.

## Other project conventions

- `examples/cooper_etal_2021/` keeps its 2021-era style frozen; style rules
  here apply to `src/` and root scripts only.
- Third-party files (`src/derivative/derivative.m`) keep their upstream style.
