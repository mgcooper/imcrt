# Changelog

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Releases before 2.0.0 have no entry here. See the
[releases page](https://github.com/mgcooper/imcrt/releases) for 1.0, 1.1, and
1.2.0.

## [2.0.0] - 2026-09-10

Five corrections change the numerical results of `src/mcrt.m`. Results from
2.x and 1.x are not interchangeable.

Note: the results in Cooper et al. 2021 came from the frozen code under
`examples/cooper_etal_2021/c_model`, not from `src/mcrt.m`. No correction
below affects the published numbers.
[`tests/reports/Post-publication corrections.md`](https://github.com/mgcooper/imcrt/blob/v2.0.0/tests/reports/Post-publication%20corrections.md)
documents each correction and its impact, and
`tests/reports/impact-report.md` quantifies them.

### Fixed

- The photon direction update used the already-updated x cosine when forming
  the y cosine, which distorted the azimuthal spread. Only radially resolved
  tallies change, by 10 to 40 percent in the half-weight radii. Tagged
  `fix-B-aliasing`.
- A scattered packet exiting within one bin of normal (grazing) was counted
  as direct. Direct now means unscattered. The diffuse hemispherical sums
  change by 0.01 to 0.1 percent, a spurious direct reflectance of order 1e-6
  goes to zero, and the direct transmittance of a thick scattering slab changes
  by up to 60 percent. Tagged `fix-A-direct-tally`.
- The direct beam's absorption was added to every radial column of the
  resolved fluence instead of the on-axis column. Only `phi_rz` changes; the
  depth profile `phi_z` and every absorption total are unchanged. Tagged
  `fix-C-fluence`.
- A Russian-roulette survivor below the `wmin` threshold was dropped instead of
  being played again. Nothing changes at albedo 0.1 and above. The measured
  bias at albedo 0.05 is about 6e-6 per packet. Tagged `fix-N-roulette`.
- The solid angle and annulus area of each bin came from the optimized plotting
  coordinates rather than the bin edges, which overstated the first bin's
  measure by 18.5 percent. The first angular and radial densities increase 18.5
  percent and the second decrease 3 percent. Every hemispherical sum is
  unchanged. Tagged `fix-S-bin-measures`.
- `mcrt` validates its inputs and reconciles the grid lengths with the tally
  sizes, which catches non-integral `Z/dz`.

### Added

- A fast `matlab.unittest` suite under `tests/`, which runs in about three
  seconds with no toolbox.
- A seeded golden-digest regression at `tests/golden/mcrt_golden.txt`.
- A performance test under `tests/perf` with a recorded baseline.
- `vdhverify`, a checkpointed driver that runs both verification cases at full
  size and writes a dated PASS or FAIL report. The `reflect` case compares
  against van de Hulst's tabulated solutions (Vol. 2, Table 35). The `fluence`
  case recreates Fig. 4 of Wang et al. 1995 (the repository holds no numeric
  fluence reference).
- `RT.se`, per-run standard errors for the reflectance and transmittance outputs.
  Absorption and fluence errors come from the spread over multiple seeded runs.
- `RT.N`, the packet count each result came from.
- Name-value options on `mcrt` for the roulette threshold and the grid.
- `tests/oracle/RotationMatrix.m` as the direction-update oracle.

### Changed

- `mcrt_verify` is a function that returns its verdict, and its cases come
  from `tests/verify/verifycases.m`.
- `buildgrid` returns a grid struct, and the tally helpers are named
  `computeReflectance`, `computeTransmittance`, and `computeAbsorption`.
- The MATLAB Project definition under `resources/project` is tracked, so a
  fresh clone opens as a Project.

### Removed

- The stale backup `src/mcrt.bk`, orphaned helpers, a dead clamp, and
  duplicated output-struct assignments.

[2.0.0]: https://github.com/mgcooper/imcrt/releases/tag/v2.0.0
