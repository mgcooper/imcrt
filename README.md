# `imcrt`

[![DOI](https://zenodo.org/badge/344242726.svg)](https://zenodo.org/badge/latestdoi/344242726)

ice Monte Carlo Radiative Transfer.  

## Dependencies

A working installation of MATLAB. The test suite runs on R2025b with no toolboxes. Developed on R2020b and Octave 7.2.0. Older releases and GNU Octave are nominally supported but untested. The project has no dependencies other than the code included here.

## Install

Run `Setup.m`. If running in Octave, check `.octaverc`.

## Usage

Run `mcrt_verify` to verify model accuracy. There are two simulations that compare model output with van De Hulst's tabulated solutions to the transfer equation. It could easily be modified for a different problem by setting the inherent optical properties and geometry to new values.

The `examples` directory includes code needed to reproduce the detector interference simulations reported in the paper below. If you wanted to investigate the influence of an instrument on optical measurements, that code would be a good place to start (e.g. see `rodintersect.m`).

For general use, there is a library of "inherent optical properties" (scattering and absorption coefficients) for water ice Ih saved in `dat/mie_iops_dE.mat`. The prefix `mie_` refers to the Mie scattering formulas used to compute the scattering coefficients. The suffix `_dE` refers to the "delta Eddington" approximation used to compute the extinction coefficients: $c=\sqrt{3ab_e}$ with extinction coefficient $c$, absorption coefficient $a$, and effective (or 'reduced') scattering coefficient $b_e$. Note that this equation is identical to the diffusion approximation, where $c$ is sometimes called the "propagation coefficient". It's inverse $1/c$ is the "transport length". See `doc/tc-2020-53-supplement.pdf` for more details on how these values enter into the Monte Carlo model. In addition, `dat/ssa_iops_dE.mat` contains the same values of absorption coefficient, but values of scattering and extinction coefficient computed with the "specific surface area" approximation, which is also called the "geometric optics" approximation. This approximation is valid for scatterers about the same size or slightly larger than the interacting wavelengths. The effective particle radii are saved in the libary as well.

## Verification

`mcrt_verify` runs one case at a time as a set of seeded simulations (8 for `reflect`, 4 for `fluence`) and prints a PASS/FAIL table. The `reflect` case compares hemispherical and angular reflectance and transmittance with van de Hulst's tabulated solutions (Vol. 2, Table 35). The `fluence` case runs the internal-fluence problem from Wang et al. (1995) with self-consistency checks, since no numeric reference is in the repository. Call `mcrt_verify('fluence')` to pick a case; the default is `reflect`, and the function returns the verdict table and the runs. Every run carries its own standard errors for the reflectance and transmittance outputs in `RT.se`; with fewer than four runs the verdict uses those instead of the spread over runs.

For the full multi-run check, run the driver from the repo root (MATLAB only: it uses `datetime` and the JVM's SHA-256). Run the dry run at `scale = 1e-2` first. The driver checkpoints every run and writes a dated PASS/FAIL report and one figure per case under `verify_out/`:

    Setup; addpath('tests/verify'); vdhverify('verify_out', 1e-2)
    Setup; addpath('tests/verify'); vdhverify('verify_out')

The cases are in `tests/verify/verifycases.m`. To add a case, add an entry there, add a verdict function for its reference, and add a branch to `mcrt_verify`. The fast test suite is `runtests('tests')`; `tests/README.md` describes it and the on-demand perf suite.

## Post-publication corrections

The results in the paper below came from the frozen kernel under `examples/cooper_etal_2021/c_model`, not from `src/mcrt.m`. That kernel saved its resolved tallies per run, but the final dataset behind the paper (`e_postprocess/mcrt_save.m`) keeps only the total transmittance, a hemispherical sum. The corrections below were found in `src/mcrt.m` after publication. `docs/impact-report.md` quantifies each one. It also runs a copy of the frozen 2021 verification script, carrying the production kernel's direct clause, against the corrected kernel on the same seeds: the two agree on every hemispherical sum run for run. A second copy with the production kernel's random-number use, on its own seeds, differs from the first by under one standard error. No correction reaches the published numbers. The detector-rod code under `examples` was not part of that audit.

The commit before the first correction carries the tag `cooper2021-as-published`. Each correction carries a tag. Where a correction changed the golden digest under `tests/golden`, the digest was re-baselined in the same commit with the change stated. `fix-N-roulette` changed no golden line:

- `fix-B-aliasing`: the direction update read the already-updated x cosine when forming the y cosine, which distorted the azimuthal spread. Only radially resolved tallies change (half-weight radii by 10 to 40 percent). The 2021 kernel called `chgdir` and never had this defect.
- `fix-A-direct-tally`: a scattered packet leaving within a bin of grazing was counted as direct. The diffuse hemispherical sums change by 0.01 to 0.1 percent. A spurious direct reflectance of order 1e-6 drops to zero. The direct transmittance of a thick scattering slab, where it is tiny, changes by up to 60 percent of itself. The 2021 production kernel counted only unscattered packets as direct.
- `fix-C-fluence`: the direct beam's absorption was added to every radial column of the resolved fluence instead of the on-axis column. Only `phi_rz` changes; the depth profile `phi_z` and every absorption total are unchanged. The 2021 runs had the fluence tallies off.
- `fix-N-roulette`: a roulette survivor still below the threshold was dropped instead of playing again. Nothing changes at albedo 0.1 and above. Below it, an affected survivor loses its remaining weight, which is less than the threshold. The measured bias at albedo 0.05 is about 6e-6 per launched packet.
- `fix-S-bin-measures`: the solid angle and annulus area of each bin came from the shifted reporting coordinates. That overstated the first bin's measure by 18.5 percent. The first angular and radial densities rise by 18.5 percent and the second fall by 3 percent; every hemispherical sum is unchanged. This one predates the paper, which reported no resolved densities.

Known limitations and deferred items are listed with their reasons in `docs/dispositions.md`.

## How do I cite this?

If you find this model useful, please consider citing the software release (see `CITATION.cff`), and/or the following paper:

Cooper, M.G., Smith, L.C., Rennermalm, A.K., Tedesco, M., Muthyala, R., Leidman, S.Z., Moustafa, S.E., Fayne, J.V., 2021. Spectral attenuation coefficients from measurements of light transmission in bare ice on the Greenland Ice Sheet. The Cryosphere 15, 1931–1953. https://doi.org/10.5194/tc-15-1931-2021

    @article{cooper_2021_TC,
      title = {Spectral Attenuation Coefficients from Measurements of Light Transmission in Bare Ice on the {{Greenland Ice Sheet}}},
      author = {Cooper, M. G. and Smith, Laurence C. and Rennermalm, {\AA}sa K. and Tedesco, Marco and Muthyala, Rohi and Leidman, Sasha Z. and Moustafa, Samiah E. and Fayne, Jessica V.},
      year = {2021},
      month = apr,
      journal = {The Cryosphere},
      volume = {15},
      number = {4},
      pages = {1931--1953},
      publisher = {{Copernicus GmbH}},
      issn = {1994-0416},
      doi = {10.5194/tc-15-1931-2021},
      langid = {english}
    }

## More details

The model roughly implements the one described in Wang et al. (1995). It is not meant to be an exact implementation. The model simulates transfer through a uniform slab. See `doc/tc-2020-53-supplement.pdf` and the two references for technical descriptions.

Wang, L., Jacques, S. L. and Zheng, L.: MCML?Monte Carlo modeling of light transport in multi layered tissues, Computer Methods and Programs
in Biomedicine, 47(2), 131?146, https://doi.org/10.1016/0169 2607(95)01640 F, 1995.