# `imcrt`

[![DOI](https://zenodo.org/badge/344242726.svg)](https://zenodo.org/badge/latestdoi/344242726)

ice Monte Carlo Radiative Transfer.  

## Dependencies

A working installation of Matlab or GNU Octave. Confirmed to run on core Matlab R2020b and core Octave 7.2.0 (no toolboxes or packages required), but it should run on any recent release. The project has no dependencies other than the code included here.

## Install

Run `Setup.m`. If running in Octave, check `.octaverc`.

## Usage

Run `mcrt_verify.m` to verify model accuracy. It runs one case and prints a PASS/FAIL table. The `reflect` case compares hemispherical and angular reflectance and transmittance with van de Hulst's tabulated solutions (Vol. 2, Table 35). The `fluence` case runs the internal-fluence problem of Wang et al. (1995) with self-consistency checks. Set `casename` before running to pick the case. For the full multi-run check, run the driver from the repo root (MATLAB only: it uses `datetime` and the JVM's SHA-256). Run the dry run at `scale = 1e-2` first. The driver checkpoints every run and writes a dated PASS/FAIL report under `verify_out/`:

    Setup; addpath('tests/verify'); vdhovernight('verify_out', 1e-2)
    Setup; addpath('tests/verify'); vdhovernight('verify_out')

The cases live in `tests/verify/verifycases.m`. To add a case, add an entry there, add a verdict function for its reference, and add a branch to the script.

The `examples` directory includes code needed to reproduce the detector interference simulations reported in the paper below. If you wanted to investigate the influence of an instrument on optical measurements, that code would be a good place to start (e.g. see `rodintersect.m`).

For general use, there is a library of "inherent optical properties" (scattering and absorption coefficients) for water ice Ih saved in `dat/mie_iops_dE.mat`. The prefix `mie_` refers to the Mie scattering formulas used to compute the scattering coefficients. The suffix `_dE` refers to the "delta Eddington" approximation used to compute the extinction coefficients: $c=\sqrt{3ab_e}$ with extinction coefficient $c$, absorption coefficient $a$, and effective (or 'reduced') scattering coefficient $b_e$. Note that this equation is identical to the diffusion approximation, where $c$ is sometimes called the "propagation coefficient". It's inverse $1/c$ is the "transport length". See `doc/tc-2020-53-supplement.pdf` for more details on how these values enter into the Monte Carlo model. In addition, `dat/ssa_iops_dE.mat` contains the same values of absorption coefficient, but values of scattering and extinction coefficient computed with the "specific surface area" approximation, which is also called the "geometric optics" approximation. This approximation is valid for scatterers about the same size or slightly larger than the interacting wavelengths. The effective particle radii are saved in the libary as well.

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