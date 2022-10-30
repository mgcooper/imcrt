# `imcrt`

[![DOI](https://zenodo.org/badge/344242726.svg)](https://zenodo.org/badge/latestdoi/344242726)

ice Monte Carlo Radiative Transfer.  

Runs on Matlab and Octave.

## Install

Run `Setup.m`. If running in Octave, check `.octaverc`.

## Usage

Run `mcrt_verify.m` to verify model accuracy. There are three simulations that compare model output with van De Hulst's tabulated solutions to the transfer equation. It could easily be modified for a different problem by setting the inherent optical properties and geometry to new values.

The `examples` directory includes code needed to reproduce the detector interference simulations reported in the paper below. If you wanted to investigate the influence of an instrument on optical measurements, that code would be a good place to start (e.g. see rodintersect.m)

## Reference

If you find this model useful, please consider citing the following paper:

Cooper, M.G., Smith, L.C., Rennermalm, A.K., Tedesco, M., Muthyala, R., Leidman, S.Z., Moustafa, S.E., Fayne, J.V., 2021. Spectral attenuation coefficients from measurements of light transmission in bare ice on the Greenland Ice Sheet. The Cryosphere 15, 1931–1953. https://doi.org/10.5194/tc-15-1931-2021
## More

The model roughly implements the one described in Wang et al. (1995). It is not meant to be an exact implementation. The model simulates transfer through a uniform slab. See doc/supplement.pdf and the two references for technical descriptions.

Wang, L., Jacques, S. L. and Zheng, L.: MCML?Monte Carlo modeling of light transport in multi layered tissues, Computer Methods and Programs
in Biomedicine, 47(2), 131?146, https://doi.org/10.1016/0169 2607(95)01640 F, 1995.