# imcrt

[![DOI](https://zenodo.org/badge/344242726.svg)](https://zenodo.org/badge/latestdoi/344242726)

This repository contains source code for the Monte Carlo radiative transfer model used in the following publication:

Cooper, M.G., Smith, L.C., Rennermalm, A.K., Tedesco, M., Muthyala, R., Leidman, S.Z., Moustafa, S.E., Fayne, J.V., 2021. Spectral attenuation coefficients from measurements of light transmission in bare ice on the Greenland Ice Sheet. The Cryosphere 15, 1931–1953. https://doi.org/10.5194/tc-15-1931-2021

If you use the model, the starting point is mcrt_verify.m. That code has options to perform three simulations that reproduce benchmark solutions to the radiative transfer equation, and therefore verifies model accuracy, but could easily be modified for a different problem by setting the inherent optical properties and geometry to new values. The remainder of the repository includes code needed to reproduce the detector interference simulations reported in the paper above. If you wanted to investigate the influence of an instrument on optical measurements, that code would be a good place to start (e.g. see rodintersect.m)
