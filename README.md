# imcrt

[![DOI](https://zenodo.org/badge/344242726.svg)](https://zenodo.org/badge/latestdoi/344242726)

This repository contains source code for the Monte Carlo radiative transfer model used in the following publication:

Cooper, M.G., Smith, L.C., Rennermalm, A.K., Tedesco, M., Muthyala, R., Leidman, S.Z., Moustafa, S.E., Fayne, J.V., 2021. Spectral attenuation coefficients from measurements of light transmission in bare ice on the Greenland Ice Sheet. The Cryosphere 15, 1931–1953. https://doi.org/10.5194/tc-15-1931-2021

Those interested in using the model are encouraged to use the code mcrt_verify.m, which has options to perform three simulations that reproduce benchmark solutions to the radiative transfer equation, and therefore verifies model accuracy. The remainder of the repository is included to reproduce the results presented in the paper given above, and contains lots of ancillary code that is only useful in the context of repeating the detector interference simulations reported in that paper.
