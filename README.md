# imcrt

[![DOI](https://zenodo.org/badge/344242726.svg)](https://zenodo.org/badge/latestdoi/344242726)

This repository contains source code for the Monte Carlo radiative transfer model used in the following publication, which has been accepted for publication in The Cryosphere:

Cooper, M. G., Smith, L. C., Rennermalm, A. K., Tedesco, M., Muthyala, R., Leidman, S. Z., Moustafa, S. E., and Fayne, J. V.: First spectral measurements of light attenuation in Greenland Ice Sheet bare ice suggest shallower subsurface radiative heating and ICESat-2 penetration depth in the ablation zone, The Cryosphere Discuss. [preprint], https://doi.org/10.5194/tc-2020-53, in review, 2020.

Those interested in using the model are encouraged to use the code mcrt_verify.m, which has options to perform three simulations that reproduce benchmark solutions to the radiative transfer equation, and therefore verifies model accuracy. The remainder of the repository is included to reproduce the results presented in the paper given above, and contains lots of ancillary code that is only useful in the context of repeating the detector interference simulations reported in that paper.
