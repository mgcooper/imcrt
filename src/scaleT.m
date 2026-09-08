function [Tdf_ra,Tdf_r,Tdf_a,Tdf,Tdr,Tt] = scaleT(Tdf_ra,Tdr,grid,N)
   % scaleT scale transmittance (T)
   % Matt Cooper, guycooper@ucla.edu, Dec 2020

   % Inputs:
   %   Tdf_ra  = diffuse photon density in spherical coord's       [-]
   %   Tdr     = unscattered photon density along the z-axis       [-]
   %   grid    = the mcrt grid: centers ai and the bin measures dA and dsr
   %             computed from the bin edges (see mcrt)
   %   N       = number of photons

   % Outputs
   %   Tdf_ra  = transmitted diffuse radiance, per incident power [1/cm2/sr]
   %   Tdf_r   = transmitted diffuse irradiance, per incident power    [1/cm2]
   %   Tdf_a   = transmitted diffuse radiant intensity, per incident   [1/sr]
   %   Tdf     = transmitted diffuse fraction of the incident power     [-]
   %   Tdr     = transmitted direct fraction of the incident power      [-]
   %   Tt      = transmitted direct+diffuse fraction                    [-]

   % differential surface area/steradians of annular rings, projection factor
   dsr = grid.dsr;      % solid angle of each angular bin        [sr]
   dA = grid.dA;        % area of each annular ring              [cm^2]
   cosa = cos(grid.ai); % projection factor

   % sum the 2-d arrays into 1-d and 0-d arrays
   Tdf_r = sum(Tdf_ra,1); % Eq. 4.5
   Tdf_a = sum(Tdf_ra,2); % Eq. 4.6
   Tdf = sum(Tdf_r);      % Eq. 4.8

   % convert raw photon weights into fractions of the incident power, per
   % cm2, per sr, and per cm2 per sr
   Tdf_ra = Tdf_ra./(dA.*dsr.*cosa.*N); % Eq. 4.10
   Tdf_r = Tdf_r./(dA.*N);              % Eq. 4.14
   Tdf_a = Tdf_a./(dsr.*N);             % Eq. 4.16
   Tdr = Tdr/N;
   Tdf = Tdf/N;                         % Eq. 4.18
   Tt = Tdf+Tdr;
end
