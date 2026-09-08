function [Rdf_ra,Rdf_r,Rdf_a,Rdf,Rdr,Rt] = scaleR(Rdf_ra,Rdr,grid,N)
   % scaleR scale reflectance (R)
   % Matt Cooper, guycooper@ucla.edu, Dec 2020

   % Inputs:
   %   Rdf_ra  = diffuse photon density in spherical coord's       [-]
   %   Rdr     = unscattered photon density along the z-axis       [-]
   %   grid    = the mcrt grid: centers ai and the bin measures dA and dsr
   %             computed from the bin edges (see mcrt)
   %   N       = number of photons

   % Outputs
   %   Rdf_ra  = reflected diffuse radiance, per unit incident power [1/cm2/sr]
   %   Rdf_r   = reflected diffuse irradiance, per incident power    [1/cm2]
   %   Rdf_a   = reflected diffuse radiant intensity, per incident   [1/sr]
   %   Rdf     = reflected diffuse fraction of the incident power     [-]
   %   Rdr     = reflected direct fraction of the incident power      [-]
   %   Rt      = reflected direct+diffuse fraction                    [-]

   % differential surface area/steradians of annular rings, projection factor
   dsr = grid.dsr;      % solid angle of each angular bin        [sr]
   dA = grid.dA;        % area of each annular ring              [cm^2]
   cosa = cos(grid.ai); % projection factor

   % sum the 2-d arrays into 1-d and 0-d arrays
   Rdf_r = sum(Rdf_ra,1); % Eq. 4.3
   Rdf_a = sum(Rdf_ra,2); % Eq. 4.4
   Rdf = sum(Rdf_r);      % Eq. 4.7

   % convert raw photon weights into fractions of the incident power, per
   % cm2, per sr, and per cm2 per sr
   Rdf_ra = Rdf_ra./(dA.*dsr.*cosa.*N); % Eq. 4.9
   Rdf_r = Rdf_r./(dA.*N);              % Eq. 4.13
   Rdf_a = Rdf_a./(dsr.*N);             % Eq. 4.15
   Rdr = Rdr/N;                         % + Rsp/N
   Rdf = Rdf/N;                         % Eq. 4.17
   Rt = Rdf+Rdr;
end
