function [Rdf_ra,Rdf_r,Rdf_a,Rdf,Rdr,Rt,se] = ...
      scaleR(Rdf_ra,Rdr,Rdf_ss,Rdr_ss,grid,N)
   % scaleR scale reflectance (R)
   % Matt Cooper, guycooper@ucla.edu, Dec 2020

   % Inputs:
   %   Rdf_ra  = diffuse photon density in spherical coord's       [-]
   %   Rdr     = unscattered photon density along the z-axis       [-]
   %   Rdf_ss  = sum of the squared weights behind Rdf_ra           [-]
   %   Rdr_ss  = sum of the squared weights behind Rdr              [-]
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
   %   se      = standard error of each output (fields ra, r, a, df, dr,
   %             t), from the sums of squares through tallyse

   % differential surface area/steradians of annular rings, projection factor
   dsr = grid.dsr;      % solid angle of each angular bin        [sr]
   dA = grid.dA;        % area of each annular ring              [cm^2]
   cosa = cos(grid.ai); % projection factor

   % sum the 2-d arrays into 1-d and 0-d arrays; the squares sum the
   % same way because one packet lands in one bin
   Rdf_r = sum(Rdf_ra,1); % Eq. 4.3
   Rdf_a = sum(Rdf_ra,2); % Eq. 4.4
   Rdf = sum(Rdf_r);      % Eq. 4.7
   Rdf_rss = sum(Rdf_ss,1);
   Rdf_ass = sum(Rdf_ss,2);
   Rdf_ssum = sum(Rdf_rss);

   % standard errors of the raw sums, then scaled like the outputs below
   se.ra = tallyse(Rdf_ss,Rdf_ra,N)./(dA.*dsr.*cosa);
   se.r = tallyse(Rdf_rss,Rdf_r,N)./dA;
   se.a = tallyse(Rdf_ass,Rdf_a,N)./dsr;
   se.df = tallyse(Rdf_ssum,Rdf,N);
   se.dr = tallyse(Rdr_ss,Rdr,N);
   se.t = tallyse(Rdf_ssum+Rdr_ss,Rdf+Rdr,N);


   % convert raw photon weights into fractions of the incident power, per
   % cm2, per sr, and per cm2 per sr
   Rdf_ra = Rdf_ra./(dA.*dsr.*cosa.*N); % Eq. 4.9
   Rdf_r = Rdf_r./(dA.*N);              % Eq. 4.13
   Rdf_a = Rdf_a./(dsr.*N);             % Eq. 4.15
   Rdr = Rdr/N;                         % + Rsp/N
   Rdf = Rdf/N;                         % Eq. 4.17
   Rt = Rdf+Rdr;
end
