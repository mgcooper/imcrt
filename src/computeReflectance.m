function [Rdf_ra, Rdf_r, Rdf_a, Rdf, Rdr, Rt, se] = computeReflectance(...
      Rdf_ra, Rdr, Rdf_ss, Rdr_ss, N, grid)
   %COMPUTEREFLECTANCE compute reflectance (R)
   %
   % Inputs:
   %   Rdf_ra  = diffuse photon density in spherical coord's         [-]
   %   Rdr     = unscattered photon density along the z-axis         [-]
   %   Rdf_ss  = sum of the squared weights of Rdf_ra                [-]
   %   Rdr_ss  = sum of the squared weights of Rdr                   [-]
   %   N       = number of photons
   %   grid    = the mcrt grid: centers ai and the bin measures dA
   %             and dsr computed from the bin edges (see mcrt)
   %
   % Outputs
   %   Rdf_ra  = reflected diffuse radiance, per unit incident power [1/cm2/sr]
   %   Rdf_r   = reflected diffuse irradiance, per incident power    [1/cm2]
   %   Rdf_a   = reflected diffuse radiant intensity, per incident   [1/sr]
   %   Rdf     = reflected diffuse fraction of the incident power    [-]
   %   Rdr     = reflected direct fraction of the incident power     [-]
   %   Rt      = reflected direct+diffuse fraction                   [-]
   %   se      = standard error of each output (fields ra, r, a, df,
   %             dr, t), from the sums of squares through mcstderr
   %
   % Matt Cooper, guycooper@ucla.edu, Dec 2020
   %
   % See also:

   % differential surface area/steradians of annular rings, projection factor
   dA = grid.dA;              % area of each annular ring              [cm^2]
   dsr = grid.dsr;            % solid angle of each angular bin        [sr]
   cosa = cos(grid.ai);       % projection factor                      [1]

   % sum the 2-d arrays into 1-d and 0-d arrays
   Rdf_r = sum(Rdf_ra, 1);    % Eq. 4.3
   Rdf_a = sum(Rdf_ra, 2);    % Eq. 4.4
   Rdf = sum(Rdf_r);          % Eq. 4.7
   Rdf_rss = sum(Rdf_ss, 1);
   Rdf_ass = sum(Rdf_ss, 2);
   Rdf_ssum = sum(Rdf_rss);

   % standard errors of the raw sums, scaled like the outputs below
   se.r  = mcstderr(Rdf_rss, Rdf_r, N) ./ dA;
   se.a  = mcstderr(Rdf_ass, Rdf_a, N) ./ dsr;
   se.t  = mcstderr(Rdf_ssum + Rdr_ss, Rdf + Rdr, N);
   se.ra = mcstderr(Rdf_ss, Rdf_ra ,N) ./ (dA .* dsr .* cosa);
   se.df = mcstderr(Rdf_ssum, Rdf, N);
   se.dr = mcstderr(Rdr_ss, Rdr, N);

   % convert raw photon weights into fractions of the incident power, per
   % cm2, per sr, and per cm2 per sr
   Rdf_ra   = Rdf_ra ./ (dA .* dsr .* cosa .* N);     % Eq. 4.9
   Rdf_r    = Rdf_r ./ (dA .* N);                     % Eq. 4.13
   Rdf_a    = Rdf_a ./ (dsr .* N);                    % Eq. 4.15
   Rdr      = Rdr / N;                                % + Rsp/N
   Rdf      = Rdf / N;                                % Eq. 4.17
   Rt       = Rdf + Rdr;
end
