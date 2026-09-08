function [Tdf_ra, Tdf_r, Tdf_a, Tdf, Tdr, Tt, se] = computeTransmittance(...
      Tdf_ra, Tdr, Tdf_ss, Tdr_ss, grid)
   %COMPUTETRANSMITTANCE compute transmittance (T)
   %
   % Inputs:
   %   Tdf_ra  = diffuse photon density in spherical coord's       [-]
   %   Tdr     = unscattered photon density along the z-axis       [-]
   %   Tdf_ss  = sum of the squared weights of Tdf_ra           [-]
   %   Tdr_ss  = sum of the squared weights of Tdr              [-]
   %   grid    = the mcrt grid: centers ai, the bin measures dA and dsr
   %             computed from the bin edges, and the photon count N
   %
   % Outputs
   %   Tdf_ra  = transmitted diffuse radiance, per incident power [1/cm2/sr]
   %   Tdf_r   = transmitted diffuse irradiance, per incident power    [1/cm2]
   %   Tdf_a   = transmitted diffuse radiant intensity, per incident   [1/sr]
   %   Tdf     = transmitted diffuse fraction of the incident power     [-]
   %   Tdr     = transmitted direct fraction of the incident power      [-]
   %   Tt      = transmitted direct+diffuse fraction                    [-]
   %   se      = standard error of each output (fields ra, r, a, df, dr,
   %             t), from the sums of squares through mcstderr
   %
   % Matt Cooper, guycooper@ucla.edu, Dec 2020
   %
   % See also:

   % differential surface area/steradians of annular rings, projection factor
   dA = grid.dA;              % area of each annular ring              [cm^2]
   dsr = grid.dsr;            % solid angle of each angular bin        [sr]
   cosa = cos(grid.ai);       % projection factor                      [1]
   N = grid.N;                % number of photons                      [-]

   % sum the 2-d arrays into 1-d and 0-d arrays
   Tdf_r = sum(Tdf_ra, 1);    % Eq. 4.5
   Tdf_a = sum(Tdf_ra, 2);    % Eq. 4.6
   Tdf = sum(Tdf_r);          % Eq. 4.8
   Tdf_rss = sum(Tdf_ss, 1);
   Tdf_ass = sum(Tdf_ss, 2);
   Tdf_ssum = sum(Tdf_rss);

   % standard errors of the raw sums, scaled like the outputs below
   se.r  = mcstderr(Tdf_rss, Tdf_r, N) ./ dA;
   se.a  = mcstderr(Tdf_ass, Tdf_a, N) ./ dsr;
   se.t  = mcstderr(Tdf_ssum + Tdr_ss, Tdf + Tdr, N);
   se.ra = mcstderr(Tdf_ss, Tdf_ra ,N) ./ (dA .* dsr .* cosa);
   se.df = mcstderr(Tdf_ssum, Tdf, N);
   se.dr = mcstderr(Tdr_ss, Tdr, N);

   % convert raw photon weights into fractions of the incident power, per
   % cm2, per sr, and per cm2 per sr
   Tdf_ra   = Tdf_ra ./ (dA .* dsr .* cosa .* N);     % Eq. 4.10
   Tdf_r    = Tdf_r ./ (dA .* N);                     % Eq. 4.14
   Tdf_a    = Tdf_a ./ (dsr .* N);                    % Eq. 4.16
   Tdr      = Tdr / N;
   Tdf      = Tdf / N;                                % Eq. 4.18
   Tt       = Tdf + Tdr;
end
