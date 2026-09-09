function [Adf_rz, Adf_z, Adf, Adr_z, phi_rz, phi_z] = computeAbsorption(...
      Adf_rz, Adr_z, kabs, N, grid)
   %COMPUTEABSORPTION compute absorption (A), also compute fluence
   %
   % Inputs:
   %   Adf_rz  = diffuse photon absorptance in cylindrical shells      [-]
   %   Adr_z   = unscattered photon absorptance along the z-axis       [-]
   %   kabs    = absorption coefficient                                [cm-1]
   %   N       = number of photons
   %   grid    = the mcrt grid: bin widths dz and the annulus areas
   %             dA computed from the bin edges (see mcrt)
   %
   % Outputs:
   %   Adf_rz  = diffuse absorptance (abs. probability / unit volume)  [1/cm3]
   %   Adf_z   = diffuse absorptance per unit depth                    [1/cm]
   %   Adf     = diffuse absorptance (total abs. probability)          [-]
   %   Adr_z   = unscattered absorptance per unit depth                [1/cm]
   %   phi_rz  = internal fluence per unit incident power              [1/cm2]
   %   phi_z   = fluence along the z-axis per unit incident power      [-]
   %
   % Matt Cooper, guycooper@ucla.edu, Dec 2020
   %
   % See also:

   % differential surface area and volume of annular rings
   dA = grid.dA;              % area of each annular ring              [cm^2]
   dz = grid.dz;              % vertical bin widths                    [cm]
   dV = dA .* dz;             % volume per (z, r) bin, nz+1 x nr+1     [cm^3]

   % sum the 2-d arrays into 1-d and 0-d arrays
   Adf_z = sum(Adf_rz, 2);    % Eq. 4.20
   Adf = sum(Adf_z);          % Eq. 4.22

   % absorption (Wang et al. 1995, Sect. 4): divide raw weight by the photon
   % count to get the absorbed fraction of incident power and by each elemental
   % control volume to get a fraction per cm3, and by the depth of each layer to
   % get a fraction per cm.
   %
   % Adf_rz becomes a volume density; Adr_z and Adf_z become per unit depth,
   % integrated over the plane. The direct beam is a pencil at r = 0, so its
   % absorption is also a volume density in the first radial bin only (Adr_rz).

   % convert raw photon weights into fractions of the incident power and
   % then into fractions per cm3 and per cm
   Adf_rz   = Adf_rz ./ dV ./ N;                      % Eq. 4.23      [1/cm3]
   Adr_rz   = Adr_z ./ dV(:,1) ./ N;                  % direct pencil [1/cm3]
   Adr_z    = Adr_z ./ dz ./ N;                       % Eq. 4.24      [1/cm]
   Adf_z    = Adf_z ./ dz ./ N;                       % Eq. 4.25      [1/cm]
   Adf      = Adf ./ N;                               % Eq. 4.27      [-]

   % compute internal fluence (phi_rz) and fluence along the z-axis (phi_z).
   % fluence per incident packet is total absorption over ka. ka times the
   % volume integral of phi_rz returns Adf + Adr. ka times the depth integral
   % of phi_z returns the same.
   phi_rz      = Adf_rz ./ kabs;                      % Eq. 4.28      [1/cm2]
   phi_rz(:,1) = phi_rz(:,1) + Adr_rz ./ kabs;        % Eq. 4.28, direct term
   phi_z       = (Adf_z + Adr_z) ./ kabs;             % Eq. 4.29      [-]
end
