function RT = mcrt(ka,ks,g,Z,dz,N)
   %MCRT Monte Carlo radiative transfer through a plane-parallel slab.
   %  RT = mcrt(ka, ks, g, Z, dz, N) launches N photon packets vertically
   %  into a slab of thickness Z with absorption coefficient ka, scattering
   %  coefficient ks (inverse of the length unit of Z), and Henyey-Greenstein
   %  asymmetry g. RT holds the reflectance (R), transmittance (T),
   %  absorption (A), and fluence (phi) tallies, diffuse (df) and direct
   %  (dr), resolved by radius (r), angle (a), and depth (z) on a grid with
   %  vertical spacing dz, plus that grid in RT.grid.

   % optical coefficients
   w = ks/(ka+ks); % single-scattering albedo          [-]
   a = 1-w;        % co-albedo                         [-]
   c = 1/(ka+ks);  % extinction path length            [cm]

   % henyey-greenstein terms (pre-computed so it runs fast)
   hg1 = 1/(2*g);
   hg2 = (1+g^2);
   hg3 = (1-g^2);
   hg4 = 1+g;
   hg5 = -2*g;
   two_pi = 2*pi;

   % default settings
   wmin = 1e-4; % photon weight below which russian roulette plays
   wrr = 10;    % 1/wrr photons are reinjected (russian roulette)

   % grid settings
   R = 2;            % cylindrical detection radius          [cm]
   A = pi/2;         % angular detection radius              [rad]
   dr = 0.001;       % radial bin width                      [cm]
   da = A/30;        % angular bin width                     [rad]
   nr = round(R/dr); % radial
   na = round(A/da); % angular
   nz = round(Z/dz); % vertical

   % initialize output grids with +1 for overflow
   Adf_rz = zeros(nz+1,nr+1); % absorption, diffuse
   Adr_z = zeros(nz+1,1);     % absorption, direct
   Tdf_ra = zeros(na,nr+1);   % transmittance, diffuse
   Rdf_ra = zeros(na,nr+1);   % reflectance, diffuse
   Tdr = 0;                   % transmittance, direct (unscattered)
   Rdr = 0;                   % reflectance, direct (unscattered)

   % monte carlo
   for n = 1:N
      wt = 1; % new photon, weight = 1
      x = 0;  % generated at position (0,0,0):
      y = 0;
      z = 0;
      ux = 0; % with trajectory 0,0,1 (mu_x, mu_y, mu_z)
      uy = 0;
      uz = 1; % = 1 for vertical, = 1-2*rand for isotropic
      ns = 0; % number of scattering events

      % Propagate the packet until roulette kills it or it exits the slab.
      while wt > 0
         l = -c*log(rand); % path length
         x = x+ux*l;       % new x-position
         y = y+uy*l;       % new y-position
         z = z+uz*l;       % new z-position

         % grid indices. The lower clamps catch r = 0 and z = 0 exactly,
         % which ceil maps to bin 0 (defect I).
         ir = ceil(sqrt(x*x+y*y)/dr); % radial index
         iz = ceil(z/dz);             % vertical index
         if ir<1; ir = 1; end         % radial on-axis
         if ir>nr; ir = nr+1; end     % radial overflow
         if iz<1; iz = 1; end         % vertical at the surface
         if iz>nz; iz = nz+1; end     % vertical overflow

         % score transmittance / reflectance. Direct means unscattered
         % (ns==0) only: a scattered grazing exit is diffuse (defect A). The
         % angular index is clamped because acos(1) = 0 gives bin 0 on axis
         % and rounding at grazing can give bin na+1 (defect J).
         if z>Z                       % transmittance
            iu = ceil(acos(uz)/da);   % angular index
            if iu<1; iu = 1; end      % on axis
            if iu>na; iu = na; end    % grazing
            if ns==0
               Tdr = Tdr+wt; % direct
            else
               Tdf_ra(iu,ir) = Tdf_ra(iu,ir)+wt; % diffuse
            end
            break % photon escapes
         end
         if z<0                       % reflection
            iu = ceil(acos(-uz)/da);  % angular index
            if iu<1; iu = 1; end      % on axis
            if iu>na; iu = na; end    % grazing
            if ns==0
               % unreachable for the vertical source (no specular term);
               % kept for the isotropic source noted at the uz = 1 line
               Rdr = Rdr+wt; % direct
            else
               Rdf_ra(iu,ir) = Rdf_ra(iu,ir)+wt; % diffuse
            end
            break % photon escapes
         end

         % score absorption/fluence
         awt = a*wt;
         if ns==0
            Adr_z(iz) = Adr_z(iz)+awt; % direct
         else
            Adf_rz(iz,ir) = Adf_rz(iz,ir)+awt; % diffuse
         end

         % absorption and scattering by ice (henyey-greenstein with azimuthal
         % symmetry)
         wt = wt*w; ns = ns+1;
         if g == 0
            us = 1-2*rand;
         else
            us = hg1*(hg2-(hg3/(hg4+hg5*rand))^2);
         end
         ps = two_pi*rand; % azimuth angle, phi_s
         % new direction cosines: an inline copy of chgdir, because a function
         % call in this loop costs more than the scattering itself. The
         % temporaries uxn and uyn keep uy from reading the updated ux (B).
         % tests/testChgdirOracle.m evaluates this block against chgdir.
         sinth = sqrt(1-uz*uz);
         sinths = sqrt(1-us*us);
         cps = cos(ps);
         sps = sin(ps);
         if sinth < 1e-12
            % initial direction straight up or down (sin(theta)=0)
            ux = sinths*cps;
            uy = sinths*sps;
            uz = sign(uz)*us;
         else
            % if initial direction not straight up or straight down
            uxn = sinths/sinth*(ux*uz*cps-uy*sps)+ux*us;
            uyn = sinths/sinth*(uy*uz*cps+ux*sps)+uy*us;
            uz = -sinths*sinth*cps+uz*us;
            ux = uxn;
            uy = uyn;
         end
         % russian roulette (Wang et al. 1995, Sect. 3.9): a packet below wmin
         % survives with probability 1/wrr carrying wrr times its weight, or
         % dies. A survivor can still sit below wmin, so the loop runs on
         % wt > 0 and it plays again next step instead of being dropped (N).
         if wt < wmin
            if rand < 1/wrr
               wt = wt*wrr;
            else
               wt = 0;
            end
         end
      end
   end

   % build a grid to calculate observable quantities (eq. 4.1/4.2 Wang)
   [ri,ai,zi,dr,da,dz] = buildgrid(R,A,Z,dr,da,dz);

   dsr = 2*pi.*sin(ai).*da; % solid angle per angular bin       [sr]
   dA = 2*pi.*ri.*dr;       % annulus area per radial bin       [cm^2]
   dV = dA.*dz;             % volume per (z, r) bin, nz+1 x nr+1 [cm^3]
   cosa = cos(ai);

   % sum the 2-d arrays into 1-d and 0-d arrays (R=reflection, T=transmission)
   Rdf_r = sum(Rdf_ra,1); % Eq. 4.3
   Rdf_a = sum(Rdf_ra,2); % Eq. 4.4
   Tdf_r = sum(Tdf_ra,1); % Eq. 4.5
   Tdf_a = sum(Tdf_ra,2); % Eq. 4.6
   Rdf = sum(Rdf_r);      % Eq. 4.7
   Tdf = sum(Tdf_r);      % Eq. 4.8

   % convert the photon counts to SI units
   Rdf_ra = Rdf_ra./(dA.*dsr.*cosa.*N); % Eq. 4.9
   Tdf_ra = Tdf_ra./(dA.*dsr.*cosa.*N); % Eq. 4.10
   Rdf_r = Rdf_r./(dA.*N);              % Eq. 4.13
   Tdf_r = Tdf_r./(dA.*N);              % Eq. 4.14
   Rdf_a = Rdf_a./(dsr.*N);             % Eq. 4.15
   Tdf_a = Tdf_a./(dsr.*N);             % Eq. 4.16
   Rdr = Rdr/N;                         % + Rsp/N
   Tdr = Tdr/N;
   Rdf = Rdf/N; % Eq. 4.17
   Tdf = Tdf/N; % Eq. 4.18
   Tt = Tdf+Tdr;

   % Absorption (Wang et al. 1995, Sect. 4). Adf_rz becomes a volume
   % density; Adr_z and Adf_z become per unit depth, integrated over the
   % plane. The direct beam is a pencil at r = 0, so its absorption is also
   % a volume density in the first radial bin only (Adr_rz).
   Adf_z = sum(Adf_rz,2);      % Eq. 4.20
   Adf = sum(Adf_z);           % Eq. 4.22
   Adf_rz = Adf_rz./dV./N;     % Eq. 4.23                       [1/cm^3]
   Adr_rz = Adr_z./dV(:,1)./N; % direct pencil in radial bin 1  [1/cm^3]
   Adr_z = Adr_z./dz./N;       % Eq. 4.24                       [1/cm]
   Adf_z = Adf_z./dz./N;       % Eq. 4.25                       [1/cm]
   Adf = Adf./N;               % Eq. 4.27                       [-]
   % phi_rz    = Adf_rz./ka;                           % Eq. 4.28
   % phi_z     = Adf_z./ka;                            % Eq. 4.29

   % fluence per incident packet is total absorption over ka. Each sum adds
   % like units: phi_rz from the two volume densities, phi_z from the two
   % per-depth totals. ka times the volume integral of phi_rz returns
   % Adf + Adr. ka times the depth integral of phi_z returns the same.
   phi_rz = Adf_rz./ka;                    % Eq. 4.28           [1/cm^2]
   phi_rz(:,1) = phi_rz(:,1)+Adr_rz./ka;   % Eq. 4.28, direct term
   phi_z = (Adf_z+Adr_z)./ka;              % Eq. 4.29           [-]

   % arrange the output
   RT.Rdf_ra = Rdf_ra;
   RT.Rdf_r = Rdf_r;
   RT.Rdf_a = Rdf_a;
   RT.Rdf = Rdf;
   RT.Rdr = Rdr;
   RT.Tdf_ra = Tdf_ra;
   RT.Tdf_r = Tdf_r;
   RT.Tdf_a = Tdf_a;
   RT.Tdf = Tdf;
   RT.Tdr = Tdr;
   RT.Tt = Tt;

   RT.Adf_z = Adf_z;
   RT.Adf = Adf;
   RT.Adf_rz = Adf_rz;
   RT.Adr_z = Adr_z;
   RT.phi_rz = phi_rz;
   RT.phi_z = phi_z;


   % return the grid
   RT.grid.ri = ri;
   RT.grid.ai = ai;
   RT.grid.zi = zi;
   RT.grid.dr = dr;
   RT.grid.da = da;
   RT.grid.dz = dz;

end
