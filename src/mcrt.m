function RT = mcrt(ka, ks, g, Z, dz, N)
   %MCRT Monte Carlo radiative transfer through a plane-parallel slab.
   %
   %  RT = mcrt(ka, ks, g, Z, dz, N) launches N photon packets vertically into a
   %  slab of thickness Z with absorption coefficient ka, scattering coefficient
   %  ks, and Henyey-Greenstein asymmetry parameter g.
   %
   %  RT holds the diffuse (df) and direct (dr) reflectance (R), transmittance
   %  (T), absorption (A), and fluence (phi), resolved by radius (r), angle (a),
   %  and depth (z) on a grid with vertical spacing dz. RT.grid holds the grid
   %  and its bin geometry, and RT.grid.N the photon packet count. RT.se holds
   %  the standard error of every reflectance and transmittance output from
   %  one run (sum of squared weights, see mcstderr, which also states the
   %  covariance between exclusive bins). Absorption and fluence carry no
   %  per-run error because a packet deposits in one bin many times; use the
   %  spread over seeded runs for those.
   %
   %  ka is a finite positive scalar (the fluence is absorption over ka).
   %  ks a finite nonnegative scalar with a finite sum.
   %  g lies in [-1, 1], with -1 backscattering and 1 forward scattering.
   %  Z and dz are positive and Z/dz must be a whole number.
   %  N is a positive integer representing the number of photon packets.
   %  All are double scalars. A bad input raises mcrt:input, or
   %  buildgrid:nonintegral for Z/dz and buildgrid:width for a dz whose
   %  half underflows, before any packet runs.
   %
   % Matt Cooper, guycooper@ucla.edu, Dec 2020
   %
   % See also: buildgrid, binindex, computeReflectance, computeTransmittance,
   % computeAbsorption, mcstderr

   % Check inputs to prevent index errors or silent NaN's inside the main loop.
   mustbe(isnumber(ka) && ka > 0, 'ka', 'a finite positive scalar');
   mustbe(isnumber(ks) && ks >= 0, 'ks', 'a finite nonnegative scalar');
   mustbe(isfinite(ka + ks), 'ka + ks', 'finite');
   mustbe(isnumber(g) && g >= -1 && g <= 1, 'g', 'a scalar in [-1, 1]');
   mustbe(isnumber(Z) && Z > 0, 'Z', 'a finite positive scalar');
   mustbe(isnumber(dz) && dz > 0, 'dz', 'a finite positive scalar');
   mustbe(isnumber(N) && N >= 1 && N == fix(N), 'N', 'a positive integer');

   % compute optical coefficients
   w = ks/(ka+ks);   % single-scattering albedo          [-]
   a = 1-w;          % co-albedo                         [-]
   c = 1/(ka+ks);    % extinction path length            [cm]

   % precompute 2*pi
   two_pi = 2*pi;

   % default settings
   wmin = 1e-4;      % photon weight below which russian roulette plays
   wrr = 10;         % 1/wrr photons are reinjected (russian roulette)

   % grid settings
   R = 2;            % cylindrical detection radius          [cm]
   A = pi/2;         % angular detection radius              [rad]
   dr = 0.001;       % radial bin width                      [cm]
   da = A/30;        % angular bin width                     [rad]
   nr = round(R/dr); % radial
   na = round(A/da); % angular
   nz = round(Z/dz); % vertical

   % Build a grid to calculate observable quantities (eq. 4.1/4.2 Wang). It
   % is built before the loop so a fractional Z/dz fails before any packet
   % runs and the random stream is untouched. The packet count joins the
   % grid so every consumer of the tallies reads one struct.
   grid = buildgrid(R, A, Z, dr, da, dz);
   grid.N = N;

   % initialize output grids with +1 for overflow
   Adf_rz = zeros(nz+1,nr+1); % absorption, diffuse
   Adr_z = zeros(nz+1,1);     % absorption, direct
   Tdf_ra = zeros(na,nr+1);   % transmittance, diffuse
   Rdf_ra = zeros(na,nr+1);   % reflectance, diffuse
   Tdr = 0;                   % transmittance, direct (unscattered)
   Rdr = 0;                   % reflectance, direct (unscattered)

   % squared weights of the exit tallies, for their standard errors
   Tdf_ss = zeros(na,nr+1);
   Rdf_ss = zeros(na,nr+1);
   Tdr_ss = 0;
   Rdr_ss = 0;

   % monte carlo
   for n = 1:N
      wt = 1;  % new photon, weight = 1
      x = 0;   % generated at position (0,0,0):
      y = 0;
      z = 0;
      ux = 0;  % with trajectory 0,0,1 (mu_x, mu_y, mu_z)
      uy = 0;
      uz = 1;  % = 1 for vertical, = 1-2*rand for isotropic
      ns = 0;  % number of scattering events

      % Propagate the packet until roulette kills it or it exits the slab.
      while wt > 0
         l = -c*log(rand); % path length
         x = x+ux*l;       % new x-position
         y = y+uy*l;       % new y-position
         z = z+uz*l;       % new z-position

         % grid indices
         ir = ceil(sqrt(x*x+y*y)/dr); % radial index
         iz = ceil(z/dz);             % vertical index

         % clamps catch r=0 and z=0, which ceil maps to bin 0. binindex does the
         % same but two calls per step measured 10-13% of the loop (tests/perf),
         % so these two are inlined and binindex is used for the exit angle
         % below, which clamps the angular index to acos(1) = 0 on axis.
         if ir<1; ir = 1; end         % radial on-axis
         if iz<1; iz = 1; end         % vertical at the surface

         % overflow catches photons beyond the grid extent.
         if ir>nr; ir = nr+1; end     % radial overflow
         if iz>nz; iz = nz+1; end     % vertical overflow

         % score transmittance / reflectance.
         if z>Z                                       % transmittance
            iu = binindex(acos(uz), da, na);          % angular index
            if ns==0
               Tdr = Tdr+wt;                          % direct
               Tdr_ss = Tdr_ss+wt*wt;
            else
               Tdf_ra(iu,ir) = Tdf_ra(iu,ir)+wt;      % diffuse
               Tdf_ss(iu,ir) = Tdf_ss(iu,ir)+wt*wt;
            end
            break % photon escapes
         end
         if z<0                                       % reflection
            iu = binindex(acos(-uz), da, na);         % angular index
            if ns==0
               % unreachable for the vertical source
               % (no specular term); kept for the
               % isotropic source noted at the
               % uz = 1 line
               Rdr = Rdr+wt;                          % direct
               Rdr_ss = Rdr_ss+wt*wt;
            else
               Rdf_ra(iu,ir) = Rdf_ra(iu,ir)+wt;      % diffuse
               Rdf_ss(iu,ir) = Rdf_ss(iu,ir)+wt*wt;
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

         % absorption and scattering (henyey-greenstein with azimuthal symmetry)
         wt = wt*w; ns = ns+1;
         us = hgcos(g);    % polar scattering cosine
         ps = two_pi*rand; % azimuth angle, phi_s

         % new direction cosines from the polar cosine and the azimuth
         [ux,uy,uz] = chgdir(ux,uy,uz,us,ps);

         % russian roulette below wmin: the packet survives with wrr times
         % its weight or dies. A survivor can still sit below wmin, so the
         % loop runs on wt > 0 and it plays again on the next step.
         if wt < wmin
            wt = roulette(wt, wrr);
         end
      end
   end

   % Convert the photon counts to SI units (Wang et al. 1995, Sect. 4):
   % computeReflectance and computeTransmittance sum the tallies and divide by
   % the bin measures and N; computeAbsorption does the same for absorption
   % and computes fluence. The squared weights give each exit output its
   % standard error.
   [Rdf_ra, Rdf_r, Rdf_a, Rdf, Rdr, Rt, seR] = computeReflectance( ...
      Rdf_ra, Rdr, Rdf_ss, Rdr_ss, grid);

   [Tdf_ra, Tdf_r, Tdf_a, Tdf, Tdr, Tt, seT] = computeTransmittance( ...
      Tdf_ra, Tdr, Tdf_ss, Tdr_ss, grid);

   [Adf_rz, Adf_z, Adf, Adr_z, phi_rz, phi_z] = computeAbsorption( ...
      Adf_rz, Adr_z, ka, grid);

   % arrange the output
   RT.Rdf_ra = Rdf_ra;
   RT.Rdf_r = Rdf_r;
   RT.Rdf_a = Rdf_a;
   RT.Rdf = Rdf;
   RT.Rdr = Rdr;
   RT.Rt = Rt;
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

   % standard errors of the exit outputs from this run
   RT.se = struct('Rdf_ra', seR.ra, 'Rdf_r', seR.r, 'Rdf_a', seR.a, ...
      'Rdf', seR.df, 'Rdr', seR.dr, 'Rt', seR.t, 'Tdf_ra', seT.ra, ...
      'Tdf_r', seT.r, 'Tdf_a', seT.a, 'Tdf', seT.df, 'Tdr', seT.dr, ...
      'Tt', seT.t);

   % return the grid
   RT.grid = grid;

end

function tf = isnumber(x)
   % A finite real double scalar. Integer classes are refused because
   % their arithmetic saturates; single is refused because buildgrid's
   % whole-bin tolerance is set for double precision.
   tf = isa(x, 'double') && isscalar(x) && isreal(x) && isfinite(x);
end

function mustbe(cond, name, what)
   % Raise mcrt:input naming the argument when its check fails.
   assert(cond, 'mcrt:input', '%s must be %s', name, what);
end
