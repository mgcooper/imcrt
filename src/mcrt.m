function RT = mcrt(ka, ks, g, Z, dz, N, varargin)
   %MCRT Monte Carlo radiative transfer through a plane-parallel slab.
   %
   %  RT = mcrt(ka, ks, g, Z, dz, N) launches N photon packets vertically into a
   %  slab of thickness Z with absorption coefficient ka, scattering coefficient
   %  ks, and Henyey-Greenstein asymmetry parameter g.
   %
   %  RT holds the diffuse (df) and direct (dr) reflectance (R), transmittance
   %  (T), absorption (A), and fluence (phi), resolved by radius (r), angle (a),
   %  and depth (z) on a grid with vertical spacing dz. RT.grid holds the grid
   %  and its bin geometry, and RT.N the photon packet count. RT.se holds
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
   %  RT = mcrt(..., Name, Value) sets an option. wmin (1e-4) is the packet
   %  weight below which Russian roulette plays and wrr (10) the survival
   %  odds, 1 in wrr, with the survivor's weight multiplied by wrr. R (2 cm)
   %  is the radial extent of the tallies, dr (0.001 cm) the radial bin
   %  width, and da (pi/60 rad) the angular bin width over the hemisphere;
   %  R/dr and (pi/2)/da must be whole numbers. Every option is a double
   %  scalar; a bad one or an unknown name raises mcrt:input.
   %
   %  The model is a plane-parallel homogeneous slab with a vertical pencil
   %  source and no refractive-index boundary, so it has no Fresnel term and
   %  no refraction. See the Limitations section of README.md for the full
   %  statement and tests/reports/Post-publication corrections.md for the
   %  corrections made after publication.
   %
   % Matt Cooper, guycooper@ucla.edu, Dec 2020
   %
   % See also: buildgrid, binindex, computeReflectance, computeTransmittance,
   % computeAbsorption, mcstderr

   % parse and validate inputs
   opts = parseinputs(ka, ks, g, Z, dz, N, varargin{:});

   % compute optical coefficients
   w = ks/(ka+ks);   % single-scattering albedo          [-]
   a = 1-w;          % co-albedo                         [-]
   c = 1/(ka+ks);    % extinction path length            [cm]

   % precompute 2*pi
   two_pi = 2*pi;

   % roulette settings
   wmin = opts.wmin; % photon weight below which russian roulette plays
   wrr = opts.wrr;   % 1/wrr photons are reinjected (russian roulette)

   % grid settings
   R = opts.R;       % cylindrical detection radius          [cm]
   A = pi/2;         % angular detection radius              [rad]
   dr = opts.dr;     % radial bin width                      [cm]
   da = opts.da;     % angular bin width                     [rad]
   nr = opts.nr;     % radial
   na = opts.na;     % angular
   nz = opts.nz;     % vertical

   % build a grid to calculate observable quantities (eq. 4.1/4.2 Wang)
   grid = buildgrid(R, A, Z, dr, da, dz);

   % initialize output grids with +1 for overflow
   Rdr = 0;                   % reflectance, direct (unscattered)
   Tdr = 0;                   % transmittance, direct (unscattered)
   Adr_z  = zeros(nz+1,1);    % absorption, direct
   Rdf_ra = zeros(na,nr+1);   % reflectance, diffuse
   Tdf_ra = zeros(na,nr+1);   % transmittance, diffuse
   Adf_rz = zeros(nz+1,nr+1); % absorption, diffuse

   % squared weights of the exit tallies, for their standard errors
   Rdf_ss = zeros(na,nr+1);
   Tdf_ss = zeros(na,nr+1);
   Rdr_ss = 0;
   Tdr_ss = 0;

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
         % below, which clamps the angular index to acos(1) = 0 on-axis.
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
         % its weight or dies. A packet can survive below wmin, since the
         % loop runs on wt > 0 it plays again on the next step.
         if wt < wmin
            wt = roulette(wt, wrr);
         end
      end
   end

   % convert the photon counts to SI units (Wang et al. 1995, Sect. 4).
   % computeReflectance and computeTransmittance sum the tallies and divide by
   % the bin measures and N; computeAbsorption does the same and computes
   % fluence; standard errors are computed from the squared weights.
   [Rdf_ra, Rdf_r, Rdf_a, Rdf, Rdr, Rt, seR] = computeReflectance( ...
      Rdf_ra, Rdr, Rdf_ss, Rdr_ss, N, grid);

   [Tdf_ra, Tdf_r, Tdf_a, Tdf, Tdr, Tt, seT] = computeTransmittance( ...
      Tdf_ra, Tdr, Tdf_ss, Tdr_ss, N, grid);

   [Adf_rz, Adf_z, Adf, Adr_z, phi_rz, phi_z] = computeAbsorption( ...
      Adf_rz, Adr_z, ka, N, grid);

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

   % assign standard errors and the packet count N
   RT.se = struct('Rdf_ra', seR.ra, 'Rdf_r', seR.r, 'Rdf_a', seR.a, ...
      'Rdf', seR.df, 'Rdr', seR.dr, 'Rt', seR.t, 'Tdf_ra', seT.ra, ...
      'Tdf_r', seT.r, 'Tdf_a', seT.a, 'Tdf', seT.df, 'Tdr', seT.dr, ...
      'Tt', seT.t);
   RT.N = N;

   % return the grid
   RT.grid = grid;

end

function opts = parseinputs(ka, ks, g, Z, dz, N, varargin)
   %PARSEINPUTS Validate the required inputs and parse the options.
   %
   %  opts = parseinputs(ka, ks, g, Z, dz, N, Name, Value, ...) checks the
   %  required inputs and the options wmin, wrr, R, dr, and da through an
   %  inputParser whose validators are validateattributes calls, fills the
   %  options with their defaults, and returns them in opts with the whole
   %  bin counts nr, na, and nz. Every failure, including an unknown option
   %  name, raises mcrt:input with the parser's message; a fractional bin
   %  count raises buildgrid:nonintegral from wholebins, the same rule
   %  buildgrid applies.

   double = {'double'};
   scalar = {'scalar', 'real', 'finite'};

   p = inputParser;
   p.FunctionName = 'mcrt';

   % a misspelled option must fail, not match a unique prefix
   p.PartialMatching = false;

   % add required arguments
   addRequired(p, 'ka', @(x) validateattributes(x, double, ...
      [scalar, {'positive'}], 'mcrt', 'ka'));
   addRequired(p, 'ks', @(x) validateattributes(x, double, ...
      [scalar, {'nonnegative'}], 'mcrt', 'ks'));
   addRequired(p, 'g', @(x) validateattributes(x, double, ...
      [scalar, {'>=', -1, '<=', 1}], 'mcrt', 'g'));
   addRequired(p, 'Z', @(x) validateattributes(x, double, ...
      [scalar, {'positive'}], 'mcrt', 'Z'));
   addRequired(p, 'dz', @(x) validateattributes(x, double, ...
      [scalar, {'positive'}], 'mcrt', 'dz'));
   addRequired(p, 'N', @(x) validateattributes(x, double, ...
      [scalar, {'integer', 'positive'}], 'mcrt', 'N'));

   % add optional name-value arguments
   addParameter(p, 'wmin', 1e-4, @(x) validateattributes(x, double, ...
      [scalar, {'>', 0, '<', 1}], 'mcrt', 'wmin'));
   addParameter(p, 'wrr', 10, @(x) validateattributes(x, double, ...
      [scalar, {'>', 1}], 'mcrt', 'wrr'));
   addParameter(p, 'R', 2, @(x) validateattributes(x, double, ...
      [scalar, {'positive'}], 'mcrt', 'R'));
   addParameter(p, 'dr', 0.001, @(x) validateattributes(x, double, ...
      [scalar, {'positive'}], 'mcrt', 'dr'));
   addParameter(p, 'da', pi/60, @(x) validateattributes(x, double, ...
      [scalar, {'>', 0, '<=', pi/2}], 'mcrt', 'da'));

   % one identifier for every bad input; the parser's message names it
   try
      parse(p, ka, ks, g, Z, dz, N, varargin{:});
   catch cause
      error('mcrt:input', '%s', cause.message);
   end
   opts = p.Results;

   % two finite coefficients can overflow when added
   if ~isfinite(ka + ks)
      error('mcrt:input', 'ka + ks must be finite');
   end

   % every extent must hold a whole number of bins, checked here so the
   % failure comes before any packet runs
   opts.nr = wholebins(opts.R, opts.dr, 'R/dr');
   opts.na = wholebins(pi/2, opts.da, 'A/da');
   opts.nz = wholebins(Z, dz, 'Z/dz');
end
