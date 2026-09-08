function [V, RTs] = mcrt_verify(casename)
   %MCRT_VERIFY Verify src/mcrt.m against its published references and print
   % a tolerance-based PASS/FAIL table.
   %
   % [V, RTs] = mcrt_verify(casename) runs one case and returns its verdict
   % table V and the M runs RTs. casename defaults to 'reflect':
   %   'reflect'  van de Hulst (1980) Vol. 2, p. 435, Table 35. It checks the
   %              hemispherical and angular reflectance and transmittance of
   %              a tau = 2 slab (albedo 0.9, g = 0.75). Every row must pass.
   %   'fluence'  Wang et al. (1995) Fig. 4: internal fluence at albedo 0.999.
   %              No numeric reference exists in this repository, so the table
   %              holds self-consistency checks only.
   % Each case runs standalone. Its run count, packets, and seeds come from
   % tests/verify/verifycases.m and are sized for seconds on a laptop. Larger
   % runs use the same case table and verdict functions. To add a case, add
   % an entry to verifycases.m. Add a verdict function for its reference.
   % Add a branch to each switch block below.
   if nargin < 1 || isempty(casename)
      casename = 'reflect';
   end

   % Resolve paths from this file so the function runs from any directory.
   here = fileparts(mfilename('fullpath'));
   run(fullfile(here, 'Setup.m'));
   addpath(fullfile(here, 'tests', 'verify'));

   % Run the M seeded, independent runs; the spread over runs is the error.
   c = verifycases(casename);
   RTs = cell(1, c.M);
   for n = 1:c.M
      rng(c.seeds(n), 'twister');
      RTs{n} = mcrt(c.ka, c.ks, c.g, c.Z, c.dz, c.N);
   end

   % Verdict table. A statistical row passes when |t| <= tmax with M - 1
   % degrees of freedom. The Rdr, energy, and surface-fluence rows use the
   % exact rules stated in the verdict functions, and their t is NaN.
   switch casename
      case 'reflect'
         V = vdhverdict(RTs, vdhtable35(), c.tmax);
      case 'fluence'
         V = fluenceverdict(RTs, c, c.tmax);
   end
   npass = printverdict(V);
   nrow = numel(V.verdict);
   if npass == nrow
      overall = 'PASS';
   else
      overall = 'FAIL';
   end
   fprintf(['VERDICT %s: %s (%d of %d rows pass; M=%d runs of N=%g, ' ...
      'tmax=%g)\n'], casename, overall, npass, nrow, c.M, c.N, c.tmax);

   % Plots: angular tables against the reference, or fluence against depth.
   verifyplot(RTs{1}, casename);
end
