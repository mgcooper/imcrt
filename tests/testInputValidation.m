function tests = testInputValidation
   % Tests for the input checks of src/mcrt.m: every bad argument raises
   % mcrt:input before the loop, a fractional Z/dz raises
   % buildgrid:nonintegral, the accepted edge values run, and the options
   % take effect or are refused by name.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % The fixture adds src and restores path and rng.
   kernelfixture(testCase);
end

function testBadArgumentsRaise(testCase)
   % One bad value per argument, plus the coupled ones: zero absorption
   % (fluence divides by ka), an overflowing sum, g outside [-1, 1], a
   % fractional or complex N, a vector ka, an integer class, single.
   good = {10, 90, 0.75, 0.02, 0.001, 100};
   bad = {{-1, 90, 0.75, 0.02, 0.001, 100}, ...
      {10, -1, 0.75, 0.02, 0.001, 100}, ...
      {0, 90, 0.75, 0.02, 0.001, 100}, ...
      {realmax, realmax, 0.75, 0.02, 0.001, 100}, ...
      {uint8(10), 90, 0.75, 0.02, 0.001, 100}, ...
      {single(10), 90, 0.75, 0.02, 0.001, 100}, ...
      {10, 90, 1.5, 0.02, 0.001, 100}, ...
      {10, 90, -1.5, 0.02, 0.001, 100}, ...
      {10, 90, 0.75, 0, 0.001, 100}, ...
      {10, 90, 0.75, 0.02, -0.001, 100}, ...
      {10, 90, 0.75, 0.02, 0.001, 0}, ...
      {10, 90, 0.75, 0.02, 0.001, 2.5}, ...
      {10, 90, 0.75, 0.02, 0.001, 2 + 1i}, ...
      {[10 10], 90, 0.75, 0.02, 0.001, 100}, ...
      {Inf, 90, 0.75, 0.02, 0.001, 100}, ...
      {10, NaN, 0.75, 0.02, 0.001, 100}};
   for n = 1:numel(bad)
      args = bad{n};
      testCase.verifyError(@() mcrt(args{:}), 'mcrt:input', ...
         sprintf('bad set %d', n));
   end
   rng(1, 'twister');
   RT = mcrt(good{:});
   testCase.verifyEqual(RT.Rdf + RT.Rdr + RT.Tdf + RT.Tdr + RT.Adf ...
      + sum(RT.Adr_z.*RT.grid.dz), 1, 'AbsTol', 1e-6);
end

function testFractionalDepthCountRaises(testCase)
   % Z/dz = 2.5 reaches buildgrid, which names the ratio, before any
   % packet runs: the random stream is where it was.
   rng(5, 'twister');
   before = rand;
   rng(5, 'twister');
   testCase.verifyError(@() mcrt(10, 90, 0.75, 0.025, 0.01, 1e6), ...
      'buildgrid:nonintegral');
   returned = rand;
   expected = before;
   testCase.verifyEqual(returned, expected);
end

function testEdgeValuesRun(testCase)
   % g = -1 (reversal), g = 1 (forward), ks = 0 (no scattering), an
   % extent of exactly one bin, and a Z/dz a rounding error below a whole
   % number all run and conserve weight.
   sets = {{10, 90, -1, 0.02, 0.001, 100}, {10, 90, 1, 0.02, 0.001, 100}, ...
      {10, 0, 0.5, 0.02, 0.001, 100}, {10, 90, 0.5, 0.01, 0.01, 100}, ...
      {10, 90, 0.5, 0.02*(1 - 1e-10), 0.001, 100}};
   for n = 1:numel(sets)
      args = sets{n};
      rng(n, 'twister');
      RT = mcrt(args{:});
      returned = RT.Rdf + RT.Rdr + RT.Tdf + RT.Tdr + RT.Adf ...
         + sum(RT.Adr_z.*RT.grid.dz);
      testCase.verifyEqual(returned, 1, 'AbsTol', 1e-6, sprintf('set %d', n));
   end
end

function testOptions(testCase)
   % Each option changes what it names: da sets the angular bin count, R
   % and dr the radial count, wmin and wrr the roulette. A bad option value
   % raises mcrt:input by name, a fractional R/dr raises
   % buildgrid:nonintegral, and an unknown name is an error.
   rng(2, 'twister');
   RT = mcrt(10, 90, 0.75, 0.02, 0.001, 100, 'da', pi/30, 'R', 1, ...
      'dr', 0.01, 'wmin', 1e-3, 'wrr', 5);
   returned = {numel(RT.grid.ai), numel(RT.grid.ri), RT.grid.N, ...
      RT.Rdf + RT.Rdr + RT.Tdf + RT.Tdr + RT.Adf + sum(RT.Adr_z.*RT.grid.dz)};
   expected = {15, 101, 100, 1};
   testCase.verifyEqual(returned, expected, 'AbsTol', 1e-6);
   bad = {{'wmin', 0}, {'wmin', 1}, {'wrr', 1}, {'R', -1}, {'dr', 0}, ...
      {'da', pi}, {'da', single(0.1)}};
   for n = 1:numel(bad)
      args = bad{n};
      testCase.verifyError(@() mcrt(10, 90, 0.75, 0.02, 0.001, 100, ...
         args{:}), 'mcrt:input', sprintf('bad option %d', n));
   end
   testCase.verifyError(@() mcrt(10, 90, 0.75, 0.02, 0.001, 100, ...
      'R', 1, 'dr', 0.3), 'buildgrid:nonintegral');
   testCase.verifyError(@() mcrt(10, 90, 0.75, 0.02, 0.001, 100, ...
      'nosuch', 1), 'mcrt:input');
   % a unique prefix of an option name is not the option
   testCase.verifyError(@() mcrt(10, 90, 0.75, 0.02, 0.001, 100, ...
      'wmi', 1e-3), 'mcrt:input');
   % The roulette options take effect: in a deep absorbing slab at albedo
   % 0.05 a packet drops below 1e-2 after two scatters and below 1e-4
   % after four, so wmin = 1e-2 draws its roulette numbers at different
   % steps from the default and the seeded tallies diverge; a different
   % wrr at the same threshold diverges again.
   rng(3, 'twister');
   base = mcrt(95, 5, 0.5, 1, 0.05, 2e3);
   rng(3, 'twister');
   early = mcrt(95, 5, 0.5, 1, 0.05, 2e3, 'wmin', 1e-2);
   rng(3, 'twister');
   odds = mcrt(95, 5, 0.5, 1, 0.05, 2e3, 'wmin', 1e-2, 'wrr', 5);
   returned = [isequal(base.Adf, early.Adf), isequal(early.Adf, odds.Adf)];
   expected = [false, false];
   testCase.verifyEqual(returned, expected);
end
