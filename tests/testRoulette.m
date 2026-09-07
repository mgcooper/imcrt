function tests = testRoulette
   % Tests for Russian roulette (defect N) in src/mcrt.m. A packet below
   % wmin must either die or survive with wrr times its weight, and a
   % survivor that is still below wmin must keep playing. Discarding it
   % loses weight without compensation.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % The fixture adds src and restores path and rng.
   kernelfixture(testCase);
end

function testDeepAbsorbingSlabConservesWeight(testCase)
   % Albedo 0.05 in a slab of optical depth 100: nearly every packet is
   % absorbed down to roulette, and every boosted survivor (wt*wrr) still
   % sits below wmin. The roulette noise on the balance is about 2e-7
   % relative at N=1e4, and dropping those survivors biases it by about
   % 6e-6, so 1e-6 (5 sigma) separates the two.
   rng(11, 'twister');
   RT = mcrt(95, 5, 0.5, 1, 0.05, 1e4);
   returned = RT.Rdf + RT.Rdr + RT.Tdf + RT.Tdr + RT.Adf ...
      + sum(RT.Adr_z .* RT.grid.dz);
   expected = 1;
   testCase.verifyEqual(returned, expected, 'AbsTol', 1e-6);
end

function testRouletteBlockIsTerminateOrBoost(testCase)
   % Evaluate the kernel's roulette block, read from src/mcrt.m, on a
   % packet below wmin. Every outcome must be 0 or wt*wrr, and the mean
   % must equal the input weight within 5 standard errors (unbiased).
   block = kernellines('% russian roulette', '');
   v = struct('wmin', 1e-4, 'wrr', 10, 'wt', 5e-5);
   rng(11, 'twister');
   ndraw = 1e5;
   returned = zeros(1, ndraw);
   for n = 1:ndraw
      w = runblock(block, v);
      returned(n) = w.wt;
   end
   testCase.verifyTrue(all(returned == 0 | returned == v.wt*v.wrr));
   expected = v.wt;
   tol = 5*std(returned)/sqrt(ndraw);
   testCase.verifyEqual(mean(returned), expected, 'AbsTol', tol);
end
