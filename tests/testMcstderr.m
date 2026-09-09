function tests = testMcstderr
   % Tests for the per-run standard errors: mcstderr on known contributions,
   % and the kernel's RT.se against the spread of many seeded runs, which
   % is the independent estimate of the same error.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % The fixture adds src and restores path and rng.
   kernelfixture(testCase);
end

function testFormulaOnKnownContributions(testCase)
   % Four packets, two of which land with weight 1: the contributions are
   % [1 0 0 1], mean 0.5, unbiased variance 1/3, so the standard error of
   % the mean is sqrt(1/12). Weighted hits use their squares. One packet
   % gives NaN, and a negative round-off inside the root is clamped to
   % zero (every packet in the bin with equal weight, so no spread).
   returned = [mcstderr(2, 2, 4), mcstderr(0.5^2 + 0.25^2, 0.75, 4), ...
      mcstderr([2 4], [2 4], 4)];
   expected = [sqrt(1/12), sqrt((0.3125 - 0.75^2/4)/12), ...
      [sqrt(1/12), 0]];
   testCase.verifyEqual(returned, expected, 'AbsTol', 1e-15);
   testCase.verifyTrue(isnan(mcstderr(1, 1, 1)));
end

function testPerRunErrorMatchesRunSpread(testCase)
   % One hundred seeded runs of the van de Hulst case at N = 2e4: the mean
   % of the per-run standard errors of Rdf, Tt, Tdr, and the first three
   % angular bins must match the standard deviation of those values over
   % the runs within 15%, the sampling error of a spread from 100 runs
   % being about 7%.
   % The same runs check the few-run angular path: the per-run error
   % propagated through the interpolation weights, averaged over runs,
   % must match the spread of the interpolated values within 15%.
   c = mcrtcases();
   c = c(strcmp({c.name}, 'vdh_reflectance'));
   ref = vdhtable35();
   M = 100;
   N = 2e4;
   vals = zeros(M, 6);
   ses = zeros(M, 6);
   RTs = cell(1, M);
   angse = zeros(M, 12);
   for n = 1:M
      rng(1000 + n, 'twister');
      RTs{n} = mcrt(c.ka, c.ks, c.g, c.Z, c.dz, N);
      RT = RTs{n};
      vals(n, :) = [RT.Rdf, RT.Tt, RT.Tdr, RT.Tdf_a(1:3)'];
      ses(n, :) = [RT.se.Rdf, RT.se.Tt, RT.se.Tdr, RT.se.Tdf_a(1:3)'];
      [~, ~, angse(n, :)] = vdhangular(RTs(n), ref, true);
   end
   returned = mean(ses, 1)./std(vals, 0, 1);
   expected = ones(1, 6);
   testCase.verifyEqual(returned, expected, 'AbsTol', 0.15, ...
      mat2str(returned, 3));
   [~, ~, spread] = vdhangular(RTs, ref);
   returned = mean(angse, 1)'./(spread*sqrt(M));
   expected = ones(12, 1);
   testCase.verifyEqual(returned, expected, 'AbsTol', 0.15, ...
      mat2str(returned', 3));
end
