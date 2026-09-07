function tests = testBinMeasures
   % Tests for the bin measures (defect S) in src/mcrt.m. The solid angles
   % and annulus areas that normalize the resolved tallies must be the
   % exact measures of the bins. The diffuse angular tallies must then
   % match the van de Hulst table at every nonzero mu.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % The fixture adds src, tests/oracle, and tests/verify.
   kernelfixture(testCase);
end

function testMeasuresMatchGeometry(testCase)
   % The angular bins tile the hemisphere and the radial bins inside R tile
   % the disc, so their measures sum to 2*pi and pi*R^2. The first bins are
   % a cap and a disc, where shifted-center formulas are 18% off.
   c = mcrtcases();
   c = c(1);
   rng(c.seed, 'twister');
   RT = mcrt(c.ka, c.ks, c.g, c.Z, c.dz, 1e3);
   na = numel(RT.grid.ai);
   nr = numel(RT.grid.ri) - 1;
   da = (pi/2)/na;
   dr = RT.grid.dr(1);
   returned = [sum(RT.grid.dsr), sum(RT.grid.dA(1:nr)), ...
      RT.grid.dsr(1), RT.grid.dA(1)];
   expected = [2*pi, pi*(nr*dr)^2, 2*pi*(1 - cos(da)), pi*dr^2];
   testCase.verifyEqual(returned, expected, 'RelTol', 1e-12);
end

function testResolvedTalliesSumToHemispherical(testCase)
   % Multiplying each resolved tally by its measure and summing must give
   % back the hemispherical fraction, which ties the measures to the
   % normalization the kernel used.
   c = mcrtcases();
   c = c(1);
   rng(c.seed, 'twister');
   RT = mcrt(c.ka, c.ks, c.g, c.Z, c.dz, c.N);
   returned = [sum(RT.Rdf_a.*RT.grid.dsr), sum(RT.Tdf_a.*RT.grid.dsr), ...
      sum(RT.Rdf_r.*RT.grid.dA), sum(RT.Tdf_r.*RT.grid.dA)];
   expected = [RT.Rdf, RT.Tdf, RT.Rdf, RT.Tdf];
   testCase.verifyEqual(returned, expected, 'RelTol', 1e-12);
end

function testVdhAngularTableWithinRunSpread(testCase)
   % The van de Hulst Table 35 sentinel. Twenty independent runs of N=5e4
   % (seeds 1 to 20) give a run-spread standard error with 19 degrees of
   % freedom; |t| <= 3.5 has a 0.24% chance per row under agreement. With
   % shifted-center measures the mu = 1 transmittance entry sits about 13
   % standard errors low.
   ref = vdhtable35();
   M = 20;
   RTs = cell(1, M);
   for m = 1:M
      rng(m, 'twister');
      RTs{m} = mcrt(ref.ka, ref.ks, ref.g, ref.Z, ref.dz, 5e4);
   end
   z = vdhangular(RTs, ref);
   returned = abs(z) <= 3.5;
   expected = true(12, 1);
   testCase.verifyEqual(returned, expected, sprintf('z = %s', mat2str(z', 3)));
end
