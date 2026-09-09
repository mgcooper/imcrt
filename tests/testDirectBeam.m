function tests = testDirectBeam
   % Tests for the direct tally (defect A) in src/mcrt.m and the index
   % clamps (I, J) in src/binindex.m. Direct means unscattered only, so
   % direct transmittance is Beer-Lambert and direct reflectance is
   % exactly zero. Scattered packets that leave exactly along the axis or
   % exactly at grazing must land in a valid angular and radial bin
   % instead of raising an index error.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % The fixture adds src and restores path and rng.
   kernelfixture(testCase);
end

function testDirectTransmittanceIsBeerLambert(testCase)
   % An unscattered packet exits iff its first free path exceeds Z, so the
   % count of direct exits is Binomial(N, exp(-tau)). Five binomial standard
   % errors at N=1e5 bound the estimate (sigma about 1.1e-3 for tau=2).
   rng(7, 'twister');
   N = 1e5;
   RT = mcrt(10, 90, 0.75, 0.02, 0.001, N);
   returned = RT.Tdr;
   expected = exp(-2);
   tol = 5*sqrt(expected*(1 - expected)/N);
   testCase.verifyEqual(returned, expected, 'AbsTol', tol);
end

function testDirectReflectanceIsZero(testCase)
   % The source points straight down and the model has no specular term,
   % so no unscattered packet can cross z < 0 on any case.
   cases = mcrtcases();
   for n = 1:numel(cases)
      c = cases(n);
      rng(c.seed, 'twister');
      RT = mcrt(c.ka, c.ks, c.g, c.Z, c.dz, c.N);
      returned = RT.Rdr;
      expected = 0;
      testCase.verifyEqual(returned, expected, c.name);
   end
end

function testOnAxisExitsLandInFirstBin(testCase)
   % g = 1 scatters every packet exactly forward, so scattered packets exit
   % with uz = 1 and r = 0: acos(1) = 0 and ceil(0/dr) = 0 would index bin
   % zero. The clamps put them in angular bin 1 and radial bin 1.
   rng(7, 'twister');
   RT = mcrt(10, 90, 1, 0.02, 0.001, 1e3);
   returned = [RT.Tdf_a(1) > 0, all(RT.Tdf_a(2:end) == 0), ...
      RT.Tdf_r(1) > 0, all(RT.Tdf_r(2:end) == 0)];
   expected = true(1, 4);
   testCase.verifyEqual(returned, expected);
end

function testOnAxisReflectionLandsInFirstBin(testCase)
   % g = -1 reverses every packet, so reflected packets leave with uz = -1
   % exactly and r = 0. The same clamps apply to the reflection branch.
   rng(7, 'twister');
   RT = mcrt(10, 90, -1, 0.02, 0.001, 1e3);
   returned = [RT.Rdf_a(1) > 0, all(RT.Rdf_a(2:end) == 0), RT.Rdr];
   expected = [true, true, 0];
   testCase.verifyEqual(returned, expected);
end

function testClampsOnCraftedIndices(testCase)
   % The upper angular clamp fires only when acos(uz)/da rounds above na,
   % and the lower clamps only when a step lands on r = 0 or z = 0
   % exactly. Neither is reachable through the public API, so call the
   % kernel's index function on crafted values. Each case first proves
   % that the unclamped index is out of range.
   na = 30;
   % Grazing exit: da is a hair below (pi/2)/30 so acos(0)/da rounds above
   % 30 and ceil gives 31; the clamp keeps the last bin.
   da = (pi/2)/30*(1 - 4*eps);
   returned = [ceil(acos(0)/da), binindex(acos(0), da, na)];
   expected = [na + 1, na];
   testCase.verifyEqual(returned, expected);
   % On the axis and on the surface: ceil(0) is 0, the clamp gives bin 1.
   returned = [ceil(0/0.001), binindex(0, 0.001, 21), binindex(0, 0.001, 21)];
   expected = [0, 1, 1];
   testCase.verifyEqual(returned, expected);
   % Radial overflow: past R the index is the overflow bin nr+1.
   returned = binindex(0.05, 0.001, 21);
   expected = 21;
   testCase.verifyEqual(returned, expected);
   % Interior: an ordinary coordinate is unchanged by the clamp.
   returned = binindex(0.0045, 0.001, 21);
   expected = 5;
   testCase.verifyEqual(returned, expected);
end

function testGrazingExitsLandInLastBin(testCase)
   % Isotropic scattering (g = 0) produces near-grazing scattered exits.
   % Direct means ns == 0 only, so they are diffuse and belong in the last
   % angular bin.
   rng(7, 'twister');
   RT = mcrt(10, 90, 0, 0.02, 0.001, 1e4);
   returned = [RT.Tdf_a(end) > 0, RT.Rdf_a(end) > 0];
   expected = [true, true];
   testCase.verifyEqual(returned, expected);
end
