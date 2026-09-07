function tests = testDirectBeam
   % Tests for the direct tally (defect A) and the index clamps (I, J) in
   % src/mcrt.m. Direct means unscattered only, so direct transmittance is
   % Beer-Lambert and direct reflectance is exactly zero. Scattered packets
   % that leave exactly along the axis or exactly at grazing must land in
   % a valid angular and radial bin instead of raising an index error.
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
   % and the vertical clamp only when a step lands on z = 0 exactly. Neither
   % is reachable through the public API, so evaluate the kernel's scoring
   % block, read from src/mcrt.m between its markers, on crafted values.
   % Each case first proves that the unclamped index is out of range.
   block = kernellines('% grid indices', '% absorption and scattering by ice');
   v = struct('dr', 0.001, 'dz', 0.001, 'nr', 20, 'nz', 20, 'na', 30, ...
      'Z', 0.02, 'a', 0.1, 'wt', 1, 'Tdr', 0, 'Rdr', 0, ...
      'Tdf_ra', zeros(30, 21), 'Rdf_ra', zeros(30, 21), ...
      'Adr_z', zeros(21, 1), 'Adf_rz', zeros(21, 21));
   % Grazing transmittance: da is a hair below (pi/2)/30 so acos(0)/da
   % rounds above 30 and ceil gives 31.
   v.da = (pi/2)/30*(1 - 4*eps);
   v.x = 0.005;
   v.y = 0;
   v.z = v.Z + v.dz;
   v.uz = 0;
   v.ns = 1;
   returned = ceil(acos(v.uz)/v.da);
   expected = v.na + 1;
   testCase.verifyEqual(returned, expected);
   w = runblock(block, v);
   returned = w.Tdf_ra(v.na, 5);
   expected = v.wt;
   testCase.verifyEqual(returned, expected);
   % The mirrored grazing reflection: z < 0 with -uz = 0 hits the same
   % rounding in the reflection branch and must land in Rdf_ra(na, :).
   v.z = -v.dz;
   w = runblock(block, v);
   returned = w.Rdf_ra(v.na, 5);
   expected = v.wt;
   testCase.verifyEqual(returned, expected);
   % Unscattered packet landing on z = 0 exactly with r = 0: both ceil
   % results are 0, and the direct absorption must go to depth bin 1.
   v.x = 0;
   v.y = 0;
   v.z = 0;
   v.uz = 0.5;
   v.ns = 0;
   returned = [ceil(sqrt(v.x^2 + v.y^2)/v.dr), ceil(v.z/v.dz)];
   expected = [0, 0];
   testCase.verifyEqual(returned, expected);
   w = runblock(block, v);
   returned = w.Adr_z(1);
   expected = v.a*v.wt;
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
