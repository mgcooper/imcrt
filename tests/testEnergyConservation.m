function tests = testEnergyConservation
   % Every packet's weight must end in reflectance, transmittance, or
   % absorption. The 1e-2 tolerance allows the Russian-roulette weight leak
   % (defect N).
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % The kernel must be on the path; the fixture restores the path afterward.
   kernelfixture(testCase);
end

function testWeightIsConserved(testCase)
   % RT.Adf is already a fraction of N. RT.Adr_z is per unit depth. Weight
   % Adr_z by the grid bin widths before adding it to the balance.
   cases = mcrtcases();
   for k = 1:numel(cases)
      c = cases(k);
      rng(c.seed, 'twister');
      RT = mcrt(c.ka, c.ks, c.g, c.Z, c.dz, c.N);
      returned = RT.Rdf + RT.Rdr + RT.Tdf + RT.Tdr + RT.Adf ...
         + sum(RT.Adr_z .* RT.grid.dz);
      expected = 1;
      testCase.verifyEqual(returned, expected, 'AbsTol', 1e-2, c.name);
   end
end
