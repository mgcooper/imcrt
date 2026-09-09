function tests = testSmoke
   % Smoke tests for mcrt on each shared case at its own N. The output struct
   % must hold every field, shapes tied to the grid, and finite non-negative
   % values.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % The kernel must be on the path; the fixture restores the path afterward.
   kernelfixture(testCase);
end

function testFieldsAndShapes(testCase)
   % Every case returns the same field set. Shapes follow the grid vectors so
   % the test survives a normalization fix that keeps the tally layout.
   cases = mcrtcases();
   for k = 1:numel(cases)
      c = cases(k);
      rng(c.seed, 'twister');
      RT = mcrt(c.ka, c.ks, c.g, c.Z, c.dz, c.N);
      returned = sort(fieldnames(RT));
      expected = sort({'Rdf_ra'; 'Rdf_r'; 'Rdf_a'; 'Rdf'; 'Rdr'; 'Tdf_ra'; ...
         'Tdf_r'; 'Tdf_a'; 'Tdf'; 'Tdr'; 'Tt'; 'Adf_z'; 'Adf'; 'Adf_rz'; ...
         'Adr_z'; 'phi_rz'; 'phi_z'; 'grid'; 'se'; 'Rt'; 'N'});
      testCase.verifyEqual(returned, expected, c.name);
      na = numel(RT.grid.ai);
      nr = numel(RT.grid.ri);
      nz = numel(RT.grid.zi);
      returned = {size(RT.Rdf_ra), size(RT.Rdf_r), size(RT.Rdf_a), ...
         size(RT.Tdf_ra), size(RT.Tdf_r), size(RT.Tdf_a), ...
         size(RT.Adf_rz), size(RT.Adr_z), size(RT.Adf_z), ...
         size(RT.phi_rz), size(RT.phi_z), ...
         size([RT.Rdf RT.Rdr RT.Tdf RT.Tdr RT.Tt RT.Adf])};
      expected = {[na nr], [1 nr], [na 1], [na nr], [1 nr], [na 1], ...
         [nz nr], [nz 1], [nz 1], [nz nr], [nz 1], [1 6]};
      testCase.verifyEqual(returned, expected, c.name);
   end
end

function testValuesFiniteNonnegative(testCase)
   % Resolved tallies are weights divided by bin areas and widths. A zero-width
   % bin shows up there as Inf or NaN. Scalars are weights divided by N. A sign
   % error shows up in any field as a negative value.
   cases = mcrtcases();
   for k = 1:numel(cases)
      c = cases(k);
      rng(c.seed, 'twister');
      RT = mcrt(c.ka, c.ks, c.g, c.Z, c.dz, c.N);
      names = fieldnames(RT);
      returned = true(numel(names), 1);
      for i = 1:numel(names)
         v = RT.(names{i});
         if isstruct(v)
            continue
         end
         returned(i) = all(isfinite(v(:)) & v(:) >= 0);
      end
      expected = true(numel(names), 1);
      testCase.verifyEqual(returned, expected, c.name);
   end
end
