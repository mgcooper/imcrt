function tests = testVerify
   % Tests for the verification layer: the case table, the van de Hulst
   % and fluence verdicts with passing and failing inputs, and the printer.
   % mcrt_verify must also run for every case and default to reflect.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % The fixture adds src and tests/verify, and the repo root holds
   % mcrt_verify. Figures the function opens are closed afterward so
   % headless runs leave nothing behind.
   kernelfixture(testCase);
   testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
      fileparts(fileparts(mfilename('fullpath')))));
   testCase.addTeardown(@() close('all'));
end

function testCasesCarryRunParameters(testCase)
   % Every case names its kernel inputs, run count, packets, seeds, and
   % cutoff. An unknown name is an error rather than a silent default.
   names = {'reflect', 'fluence'};
   for n = 1:numel(names)
      c = verifycases(names{n});
      returned = sort(fieldnames(c));
      expected = sort({'name'; 'ka'; 'ks'; 'g'; 'Z'; 'dz'; 'M'; 'N'; ...
         'seeds'; 'tmax'; 'full'});
      testCase.verifyEqual(returned, expected, names{n});
      returned = [numel(c.seeds), numel(c.full.seeds)];
      expected = [c.M, c.full.M];
      testCase.verifyEqual(returned, expected, names{n});
   end
   testCase.verifyError(@() verifycases('nope'), 'verifycases:unknownCase');
end

function testVdhVerdictRules(testCase)
   % Four short runs give the 16-row verdict. Rdr passes by the exact-zero
   % rule and fails once any run reports a nonzero Rdr. Scaling the
   % reference reflectance by 1.5 flips Rd to FAIL and nothing else. The
   % absolute allowance applies to Rd and Tt only.
   ref = vdhtable35();
   RTs = runcase(ref, 4, 2e4, 1:4);
   V = vdhverdict(RTs, ref, 4);
   returned = {numel(V.verdict), V.metric{1}, V.verdict{4}, isnan(V.t(4))};
   expected = {16, 'Rd', 'PASS', true};
   testCase.verifyEqual(returned, expected);
   bad = ref;
   bad.Rd = 1.5*ref.Rd;
   V = vdhverdict(RTs, bad, 4);
   returned = {V.verdict{1}, nnz(strcmp(V.verdict, 'FAIL'))};
   expected = {'FAIL', 1};
   testCase.verifyEqual(returned, expected);
   leaky = RTs;
   leaky{2}.Rdr = 1e-3;
   V = vdhverdict(leaky, ref, 4);
   returned = V.verdict{4};
   expected = 'FAIL';
   testCase.verifyEqual(returned, expected);
   % An absolute allowance rescues a hemispherical row that misses by less
   % than the allowance and marks it in byallowance; without it the row fails.
   % Tt may be rescued too at this tiny tmax; Tdr and the angular rows
   % never are.
   V = vdhverdict(RTs, ref, 4);
   near = ref;
   near.Rd = V.model(1) + 3e-4;
   V = vdhverdict(RTs, near, 0.1, 5e-4);
   returned = {V.verdict{1}, V.byallowance(1), any(V.byallowance(3:end))};
   expected = {'PASS', true, false};
   testCase.verifyEqual(returned, expected);
   V = vdhverdict(RTs, near, 0.1);
   returned = V.verdict{1};
   expected = 'FAIL';
   testCase.verifyEqual(returned, expected);
   % The relative allowance rescues an angular row (row 11, T at mu = 0.1)
   % that misses by 1% when reltol is 2%; without reltol it fails. Rows 1
   % to 4 never use the relative allowance, so none of them is marked.
   V = vdhverdict(RTs, ref, 4);
   skew = ref;
   skew.T_sr(2) = V.model(11)/1.01;
   V = vdhverdict(RTs, skew, 0.1, 0, 0.02);
   returned = {V.verdict{11}, V.byallowance(11), any(V.byallowance(1:4))};
   expected = {'PASS', true, false};
   testCase.verifyEqual(returned, expected);
   V = vdhverdict(RTs, skew, 0.1);
   returned = V.verdict{11};
   expected = 'FAIL';
   testCase.verifyEqual(returned, expected);
end

function testVdhVerdictWithFewRuns(testCase)
   % Below four runs the verdict takes its standard errors from RT.se
   % instead of the run spread. Two runs at N = 2e4 give finite, positive
   % errors on every statistical row and a passing table at tmax 4; the
   % Rdr row keeps its exact rule.
   ref = vdhtable35();
   RTs = runcase(ref, 2, 2e4, 1:2);
   V = vdhverdict(RTs, ref, 4, 5e-4, 0.02);
   rows = [1:3, 5:16];
   returned = {all(isfinite(V.se(rows)) & V.se(rows) > 0), isnan(V.t(4)), ...
      nnz(strcmp(V.verdict, 'PASS'))};
   expected = {true, true, 16};
   testCase.verifyEqual(returned, expected, mat2str(V.t', 3));
   % A run with no diffuse exit on one side (no scattering) gives zero
   % angular errors on that side, not NaN; a one-packet run, whose errors
   % are undefined, keeps NaN.
   rng(3, 'twister');
   RT = mcrt(10, 0, 0.5, 0.02, 0.001, 200);
   [~, ~, se] = vdhangular({RT}, ref, true);
   RT1 = mcrt(10, 0, 0.5, 0.02, 0.001, 1);
   [~, ~, se1] = vdhangular({RT1}, ref, true);
   returned = {se(1:6)', all(isnan(se1(1:6)))};
   expected = {zeros(1, 6), true};
   testCase.verifyEqual(returned, expected);
end

function testFluenceVerdictRules(testCase)
   % The configured fluence case must pass all 12 rows. Each rule also
   % gets a failing input. An inflated absorbed weight fails the energy
   % row. A doubled first direct-beam bin fails its t row. A surface
   % fluence pushed below 1 fails the last row.
   c = verifycases('fluence');
   RTs = runcase(c, c.M, c.N, c.seeds);
   V = fluenceverdict(RTs, c, c.tmax);
   returned = V.verdict';
   expected = repmat({'PASS'}, 1, 12);
   testCase.verifyEqual(returned, expected, mat2str(V.t', 3));
   heavy = RTs;
   heavy{1}.Adf = heavy{1}.Adf + 0.01;
   V = fluenceverdict(heavy, c, c.tmax);
   returned = V.verdict{1};
   expected = 'FAIL';
   testCase.verifyEqual(returned, expected);
   skewed = RTs;
   for n = 1:numel(skewed)
      skewed{n}.Adr_z(1) = 2*skewed{n}.Adr_z(1);
   end
   V = fluenceverdict(skewed, c, c.tmax);
   returned = V.verdict{2};
   expected = 'FAIL';
   testCase.verifyEqual(returned, expected);
   dim = RTs;
   for n = 1:numel(dim)
      dim{n}.phi_z(1) = 0.5;
   end
   V = fluenceverdict(dim, c, c.tmax);
   returned = V.verdict{end};
   expected = 'FAIL';
   testCase.verifyEqual(returned, expected);
end

function testVarianceCheckThreshold(testCase)
   % Synthetic runs with a known spread and known per-run errors: a ratio
   % of spread to error inside [0.5, 2] passes, and ratios above 2 and
   % below 0.5 fail, for Rd and Tt independently.
   p = 0.1;
   sigma = 1e-3;
   RTs = cell(1, 4);
   for n = 1:4
      RTs{n} = struct('Rdf', p + sigma*(n - 2.5), ...
         'Tt', p + 3*sigma*(n - 2.5), 'se', struct('Rdf', sigma, 'Tt', sigma));
   end
   W = variancecheck(RTs);
   returned = {W.verdict{1}, W.verdict{2}, W.model(1), W.model(2) > 2};
   expected = {'PASS', 'FAIL', std([-1.5 -0.5 0.5 1.5]), true};
   testCase.verifyEqual(returned, expected, 'AbsTol', 1e-12);
   for n = 1:4
      RTs{n}.se.Rdf = 10*sigma;
   end
   W = variancecheck(RTs);
   returned = {W.verdict{1}, W.model(1) < 0.5};
   expected = {'FAIL', true};
   testCase.verifyEqual(returned, expected);
end

function testAngularUsesPchip(testCase)
   % Synthetic angular tallies that curve like the transmittance near
   % grazing: the helper must return the pchip interpolant, not the linear
   % one, at the table's mu values.
   ref = vdhtable35();
   RT = struct();
   RT.grid.ai = ((0.5:29.5)'*pi/60);
   RT.Rdf_a = exp(-3*RT.grid.ai);
   RT.Tdf_a = cos(RT.grid.ai).^4;
   [~, model] = vdhangular({RT, RT}, ref);
   theta = acos(ref.mu(2:end));
   returned = model;
   expected = [interp1(RT.grid.ai, RT.Rdf_a, theta, 'pchip', RT.Rdf_a(1)), ...
      interp1(RT.grid.ai, RT.Tdf_a, theta, 'pchip', RT.Tdf_a(1))]';
   testCase.verifyEqual(returned, expected, 'AbsTol', 1e-15);
   linear = [interp1(RT.grid.ai, RT.Rdf_a, theta, 'linear', RT.Rdf_a(1)), ...
      interp1(RT.grid.ai, RT.Tdf_a, theta, 'linear', RT.Tdf_a(1))]';
   testCase.verifyNotEqual(returned, linear);
end

function testPrintVerdictCountsPasses(testCase)
   % The printer returns the number of PASS rows and prints one line per
   % row plus the header.
   V = struct('metric', {{'a'; 'b'; 'c'}}, 'model', [1; 2; 3], ...
      'reference', [1; 2; 4], 'se', [0.1; 0.1; 0.1], 't', [0; 0; -10], ...
      'verdict', {{'PASS'; 'PASS'; 'FAIL'}});
   out = evalc('npass = printverdict(V);');
   returned = {npass, nnz(out == newline), numel(V.verdict)};
   expected = {2, 4, 3};
   testCase.verifyEqual(returned, expected);
end

function testFunctionRunsEachCase(testCase)
   % mcrt_verify runs each case at its interactive size, prints a passing
   % VERDICT line, and returns a verdict table with every row PASS and
   % the M runs. With no argument it runs reflect.
   [out, V, RTs] = evalc('mcrt_verify(''reflect'')');
   returned = {contains(out, 'VERDICT reflect: PASS'), ...
      all(strcmp(V.verdict, 'PASS')), numel(RTs)};
   expected = {true, true, verifycases('reflect').M};
   testCase.verifyEqual(returned, expected, out(max(1, end-300):end));
   [out, V, RTs] = evalc('mcrt_verify(''fluence'')');
   returned = {contains(out, 'VERDICT fluence: PASS'), ...
      all(strcmp(V.verdict, 'PASS')), numel(RTs)};
   expected = {true, true, verifycases('fluence').M};
   testCase.verifyEqual(returned, expected, out(max(1, end-300):end));
   out = evalc('mcrt_verify()');
   testCase.verifyTrue(contains(out, 'VERDICT reflect: PASS'), ...
      out(max(1, end-300):end));
end

function RTs = runcase(c, M, N, seeds)
   % M seeded runs of one case at N packets.
   RTs = cell(1, M);
   for n = 1:M
      rng(seeds(n), 'twister');
      RTs{n} = mcrt(c.ka, c.ks, c.g, c.Z, c.dz, N);
   end
end
