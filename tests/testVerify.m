function tests = testVerify
   % Tests for the verification layer: the case table, the van de Hulst
   % and fluence verdicts with passing and failing inputs, and the printer.
   % mcrt_verify.m must also run standalone for every case.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % The fixture adds src and tests/verify. Figures the script opens are
   % closed afterward so headless runs leave nothing behind.
   kernelfixture(testCase);
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
         'seeds'; 'tmax'});
      testCase.verifyEqual(returned, expected, names{n});
      returned = numel(c.seeds);
      expected = c.M;
      testCase.verifyEqual(returned, expected, names{n});
   end
   testCase.verifyError(@() verifycases('nope'), 'verifycases:unknownCase');
end

function testVdhVerdictRules(testCase)
   % Four short runs give the 16-row verdict. Rdr passes by the exact-zero
   % rule and fails once any run reports a nonzero Rdr. Scaling the
   % reference reflectance by 1.5 flips Rd to FAIL and nothing else.
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

function testScriptRunsEachCase(testCase)
   % mcrt_verify.m honors a preset casename, runs at its interactive size,
   % and prints a passing VERDICT line for both cases. The runner makes
   % tests/ the current folder, so the script runs by path. The script
   % leaves c in this workspace, which proves it read casename.
   script = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
      'mcrt_verify.m');
   command = ['run(''' script ''')'];
   casename = 'reflect';
   out = evalc(command);
   returned = {contains(out, 'VERDICT reflect: PASS'), c.name};
   expected = {true, casename};
   testCase.verifyEqual(returned, expected, out(max(1, end-300):end));
   casename = 'fluence';
   out = evalc(command);
   returned = {contains(out, 'VERDICT fluence: PASS'), c.name};
   expected = {true, casename};
   testCase.verifyEqual(returned, expected, out(max(1, end-300):end));
end

function RTs = runcase(c, M, N, seeds)
   % M seeded runs of one case at N packets.
   RTs = cell(1, M);
   for n = 1:M
      rng(seeds(n), 'twister');
      RTs{n} = mcrt(c.ka, c.ks, c.g, c.Z, c.dz, N);
   end
end
