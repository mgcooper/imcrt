function tests = testPerf
   % Performance regression tests, run on demand with runtests('tests/perf')
   % and kept out of the deterministic fast suite because timing on a
   % loaded laptop is not repeatable. The kernel's time on the perf cases,
   % divided by the reference loop timed in the same run, must stay within
   % the tolerance of the checked-in baseline, measured on one host and
   % one MATLAB release; on any other host or release the comparison is
   % filtered, not failed. The overhead microbenchmark must return finite
   % positive times.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % tests/ holds kernelfixture and mcrtcases; the fixture then adds src
   % and tests/perf and restores path and rng.
   testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
      fileparts(fileparts(mfilename('fullpath')))));
   kernelfixture(testCase);
end

function testKernelWithinBaseline(testCase)
   % A regression is a normalized time (case over its paired reference)
   % above the baseline's by more than the tolerance. The reference
   % cancels most of the machine's load and thermal drift, perfbench
   % keeps the least disturbed of three paired timings, and the tolerance
   % is still a factor of two: the test catches a tripled loop, not a
   % ten-percent change.
   tol = 1.0;
   b = perfbench('read');
   testCase.assumeEqual(b.computer, computer, ...
      'baseline is from another architecture');
   testCase.assumeEqual(b.host, perfbench('host'), ...
      'baseline is from another host');
   testCase.assumeEqual(b.version, version, ...
      'baseline is from another MATLAB release');
   T = perfbench();
   for k = 1:numel(T)
      bc = b.cases(strcmp({b.cases.name}, T(k).name));
      base = bc.seconds/bc.reference;
      returned = T(k).seconds/T(k).reference;
      testCase.verifyLessThanOrEqual(returned, base*(1 + tol), ...
         sprintf('%s: %.3f of the reference against baseline %.3f', ...
         T(k).name, returned, base));
   end
end

function testBaselineWriteRead(testCase)
   % A written baseline reads back with this architecture, host, release,
   % and the measured times; the shipped baseline has the perf cases and
   % sits at the path the harness reports.
   fixture = testCase.applyFixture( ...
      matlab.unittest.fixtures.TemporaryFolderFixture);
   file = fullfile(fixture.Folder, 'baseline.txt');
   T = perfbench('write', file);
   b = perfbench('read', file);
   shipped = perfbench('read');
   root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
   returned = {b.computer, b.host, b.version, {b.cases.name}, ...
      [b.cases.seconds], [b.cases.reference], {shipped.cases.name}, ...
      perfbench('path')};
   expected = {computer, perfbench('host'), version, {T.name}, ...
      [T.seconds], [T.reference], {perfcases().name}, ...
      fullfile(root, 'tests', 'perf', 'baseline.txt')};
   testCase.verifyEqual(returned, expected, 'AbsTol', 1e-6);
end

function testOverheadNumbersAreSane(testCase)
   % Every timed loop, at the default length and at a given length, costs
   % a finite positive time per iteration, and the chgdir loop costs at
   % most a few times the inline loop.
   for o = [perfoverhead(), perfoverhead(2e4)]
      t = cell2mat(struct2cell(rmfield(o, 'ratio')));
      returned = [all(isfinite(t) & t > 0), o.ratio > 0.5, o.ratio < 10];
      expected = true(1, 3);
      testCase.verifyEqual(returned, expected, sprintf( ...
         'call %.3g s, inline %.3g s, chgdir %.3g s', o.call, o.inline, ...
         o.chgdir));
   end
end
