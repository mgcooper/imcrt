function tests = testOvernight
   % Tests for vdhovernight at a tiny scale. It must create the output
   % folder, write one checkpoint per run and a dated report, resume from
   % matching checkpoints, recompute stale or incomplete ones, keep every
   % report, and reach OVERALL FAIL when the t cutoff is zero.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % The fixture adds src and tests/verify.
   kernelfixture(testCase);
end

function testReportNameSuffixOnCollision(testCase)
   % An existing report with the same stamp gets a numeric suffix, and the
   % suffix counts up while names are taken.
   folder = testCase.applyFixture( ...
      matlab.unittest.fixtures.TemporaryFolderFixture);
   stamp = '20260101_000000';
   first = reportname(folder.Folder, stamp);
   fclose(fopen(first, 'w'));
   second = reportname(folder.Folder, stamp);
   fclose(fopen(second, 'w'));
   third = reportname(folder.Folder, stamp);
   returned = {first, second, third};
   expected = {fullfile(folder.Folder, 'vdh_report_20260101_000000.txt'), ...
      fullfile(folder.Folder, 'vdh_report_20260101_000000_2.txt'), ...
      fullfile(folder.Folder, 'vdh_report_20260101_000000_3.txt')};
   testCase.verifyEqual(returned, expected);
end

function testReportFailurePaths(testCase)
   % A read-only output folder makes the report unwritable, which is an
   % error with the driver's own identifier. A checkpoint whose RT is not
   % a struct passes the freshness check and then breaks the verdict; the
   % cleanup must close the report so no file identifier leaks.
   folder = testCase.applyFixture( ...
      matlab.unittest.fixtures.TemporaryFolderFixture);
   locked = fullfile(folder.Folder, 'locked');
   mkdir(locked);
   fileattrib(locked, '-w');
   testCase.addTeardown(@() fileattrib(locked, '+w'));
   testCase.verifyError(@() vdhovernight(locked, 1e-3), ...
      'vdhovernight:report');
   outdir = fullfile(folder.Folder, 'broken');
   evalc('vdhovernight(outdir, 1e-3)');
   files = dir(fullfile(outdir, 'checkpoint_reflect_*.mat'));
   ck = fullfile(outdir, files(1).name);
   S = load(ck);
   RT = 0;
   meta = S.meta;
   save(ck, 'RT', 'meta');
   nopen = numel(openedFiles);
   testCase.verifyError(@() evalc('vdhovernight(outdir, 1e-3)'), ...
      ?MException);
   returned = numel(openedFiles);
   expected = nopen;
   testCase.verifyEqual(returned, expected);
end

function testCheckpointsResumeAndVerdicts(testCase)
   % Scale 1e-3 keeps the two cases to about a second per call. The output
   % folder does not exist yet, so the first call creates it.
   folder = testCase.applyFixture( ...
      matlab.unittest.fixtures.TemporaryFolderFixture);
   outdir = fullfile(folder.Folder, 'nightly');
   [out, reportfile] = evalc('vdhovernight(outdir, 1e-3)');
   files = dir(fullfile(outdir, 'checkpoint_*.mat'));
   tmpfiles = dir(fullfile(outdir, '*_tmp.mat'));
   text = fileread(reportfile);
   returned = {numel(files), numel(tmpfiles), exist(reportfile, 'file'), ...
      contains(out, 'OVERALL'), contains(text, 'var(Rd)/binomial'), ...
      contains(text, 'var(Tt)/binomial'), contains(text, 'sha256'), ...
      contains(text, 'phi_z profile')};
   expected = {36, 0, 2, true, true, true, true, true};
   testCase.verifyEqual(returned, expected, out(max(1, end-400):end));
   % A second call loads every checkpoint and still reports run times.
   [out, reportfile] = evalc('vdhovernight(outdir, 1e-3)');
   returned = {numel(strfind(out, 'loaded')), ...
      isempty(regexp(fileread(reportfile), 'over all runs at scale', 'once'))};
   expected = {36, false};
   testCase.verifyEqual(returned, expected);
   % A checkpoint whose kernel hash differs is stale and is recomputed.
   ck = fullfile(outdir, files(1).name);
   S = load(ck);
   meta = S.meta;
   meta.kernel = 'stale';
   RT = S.RT;
   save(ck, 'RT', 'meta');
   out = evalc('vdhovernight(outdir, 1e-3)');
   returned = {numel(strfind(out, 'stale checkpoint')), ...
      numel(strfind(out, 'loaded'))};
   expected = {1, 35};
   testCase.verifyEqual(returned, expected);
   % Every freshness discriminator makes a checkpoint stale on its own, not
   % an error: a foreign version, changed inputs, a missing RT, a missing
   % meta, a meta without its fields, a nonnumeric or negative run time,
   % and a cell-valued kernel hash. A string output path works too.
   good = load(ck);
   broken = {@(T) setfield(T, 'meta', setfield(T.meta, 'version', 'x')), ...
      @(T) setfield(T, 'meta', ...
         setfield(T.meta, 'inputs', T.meta.inputs + 1)), ...
      @(T) rmfield(T, 'RT'), @(T) rmfield(T, 'meta'), ...
      @(T) setfield(T, 'meta', struct('kernel', T.meta.kernel)), ...
      @(T) setfield(T, 'meta', setfield(T.meta, 'elapsed', 'x')), ...
      @(T) setfield(T, 'meta', setfield(T.meta, 'elapsed', -1)), ...
      @(T) setfield(T, 'meta', setfield(T.meta, 'kernel', {'a', 'b'}))};
   for n = 1:numel(broken)
      T = broken{n}(good);
      save(ck, '-struct', 'T');
      out = evalc('vdhovernight(string(outdir), 1e-3)');
      returned = numel(strfind(out, 'stale checkpoint'));
      expected = 1;
      testCase.verifyEqual(returned, expected, sprintf('corruption %d', n));
   end
   % An unreadable checkpoint file is stale, not an error.
   fid = fopen(ck, 'w');
   fprintf(fid, 'not a MAT-file');
   fclose(fid);
   out = evalc('vdhovernight(outdir, 1e-3)');
   returned = numel(strfind(out, 'stale checkpoint'));
   expected = 1;
   testCase.verifyEqual(returned, expected);
   % tmax = 0 fails every statistical row without an allowance, so the
   % report ends in OVERALL FAIL. Two reports in one second get distinct
   % names.
   [out, first] = evalc('vdhovernight(outdir, 1e-3, 0)');
   [~, second] = evalc('vdhovernight(outdir, 1e-3, 0)');
   returned = {contains(out, 'OVERALL FAIL'), strcmp(first, second), ...
      exist(first, 'file') == 2 && exist(second, 'file') == 2};
   expected = {true, false, true};
   testCase.verifyEqual(returned, expected);
end
