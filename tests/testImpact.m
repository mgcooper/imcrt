function tests = testImpact
   % Tests for the impact-report tooling: the attribution rule, the
   % quantity summary, the change text, the source check, and the report
   % generator end to end at a tiny size.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % The fixture adds src and tests/verify.
   kernelfixture(testCase);
end

function testAttributionNamesChangedVersions(testCase)
   % A row changed at the second and fourth versions names both fixes; a
   % row that only moves by rounding names none; a row that becomes
   % defined and then undefined names both transitions.
   values = [1, 1, 1.5, 1.5, 3; 2, 2 + 1e-9, 2, 2, 2; NaN, NaN, 1, 1, NaN];
   returned = impactattribution(values, {'B', 'A', 'C', 'S'});
   expected = {'A, S'; 'none'; 'A, S'};
   testCase.verifyEqual(returned, expected);
end

function testQuantitiesFromOneRun(testCase)
   % Thirteen named quantities with notes, also from the no-argument
   % call; the half-weight radii sit inside the detection radius; the
   % hemispherical entry matches the struct.
   c = mcrtcases();
   c = c(1);
   rng(c.seed, 'twister');
   RT = mcrt(c.ka, c.ks, c.g, c.Z, c.dz, 2e3);
   [names, values, notes] = impactquantities(RT);
   [names0, values0] = impactquantities();
   returned = {numel(names), numel(values), numel(notes), names0, ...
      values0, values(1), values(10) > 0 && values(10) < 2, ...
      values(11) > 0 && values(11) < 2};
   expected = {13, 13, 13, names, [], RT.Rdf, true, true};
   testCase.verifyEqual(returned, expected);
end

function testDefaultsAndTextHelpers(testCase)
   % impactargs fills missing or empty arguments with the documented
   % defaults, passes given values through, and rejects a single run
   % and fractional counts.
   % srctext and agreetext are the two report sentences whose branches a
   % report run cannot reach on demand; paired NaN entries agree, and a
   % one-sided NaN is excluded from the largest difference.
   [f1, n1, m1] = impactargs();
   [f2, n2, m2] = impactargs('r.md', [], 3);
   testCase.verifyError(@() impactargs('r.md', 5, 1), 'impactargs:runs');
   testCase.verifyError(@() impactargs('r.md', 5, 2.5), 'impactargs:runs');
   testCase.verifyError(@() impactargs('r.md', 2.5, 2), ...
      'impactargs:packets');
   testCase.verifyError(@() impactargs('r.md', 2 + 1i, 2), ...
      'impactargs:packets');
   testCase.verifyError(@() impactargs('r.md', 5, 2 + 1i), 'impactargs:runs');
   returned = {f1, n1, m1, f2, n2, m2, ...
      srctext(true, true), srctext(true, false), ...
      agreetext([0 NaN; 1 2], [0 NaN; 1 2]), ...
      agreetext([0 NaN; 1 2], [0 3; 1.5 2])};
   expected = {'docs/impact-report.md', 1e5, 8, 'r.md', 1e5, 3, ...
      ['src/ is unchanged since tag fix-S-bin-measures, so the ' ...
      'post-fix column is the shipped kernel.'], ...
      ['src/ differs from tag fix-S-bin-measures or has uncommitted ' ...
      'edits, so the retrospective runs the worktree kernel.'], ...
      'are identical run for run', ...
      'differ by at most 0.5 where both are defined'};
   testCase.verifyEqual(returned, expected);
end

function testPairedRatioText(testCase)
   % The three outcomes of both ratio texts: a ratio, 'none' for
   % identical samples, and 'n/a' for a constant nonzero difference with
   % no spread. The unpaired ratio scales by the combined spread of two
   % samples of equal size.
   returned = {pairedratio([1 2 3 4]), pairedratio([0 0 0]), ...
      pairedratio([2 2 2]), unpairedratio([1 2 3 4], [1 2 3 5]), ...
      unpairedratio([2 2], [2 2]), unpairedratio([2 2], [3 3])};
   expected = {'3.87', 'none', 'n/a', '-0.23', 'none', 'n/a'};
   testCase.verifyEqual(returned, expected);
end

function testRelativeChangeText(testCase)
   % Signed percent text, and n/a for a zero or undefined baseline.
   returned = {relchange(2, 3), relchange(0, 1), relchange(NaN, 1), ...
      relchange(1, NaN)};
   expected = {'+50.00%', 'n/a', 'n/a', 'n/a'};
   testCase.verifyEqual(returned, expected);
end

function testSourceMatchesScratchRepository(testCase)
   % A scratch repository holds the three checked files at a tag. The
   % check reports the tag match and a clean tree, still reports clean
   % when the worktree copy only changes LF to CRLF, reports dirty after
   % a content edit, and reports a tree mismatch once that edit is
   % committed. A scratch repository keeps the suite independent of the
   % live worktree, which may hold legitimate uncommitted edits.
   % Assign the fixture first: R2020b cannot dot-index a call result.
   fixture = testCase.applyFixture( ...
      matlab.unittest.fixtures.TemporaryFolderFixture);
   folder = fixture.Folder;
   git = @(args) system(sprintf( ...
      ['git -C "%s" -c user.name=t -c user.email=t@t ' ...
      '-c commit.gpgsign=false --no-pager %s'], folder, args));
   mkdir(fullfile(folder, 'src', 'derivative'));
   files = {fullfile('src', 'mcrt.m'), fullfile('src', 'buildgrid.m'), ...
      fullfile('src', 'derivative', 'derivative.m')};
   for n = 1:numel(files)
      writetext(fullfile(folder, files{n}), ...
         sprintf('function x = f\n   x = 1;\nend\n'));
   end
   status = [git('init -q'), git('add -A'), git('commit -q -m base'), ...
      git('tag t1')];
   testCase.assertEqual(status, zeros(1, 4), 'git setup failed');
   [same1, clean1] = srcmatches(folder, 't1');
   % LF to CRLF only: the clean check normalizes line endings.
   writetext(fullfile(folder, files{1}), ...
      sprintf('function x = f\r\n   x = 1;\r\nend\r\n'));
   [same2, clean2] = srcmatches(folder, 't1');
   % A content edit is dirty until committed; then the tree differs.
   writetext(fullfile(folder, files{1}), ...
      sprintf('function x = f\n   x = 2;\nend\n'));
   [same3, clean3] = srcmatches(folder, 't1');
   status = [git('add -A'), git('commit -q -m edit')];
   testCase.assertEqual(status, zeros(1, 2), 'git commit failed');
   [same4, clean4] = srcmatches(folder, 't1');
   returned = {same1, clean1, same2, clean2, same3, clean3, same4, clean4};
   expected = {true, true, true, true, true, false, false, true};
   testCase.verifyEqual(returned, expected);
end

function testUndefinedHalfWeightRadius(testCase)
   % With no scattering nothing is diffuse inside R, so both r50 entries
   % are NaN rather than the first bin.
   rng(1, 'twister');
   RT = mcrt(10, 0, 0.5, 0.02, 0.001, 200);
   [~, values] = impactquantities(RT);
   returned = isnan(values(10:11))';
   expected = [true, true];
   testCase.verifyEqual(returned, expected);
end

function testReportEndToEnd(testCase)
   % A tiny report (N = 200, M = 2) must hold both sections, one table per
   % case including paper_radial, the retrospective table, and the bounding
   % paragraph. A second call leaves M at its default of 8 runs.
   folder = testCase.applyFixture( ...
      matlab.unittest.fixtures.TemporaryFolderFixture);
   file = fullfile(folder.Folder, 'impact.md');
   evalc('impactreport(file, 200, 2)');
   text = fileread(file);
   returned = {contains(text, '## 1. Impact of each fix'), ...
      contains(text, '## 2. Paper-kernel retrospective'), ...
      numel(strfind(text, '### ')), ...
      numel(strfind(text, '| Rdf |')), ...
      contains(text, 'paper_radial'), ...
      contains(text, 'How the published numbers relate')};
   expected = {true, true, 5, 5, true, true};
   testCase.verifyEqual(returned, expected, text(1:min(400, end)));
   evalc('impactreport(file, 200)');
   returned = contains(fileread(file), 'M = 8 runs');
   expected = true;
   testCase.verifyEqual(returned, expected);
end

function writetext(path, text)
   % Write text verbatim with fopen and fprintf, which R2020b has;
   % writelines arrived in R2022a.
   fid = fopen(path, 'w');
   fprintf(fid, '%s', text);
   fclose(fid);
end
