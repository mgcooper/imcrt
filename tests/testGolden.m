function tests = testGolden
   % Golden-digest regression: the seeded cases must reproduce the tracked
   % baseline line for line. The baseline matches the kernel at HEAD, so a
   % behavior-neutral commit must leave it bit-identical. A physics fix
   % re-baselines in the same commit with mcrtgolden('write') and states
   % the delta in its commit message.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % The kernel must be on the path; the fixture restores the path afterward.
   kernelfixture(testCase);
end

function testDigestMatchesBaseline(testCase)
   % Compare line cells, not one string, so a mismatch names the first
   % differing line in the diagnostic.
   returned = mcrtgolden();
   expected = readlines(mcrtgolden('path'));
   testCase.verifyEqual(returned, expected);
end

function testWriteRoundTrips(testCase)
   % The re-baseline command must write exactly the digest lines. Write into
   % a temporary folder that the fixture deletes, never the tracked file.
   folder = testCase.applyFixture( ...
      matlab.unittest.fixtures.TemporaryFolderFixture);
   file = fullfile(folder.Folder, 'digest.txt');
   expected = mcrtgolden('write', file);
   returned = readlines(file);
   testCase.verifyEqual(returned, expected);
end

function lines = readlines(file)
   % Split on either line ending so a CRLF checkout compares the same lines.
   lines = regexp(strtrim(fileread(file)), '\r?\n', 'split');
end
