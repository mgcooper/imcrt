function tests = testSetup
   % Tests for Setup.m path scoping (defect G): Setup must add exactly the repo
   % root, src, and src/derivative, and must drop stale example dirs quietly.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % Snapshot the path, clear the repo dirs so each addition is observable, bind
   % a handle to Setup while it is visible, then leave the repo so Setup cannot
   % rely on the working directory.
   testCase.TestData.path = path;
   root = fileparts(fileparts(mfilename('fullpath')));
   testCase.TestData.root = root;
   testCase.TestData.repodirs = {root, fullfile(root, 'src'), ...
      fullfile(root, 'src', 'derivative')};
   onpath = ismember(testCase.TestData.repodirs, strsplit(path, pathsep));
   if any(onpath)
      rmpath(testCase.TestData.repodirs{onpath});
   end
   testCase.TestData.startDir = cd(root);
   testCase.TestData.Setup = @Setup;
   cd(tempdir);
end

function teardownOnce(testCase)
   % Undo the process-global path and working-directory changes so they do not
   % leak into later suites or the caller's session.
   path(testCase.TestData.path);
   cd(testCase.TestData.startDir);
end

function testAddsOnlyRepoDirs(testCase)
   % Setup adds exactly root, src, and src/derivative, so no sandbox or archive
   % entry can shadow src/mcrt.m.
   before = strsplit(path, pathsep);
   testCase.TestData.Setup();
   returned = setdiff(strsplit(path, pathsep), before);
   expected = sort(testCase.TestData.repodirs);
   testCase.verifyEqual(returned, expected);
   returned = which('mcrt');
   expected = fullfile(testCase.TestData.root, 'src', 'mcrt.m');
   testCase.verifyEqual(returned, expected);
end

function testExamplesRemovedQuietly(testCase)
   % A stale session may hold the frozen example kernel on its path. Setup must
   % remove it and raise no warning while doing so.
   exdir = fullfile(testCase.TestData.root, 'examples', ...
      'cooper_etal_2021', 'c_model');
   addpath(exdir);
   lastwarn('');
   testCase.TestData.Setup();
   returned = ismember(exdir, strsplit(path, pathsep));
   expected = false;
   testCase.verifyEqual(returned, expected);
   returned = lastwarn;
   expected = '';
   testCase.verifyEqual(returned, expected);
end
