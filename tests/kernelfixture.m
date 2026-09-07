function kernelfixture(testCase)
   % Put src and src/derivative on the path for the life of one test file.
   % The fixture restores the path afterward, so runtests works without Setup
   % and the Setup.m path tests stay independent. The teardown also restores
   % the caller's global random stream, which the seeded tests replace.
   root = fileparts(fileparts(mfilename('fullpath')));
   testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
      {fullfile(root, 'src'), fullfile(root, 'src', 'derivative')}));
   state = rng;
   testCase.addTeardown(@() rng(state));
end
