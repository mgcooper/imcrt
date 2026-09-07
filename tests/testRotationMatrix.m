function tests = testRotationMatrix
   % Tests for the tests/oracle/RotationMatrix.m oracle. It must be the
   % minimal proper rotation that maps a onto b for generic, parallel,
   % antiparallel, and near-antiparallel pairs. It must also resolve from
   % the suite's path.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % The fixture puts tests/oracle on the path and restores it afterward.
   kernelfixture(testCase);
end

function testResolvesFromOracleFolder(testCase)
   % The oracle lives outside src so it never shadows or joins production.
   returned = which('RotationMatrix');
   expected = fullfile(fileparts(mfilename('fullpath')), 'oracle', ...
      'RotationMatrix.m');
   testCase.verifyEqual(returned, expected);
end

function testAlignsAToB(testCase)
   % Random pairs: M*a must land on b to round-off, and M must be a proper
   % rotation (orthonormal with determinant +1).
   rng(7, 'twister');
   for n = 1:100
      verifyMapsAToB(testCase, randn(3, 1), randn(3, 1));
   end
end

function testParallelIsIdentity(testCase)
   % cross(a, a) is exactly zero, so the Rodrigues form must return I with
   % no division by zero; the kernel's initial direction is exactly +z.
   returned = RotationMatrix([0 0 1], [0 0 1]);
   expected = eye(3);
   testCase.verifyEqual(returned, expected);
end

function testAntiparallelIsHalfTurn(testCase)
   % The half-turn branch must map a to -a as a proper rotation, on axis
   % and off axis. [1 3 10] is a pair whose normalized dot product need not
   % round to exactly -1.
   pairs = {[0 0 1], [0 0 -1]; [1 2 3], [-1 -2 -3]; [1 0 0], [-2 0 0]; ...
      [1 3 10], [-1 -3 -10]};
   for n = 1:size(pairs, 1)
      verifyMapsAToB(testCase, pairs{n, 1}(:), pairs{n, 2}(:));
   end
end

function testNearAntiparallelMapsAToB(testCase)
   % Within about sqrt(eps) of antiparallel the dot product rounds to -1
   % while the cross product does not vanish. The oracle must still land on
   % b, not on -a.
   pairs = {[1 0 0], [-cos(1e-8) sin(1e-8) 0]; [1 0 0], [-1 1e-8 0]; ...
      [0 0 1], [1e-9 -1e-9 -1]};
   for n = 1:size(pairs, 1)
      verifyMapsAToB(testCase, pairs{n, 1}(:), pairs{n, 2}(:));
   end
end

function testObtusePairKeepsAxisFixed(testCase)
   % An obtuse pair that is not antiparallel must get the minimal rotation
   % about cross(a, b), not a rotation with an added twist. The axis itself
   % must therefore stay fixed.
   pairs = {[1 0 0], [-1 1 0]/sqrt(2); [1 0 0], [-0.5 sqrt(3)/2 0]; ...
      [1 2 3], [-3 -2 1]};
   for n = 1:size(pairs, 1)
      a = pairs{n, 1}(:);
      b = pairs{n, 2}(:);
      M = RotationMatrix(a, b);
      returned = M * cross(a, b);
      expected = cross(a, b);
      testCase.verifyEqual(returned, expected, 'AbsTol', 1e-12);
   end
end

function testRejectsBadInput(testCase)
   % Both inputs must be nonzero 3-vectors; the asserts stop anything else.
   testCase.verifyError(@() RotationMatrix([1 2], [1 2 3]), ...
      'MATLAB:assertion:failed');
   testCase.verifyError(@() RotationMatrix([1 2 3], eye(3)), ...
      'MATLAB:assertion:failed');
   testCase.verifyError(@() RotationMatrix([0 0 0], [1 2 3]), ...
      'MATLAB:assertion:failed');
   testCase.verifyError(@() RotationMatrix([1 2 3], [0 0 0]), ...
      'MATLAB:assertion:failed');
end

function verifyMapsAToB(testCase, a, b)
   % Shared check: M*a/|a| equals b/|b| to round-off and M is a proper
   % rotation. The diagnostic names the pair so a failing case is visible.
   M = RotationMatrix(a, b);
   label = sprintf('a=%s b=%s', mat2str(a', 3), mat2str(b', 3));
   returned = M * a / norm(a);
   expected = b / norm(b);
   testCase.verifyEqual(returned, expected, 'AbsTol', 1e-12, label);
   returned = {M' * M, det(M)};
   expected = {eye(3), 1};
   testCase.verifyEqual(returned, expected, 'AbsTol', 1e-12, label);
   % The axis is conditioned as eps/norm(cross(a, b)), so check that it
   % stays fixed only when it is well determined.
   v = cross(a / norm(a), b / norm(b));
   if norm(v) > 1e-6
      returned = M * v;
      expected = v;
      testCase.verifyEqual(returned, expected, 'AbsTol', 1e-12, label);
   end
end
