function tests = testBuildgrid
   % Unit tests for buildgrid: grid vector lengths, orientation, optimized
   % center coordinates, and bin widths for integral R/dr, A/da, and Z/dz.
   % The non-integral test records current behavior only. Bead .15 (defect
   % K) decides whether buildgrid rejects non-integral inputs or matches the
   % kernel's round(), and updates that test.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % buildgrid calls derivative from src/derivative; the fixture adds both.
   kernelfixture(testCase);
end

function testIntegralLengthsAndOrientation(testCase)
   % The kernel default grid: 2000 radial bins, 30 angular, 20 vertical.
   % r and z carry one overflow bin and a does not. r is a row; a and z
   % are columns, which the kernel's broadcasting relies on.
   [ri, ai, zi, dr, da, dz] = buildgrid(2, pi/2, 0.02, 0.001, pi/60, 0.001);
   returned = {size(ri), size(ai), size(zi), size(dr), size(da), size(dz)};
   expected = {[1 2001], [30 1], [21 1], [1 2001], [30 1], [21 1]};
   testCase.verifyEqual(returned, expected);
end

function testOptimizedCenters(testCase)
   % r and a centers carry the bin-shift corrections cited in buildgrid
   % (Eqs. 8 and 14). z centers are plain bin centers with the overflow
   % center one bin above the slab.
   R = 2;
   A = pi/2;
   Z = 0.02;
   dr = 0.001;
   da = A/30;
   dz = 0.001;
   [ri, ai, zi] = buildgrid(R, A, Z, dr, da, dz);
   rc = [dr/2:dr:R-dr/2, R+dr/2];
   ac = (da/2:da:A-da/2)';
   returned = ri;
   expected = rc + (dr*dr/12)./rc;
   testCase.verifyEqual(returned, expected, 'AbsTol', 1e-15);
   returned = ai;
   expected = ac + cot(ac).*(1-da/2*cot(da/2));
   testCase.verifyEqual(returned, expected, 'AbsTol', 1e-15);
   returned = zi;
   expected = [dz/2:dz:Z-dz/2, Z+dz/2]';
   testCase.verifyEqual(returned, expected, 'AbsTol', 1e-15);
end

function testWidthsFollowCenters(testCase)
   % The returned widths are derivative() of each center vector, not the
   % input bin sizes, so the shifted r and a centers get shifted widths.
   % The z widths sum back to the slab thickness because z is unshifted.
   Z = 0.02;
   [ri, ai, zi, dr, da, dz] = buildgrid(2, pi/2, Z, 0.001, pi/60, 0.001);
   returned = {dr, da, dz};
   expected = {derivative(ri), derivative(ai), derivative(zi)};
   testCase.verifyEqual(returned, expected);
   returned = sum(dz(1:end-1));
   expected = Z;
   testCase.verifyEqual(returned, expected, 'RelTol', 1e-12);
end

function testNonIntegralLengthsCurrentBehavior(testCase)
   % Characterization of current behavior. With R/dr = Z/dz = 2.5 and A/da = 3.5
   % the colon operator yields floor(n) centers, while mcrt sizes its tallies
   % with round(n) (defect K). The two disagree at .5 and above.
   R = 0.025;
   A = pi/2;
   Z = 0.025;
   dr = 0.01;
   da = pi/7;
   dz = 0.01;
   [ri, ai, zi] = buildgrid(R, A, Z, dr, da, dz);
   returned = {numel(ri), numel(ai), numel(zi)};
   expected = {floor(R/dr) + 1, floor(A/da), floor(Z/dz) + 1};
   testCase.verifyEqual(returned, expected);
   returned = {round(R/dr) + 1, round(A/da), round(Z/dz) + 1};
   testCase.verifyNotEqual(returned, expected);
end
