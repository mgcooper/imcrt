function tests = testBuildgrid
   % Unit tests for buildgrid: grid vector lengths, orientation, optimized
   % center coordinates, bin widths and measures for integral R/dr, A/da,
   % and Z/dz, the rejection of a fractional bin count (defect K), and the
   % comparison plot.
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
   [~, ri, ai, zi, dr, da, dz] = buildgrid(2, pi/2, 0.02, 0.001, pi/60, ...
      0.001);
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
   [~, ri, ai, zi] = buildgrid(R, A, Z, dr, da, dz);
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

function testWidthsAndMeasures(testCase)
   % The returned widths are the nominal bin widths, one per bin with the
   % overflow bins included, and the grid struct carries them with the
   % annulus areas and solid angles computed from the bin edges: the areas
   % sum to pi*R^2 inside R, and the solid angles to 2*pi over the
   % hemisphere.
   R = 2;
   Z = 0.02;
   [grid, ~, ~, ~, dr, da, dz] = buildgrid(R, pi/2, Z, 0.001, pi/60, 0.001);
   returned = {dr, da, dz, sum(grid.dA(1:end-1)), sum(grid.dsr), ...
      sum(dz(1:end-1)), sort(fieldnames(grid))};
   expected = {0.001*ones(1, 2001), pi/60*ones(30, 1), 0.001*ones(21, 1), ...
      pi*R^2, 2*pi, Z, sort({'ri'; 'ai'; 'zi'; 'dr'; 'da'; 'dz'; 'dA'; 'dsr'})};
   testCase.verifyEqual(returned, expected, 'RelTol', 1e-12);
end

function testFractionalCountsAreRejected(testCase)
   % With R/dr = Z/dz = 2.5 and A/da = 3.5 the colon operator would build
   % floor(n) centers while mcrt sizes its tallies with round(n), so
   % buildgrid refuses each fractional count by name.
   R = 0.025;
   A = pi/2;
   Z = 0.025;
   dr = 0.01;
   da = pi/7;
   dz = 0.01;
   testCase.verifyError(@() buildgrid(R, A, 0.02, dr, pi/6, 0.01), ...
      'buildgrid:nonintegral');
   testCase.verifyError(@() buildgrid(0.02, A, 0.02, dr, da, 0.01), ...
      'buildgrid:nonintegral');
   testCase.verifyError(@() buildgrid(0.02, A, Z, dr, pi/6, dz), ...
      'buildgrid:nonintegral');
   % An extent below one bin is a fractional count too, and a ratio that
   % rounds to zero within the tolerance fails the one-or-more bound.
   testCase.verifyError(@() buildgrid(0.02, A, 0.005, dr, pi/6, 0.01), ...
      'buildgrid:nonintegral');
   testCase.verifyError(@() buildgrid(0.02, A, 1e-12, dr, pi/6, 1), ...
      'buildgrid:nonintegral');
   % A width whose half underflows cannot start a grid.
   testCase.verifyError(@() buildgrid(0.02, A, eps(0), dr, pi/6, eps(0)), ...
      'buildgrid:width');
end

function testWholeCountsWithRoundoff(testCase)
   % Ratios that are whole on paper but not in floating point, such as
   % 0.02/0.001 and 0.1/0.005, and a ratio 1e-10 relative below a whole
   % number, pass the check and give round(n) bins, the sizes mcrt uses
   % for its tallies.
   [~, ri, ~, zi] = buildgrid(2, pi/2, 0.02, 0.001, pi/60, 0.001);
   [~, ri2, ~, zi2] = buildgrid(0.1, pi/2, 0.1, 0.005, pi/60, 0.005);
   [~, ~, ~, zi3] = buildgrid(0.1, pi/2, 20*(1 - 1e-10), 0.005, pi/60, 1);
   % one bin a rounding error short: the colon is empty and the center
   % is built directly
   [~, ~, ~, zi4] = buildgrid(0.1, pi/2, 1 - 1e-10, 0.005, pi/60, 1);
   returned = {numel(ri), numel(zi), numel(ri2), numel(zi2), numel(zi3), ...
      zi4};
   expected = {2001, 21, 21, 21, 21, [0.5; 1.5]};
   testCase.verifyEqual(returned, expected, 'AbsTol', 1e-15);
end

function testOneAngularBin(testCase)
   % A single angular bin (A = da) is allowed: its center is the shifted
   % half-width and its width is the input width.
   [~, ~, ai, ~, ~, da] = buildgrid(0.02, pi/2, 0.02, 0.001, pi/2, 0.001);
   ac = pi/4;
   returned = {numel(ai), ai, da};
   expected = {1, ac + cot(ac)*(1 - pi/4*cot(pi/4)), pi/2};
   testCase.verifyEqual(returned, expected, 'AbsTol', 1e-15);
end

function testPlotOption(testCase)
   % The seventh argument draws the comparison of bin centers and optimized
   % coordinates: one new figure with four axes, and only that figure is
   % closed afterward.
   before = findall(0, 'Type', 'figure');
   buildgrid(0.02, pi/2, 0.02, 0.001, pi/60, 0.001, true);
   after = findall(0, 'Type', 'figure');
   new = setdiff(after, before);
   testCase.addTeardown(@() close(new));
   returned = {numel(new), numel(findall(new, 'Type', 'axes'))};
   expected = {1, 4};
   testCase.verifyEqual(returned, expected);
end
