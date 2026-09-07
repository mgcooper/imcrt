function tests = testChgdirOracle
   % Tests for the direction-cosine update (defect B) in src/chgdir.m and in
   % the inline copy inside the src/mcrt.m photon loop. The update must keep
   % the unit norm and deflect by exactly the sampled cosine. It must spread
   % azimuths uniformly per the RotationMatrix oracle and agree between the
   % two copies. The inline copy is read from the source between its marker
   % comments and evaluated, so the hot loop keeps no function call.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % The fixture adds src and tests/oracle and restores path and rng.
   kernelfixture(testCase);
end

function testUnitNormAndDeflection(testCase)
   % A rotation keeps the norm, and the cosine between the old and the new
   % direction is the sampled polar cosine us. Both fail by up to 40% when
   % the uy expression reads the already-updated ux. The norm tolerance is
   % 1e-10, not 1e-12. The formula divides by sin(theta), so a direction
   % within about 1e-2 of the polar axis loses eps/sin(theta)^2 in the norm.
   rng(11, 'twister');
   ndraw = 1e5;
   u = randn(3, ndraw);
   u = u ./ sqrt(sum(u.^2, 1));
   us = 1 - 2*rand(1, ndraw);
   ps = 2*pi*rand(1, ndraw);
   v = zeros(3, ndraw);
   for n = 1:ndraw
      [v(1, n), v(2, n), v(3, n)] = chgdir(u(1, n), u(2, n), u(3, n), ...
         us(n), ps(n));
   end
   returned = sqrt(sum(v.^2, 1));
   expected = ones(1, ndraw);
   testCase.verifyEqual(returned, expected, 'AbsTol', 1e-10);
   returned = sum(u .* v, 1);
   expected = us;
   testCase.verifyEqual(returned, expected, 'AbsTol', 1e-12);
end

function testOnAxisBranch(testCase)
   % From exactly +z or -z the update takes the sin(theta) = 0 branch. The
   % new direction must still sit on the cone about the old one.
   for uz = [1, -1]
      [vx, vy, vz] = chgdir(0, 0, uz, 0.3, 1.1);
      returned = [vx^2 + vy^2 + vz^2, uz*vz];
      expected = [1, 0.3];
      testCase.verifyEqual(returned, expected, 'AbsTol', 1e-12);
   end
end

function testAzimuthUniformVsOracle(testCase)
   % Sweep the azimuth in K equal steps about a fixed off-axis direction.
   % In the oracle frame that carries +z onto that direction, the new
   % directions must sit on the cone at height us. Their azimuths must be
   % equally spaced by 2*pi/K, up to the sign convention of phi.
   u = [0.6; 0; 0.8];
   us = 0.3;
   K = 360;
   phis = 2*pi*(0:K-1)/K;
   v = zeros(3, K);
   for n = 1:K
      [v(1, n), v(2, n), v(3, n)] = chgdir(u(1), u(2), u(3), us, phis(n));
   end
   local = RotationMatrix([0 0 1], u)' * v;
   returned = local(3, :);
   expected = us*ones(1, K);
   testCase.verifyEqual(returned, expected, 'AbsTol', 1e-12);
   returned = abs(diff(unwrap(atan2(local(2, :), local(1, :)))));
   expected = 2*pi/K*ones(1, K-1);
   testCase.verifyEqual(returned, expected, 'AbsTol', 1e-9);
end

function testInlineBlockMatchesChgdir(testCase)
   % The kernel keeps an inline copy of chgdir for speed. Evaluate that
   % copy, read from src/mcrt.m between its marker comments, together with
   % the kernel's Henyey-Greenstein sampler. Check it against chgdir on the
   % same draws. Check that one scatter leaves weight w and count 1. Also
   % check the sampler's first moment, which is g.
   precompute = kernellines('% henyey-greenstein terms', '');
   scatter = kernellines('% absorption and scattering by ice', ...
      '% russian roulette');
   g = 0.75;
   w = 0.9;
   eval(precompute);
   rng(12, 'twister');
   ndraw = 2e4;
   u = randn(3, ndraw);
   u = u ./ sqrt(sum(u.^2, 1));
   u(:, 1) = [0; 0; 1];
   u(:, 2) = [0; 0; -1];
   returned = zeros(6, ndraw);
   expected = zeros(6, ndraw);
   for n = 1:ndraw
      ux = u(1, n);
      uy = u(2, n);
      uz = u(3, n);
      wt = 1;
      ns = 0;
      eval(scatter);
      [cx, cy, cz] = chgdir(u(1, n), u(2, n), u(3, n), us, ps);
      returned(:, n) = [ux; uy; uz; ux*u(1, n) + uy*u(2, n) + uz*u(3, n); ...
         wt; ns];
      expected(:, n) = [cx; cy; cz; us; w; 1];
   end
   testCase.verifyEqual(returned, expected, 'AbsTol', 1e-12);
   % Five standard errors of the sample mean bound the first moment.
   tol = 5*std(expected(4, :))/sqrt(ndraw);
   returned = mean(expected(4, :));
   expected = g;
   testCase.verifyEqual(returned, expected, 'AbsTol', tol);
end
