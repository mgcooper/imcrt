function tests = testChgdirOracle
   % Tests for the direction-cosine update (defect B) in src/chgdir.m, which
   % the src/mcrt.m photon loop calls. The update must keep the unit norm
   % and deflect by exactly the sampled cosine. It must spread azimuths
   % uniformly per the RotationMatrix oracle.
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
