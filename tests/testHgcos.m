function tests = testHgcos
   % Tests for the Henyey-Greenstein polar-cosine sampler src/hgcos.m,
   % which the src/mcrt.m photon loop calls. Every draw is a cosine, the
   % sample mean is the asymmetry g, and g = 0 takes the isotropic branch.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % The fixture adds src and restores path and rng.
   kernelfixture(testCase);
end

function testFirstMomentIsAsymmetry(testCase)
   % The mean cosine of Henyey-Greenstein is g. Five standard errors of
   % the sample mean bound it, and every draw lies in [-1, 1].
   rng(12, 'twister');
   ndraw = 2e4;
   for g = [0.75, -0.5]
      us = zeros(1, ndraw);
      for n = 1:ndraw
         us(n) = hgcos(g);
      end
      testCase.verifyTrue(all(us >= -1 & us <= 1));
      tol = 5*std(us)/sqrt(ndraw);
      returned = mean(us);
      expected = g;
      testCase.verifyEqual(returned, expected, 'AbsTol', tol);
   end
end

function testIsotropicBranch(testCase)
   % g = 0 draws 1 - 2*rand: uniform on [-1, 1] with mean 0, and the same
   % value the uniform draw gives from the same stream state.
   rng(12, 'twister');
   expected = 1 - 2*rand;
   rng(12, 'twister');
   returned = hgcos(0);
   testCase.verifyEqual(returned, expected);
   us = zeros(1, 2e4);
   for n = 1:numel(us)
      us(n) = hgcos(0);
   end
   testCase.verifyEqual(mean(us), 0, 'AbsTol', 5*std(us)/sqrt(numel(us)));
end
