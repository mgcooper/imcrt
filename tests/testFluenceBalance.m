function tests = testFluenceBalance
   % Tests for the fluence normalization (defects C, D) in src/mcrt.m. The
   % fluence is total absorption over ka. ka times its volume integral must
   % return the absorbed weight. The direct beam's absorption must sit in
   % the first radial bin only. A purely absorbing slab must give the
   % Beer-Lambert depth profile.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % The fixture adds src and restores path and rng.
   kernelfixture(testCase);
end

function testAbsorptionBalance(testCase)
   % ka*sum(phi_rz.*dV) and ka*sum(phi_z.*dz) both equal the absorbed
   % weight Adf + Adr to round-off, because phi is built from the same
   % tallies. Adding the per-depth direct term to every radial column
   % (defect C) inflates the phi_rz integral by the slab's area.
   cases = mcrtcases();
   for n = 1:numel(cases)
      c = cases(n);
      rng(c.seed, 'twister');
      RT = mcrt(c.ka, c.ks, c.g, c.Z, c.dz, c.N);
      dA = 2*pi*RT.grid.ri.*RT.grid.dr;
      dV = dA.*RT.grid.dz;
      absorbed = RT.Adf + sum(RT.Adr_z.*RT.grid.dz);
      returned = [c.ka*sum(RT.phi_rz(:).*dV(:)), ...
         c.ka*sum(RT.phi_z.*RT.grid.dz)];
      expected = [absorbed, absorbed];
      testCase.verifyEqual(returned, expected, 'RelTol', 1e-9, c.name);
   end
end

function testDirectFluenceConfinedToAxis(testCase)
   % With ks = 0 nothing scatters, so all absorption is the direct pencil
   % at r = 0 and phi_rz is nonzero in the first radial column only.
   rng(5, 'twister');
   RT = mcrt(10, 0, 0.5, 0.2, 0.01, 1e4);
   returned = [any(RT.phi_rz(:, 1) > 0), any(any(RT.phi_rz(:, 2:end) ~= 0))];
   expected = [true, false];
   testCase.verifyEqual(returned, expected);
end

function testPureAbsorptionIsBeerLambert(testCase)
   % With ks = 0 the absorbed fraction in depth bin i is the Beer-Lambert
   % increment exp(-ka*z_lo) - exp(-ka*z_hi), a binomial proportion. Five
   % binomial standard errors at N=1e6 bound each bin (about 1.5% relative
   % at the surface and 4% at tau = 2). ka*phi_z*dz recovers that fraction,
   % which checks the phi_z units end to end.
   rng(5, 'twister');
   N = 1e6;
   ka = 10;
   Z = 0.2;
   dz = 0.01;
   RT = mcrt(ka, 0, 0.5, Z, dz, N);
   nz = round(Z/dz);
   edges = (0:nz)'*dz;
   expected = exp(-ka*edges(1:nz)) - exp(-ka*edges(2:nz+1));
   returned = ka*RT.phi_z(1:nz).*RT.grid.dz(1:nz);
   tol = 5*sqrt(expected.*(1 - expected)/N);
   testCase.verifyEqual(returned, expected, 'AbsTol', tol);
end
