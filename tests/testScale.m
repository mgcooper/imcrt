function tests = testScale
   % Tests for the normalization functions src/scaleR.m, src/scaleT.m, and
   % src/scaleA.m on synthetic raw tallies with a small grid: the sums are
   % the packet weights over N, the resolved tallies are the weights over
   % measure and N, the totals add direct and diffuse, and the fluence
   % closes the absorption balance.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % The fixture adds src and restores path and rng.
   kernelfixture(testCase);
end

function grid = smallgrid()
   % A 3 x 2 (angle x radius, overflow included) and 2 x 2 (depth x
   % radius) grid with measures from its edges, as mcrt builds them.
   dr = 0.5;
   da = pi/6;
   dz = 0.25;
   redge = (0:2)*dr;
   aedge = (0:3)'*da;
   grid.ri = [0.25 0.75];
   grid.ai = (aedge(1:end-1) + aedge(2:end))/2;
   grid.zi = [0.125; 0.375];
   grid.dr = dr*ones(1, 2);
   grid.da = da*ones(3, 1);
   grid.dz = dz*ones(2, 1);
   grid.dA = pi*(redge(2:end).^2 - redge(1:end-1).^2);
   grid.dsr = 2*pi*(cos(aedge(1:end-1)) - cos(aedge(2:end)));
end

function testReflectanceAndTransmittance(testCase)
   % Both functions apply the same rule: sums over N, resolved tallies
   % over measure, projection factor, and N. Rt and Tt add the direct
   % part. With unit weights the squares equal the sums, and every
   % standard error is tallyse of the raw sum scaled like its output.
   grid = smallgrid();
   N = 10;
   raw = [1 2; 0 3; 4 0];
   direct = 2;
   [ra, r, a, df, dr, t, se] = scaleR(raw, direct, raw, direct, grid, N);
   [tra, tr, ta, tdf, tdr, tt, tse] = scaleT(raw, direct, raw, direct, ...
      grid, N);
   cosa = cos(grid.ai);
   returned = {ra, r, a, df, dr, t, tra, tr, ta, tdf, tdr, tt};
   expected = {raw./(grid.dA.*grid.dsr.*cosa.*N), sum(raw, 1)./(grid.dA*N), ...
      sum(raw, 2)./(grid.dsr*N), sum(raw(:))/N, direct/N, ...
      (sum(raw(:)) + direct)/N};
   expected = [expected, expected];
   testCase.verifyEqual(returned, expected, 'AbsTol', 1e-15);
   returned = {se.ra, se.r, se.a, se.df, se.dr, se.t, tse};
   expected = {tallyse(raw, raw, N)./(grid.dA.*grid.dsr.*cosa), ...
      tallyse(sum(raw, 1), sum(raw, 1), N)./grid.dA, ...
      tallyse(sum(raw, 2), sum(raw, 2), N)./grid.dsr, ...
      tallyse(sum(raw(:)), sum(raw(:)), N), tallyse(direct, direct, N), ...
      tallyse(sum(raw(:)) + direct, sum(raw(:)) + direct, N), se};
   testCase.verifyEqual(returned, expected, 'AbsTol', 1e-15);
end

function testAbsorptionAndFluence(testCase)
   % Adf_rz is per volume, Adf_z and Adr_z per depth, Adf a fraction; the
   % fluence is total absorption over ka, the direct pencil in radial bin
   % 1 only, and ka times the depth integral of phi_z returns Adf + Adr.
   grid = smallgrid();
   N = 10;
   ka = 2;
   raw = [1 2; 3 0];
   rawdr = [4; 1];
   [rz, z, df, drz, prz, pz] = scaleA(raw, rawdr, ka, grid, N);
   dV = grid.dA.*grid.dz;
   erz = raw./dV/N;
   pencil = rawdr./dV(:, 1)/N;
   eprz = erz/ka;
   eprz(:, 1) = eprz(:, 1) + pencil/ka;
   returned = {rz, z, df, drz, prz, pz, ka*sum(pz.*grid.dz)};
   expected = {erz, sum(raw, 2)./grid.dz/N, sum(raw(:))/N, ...
      rawdr./grid.dz/N, eprz, (sum(raw, 2) + rawdr)./grid.dz/N/ka, ...
      (sum(raw(:)) + sum(rawdr))/N};
   testCase.verifyEqual(returned, expected, 'AbsTol', 1e-15);
end
