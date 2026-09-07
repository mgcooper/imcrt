function tests = testDeterminism
   % Seeding convention: rng(seed,'twister') immediately before the kernel call
   % must reproduce every output exactly, and a different seed must not.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % The kernel must be on the path; the fixture restores the path afterward.
   kernelfixture(testCase);
end

function testSameSeedSameOutput(testCase)
   % Two runs from the same seed are bit-identical, which the golden-digest
   % regression relies on.
   c = mcrtcases();
   c = c(1);
   rng(c.seed, 'twister');
   returned = mcrt(c.ka, c.ks, c.g, c.Z, c.dz, 1e3);
   rng(c.seed, 'twister');
   expected = mcrt(c.ka, c.ks, c.g, c.Z, c.dz, 1e3);
   testCase.verifyEqual(returned, expected);
end

function testDifferentSeedDifferentOutput(testCase)
   % The kernel draws from the global stream, so a new seed changes the tallies.
   c = mcrtcases();
   c = c(1);
   rng(c.seed, 'twister');
   returned = mcrt(c.ka, c.ks, c.g, c.Z, c.dz, 1e3);
   rng(c.seed + 1, 'twister');
   expected = mcrt(c.ka, c.ks, c.g, c.Z, c.dz, 1e3);
   testCase.verifyNotEqual(returned.Rdf_ra, expected.Rdf_ra);
end
