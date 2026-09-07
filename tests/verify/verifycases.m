function c = verifycases(casename)
   % Verification cases for mcrt_verify.m and for larger runs, keyed by
   % name. Each case carries the kernel inputs, the interactive run
   % count M, the packets N, the seeds, and tmax. M and N are sized for
   % seconds on a laptop. tmax is the run-spread t cutoff for a PASS. The
   % full field holds the same four values for vdhovernight: 1e8
   % packets for reflect and 1e6 for fluence. reflect uses the DesignSpec
   % cutoff of 3 standard errors over 20 runs (0.7% false-fail per row at
   % 19 degrees of freedom; Rd, Tt, and the angular rows also carry the
   % driver's allowances). fluence uses 3.5 over 16 runs because its ten
   % direct-beam rows have no allowance.
   % 'reflect' takes its inputs and references from vdhtable35. 'fluence'
   % is the internal-fluence case of Wang et al. 1995, Fig. 4 (albedo
   % 0.999, g = 0.9, 4 cm slab). Wang compared it with van de Hulst Vol. 1,
   % p. 262, Table 12 through the similarity relation. This repository
   % holds no numeric values from either source, so the fluence checks are
   % self-consistency only.
   switch casename
      case 'reflect'
         ref = vdhtable35();
         c = struct('name', 'reflect', 'ka', ref.ka, 'ks', ref.ks, ...
            'g', ref.g, 'Z', ref.Z, 'dz', ref.dz, 'M', 8, 'N', 1.25e5, ...
            'seeds', 1:8, 'tmax', 4, 'full', struct('M', 20, 'N', 5e6, ...
            'seeds', 1:20, 'tmax', 3));
      case 'fluence'
         c = struct('name', 'fluence', 'ka', 0.1, 'ks', 100, 'g', 0.9, ...
            'Z', 4, 'dz', 0.005, 'M', 4, 'N', 500, 'seeds', 21:24, ...
            'tmax', 4, 'full', struct('M', 16, 'N', 6.25e4, ...
            'seeds', 21:36, 'tmax', 3.5));
      otherwise
         error('verifycases:unknownCase', ...
            'Unknown case "%s"; use reflect or fluence.', casename);
   end
end
