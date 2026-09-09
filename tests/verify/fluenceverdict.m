function V = fluenceverdict(RTs, c, tmax)
   % Self-consistency verdict for the fluence case from M independent runs,
   % because this repository holds no numeric fluence reference. V is a
   % struct of column arrays in the layout of vdhverdict. Row 1 is the
   % energy balance: model is the largest |balance - 1| over runs and it
   % passes at or below 1e-4. Rows 2 to 11 are the direct beam's absorbed
   % fraction in the first ten depth bins. Their reference is
   % a*(exp(-mt*z_lo) - exp(-mt*z_hi)) with mt = ka + ks and a = ka/mt.
   % model is the mean over runs and se the run spread over sqrt(M); the
   % row passes when |t| <= tmax. Row 12 is the surface fluence phi_z(1).
   % Multiple scattering at albedo 0.999 lifts it above 1 (Wang et al.
   % 1995, Fig. 4). The row passes when the mean over runs exceeds 1.
   M = numel(RTs);
   nbin = 10;
   balance = zeros(M, 1);
   direct = zeros(M, nbin);
   surface = zeros(M, 1);
   for n = 1:M
      RT = RTs{n};
      balance(n) = RT.Rdf + RT.Rdr + RT.Tdf + RT.Tdr + RT.Adf ...
         + sum(RT.Adr_z.*RT.grid.dz);
      direct(n, :) = (RT.Adr_z(1:nbin).*RT.grid.dz(1:nbin))';
      surface(n) = RT.phi_z(1);
   end
   mt = c.ka + c.ks;
   edges = (0:nbin)'*c.dz;
   pdirect = c.ka/mt*(exp(-mt*edges(1:nbin)) - exp(-mt*edges(2:nbin+1)));
   bins = arrayfun(@(k) sprintf('Adr(bin %d)', k), 1:nbin, ...
      'UniformOutput', false)';
   V.metric = [{'energy'}; bins; {'phi_z(1) > 1'}];
   V.model = [max(abs(balance - 1)); mean(direct, 1)'; mean(surface)];
   V.reference = [0; pdirect; 1];
   V.se = [NaN; std(direct, 0, 1)'/sqrt(M); std(surface)/sqrt(M)];
   V.t = [NaN; (V.model(2:nbin+1) - V.reference(2:nbin+1)) ...
      ./V.se(2:nbin+1); NaN];
   V.verdict = repmat({'FAIL'}, numel(V.t), 1);
   V.verdict(abs(V.t) <= tmax) = {'PASS'};
   if V.model(1) <= 1e-4
      V.verdict{1} = 'PASS';
   end
   if V.model(end) > 1
      V.verdict{end} = 'PASS';
   end
end
