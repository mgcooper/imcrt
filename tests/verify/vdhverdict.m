function V = vdhverdict(RTs, ref, tmax)
   % Verdict for the van de Hulst Table 35 case from M independent runs
   % (cell array RTs). V is a struct of column arrays: metric and verdict
   % are cellstr; model, reference, se, and t are double. The rows are the
   % hemispherical Rd, Tt, Tdr, and Rdr, then twelve angular entries.
   % model is the mean over runs and se the run spread over sqrt(M); t is
   % the distance to the reference in standard errors. A statistical row
   % passes when |t| <= tmax. Rdr has no spread test and passes only when
   % every run returns exactly zero. The model has no specular term and
   % the source points straight down, so any other value is a defect.
   % Plain arrays keep the script runnable on core Octave, which has no
   % table type.
   M = numel(RTs);
   hemi = zeros(M, 4);
   for n = 1:M
      hemi(n, :) = [RTs{n}.Rdf, RTs{n}.Tt, RTs{n}.Tdr, RTs{n}.Rdr];
   end
   [ta, ma, sa] = vdhangular(RTs, ref);
   mu = ref.mu(2:end);
   label = @(fmt) arrayfun(@(x) sprintf(fmt, x), mu, 'UniformOutput', false)';
   V.metric = [{'Rd'; 'Tt'; 'Tdr'; 'Rdr'}; label('R(mu=%.1f)'); ...
      label('T(mu=%.1f)')];
   V.model = [mean(hemi, 1)'; ma];
   V.reference = [ref.Rd; ref.Tt; ref.Tdr; 0; ref.R_sr(2:end)'; ...
      ref.T_sr(2:end)'];
   V.se = [std(hemi, 0, 1)'/sqrt(M); sa];
   V.t = [(V.model(1:3) - V.reference(1:3))./V.se(1:3); NaN; ta];
   V.verdict = repmat({'FAIL'}, numel(V.t), 1);
   V.verdict(abs(V.t) <= tmax) = {'PASS'};
   if all(hemi(:, 4) == 0)
      V.verdict{4} = 'PASS';
   end
end
