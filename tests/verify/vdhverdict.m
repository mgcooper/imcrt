function V = vdhverdict(RTs, ref, tmax, abstol, reltol)
   % Verdict for the van de Hulst Table 35 case from M independent runs
   % (cell array RTs). V is a struct of column arrays: metric and verdict
   % are cellstr; model, reference, se, and t are double. The rows are the
   % hemispherical Rd, Tt, Tdr, and Rdr, then twelve angular entries.
   % model is the mean over runs and se the run spread over sqrt(M); t is
   % the distance to the reference in standard errors. A statistical row
   % passes when |t| <= tmax. Rdr has no spread test and passes only when
   % every run returns exactly zero. The model has no specular term and
   % the source points straight down, so any other value is a defect.
   % With a fourth argument abstol, the hemispherical Rd and Tt rows also
   % pass when |model - reference| <= abstol, which allows for the table's
   % own precision. With a fifth argument reltol, an angular row also
   % passes when |model - reference| <= reltol*|reference|. That allows
   % for the 3-degree bin averaging and the interpolation between bin
   % centers: 1.6% at the mu = 1 peak, 0.8% at mu = 0.1, under 0.5%
   % elsewhere. V.byallowance marks the rows that passed only by an
   % allowance; a zero allowance never marks a row. Plain arrays keep the
   % script runnable on core Octave, which has no table type.
   if nargin < 4
      abstol = 0;
   end
   if nargin < 5
      reltol = 0;
   end
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
   V.byallowance = false(numel(V.t), 1);
   for n = 1:numel(V.t)
      if n <= 2
         allowance = abstol;
      elseif n >= 5
         allowance = reltol*abs(V.reference(n));
      else
         allowance = 0;
      end
      if allowance > 0 && strcmp(V.verdict{n}, 'FAIL') ...
            && abs(V.model(n) - V.reference(n)) <= allowance
         V.verdict{n} = 'PASS';
         V.byallowance(n) = true;
      end
   end
end
