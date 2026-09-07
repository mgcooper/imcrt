function W = variancecheck(RTs, N)
   % Run-to-run variance of Rd and Tt against the binomial estimate
   % p(1-p)/N with p the mean over the M runs in RTs. Packet weights below
   % 1 lower the true variance, so the ratio sits at or below 1. A ratio
   % above 2 points to correlated runs (bad seeding) or a broken tally and
   % fails. W has the layout of vdhverdict.
   M = numel(RTs);
   x = zeros(M, 2);
   for n = 1:M
      x(n, :) = [RTs{n}.Rdf, RTs{n}.Tt];
   end
   p = mean(x, 1);
   ratio = (var(x, 0, 1)./(p.*(1 - p)/N))';
   W.metric = {'var(Rd)/binomial'; 'var(Tt)/binomial'};
   W.model = ratio;
   W.reference = [1; 1];
   W.se = [NaN; NaN];
   W.t = [NaN; NaN];
   W.verdict = repmat({'FAIL'}, 2, 1);
   W.verdict(ratio <= 2) = {'PASS'};
end
