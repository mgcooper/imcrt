function W = variancecheck(RTs)
   % Run-to-run spread of Rd and Tt against the per-run standard errors in
   % RT.se: the ratio of the standard deviation over the M runs in RTs to
   % the mean per-run standard error, which should be 1. The spread from M
   % runs is itself uncertain by about 1/sqrt(2(M-1)), 16% at M = 20, so
   % the row passes for a ratio between 0.5 and 2. A ratio above 2 points
   % to correlated runs (bad seeding) or a broken tally; a ratio below 0.5
   % to an inflated error estimate. W has the layout of vdhverdict.
   M = numel(RTs);
   x = zeros(M, 2);
   s = zeros(M, 2);
   for n = 1:M
      x(n, :) = [RTs{n}.Rdf, RTs{n}.Tt];
      s(n, :) = [RTs{n}.se.Rdf, RTs{n}.se.Tt];
   end
   ratio = (std(x, 0, 1)./mean(s, 1))';
   W.metric = {'spread(Rd)/se(Rd)'; 'spread(Tt)/se(Tt)'};
   W.model = ratio;
   W.reference = [1; 1];
   W.se = [NaN; NaN];
   W.t = [NaN; NaN];
   W.verdict = repmat({'FAIL'}, 2, 1);
   W.verdict(ratio >= 0.5 & ratio <= 2) = {'PASS'};
end
