function [z, model, se] = vdhangular(RTs, ref)
   % Compare the diffuse angular tallies of M independent runs (cell array
   % RTs) with the van de Hulst table at its nonzero mu. Rows 1 to 6 are
   % reflectance and rows 7 to 12 transmittance at mu = [0.1 0.3 0.5 0.7
   % 0.9 1.0]. Each run is interpolated linearly between the shifted bin
   % centers RT.grid.ai, the points where a bin average equals the point
   % value; the mu = 1 row uses the first bin because theta = 0 lies inside
   % it. The tallies hold packet weights, so no per-bin binomial error
   % exists: the standard error is the spread of the interpolated values
   % over runs divided by sqrt(M), which also covers the interpolation.
   % mu = 0 is left out because its reference is zero through the mu factor.
   mu = ref.mu(2:end);
   theta = acos(mu);
   M = numel(RTs);
   V = zeros(M, 12);
   for n = 1:M
      % Reflectance fills columns 1 to 6 and transmittance 7 to 12, one
      % row per run, so the spread over rows is the run-to-run error.
      ai = RTs{n}.grid.ai;
      V(n, 1:6) = interp1(ai, RTs{n}.Rdf_a, theta, 'linear', RTs{n}.Rdf_a(1));
      V(n, 7:12) = interp1(ai, RTs{n}.Tdf_a, theta, 'linear', RTs{n}.Tdf_a(1));
   end
   model = mean(V, 1)';
   se = std(V, 0, 1)'/sqrt(M);
   z = (model - [ref.R_sr(2:end) ref.T_sr(2:end)]')./se;
end
