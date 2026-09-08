function [z, model, se] = vdhangular(RTs, ref, perrun)
   % Compare the diffuse angular tallies of M independent runs (cell array
   % RTs) with the van de Hulst table at its nonzero mu. Rows 1 to 6 are
   % reflectance and rows 7 to 12 transmittance at mu = [0.1 0.3 0.5 0.7
   % 0.9 1.0]. Each run is interpolated with pchip between the shifted bin
   % centers RT.grid.ai, the points where a bin average equals the point
   % value. pchip halves the chord bias of linear interpolation where the
   % tallies curve near grazing. The mu = 1 row uses the first bin because
   % theta = 0 lies inside it; a 3-degree bin average sits about 1.6% below
   % the transmittance peak there. The tallies hold packet weights, so no
   % per-bin binomial error exists: the standard error is the spread of the
   % interpolated values over runs divided by sqrt(M), which also covers
   % the interpolation. With perrun true (few runs), each run's error at a
   % table angle is propagated from the per-bin errors in RT.se through
   % the pchip interpolant by the delta method: the sensitivity of the
   % interpolated value to each bin, from a finite difference, applied to
   % the covariance of the bins, whose diagonal is the bin variance and
   % whose off-diagonal, -x_i*x_j/(N-1), comes from the bins being
   % exclusive. The runs combine as sqrt(sum of squares)/M. That leaves
   % the interpolation error out. mu = 0 is left out because its
   % reference is zero through the mu factor.
   if nargin < 3
      perrun = false;
   end
   mu = ref.mu(2:end);
   theta = acos(mu);
   M = numel(RTs);
   V = zeros(M, 12);
   S = zeros(M, 12);
   for n = 1:M
      % Reflectance fills columns 1 to 6 and transmittance 7 to 12, one
      % row per run, so the spread over rows is the run-to-run error.
      ai = RTs{n}.grid.ai;
      V(n, 1:6) = interp1(ai, RTs{n}.Rdf_a, theta, 'pchip', RTs{n}.Rdf_a(1));
      V(n, 7:12) = interp1(ai, RTs{n}.Tdf_a, theta, 'pchip', RTs{n}.Tdf_a(1));
      if perrun
         S(n, 1:6) = pchipvar(ai, RTs{n}.Rdf_a, RTs{n}.se.Rdf_a, theta, ...
            RTs{n}.grid.N);
         S(n, 7:12) = pchipvar(ai, RTs{n}.Tdf_a, RTs{n}.se.Tdf_a, theta, ...
            RTs{n}.grid.N);
      end
   end
   model = mean(V, 1)';
   if perrun
      se = sqrt(sum(S, 1))'/M;
   else
      se = std(V, 0, 1)'/sqrt(M);
   end
   z = (model - [ref.R_sr(2:end) ref.T_sr(2:end)]')./se;
end

function v = pchipvar(ai, x, se, theta, N)
   % Variance of the pchip interpolant of the bin values x at the angles
   % theta, from the standard errors se of the bins and the packet count
   % N by the delta method: w' C w with w the sensitivities and C the bin
   % covariance, se.^2 on the diagonal and -x_i*x_j/(N-1) off it, the
   % covariance of exclusive bins. pchip is nonlinear in its data (its
   % slopes are harmonic means of secants), so each sensitivity is a
   % finite difference with a step small against the values. The mu = 1
   % angle sits inside bin 1 and returns bin 1 alone, which the
   % difference reproduces. A side with no exit has zero tallies: its
   % variance is zero when its errors are zero, and stays NaN when one
   % packet left them undefined.
   v = zeros(1, numel(theta));
   h = 1e-6*max(abs(x));
   if h == 0
      if any(isnan(se))
         v(:) = NaN;
      end
      return
   end
   base = interp1(ai, x, theta, 'pchip', x(1));
   W = zeros(numel(theta), numel(x));
   for i = 1:numel(x)
      xp = x;
      xp(i) = xp(i) + h;
      W(:, i) = (interp1(ai, xp, theta, 'pchip', xp(1)) - base)/h;
   end
   x = x(:);
   C = -(x*x')/(N - 1);
   C(1:numel(x)+1:end) = se(:).^2;
   v = sum((W*C).*W, 2)';
end
