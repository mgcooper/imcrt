function [grid,ri,ai,zi,dr,da,dz] = buildgrid(R,A,Z,dr,da,dz,N,makeplot)
   %BUILDGRID Build grids for scoring observable quantities.
   %
   %   [RI,AI,ZI,DR,DA,DZ] = BUILDGRID(R,A,Z,DR,DA,DZ) builds radial,
   %   angular, and vertical scoring grids from the extents R, A, and Z and
   %   nominal bin widths DR, DA, and DZ.
   %
   %   Each extent must contain a whole number of bins. R/DR, A/DA, and Z/DZ are
   %   checked to ensure they are whole numbers within a relative tolerance of
   %   1e-9. The resulting grid therefore has exactly the bin counts used to
   %   size the tallies with ROUND. Fractional bin counts raise the error
   %   buildgrid:nonintegral.
   %
   %   RI and ZI include one overflow coordinate, one bin width beyond the final
   %   regular bin. The radial overflow bin pools all values with r > R, so a
   %   density normalized by its nominal area is not meaningful when substantial
   %   light exits beyond R.
   %
   %   The radial and angular centers are shifted according to Equations 8 and
   %   14 of Wang et al. 1995 to reduce scoring error. These shifted centers are
   %   reporting coordinates, not bin measures; bin measures must be computed
   %   from the bin edges. The optimal vertical reporting coordinate is the
   %   geometric bin center.
   %
   %   Angular coordinates are expressed as theta in radians; mcrt scoring uses
   %   u = cos(theta).
   %
   %   BUILDGRID(...,makeplot=true) plots the original and optimized radial and
   %   angular grid centers.
   %
   %   Example:
   %      [ri,ai,zi,dr,da,dz] = buildgrid(200,pi/2,100,4,pi/50,5,Plot=true);
   %
   % See also:

   if nargin < 8
      makeplot = false;
   else
      assert(islogical(makeplot))
   end

   % Each extent must contain a whole number of bins so the grid and tallies
   % have the same length.
   nr = wholebins(R, dr, 'R/dr');
   na = wholebins(A, da, 'A/da');
   nz = wholebins(Z, dz, 'Z/dz');

   % Build the grid-center coordinates and complete the expected bin count.
   ri = centers(dr, R, nr);
   ai = centers(da, A, na);
   zi = centers(dz, Z, nz);

   % Keep the original bin centers for the comparison plot.
   if makeplot
      ri0 = ri;
      ai0 = ai(:);
   end

   % Add an overflow coordinate to the radial and vertical grids.
   ri(end+1) = ri(end) + dr;
   zi(end+1) = zi(end) + dz;

   % Adjust the radial and angular coordinates to minimize scoring error.
   % https://core.ac.uk/download/pdf/84314628.pdf
   ri = ri + (dr*dr/12)./ri;               % Equation 8
   ai = ai + cot(ai).*(1-da/2*cot(da/2));  % Equation 14
   % The optimal vertical coordinate is the center of each element.

   % Make the angular and vertical coordinates column vectors.
   ai = ai(:);
   zi = zi(:);

   % Compare the bin centers to the optimized coordinates.
   if makeplot
      plotgrids(ri0, ri(1:nr), ai0, ai);
   end

   % Compute the bin edges.
   redge = (0:nr+1)*dr;          % radial edges, overflow included      [cm]
   aedge = (0:na)'*da;           % angular edges                        [rad]

   % Compute the bin widths.
   dr = dr*ones(1,nr+1);         % radial widths                        [cm]
   da = da*ones(na,1);           % angular widths                       [rad]
   dz = dz*ones(nz+1,1);         % vertical widths, overflow included   [cm]

   % Compute the area of each annulus [cm^2].
   dA = pi*(redge(2:end).^2-redge(1:end-1).^2);

   % Compute the solid angle of each angular bin [sr].
   dsr = 2*pi*(cos(aedge(1:end-1))-cos(aedge(2:end)));

   % Assign the grid output struct.
   grid = struct('ri', ri, 'ai', ai, 'zi', zi, ...
      'dr', dr, ...
      'da', da, ...
      'dz', dz, ...
      'dA', dA, ...
      'dsr', dsr, ...
      'N', N);
end

function c = centers(width, extent, n)
   %CENTERS Return the bin-center coordinates.
   %   A rounding error can make the colon operator stop one center short, so
   %   the final center is appended when needed.

   c = width/2:width:extent-width/2;
   if numel(c) < n
      c(end+1) = (numel(c) + 0.5) * width;
   end
end

function plotgrids(ri0, ri, ai0, ai)
   %PLOTGRIDS Compare the original and optimized grid centers.

   % subplot rather than tiledlayout, which Octave lacks
   figure('Units', 'inches', 'Position', [2 2 10 6])

   subplot(2, 2, 1)
   plot(1:numel(ri0), ri0, 'o')
   hold on
   plot(1:numel(ri), ri, '.')
   set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 11)
   xlabel('radial bin')
   ylabel('r [cm]')
   legend('bin center', 'optimized (Eq. 8)', 'Location', 'northwest')

   subplot(2, 2, 2)
   plot(1:numel(ri0), ri - ri0, 'o')
   set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 11)
   xlabel('radial bin')
   ylabel('shift, optimized - center [cm]')

   subplot(2, 2, 3)
   plot(1:numel(ai0), ai0, 'o')
   hold on
   plot(1:numel(ai), ai, '.')
   set(gca, 'FontSize', 11)
   xlabel('angular bin')
   ylabel('\theta [rad]')
   legend('bin center', 'optimized (Eq. 14)', 'Location', 'northwest')

   subplot(2, 2, 4)
   plot(1:numel(ai0), ai - ai0, 'o')
   set(gca, 'YScale', 'log', 'FontSize', 11)
   xlabel('angular bin')
   ylabel('shift, optimized - center [rad]')
end
