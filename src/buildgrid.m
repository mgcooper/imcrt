function [ri,ai,zi,dr,da,dz] = buildgrid(R,A,Z,dr,da,dz)
   %BUILDGRID build grids for scoring observable quantities
   %
   % Each extent must be a whole number of bins: R/dr, A/da, and Z/dz are
   % checked to be integral within 1e-9 relative, and a center the colon
   % operator drops when the ratio sits a rounding error below the integer
   % is appended, so the grid has the bins mcrt sizes its tallies with. A
   % fractional count raises buildgrid:nonintegral.

   % note that U/du are in radians here i.e. theta, but for scoring within the
   % mcrt program they need to be in u = cos(theta)

   % example:
   % R           = 200;        % radius of detection [cm]
   % A           = pi/2;       % angular detection radius
   % Z           = 100;        % thickness of medium
   % nr          = 50;         % number of radial bins
   % na          = 25;         % number of angular bins

   % number of bins in each dimension
   %     nr      = roundn(R/dr,0);     % radial
   %     na      = roundn(A/da,0);     % angular
   %     nz      = roundn(Z/dz,0);     % number

   % every extent is a whole number of bins, or the grid and the tallies
   % would disagree in length
   nr = wholebins(R, dr, 'R/dr');
   na = wholebins(A, da, 'A/da');
   nz = wholebins(Z, dz, 'Z/dz');

   % grid center coordinates, completed to the whole count
   ri = centers(dr, R, nr);
   zi = centers(dz, Z, nz);
   ai = centers(da, A, na);

   % r and z need an extra overflow coordinat
   zi(end+1) = zi(end)+dz;
   ri(end+1) = ri(end)+dr;

   % optimized to minimize error (https://core.ac.uk/download/pdf/84314628.pdf)
   ri = ri + (dr*dr/12)./ri;              % Eq. 8
   ai = ai + cot(ai).*(1-da/2*cot(da/2)); % Eq. 14
   % the optimal z-coordinate is the center of each element i.e. zi

   % make ai and zi columns
   ai = ai(:);
   zi = zi(:);

   % get new da/dr/dz ('derivative' can be downloaded on the file exchange).
   % derivative needs two points; a one-bin angular grid keeps its width.
   if numel(ai) > 1
      da = derivative(ai);
   end
   dr = derivative(ri);
   dz = derivative(zi);

   % see the difference
   %     figure;
   %     tiledlayout(1,2); nexttile;
   %     plot(dr/2:dr:R-dr/2,'o'); hold on; ylabel('radial grid')
   %     plot(ri,'o'); legend('grid centers','optimized');
   %     nexttile;
   %     plot((dr/2:dr:R-dr/2)-ri,'o'); legend('grid centers - optimized');
   %
   %     figure;
   %     tiledlayout(1,2); nexttile;
   %     plot(du/2:du:U-du/2,'o'); hold on; ylabel('Angular grid')
   %     plot(ui,'o'); legend('grid centers','optimized');
   %     nexttile;
   %     plot((du/2:du:U-du/2)-ui,'o'); legend('grid centers - optimized');
end

function c = centers(width, extent, n)
   % Bin centers from the colon operator, which keeps the values of the
   % original grid; when a ratio a rounding error below the integer makes
   % the colon stop one bin short, the last center is appended.
   c = width/2:width:extent-width/2;
   if numel(c) < n
      c(end+1) = (numel(c) + 0.5)*width;
   end
end

function n = wholebins(extent, width, name)
   % The extent divided by the width must be a whole number within 1e-9
   % relative, which lets 0.02/0.001 pass despite double-precision error
   % (single precision is not served: its error is about 1e-7). Returns
   % that whole number, which must be one or more.
   % a width whose half underflows to zero cannot start a center grid
   assert(width/2 > 0, 'buildgrid:width', ...
      'the bin width for %s is too small for a half-width', name);
   ratio = extent/width;
   n = round(ratio);
   assert(abs(ratio - n) <= 1e-9*max(1, abs(ratio)) && n >= 1, ...
      'buildgrid:nonintegral', ...
      '%s = %.12g must be a whole number of bins, one or more', name, ratio);
end
