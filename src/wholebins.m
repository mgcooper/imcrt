function n = wholebins(extent, width, name)
   %WHOLEBINS Return the number of whole bins in an extent.
   %   The extent-to-width ratio must be a positive whole number within a
   %   relative tolerance of 1e-9.

   % A width whose half underflows to zero cannot define the first grid center.
   assert(width/2 > 0, 'buildgrid:width', ...
      'the bin width for %s is too small for a half-width', name);

   ratio = extent/width;
   n = round(ratio);

   assert(abs(ratio - n) <= 1e-9*max(1, abs(ratio)) && n >= 1, ...
      'buildgrid:nonintegral', ...
      '%s = %.12g must be a whole number of bins, one or more', name, ratio);
end
