function i = binindex(v, d, n)
   %BININDEX Return the bin index for a coordinate on a uniform grid.
   %
   %   I = BININDEX(V,D,N) returns the bin containing coordinate V for a
   %   grid with bin width D and N bins. The result is clamped to 1:N.
   %
   %   V = 0 is assigned to the first bin. Values that round beyond the
   %   final edge are assigned to bin N. If N includes an overflow bin,
   %   overflow values are assigned to the overflow bin.

   i = ceil(v/d);

   if i < 1
      i = 1;
   elseif i > n
      i = n;
   end
end
