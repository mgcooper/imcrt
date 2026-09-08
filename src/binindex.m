function i = binindex(v, d, n)
   % Bin index of a coordinate v on a grid of width d with n bins, clamped
   % into 1:n. ceil maps v = 0 exactly (a packet on the axis or on the
   % surface) to bin 0, and rounding at the far edge can give n+1; both
   % must land in a valid bin. The caller passes n as the count including
   % any overflow bin, so the clamp above is the overflow bin.
   i = ceil(v/d);
   if i < 1
      i = 1;
   elseif i > n
      i = n;
   end
end
