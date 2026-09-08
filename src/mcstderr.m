function se = mcstderr(ss, s, N)
   %MCSTDERR Compute the standard error of a Monte Carlo tally.
   % 
   %   SE = MCSTDERR(SS,S,N) computes the standard error of the mean tally
   %   per packet from S, the sum of packet contributions, SS, the sum of
   %   squared packet contributions, and N, the number of packets:
   %
   %      SE = sqrt((SS - S.^2/N)/(N*(N-1))).
   %
   %   The estimate assumes each packet contributes to a given tally bin at
   %   most once. It therefore applies to exit tallies, where each packet
   %   lands in one bin, but not to tallies that accumulate multiple
   %   contributions from the same packet in the same bin. A single packet
   %   provides no variance estimate and returns NaN. Array inputs are handled
   %   elementwise.
   %
   %   Tallies may be summed across mutually exclusive bins by summing their
   %   corresponding squared contributions because a packet cannot contribute
   %   to more than one such bin.
   %
   %   Mutually exclusive bins are negatively correlated: a packet that lands
   %   in one cannot land in another. For the normalized means x_i = S_i/N the
   %   covariance is -x_i*x_j/(N-1), which the standard errors above leave
   %   out because each refers to one bin or to one sum over bins. A quantity
   %   that combines several bins with different weights, such as an
   %   interpolated angular value, needs that covariance; vdhangular applies
   %   it.

   se = sqrt(max(ss - s.^2/N, 0)/(N*(N-1)));
end
