function se = tallyse(ss, s, N)
   % Standard error of a tally s, the sum of packet weights that landed
   % in a bin over N packets, from ss, the sum of their squares. Each
   % packet contributes its weight once or nothing, so the per-packet
   % contributions have mean s/N and unbiased sample variance
   % (ss - s^2/N)/(N-1), and the mean over N packets has standard error
   % sqrt((ss - s^2/N)/(N*(N-1))). One packet gives no estimate (NaN). The
   % estimate holds for tallies a packet touches at most once (an exit
   % lands in one bin); a packet absorbs in the same bin many times, so
   % its per-step squares would understate that variance. Arrays are
   % handled elementwise, and a bin summed over other bins uses the
   % summed squares, which stay exact because the contributions of one
   % packet to different bins never overlap.
   se = sqrt(max(ss - s.^2/N, 0)/(N*(N-1)));
end
