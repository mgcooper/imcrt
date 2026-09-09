function wt = roulette(wt, wrr)
   % Russian roulette (Wang et al. 1995, Sect. 3.9) on a packet of weight
   % wt: it survives with probability 1/wrr carrying wrr times its weight,
   % or dies with weight 0. The expected weight is unchanged. The caller
   % plays only below its threshold wmin; a survivor can still sit below
   % it and then plays again on the next step.
   if rand < 1/wrr
      wt = wt*wrr;
   else
      wt = 0;
   end
end
