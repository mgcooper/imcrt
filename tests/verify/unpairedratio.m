function s = unpairedratio(a, b)
   % Text for the difference of two sample means, a minus b, in units of
   % its standard error for samples of equal size that are not paired:
   % the ratio, 'none' when both samples are the same constant, or 'n/a'
   % when they are different constants (no spread to scale by).
   se = sqrt((var(a) + var(b))/numel(a));
   d = mean(a) - mean(b);
   if se > 0
      s = sprintf('%.2f', d/se);
   elseif d == 0
      s = 'none';
   else
      s = 'n/a';
   end
end
