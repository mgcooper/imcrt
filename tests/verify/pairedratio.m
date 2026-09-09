function s = pairedratio(d)
   % Text for the mean of paired differences d in units of its standard
   % error: the ratio, 'none' when every difference is zero, or 'n/a' when
   % the differences are constant but not zero (no spread to scale by).
   sd = std(d)/sqrt(numel(d));
   if sd > 0
      s = sprintf('%.2f', mean(d)/sd);
   elseif mean(d) == 0
      s = 'none';
   else
      s = 'n/a';
   end
end
