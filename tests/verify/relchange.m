function s = relchange(a, b)
   % Relative change from a to b as signed percent text, or 'n/a' when a
   % is zero or either value is undefined (NaN), so a table never divides
   % by zero.
   if a == 0 || isnan(a) || isnan(b)
      s = 'n/a';
   else
      s = sprintf('%+.2f%%', (b - a)/abs(a)*100);
   end
end
