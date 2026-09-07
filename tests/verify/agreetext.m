function s = agreetext(a, b)
   % Phrase for two arrays of results: identical when every entry agrees,
   % with paired NaN entries (an undefined r50 on both sides) counted as
   % agreement, otherwise the largest difference over the entries that
   % both sides define. The report states agreement from the data, not
   % from prose.
   if isequaln(a, b)
      s = 'are identical run for run';
   else
      s = sprintf('differ by at most %.3g where both are defined', ...
         max(abs(a(:) - b(:)), [], 'omitnan'));
   end
end
