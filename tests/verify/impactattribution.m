function fixes = impactattribution(values, labels)
   % Attribute each row of values (quantities x versions, in fix order) to
   % the versions at which it changed by more than 1e-6 relative to the
   % previous version. A value that becomes defined or undefined (one side
   % NaN) is a change; two NaNs are unchanged. labels names the versions
   % from the second column on. fixes is a cellstr, one entry per row: the
   % changed labels joined by commas, or 'none' when every change is
   % rounding only.
   nq = size(values, 1);
   fixes = cell(nq, 1);
   for n = 1:nq
      prev = values(n, 1:end-1);
      next = values(n, 2:end);
      changed = abs(next - prev) > 1e-6*max(abs(prev), 1e-300) ...
         | xor(isnan(prev), isnan(next));
      if any(changed)
         fixes{n} = strjoin(labels(changed), ', ');
      else
         fixes{n} = 'none';
      end
   end
end
