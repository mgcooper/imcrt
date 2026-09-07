function npass = printverdict(V)
   % Print a verdict struct (see vdhverdict) as a fixed-width table with
   % fprintf, which runs on core Octave as well as MATLAB. Return the
   % number of rows that pass.
   fprintf('%-14s %12s %12s %11s %8s  %s\n', 'metric', 'model', ...
      'reference', 'se', 't', 'verdict');
   for n = 1:numel(V.metric)
      fprintf('%-14s %12.6g %12.6g %11.3g %8.2f  %s\n', V.metric{n}, ...
         V.model(n), V.reference(n), V.se(n), V.t(n), V.verdict{n});
   end
   npass = nnz(strcmp(V.verdict, 'PASS'));
end
