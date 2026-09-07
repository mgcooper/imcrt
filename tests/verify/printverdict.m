function npass = printverdict(V, fid)
   % Print a verdict struct (see vdhverdict) as a fixed-width table with
   % fprintf, which runs on core Octave as well as MATLAB. Return the
   % number of rows that pass. The optional fid sends the table to an open
   % file instead of the screen, for the overnight report.
   if nargin < 2
      fid = 1;
   end
   fprintf(fid, '%-18s %12s %12s %11s %8s  %s\n', 'metric', 'model', ...
      'reference', 'se', 't', 'verdict');
   for n = 1:numel(V.metric)
      fprintf(fid, '%-18s %12.6g %12.6g %11.3g %8.2f  %s\n', V.metric{n}, ...
         V.model(n), V.reference(n), V.se(n), V.t(n), V.verdict{n});
   end
   npass = nnz(strcmp(V.verdict, 'PASS'));
end
