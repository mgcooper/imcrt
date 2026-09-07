function [reportfile, N, M] = impactargs(reportfile, N, M)
   % Defaults for impactreport. A missing or empty argument takes the
   % documented value. The defaults live apart from the generator so a
   % test can check them without a report run at the default N.
   if nargin < 1 || isempty(reportfile)
      reportfile = 'docs/impact-report.md';
   end
   if nargin < 2 || isempty(N)
      N = 1e5;
   end
   if nargin < 3 || isempty(M)
      M = 8;
   end
   % A fractional N would run fix(N) packets and normalize by N, and a
   % fractional M fails only after the report is partly written. The
   % retrospective's standard errors also need a run spread.
   assert(isnumeric(N) && isreal(N) && isscalar(N) && isfinite(N) ...
      && N == fix(N) && N >= 1, 'impactargs:packets', ...
      'N must be a positive integer');
   assert(isnumeric(M) && isreal(M) && isscalar(M) && isfinite(M) ...
      && M == fix(M) && M >= 2, 'impactargs:runs', ...
      'M must be an integer of at least 2');
end
