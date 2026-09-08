function out = perfbench(mode, file)
   % Time the kernel on perfcases and keep a baseline for regression.
   %
   %  T = perfbench() returns a struct array with fields name, seconds,
   %  and reference. Each case is timed three times with timeit, and a
   %  fixed scalar reference loop is timed right after each of those;
   %  seconds is the least of the three case timings and reference is the
   %  reference timing paired with it. The least of three is the timing
   %  least disturbed by other processes, and the pairing keeps the two
   %  numbers from the same machine state. A machine under load slows the
   %  kernel and the reference alike, so the regression test compares
   %  seconds divided by reference, not raw seconds. Seeding inside the
   %  timed function keeps every repetition on the same packet paths.
   %  perfbench('path') returns the baseline file tests/perf/baseline.txt.
   %  perfbench('write') measures and writes the baseline; a second
   %  argument writes to another file. perfbench('read') or
   %  perfbench('read', file) parses a baseline into a struct with fields
   %  computer, host, version, and cases (name, seconds, reference). The
   %  baseline
   %  holds the architecture, the host name, and the MATLAB release it was
   %  measured on, because a time from another machine is not a
   %  regression bar. perfbench('host') returns this machine's host name,
   %  the first label of hostname, which is stable across networks.
   if nargin < 1
      mode = '';
   end
   here = fileparts(mfilename('fullpath'));
   if nargin < 2
      file = fullfile(here, 'baseline.txt');
   end
   switch mode
      case 'path'
         out = file;
         return
      case 'read'
         out = readbaseline(file);
         return
      case 'host'
         out = hostname();
         return
   end
   cases = perfcases();
   out = struct('name', {cases.name}, ...
      'seconds', num2cell(zeros(1, numel(cases))), ...
      'reference', num2cell(zeros(1, numel(cases))));
   for k = 1:numel(cases)
      c = cases(k);
      [out(k).seconds, out(k).reference] = leastpaired( ...
         @() seededrun(c), @() refloop(1e6), 3);
   end
   if strcmp(mode, 'write')
      fid = fopen(file, 'w');
      assert(fid > 0, 'perfbench:open', 'Cannot open %s', file);
      fprintf(fid, ['# imcrt perf baseline: computer %s, host %s, ' ...
         'MATLAB %s, %s\n'], computer, hostname(), version, ...
         char(datetime('now', 'Format', 'yyyy-MM-dd')));
      for k = 1:numel(out)
         fprintf(fid, '%s %.6f %.6f\n', out(k).name, out(k).seconds, ...
            out(k).reference);
      end
      fclose(fid);
   end
end

function [t, r] = leastpaired(f, g, k)
   % The least of k timeit measurements of f, and the measurement of g
   % taken right after that one, so the pair comes from one machine
   % state.
   t = Inf;
   r = Inf;
   for n = 1:k
      tn = timeit(f);
      rn = timeit(g);
      if tn < t
         t = tn;
         r = rn;
      end
   end
end

function s = refloop(ndraw)
   % A scalar loop of the kernel's flavor, a few flops and a draw per
   % iteration, that takes about a tenth of a second. Its time tracks the
   % machine's speed at the moment of measurement.
   s = 0;
   for n = 1:ndraw
      x = rand;
      s = s + sqrt(1 - x*x)*cos(x) - log(x);
   end
end

function RT = seededrun(c)
   % One seeded run of a case, so every timeit repetition is the same work.
   rng(c.seed, 'twister');
   RT = mcrt(c.ka, c.ks, c.g, c.Z, c.dz, c.N);
end

function h = hostname()
   % First label of the host name, so a laptop keeps one name on every
   % network. hostname is a command on every platform MATLAB runs on.
   [~, h] = system('hostname');
   h = regexprep(strtrim(h), '\..*$', '');
end

function b = readbaseline(file)
   % Parse the header line for the architecture, host, and release, then
   % one case per line as "name seconds".
   text = fileread(file);
   lines = regexp(strtrim(text), '\r?\n', 'split');
   tok = regexp(lines{1}, ['^# imcrt perf baseline: computer (\S+), ' ...
      'host (\S+), MATLAB (.*), (\S+)$'], 'tokens', 'once');
   assert(~isempty(tok), 'perfbench:header', ...
      'Bad baseline header in %s', file);
   b = struct('computer', tok{1}, 'host', tok{2}, 'version', tok{3}, ...
      'date', tok{4});
   b.cases = struct('name', {}, 'seconds', {}, 'reference', {});
   for n = 2:numel(lines)
      parts = strsplit(strtrim(lines{n}));
      b.cases(end+1) = struct('name', parts{1}, ...
         'seconds', str2double(parts{2}), 'reference', str2double(parts{3}));
   end
end
