function reportfile = vdhovernight(outdir, scale, tmax)
   %VDHOVERNIGHT Multi-run verification of mcrt with checkpoints and a
   % dated PASS/FAIL report.
   %
   %  reportfile = vdhovernight(outdir) runs case reflect (van de Hulst
   %  Table 35) and case fluence (self-consistency) at the full sizes in
   %  verifycases, M independent seeded runs each. It saves every run to
   %  outdir/checkpoint_<case>_N<N>_seed<s>.mat as it completes and writes
   %  outdir/vdh_report_<yyyyMMdd_HHmmss>.txt, whose path it returns.
   %
   %  A checkpoint holds the run's RT and a meta struct with:
   %    - the SHA-256 of every file in kernelfiles together (the kernel,
   %      the functions its loop calls, the grid builder, and the
   %      derivative; every file that shapes a
   %      result);
   %    - the MATLAB version;
   %    - the kernel inputs [ka ks g Z dz N seed];
   %    - the run time in seconds.
   %  A checkpoint whose meta matches the current kernel and inputs is
   %  loaded instead of recomputed, so an interrupted job resumes where it
   %  stopped. A checkpoint with different meta, a missing field, or a file
   %  that does not load is stale and is recomputed. Each file is written
   %  to a temporary name and renamed, so a stop during save leaves no
   %  half-written checkpoint behind.
   %
   %  vdhovernight(outdir, scale) multiplies every N by scale. Run
   %  scale = 1e-2 first as the mandatory dry run. It exercises every step
   %  in seconds and prints the runtime estimate for scale = 1. The N in a
   %  checkpoint name keeps dry-run files apart from full-run files.
   %  vdhovernight(outdir, scale, tmax) replaces every case's t cutoff;
   %  the tests use tmax = 0 to reach the OVERALL FAIL path.
   %
   %  Verdict rules. A statistical row passes when |t| <= tmax with M - 1
   %  degrees of freedom, t from the spread over runs. The hemispherical
   %  Rd and Tt rows also pass within 5e-4 of the table, which allows for
   %  the table's own precision. An angular row also passes within 2% of
   %  the table. That allows for the 3-degree bin averaging and the
   %  interpolation between bin centers, which a 1e8-packet run resolves
   %  at several standard errors. The report names every row that passed
   %  only by an allowance. Tdr has no allowance. Rdr must be exactly zero
   %  in every run. The run-to-run variance of Rd and Tt must be at most
   %  twice the binomial estimate. The fluence case has no numeric
   %  reference in this repository and reports self-consistency only.
   %  OVERALL is PASS when every row of both cases passes.
   if nargin < 2
      scale = 1;
   end
   if nargin < 3
      tmax = [];
   end
   outdir = char(outdir);
   if ~exist(outdir, 'dir')
      mkdir(outdir);
   end

   % The kernel and helpers must resolve from wherever the caller sits.
   root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
   run(fullfile(root, 'Setup.m'));
   addpath(fullfile(root, 'tests', 'verify'));

   kernel = kernelhash(fullfile(root, kernelfiles()));
   reportfile = reportname(outdir, ...
      char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')));
   fid = fopen(reportfile, 'w');
   if fid < 0
      error('vdhovernight:report', 'Cannot open %s', reportfile);
   end
   % An interrupt or an error in a run must not leave the report open and
   % unflushed; the cleanup closes it if the normal fclose is skipped.
   cleaner = onCleanup(@() closeifopen(fid));
   fprintf(fid, 'imcrt overnight verification, %s\nMATLAB %s\n', ...
      char(datetime('now')), version);
   fprintf(fid, 'kernel sha256 %s of %s\n', kernel, ...
      strjoin(kernelfiles(), '+'));
   fprintf(fid, 'scale %g (every N multiplied by scale)\n\n', scale);
   overall = true;
   tnew = 0;
   tall = 0;
   names = {'reflect', 'fluence'};
   for n = 1:numel(names)
      c = verifycases(names{n});
      f = c.full;
      if ~isempty(tmax)
         f.tmax = tmax;
      end
      N = max(1, round(f.N*scale));
      fprintf(fid, '== case %s: M=%d runs, N=%d, seeds %d:%d, tmax %g\n', ...
         c.name, f.M, N, f.seeds(1), f.seeds(end), f.tmax);
      [RTs, tcase, trecorded] = runcheckpointed(c, f, N, outdir, kernel, fid);
      [npass, nrow] = reportcase(RTs, c, N, f.tmax, fid);
      fprintf(fid, 'case %s: %d/%d pass; %.0f s new, %.0f s all\n\n', ...
         c.name, npass, nrow, tcase, trecorded);
      overall = overall && npass == nrow;
      tnew = tnew + tcase;
      tall = tall + trecorded;
   end
   fprintf(fid, 'compute %.0f s now, %.0f s over all runs at scale %g\n', ...
      tnew, tall, scale);
   fprintf(fid, 'scale-1 estimate from all runs: %.0f s\n', tall/scale);
   if overall
      word = 'PASS';
   else
      word = 'FAIL';
   end
   fprintf(fid, 'OVERALL %s\n', word);
   fclose(fid);
   clear cleaner
   fprintf('%sreport: %s\n', fileread(reportfile), reportfile);
end

function printprofile(RTs, fid)
   % Write the fluence depth profile phi_z, mean and run spread over the M
   % runs, at every 20th depth bin. It is the main fluence result (Wang et
   % al. 1995, Fig. 4) and has no numeric reference in this repository;
   % the checkpoints hold every run's full profile.
   M = numel(RTs);
   z = RTs{1}.grid.zi;
   P = zeros(numel(z), M);
   for n = 1:M
      P(:, n) = RTs{n}.phi_z;
   end
   % The last grid coordinate is the overflow bin beyond the slab; the
   % physical profile stops one bin before it.
   rows = 1:20:numel(z)-1;
   fprintf(fid, 'phi_z profile (mean and standard error over %d runs)\n', M);
   fprintf(fid, '%10s %12s %12s\n', 'z', 'phi_z', 'se');
   for n = rows
      fprintf(fid, '%10.4f %12.5g %12.3g\n', z(n), mean(P(n, :)), ...
         std(P(n, :))/sqrt(M));
   end
end

function closeifopen(fid)
   % Close the report when an error or interrupt skipped the normal fclose.
   % fclose on an identifier that is already closed raises, and that is the
   % normal case after a clean run, so the error is ignored.
   try
      fclose(fid);
   catch
   end
end

function [RTs, tnew, tall] = runcheckpointed(c, f, N, outdir, kernel, fid)
   % Run or load the M runs of one case. A checkpoint is reused only when
   % its meta matches the current kernel hash, MATLAB version, and inputs;
   % otherwise it is stale and recomputed. tnew is the compute time of
   % this invocation and tall the recorded time of every run, loaded or
   % new, so the runtime estimate survives a resume.
   RTs = cell(1, f.M);
   tnew = 0;
   tall = 0;
   for n = 1:f.M
      seed = f.seeds(n);
      inputs = [c.ka, c.ks, c.g, c.Z, c.dz, N, seed];
      ck = fullfile(outdir, sprintf('checkpoint_%s_N%d_seed%d.mat', ...
         c.name, N, seed));
      fresh = false;
      if exist(ck, 'file')
         % Validate every meta field before reading it, so a checkpoint
         % whose metadata is older, partial, unreadable, or malformed counts
         % as stale, not as an error. Any failure inside the check, such as
         % a nonscalar text field, also means stale. The RT payload itself
         % is trusted once its metadata matches.
         try
            S = load(ck);
            fresh = isfield(S, 'RT') && isfield(S, 'meta') ...
               && isstruct(S.meta) && isscalar(S.meta) ...
               && all(isfield(S.meta, ...
                  {'kernel', 'version', 'inputs', 'elapsed'})) ...
               && ischar(S.meta.kernel) && strcmp(S.meta.kernel, kernel) ...
               && ischar(S.meta.version) ...
               && strcmp(S.meta.version, version) ...
               && isequal(S.meta.inputs, inputs) ...
               && isnumeric(S.meta.elapsed) && isscalar(S.meta.elapsed) ...
               && isfinite(S.meta.elapsed) && S.meta.elapsed >= 0;
         catch
            fresh = false;
         end
      end
      if fresh
         RTs{n} = S.RT;
         tall = tall + S.meta.elapsed;
         fprintf(fid, 'run %2d seed %2d: loaded %s\n', n, seed, ck);
      else
         if exist(ck, 'file')
            fprintf(fid, 'run %2d seed %2d: stale checkpoint, recompute\n', ...
               n, seed);
         end
         rng(seed, 'twister');
         t0 = tic;
         RT = mcrt(c.ka, c.ks, c.g, c.Z, c.dz, N);
         elapsed = toc(t0);
         meta = struct('kernel', kernel, 'version', version, ...
            'inputs', inputs, 'elapsed', elapsed);
         % Write to a temporary name and rename so a stop during save never
         % leaves a half-written file under the checkpoint name.
         tmp = [ck(1:end-4) '_tmp.mat'];
         save(tmp, 'RT', 'meta');
         movefile(tmp, ck);
         RTs{n} = RT;
         tnew = tnew + elapsed;
         tall = tall + elapsed;
         fprintf(fid, 'run %2d seed %2d: %.1f s, saved %s\n', n, seed, ...
            elapsed, ck);
      end
   end
end

function [npass, nrow] = reportcase(RTs, c, N, tmax, fid)
   % Print the verdict of one case to the report and return the pass and
   % row counts. reflect adds the variance check; fluence has none.
   abstol = 5e-4;
   reltol = 0.02;
   switch c.name
      case 'reflect'
         V = vdhverdict(RTs, vdhtable35(), tmax, abstol, reltol);
         W = variancecheck(RTs, N);
      case 'fluence'
         V = fluenceverdict(RTs, c, tmax);
         W = [];
         printprofile(RTs, fid);
   end
   npass = printverdict(V, fid);
   nrow = numel(V.verdict);
   if isfield(V, 'byallowance') && any(V.byallowance)
      fprintf(fid, 'passed only by allowance (%g abs, %g rel): %s\n', ...
         abstol, reltol, strjoin(V.metric(V.byallowance), ', '));
   end
   if ~isempty(W)
      npass = npass + printverdict(W, fid);
      nrow = nrow + numel(W.verdict);
   end
end

function h = kernelhash(files)
   % SHA-256 of the concatenated bytes of the files (cellstr) as lowercase
   % hex, through the JVM that MATLAB ships. Core MATLAB has no digest
   % function; this tool is MATLAB-only, so the JVM dependency is
   % acceptable here. Every file that shapes a result belongs in the list.
   md = java.security.MessageDigest.getInstance('SHA-256');
   for n = 1:numel(files)
      md.update(uint8(fileread(files{n})));
   end
   h = lower(reshape(dec2hex(typecast(md.digest(), 'uint8'), 2)', 1, []));
end
