function reportfile = impactreport(varargin)
   %IMPACTREPORT Write the impact and paper-kernel retrospective report.
   %
   %  impactreport(reportfile, N, M) writes a Markdown report with two
   %  sections and returns its path. Defaults (impactargs): reportfile
   %  tests/reports/impact-report.md, N 1e5 packets, M 8 runs.
   %
   %  Section 1, impact. The cases are the shared cases in mcrtcases plus
   %  a paper-scale radial case (albedo 0.999, g 0.9, optical depth 10).
   %  Its lengths are scaled to the kernel's fixed 2 cm radial grid; the
   %  paper's rod geometry is not in src/mcrt.m. The kernel at the pre-fix
   %  commit and at each fix tag runs every case at N with the case's
   %  seed. One table per case lists every quantity from
   %  impactquantities. Its columns are the pre-fix value, the post-fix
   %  value, the relative change, the fixes at which it changed
   %  (impactattribution, bisected through the tags), and a physical note.
   %  The kernels come from git (src/mcrt.m, src/buildgrid.m, and
   %  src/derivative/derivative.m at each ref) into a temporary folder,
   %  as do the frozen 2021 script and its helpers at the pre-fix commit.
   %
   %  Section 2, retrospective: the van de Hulst case through the frozen
   %  2021 verification script examples/cooper_etal_2021/c_model/
   %  mcrt_verify.m and through the fixed kernel, M seeded runs of N each.
   %  The script runs from a scratch copy with six edits (scratch2021):
   %  the geometry load line is removed, N and wmin come from the
   %  workspace, roundn becomes round, the text stops before the
   %  validation section, and the direct clause is the production
   %  kernel's. The paired runs use the script's wmin = 1e-4, which the
   %  fixed kernel shares. The production kernel used 1e-5, so the copy
   %  runs the same seeds again at that threshold, and the text under
   %  the table states whether the complete tallies of the two
   %  thresholds agree. A second copy also takes the production kernel's
   %  random-number use (path length from log(1-rand), hg4 = 1-g,
   %  hg5 = 2*g), which relabels the random stream without changing any
   %  sampled distribution; its hemispherical shifts from the paired
   %  copy are listed in standard errors of the unpaired difference. The
   %  production script itself cannot run the slab case: it scores
   %  transmittance only, has no lower boundary, and carries the rod
   %  machinery. The
   %  two kernels share seeds, so their runs are paired. The table lists
   %  each quantity's mean and standard error for both kernels and the
   %  mean paired difference in units of its own standard error, which
   %  bounds how the published numbers relate to the corrected model.
   [reportfile, N, M] = impactargs(varargin{:});
   root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
   run(fullfile(root, 'Setup.m'));
   addpath(fullfile(root, 'tests'), fullfile(root, 'tests', 'verify'));

   % Every extracted kernel and the scratch copy of the 2021 script live
   % in one owned temporary root that the cleanup removes. The runs also
   % execute from that root: MATLAB resolves current-folder files before
   % the path, so a caller sitting in src/ would otherwise shadow every
   % extracted kernel with the worktree one.
   scratch = tempname;
   mkdir(scratch);
   cleaner = onCleanup(@() rmdir(scratch, 's'));

   refs = {'b9d51d5', 'pre-fix'; 'fix-B-aliasing', 'B'; ...
      'fix-A-direct-tally', 'A'; 'fix-C-fluence', 'C'; ...
      'fix-N-roulette', 'N'; 'fix-S-bin-measures', 'S'};
   labels = refs(2:end, 2)';
   [srcsame, srcclean] = srcmatches(root, 'fix-S-bin-measures');

   fid = fopen(reportfile, 'w');
   assert(fid > 0, 'impactreport:open', 'Cannot open %s', reportfile);
   closer = onCleanup(@() closeifopen(fid));
   % Move into the scratch root only after the report is open, so a
   % relative report path still resolves against the caller's folder.
   startdir = cd(scratch);
   goback = onCleanup(@() cd(startdir));
   [~, head] = system(sprintf( ...
      'git -C "%s" --no-pager rev-parse --short HEAD', root));
   fprintf(fid, '# imcrt impact report\n\n');
   fprintf(fid, 'Generated %s at commit %s by tests/verify/impactreport.m ', ...
      char(datetime('now')), strtrim(head));
   fprintf(fid, 'with N = %g packets per run and M = %d runs.\n\n', N, M);
   fprintf(fid, '%s\n\n', srctext(srcsame, srcclean));
   % Table rows come from impactquantities, so its row count sizes every
   % value array here.
   nq = numel(impactquantities());

   % Section 1: each case through each kernel version.
   fprintf(fid, '## 1. Impact of each fix\n\n');
   fprintf(fid, ['Versions in fix order: pre-fix commit %s, then tags %s. ' ...
      'Every version runs the same seed, so each table is one paired ' ...
      'realization. A fix that changes trajectories or random-number ' ...
      'use still leaves sampling error in that realization. A change ' ...
      'is attributed to a fix when it exceeds 1e-6 ' ...
      'relative. That separates a real change from summation-order ' ...
      'rounding. r50 is quantized to one radial bin. The paper_radial ' ...
      'case has the paper albedo and asymmetry (0.999, 0.9) at optical ' ...
      'depth 10. Its lengths are scaled to the kernel''s 2 cm radial ' ...
      'grid. The paper''s rod geometry is not part of src/mcrt.m.\n\n'], ...
      refs{1, 1}, strjoin(refs(2:end, 1)', ', '));
   % A helper that a ref predates (its code was inline then) is skipped by
   % extractfiles, so one list serves every ref. Historical grid builders
   % call the third-party derivative, so it joins the list here.
   kfiles = [kernelfiles(), {'src/derivative/derivative.m'}];
   cases = [mcrtcases(), struct('name', 'paper_radial', 'ka', 0.01, ...
      'ks', 9.99, 'g', 0.9, 'Z', 1, 'dz', 0.05, 'N', N, 'seed', 45)];
   for k = 1:numel(cases)
      c = cases(k);
      values = zeros(nq, size(refs, 1));
      for v = 1:size(refs, 1)
         kernel = extractfiles(root, scratch, refs{v, 1}, refs{v, 1}, ...
            kfiles);
         RT = runwith(kernel, c, N);
         [names, values(:, v), notes] = impactquantities(RT);
      end
      fixes = impactattribution(values, labels);
      fprintf(fid, '### %s (ka %g, ks %g, g %g, Z %g, dz %g, seed %d)\n\n', ...
         c.name, c.ka, c.ks, c.g, c.Z, c.dz, c.seed);
      fprintf(fid, ['| quantity | pre-fix | post-fix | change | fixes ' ...
         '| note |\n']);
      fprintf(fid, '|---|---|---|---|---|---|\n');
      for n = 1:numel(names)
         fprintf(fid, '| %s | %.6g | %.6g | %s | %s | %s |\n', names{n}, ...
            values(n, 1), values(n, end), relchange(values(n, 1), ...
            values(n, end)), fixes{n}, notes{n});
      end
      fprintf(fid, '\n');
   end

   % Section 2: the 2021 script against the fixed kernel.
   fprintf(fid, '## 2. Paper-kernel retrospective\n\n');
   ref = vdhtable35();
   % The frozen script and the two helpers it calls come from git at the
   % pre-fix commit, as does the derivative that the 2021 grid builder
   % calls. That pins the historical side to committed content instead
   % of worktree files.
   c2021 = 'examples/cooper_etal_2021/c_model/';
   frozen = extractfiles(root, scratch, refs{1, 1}, 'c2021', ...
      {[c2021 'mcrt_verify.m'], [c2021 'func/main/chgdir.m'], ...
      [c2021 'func/util/mcrt_build_grid.m']});
   pinned = extractfiles(root, scratch, refs{1, 1}, refs{1, 1}, kfiles);
   script = scratch2021(frozen, scratch, 'script');
   scriptp = scratch2021(frozen, scratch, 'production');
   seeds = 100 + (1:M);
   % The production-stream copy gets its own seeds, disjoint from the
   % first copy's for any M: with shared seeds its runs would still share
   % draws such as the azimuth and be correlated, which the unpaired
   % standard error does not carry.
   seedsx = seeds(end) + (1:M);
   Q21 = zeros(nq, M);
   Qnow = zeros(nq, M);
   Q21p = zeros(nq, M);
   Q21x = zeros(nq, M);
   T21 = cell(1, M);
   T21p = cell(1, M);
   for m = 1:M
      [Q21(:, m), T21{m}] = run2021(frozen, script, pinned, N, ...
         seeds(m), 1e-4);
      [Q21p(:, m), T21p{m}] = run2021(frozen, script, pinned, N, ...
         seeds(m), 1e-5);
      Q21x(:, m) = run2021(frozen, scriptp, pinned, N, seedsx(m), 1e-5);
      rng(seeds(m), 'twister');
      RT = mcrt(ref.ka, ref.ks, ref.g, ref.Z, ref.dz, N);
      [names, Qnow(:, m)] = impactquantities(RT);
   end
   fprintf(fid, ['The retrospective runs the van de Hulst case M = %d ' ...
      'times at N = %g, seeds %d to %d. A scratch copy of the frozen ' ...
      '2021 verification script, extracted with its helpers from commit ' ...
      '%s, and the fixed kernel each run every seed. The copy carries ' ...
      'the production kernel''s direct clause ' ...
      '(ns==0 || ia==0) in place of the script''s grazing clause. It ' ...
      'therefore scores exits the way the published runs did. It keeps ' ...
      'the out-of-place chgdir (no defect B) and the shifted bin ' ...
      'measures (defect S). It also keeps the old roulette (defect N) ' ...
      'at the script''s wmin = 1e-4, which the fixed kernel shares. ' ...
      'The production kernel used wmin = 1e-5. The text under the ' ...
      'table gives the measured effect of that threshold and of the ' ...
      'production kernel''s random-number use. se is the ' ...
      'run spread over sqrt(M). The two kernels share seeds and the same ' ...
      'random-number order, so their runs are paired. diff/se is the ' ...
      'mean paired difference (fixed minus 2021) in units of its own ' ...
      'standard error; identical runs give none.\n\n'], M, N, ...
      seeds(1), seeds(end), refs{1, 1});
   fprintf(fid, ['| quantity | 2021 mean | 2021 se | fixed mean ' ...
      '| fixed se | paired diff | diff/se |\n']);
   fprintf(fid, '|---|---|---|---|---|---|---|\n');
   for n = 1:numel(names)
      m21 = mean(Q21(n, :));
      s21 = std(Q21(n, :))/sqrt(M);
      mnow = mean(Qnow(n, :));
      snow = std(Qnow(n, :))/sqrt(M);
      d = Qnow(n, :) - Q21(n, :);
      fprintf(fid, '| %s | %.6g | %.2g | %.6g | %.2g | %.3g | %s |\n', ...
         names{n}, m21, s21, mnow, snow, mean(d), pairedratio(d));
   end
   fprintf(fid, '\nReference values: Rd %.5f, Tt %.5f, Tdr %.5f.\n\n', ...
      ref.Rd, ref.Tt, ref.Tdr);
   % Rows 1 to 5 are the hemispherical sums. The two thresholds share
   % seeds, so their runs are paired and each row gets pairedratio's
   % text, which prints none for an all-zero shift instead of 0/0.
   % agreetext compares the complete tallies, not the summaries: a
   % roulette draw at 1e-4 alone would have changed some tally.
   fprintf(fid, ['The copy''s complete tallies at wmin = 1e-5, the ' ...
      'production threshold, and at wmin = 1e-4 %s. A packet that ' ...
      'reaches wmin = 1e-4 draws a roulette number at that threshold ' ...
      'only, so identical tallies mean that no packet reached it. ' ...
      'Roulette is unbiased, so the threshold changes the variance ' ...
      'only. The hemispherical shifts in standard errors of the paired ' ...
      'difference:\n\n'], agreetext([T21{:}], [T21p{:}]));
   for n = 1:5
      fprintf(fid, '- %s: %s\n', names{n}, ...
         pairedratio(Q21p(n, :) - Q21(n, :)));
   end
   % The production-stream copy runs on its own seeds, so its rows are
   % independent of the first copy's and unpairedratio applies. It is
   % compared with the copy at the same threshold, 1e-5, so the listed
   % shifts carry the stream change alone.
   fprintf(fid, ['\nA second copy also takes the production kernel''s ' ...
      'random-number use: the path length from log(1-rand) and the ' ...
      'Henyey-Greenstein terms hg4 = 1-g, hg5 = 2*g. That relabels the ' ...
      'random stream and leaves every sampled distribution unchanged. ' ...
      'This copy runs seeds %d to %d, independent of the first copy''s. ' ...
      'Its hemispherical means at wmin = 1e-5, in standard errors of ' ...
      'the unpaired difference from the first copy at wmin = 1e-5, ' ...
      'shift by:\n\n'], seedsx(1), seedsx(end));
   for n = 1:5
      fprintf(fid, '- %s: %s\n', names{n}, ...
         unpairedratio(Q21x(n, :), Q21p(n, :)));
   end
   fprintf(fid, '\n');
   % Rows 1 to 5 are the hemispherical sums. State from the data whether
   % the two kernels agree run for run rather than asserting it in prose.
   hemitext = agreetext(Qnow(1:5, :), Q21(1:5, :));
   fprintf(fid, ['### How the published numbers relate to the corrected ' ...
      'model\n\n' ...
      'The 2021 kernel used the out-of-place chgdir, so defect B never ' ...
      'touched it. Its radial spread (r50) matches the fixed kernel. Its ' ...
      'production loop counted only unscattered packets as direct: ia is ' ...
      '0 only for an exit exactly along the normal. Defect A therefore ' ...
      'applies only to its verification script, and section 1 measures ' ...
      'that clause as the A column. Its grid builder is the shipped one, ' ...
      'so defect S applies. In any 2021 output that reports them, the ' ...
      'first angular ' ...
      'and radial bins per steradian or per area are 15.6%% below the ' ...
      'corrected value. The corrected value is 18.5%% higher. The ' ...
      'second bins are 3.3%% above it. Defect S changes the per-bin ' ...
      'densities only, not the hemispherical sums. Its ' ...
      'fluence density adds the per-depth direct absorption to every ' ...
      'radial column (defect C). The fix confines that pencil to bin 1 ' ...
      'as a volume density. Defect N drops a boosted packet that still ' ...
      'sits below wmin. That needs an albedo below 1/wrr = 0.1, so it ' ...
      'never fires at the paper albedos. The ' ...
      'hemispherical reflectance, transmittance, and diffuse absorption ' ...
      'of the fixed kernel and the paired copy %s, as the paired ' ...
      'differences ' ...
      'show. The production-stream copy differs from that copy by the ' ...
      'listed standard-error ratios; a relabeled random stream predicts ' ...
      'ratios of order one.\n'], hemitext);
   fclose(fid);
   clear closer
end

function RT = runwith(kernel, c, N)
   % Run one case with the kernel folder at the front of the path, then
   % take it off again, also on error or interrupt, so no historical
   % kernel stays active. The clear drops the cached functions so the next
   % version's files are the ones that run.
   addpath(kernel);
   restore = onCleanup(@() dropkernel(kernel));
   clear mcrt buildgrid derivative chgdir hgcos roulette binindex
   % The extracted kernel must be the one that runs, not the worktree's.
   assert(startsWith(which('mcrt'), kernel), 'impactreport:path', ...
      'mcrt resolves outside %s', kernel);
   rng(c.seed, 'twister');
   RT = mcrt(c.ka, c.ks, c.g, c.Z, c.dz, N);
   clear restore
end

function dropkernel(folder)
   % Take a kernel or helper folder off the path and drop the cached
   % functions it provided.
   rmpath(folder);
   clear mcrt buildgrid derivative chgdir hgcos roulette binindex ...
      mcrt_build_grid
end

function script = scratch2021(frozen, scratch, stream)
   % Scratch copy of the frozen 2021 verification script with six edits,
   % and two more when stream is 'production'.
   % The geometry load, which needs a data file the case does not use, is
   % removed. N and wmin come from the workspace instead of the script, so
   % the copy can run at the production threshold as well. roundn,
   % which needs the Mapping Toolbox, becomes round. The text stops before
   % the validation section, which adds the direct beam into the first
   % angular bin for its plots and would corrupt Rdf_a(1) and Tdf_a(1).
   % The two grazing direct clauses (defect A, verification only) become
   % the production kernel's clause, ns==0 || ia==0, written with the
   % script's angular index iu. That makes the copy score exits the way
   % c_model/mcrt.m did for the published runs. The production kernel
   % scores transmittance only; the reflection side gets the same clause.
   % With stream 'production' the path length draws log(1-rand) and the
   % Henyey-Greenstein terms are hg4 = 1-g, hg5 = 2*g, as in
   % c_model/mcrt.m; both use the random stream the way the published
   % runs did, which changes packet paths but no sampled distribution.
   text = fileread(fullfile(frozen, 'mcrt_verify.m'));
   text = regexprep(text, '\nload\(\[opts\.path\.data[^\n]*\n', '\n');
   text = strrep(text, 'N   = 5e5;', 'N   = Nrun;');
   assert(contains(text, 'wmin    = 1e-4;'), 'impactreport:script', ...
      'wmin line missing');
   text = strrep(text, 'wmin    = 1e-4;', 'wmin    = Wmin;');
   text = regexprep(text, 'roundn\(([^,]+),0\)', 'round($1)');
   cut = strfind(text, '%% validate by comparison');
   assert(~isempty(cut), 'impactreport:script', 'validation marker missing');
   text = text(1:cut(1)-1);
   grazing = 'ns==0 \|\| -?uz < du/2';
   nrep = numel(regexp(text, grazing, 'start'));
   assert(nrep == 2, 'impactreport:script', ...
      'expected two grazing clauses, found %d', nrep);
   text = regexprep(text, grazing, 'ns==0 || iu == 0');
   script = fullfile(scratch, 'verify2021.m');
   if strcmp(stream, 'production')
      lines = {'l = -c*log(rand);', 'hg4     = 1+g;', 'hg5     = -2*g;'};
      for n = 1:numel(lines)
         assert(contains(text, lines{n}), 'impactreport:script', ...
            'line missing: %s', lines{n});
      end
      text = strrep(text, lines{1}, 'l = -c*log(1-rand);');
      text = strrep(text, lines{2}, 'hg4     = 1-g;');
      text = strrep(text, lines{3}, 'hg5     = 2*g;');
      script = fullfile(scratch, 'verify2021p.m');
   end
   fid = fopen(script, 'w');
   fprintf(fid, '%s', text);
   fclose(fid);
end

function [q, tallies] = run2021(frozen, script, pinned, Nrun, seed, Wmin)
   % Run the scratch 2021 script at Nrun packets and roulette threshold
   % Wmin, which the script reads from this workspace. The frozen
   % helpers (chgdir, mcrt_build_grid) go first on the path with the
   % pinned derivative behind them. The quantities of impactquantities
   % and one column holding every tally the script scores come from the
   % workspace the script leaves behind.
   addpath(pinned);
   restorepin = onCleanup(@() dropkernel(pinned));
   addpath(frozen);
   restorefrozen = onCleanup(@() dropkernel(frozen));
   clear chgdir mcrt_build_grid derivative
   rng(seed, 'twister');
   evalc(['run(''' script ''')']);
   clear restorefrozen restorepin
   RT = struct('Rdf', Rdf, 'Tdf', Tdf, 'Tdr', Tdr, 'Rdr', Rdr, 'Adf', Adf, ...
      'Rdf_a', Rdf_a, 'Tdf_a', Tdf_a, 'Rdf_r', Rdf_r, 'Tdf_r', Tdf_r, ...
      'phi_rz', phi_rz, 'phi_z', phi_z, 'Adr_z', Adr_z, ...
      'Rdf_ra', Rdf_ra, 'Tdf_ra', Tdf_ra, 'Adf_rz', Adf_rz);
   RT.grid = struct('ri', ri, 'dr', dr, 'dz', dz);
   [~, q] = impactquantities(RT);
   % The tallies come out of the struct: the script's variables are not
   % visible to the analyzer, which reads indexing on them as a call.
   tallies = [RT.Rdf_ra(:); RT.Tdf_ra(:); RT.Adf_rz(:); RT.Adr_z(:); ...
      RT.Rdr; RT.Tdr];
   fprintf('2021 script: seed %d, N %g, wmin %g\n', seed, Nrun, Wmin);
end

function closeifopen(fid)
   % Close the report when an error skipped the normal fclose. After a
   % clean run the identifier is already closed and fclose raises; that
   % error is the expected case and is ignored.
   try
      fclose(fid);
   catch
   end
end
