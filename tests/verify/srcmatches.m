function [tagsame, clean] = srcmatches(root, ref)
   % Compare the kernel that MATLAB will run with a git ref. tagsame is
   % true when the src tree at HEAD equals the src tree at ref (rev-parse
   % tree hashes). Every git call runs with --no-pager: MATLAB's system
   % attaches a pseudo-terminal, so a paging command waits for a key.
   % clean is true when the worktree copies of src/mcrt.m,
   % src/buildgrid.m, and src/derivative/derivative.m equal their HEAD
   % content, so uncommitted edits to the kernel are reported instead of
   % silently run.
   [status, treenow] = system(sprintf( ...
      'git -C "%s" --no-pager rev-parse HEAD:src', root));
   assert(status == 0, 'srcmatches:git', 'git rev-parse failed: %s', treenow);
   [status, treeref] = system(sprintf( ...
      'git -C "%s" --no-pager rev-parse %s:src', root, ref));
   assert(status == 0, 'srcmatches:git', 'git rev-parse failed: %s', treeref);
   tagsame = strcmp(strtrim(treenow), strtrim(treeref));
   files = {'mcrt.m', 'buildgrid.m', 'derivative/derivative.m'};
   clean = true;
   for n = 1:numel(files)
      [status, headtext] = system(sprintf( ...
         'git -C "%s" --no-pager show HEAD:src/%s', root, files{n}));
      assert(status == 0, 'srcmatches:git', 'git show failed: %s', headtext);
      worktext = fileread(fullfile(root, 'src', files{n}));
      % Compare with line endings normalized: a checkout under text=auto
      % may hold CRLF while git show returns the LF blob.
      clean = clean && strcmp(regexprep(headtext, '\r\n', '\n'), ...
         regexprep(worktext, '\r\n', '\n'));
   end
end
