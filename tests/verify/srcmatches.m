function [tagsame, clean] = srcmatches(root, ref)
   % Compare the kernel that MATLAB will run with a git ref. tagsame is
   % true when the src tree at HEAD equals the src tree at ref (rev-parse
   % tree hashes). Every git call runs with --no-pager: MATLAB's system
   % attaches a pseudo-terminal, so a paging command waits for a key.
   % clean is true when the worktree copy of every file in kernelfiles
   % equals its HEAD content, so uncommitted edits to the kernel or its
   % helpers are reported instead of silently run. A file that HEAD does
   % not hold is an uncommitted helper, so it counts as dirty.
   [status, treenow] = system(sprintf( ...
      'git -C "%s" --no-pager rev-parse HEAD:src', root));
   assert(status == 0, 'srcmatches:git', 'git rev-parse failed: %s', treenow);
   [status, treeref] = system(sprintf( ...
      'git -C "%s" --no-pager rev-parse %s:src', root, ref));
   assert(status == 0, 'srcmatches:git', 'git rev-parse failed: %s', treeref);
   tagsame = strcmp(strtrim(treenow), strtrim(treeref));
   files = kernelfiles();
   clean = true;
   for n = 1:numel(files)
      [status, headtext] = system(sprintf( ...
         'git -C "%s" --no-pager show HEAD:%s', root, files{n}));
      if status ~= 0
         clean = false;
         continue
      end
      worktext = fileread(fullfile(root, files{n}));
      % Compare with line endings normalized: a checkout under text=auto
      % may hold CRLF while git show returns the LF blob.
      clean = clean && strcmp(regexprep(headtext, '\r\n', '\n'), ...
         regexprep(worktext, '\r\n', '\n'));
   end
end
