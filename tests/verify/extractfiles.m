function folder = extractfiles(root, scratch, ref, name, files)
   % Folder scratch/name holding the named repository files as of ref,
   % by base name, extracted once per name with git show. Historical
   % kernels come out this way so a historical buildgrid never resolves
   % the worktree derivative, and the frozen 2021 script and its helpers
   % come out this way so an uncommitted edit under the frozen example
   % cannot pass as the 2021 kernel. --no-pager keeps git from waiting on
   % the pseudo-terminal that MATLAB's system attaches.
   % A finished extraction leaves a marker holding its ref, written last,
   % so a folder from a failed run is never reused and a folder reused
   % under another ref is refused instead of handing over the wrong code.
   folder = fullfile(scratch, name);
   marker = fullfile(folder, 'extracted.ref');
   if exist(marker, 'file')
      assert(strcmp(strtrim(fileread(marker)), ref), 'extractfiles:cache', ...
         '%s holds another ref, not %s', folder, ref);
      return
   end
   if exist(folder, 'dir')
      rmdir(folder, 's');
   end
   % The ref itself must resolve first: an unknown tag would make every
   % path look absent, and a folder without a kernel would let the
   % worktree kernel run under the ref's name.
   % The revision is quoted: on Windows, system runs through cmd.exe,
   % where a bare caret is an escape character.
   [status, out] = system(sprintf( ...
      'git -C "%s" --no-pager rev-parse --verify --quiet "%s^{commit}"', ...
      root, ref));
   assert(status == 0, 'extractfiles:ref', 'ref %s does not resolve: %s', ...
      ref, out);
   mkdir(folder);
   for n = 1:numel(files)
      % A file absent at a valid ref is not an error: a helper the kernel
      % calls today was inline code at an older ref. The kernel itself is
      % never optional.
      % ref:path is quoted for the same cmd.exe caret reason as above.
      status = system(sprintf( ...
         'git -C "%s" --no-pager cat-file -e "%s:%s"', root, ref, files{n}));
      if status ~= 0
         assert(~endsWith(files{n}, 'mcrt.m'), 'extractfiles:git', ...
            '%s is absent at %s', files{n}, ref);
         continue
      end
      [~, base, ext] = fileparts(files{n});
      [status, out] = system(sprintf( ...
         'git -C "%s" --no-pager show "%s:%s" > "%s"', ...
         root, ref, files{n}, fullfile(folder, [base ext])));
      assert(status == 0, 'extractfiles:git', 'git show %s:%s failed: %s', ...
         ref, files{n}, out);
   end
   fid = fopen(marker, 'w');
   fprintf(fid, '%s\n', ref);
   fclose(fid);
end
