function s = srctext(tagsame, clean)
   % Report sentence for the source check: which kernel the post-fix
   % column and the retrospective ran. Both srcmatches flags must hold for
   % the shipped-kernel statement; otherwise the worktree kernel ran.
   if tagsame && clean
      s = ['src/ is unchanged since tag fix-S-bin-measures, so the ' ...
         'post-fix column is the shipped kernel.'];
   else
      s = ['src/ differs from tag fix-S-bin-measures or has ' ...
         'uncommitted edits, so the retrospective runs the worktree ' ...
         'kernel.'];
   end
end
