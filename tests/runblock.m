function w = runblock(block, v)
   % Evaluate one kernel block from src/mcrt.m with the fields of v as its
   % workspace and return every variable afterward as a struct. The block's
   % break lines are removed because eval parses the text outside the
   % kernel's loop; a block without break is unaffected.
   names = fieldnames(v);
   for n = 1:numel(names)
      eval([names{n} ' = v.' names{n} ';']);
   end
   eval(regexprep(block, '\n[ \t]*break\>[^\n]*', ''));
   vars = who;
   w = struct();
   for n = 1:numel(vars)
      w.(vars{n}) = eval(vars{n});
   end
end
