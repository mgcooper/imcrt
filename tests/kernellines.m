function code = kernellines(startMarker, endMarker)
   % Return the lines of src/mcrt.m from the line containing startMarker
   % (inclusive) up to the line containing endMarker (exclusive). An empty
   % endMarker stops at the next blank line or at the first line indented
   % less than the start line, whichever comes first, so a block inside
   % the photon loop ends where its enclosing loop closes. Tests evaluate
   % the returned text so the kernel's inline blocks are covered without a
   % function call in the hot loop.
   root = fileparts(fileparts(mfilename('fullpath')));
   lines = regexp(fileread(fullfile(root, 'src', 'mcrt.m')), '\r?\n', ...
      'split');
   first = find(contains(lines, startMarker), 1);
   if isempty(endMarker)
      indent = @(s) numel(s) - numel(strip(s, 'left'));
      last = first;
      while last < numel(lines) && strlength(strtrim(lines{last+1})) > 0 ...
            && indent(lines{last+1}) >= indent(lines{first})
         last = last + 1;
      end
   else
      last = find(contains(lines, endMarker), 1) - 1;
   end
   code = strjoin(lines(first:last), newline);
end
