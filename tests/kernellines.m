function code = kernellines(startMarker, endMarker)
   % Return the lines of src/mcrt.m from the line containing startMarker
   % (inclusive) up to the line containing endMarker (exclusive). An empty
   % endMarker stops at the next blank line. Tests evaluate the returned
   % text so the kernel's inline blocks are covered without a function
   % call in the hot loop.
   root = fileparts(fileparts(mfilename('fullpath')));
   lines = regexp(fileread(fullfile(root, 'src', 'mcrt.m')), '\r?\n', ...
      'split');
   first = find(contains(lines, startMarker), 1);
   if isempty(endMarker)
      last = first + find(strlength(strtrim(lines(first:end))) == 0, 1) - 2;
   else
      last = find(contains(lines, endMarker), 1) - 1;
   end
   code = strjoin(lines(first:last), newline);
end
