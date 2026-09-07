function out = mcrtgolden(mode, file)
   % Golden digest of every mcrt output for the seeded cases in mcrtcases.
   % out = mcrtgolden() returns the digest as a cell array of text lines.
   % out = mcrtgolden('path') returns the tracked baseline path.
   % mcrtgolden('write') overwrites that baseline. It is the one re-baseline
   % command. Run it only in the commit that carries a physics fix, and
   % quantify the delta in that commit message.
   % mcrtgolden('write', file) writes to file instead; the tests use this.
   % Scalars are stored exactly. Each array is stored as size, sum, sum of
   % squares, first and last element, and a polynomial hash of its bytes.
   % Every number is printed with %.17g so the text round-trips exactly.
   baseline = fullfile(fileparts(mfilename('fullpath')), 'golden', ...
      'mcrt_golden.txt');
   if nargin > 0 && strcmp(mode, 'path')
      out = baseline;
      return
   end
   cases = mcrtcases();
   parts = cell(1, numel(cases));
   for k = 1:numel(cases)
      c = cases(k);
      rng(c.seed, 'twister');
      RT = mcrt(c.ka, c.ks, c.g, c.Z, c.dz, c.N);
      parts{k} = digest(c.name, '', RT);
   end
   out = [parts{:}];
   if nargin > 0 && strcmp(mode, 'write')
      if nargin < 2
         file = baseline;
      end
      fid = fopen(file, 'w');
      fprintf(fid, '%s\n', out{:});
      fclose(fid);
   end
end

function lines = digest(name, prefix, s)
   % One line per scalar and six per array. Nested structs recurse with a
   % dotted prefix so grid fields read as grid.ri.size and so on.
   names = fieldnames(s);
   parts = cell(1, numel(names));
   for i = 1:numel(names)
      v = s.(names{i});
      key = [prefix names{i}];
      if isstruct(v)
         parts{i} = digest(name, [key '.'], v);
      elseif isscalar(v)
         parts{i} = {sprintf('%s %s %.17g', name, key, v)};
      else
         parts{i} = { ...
            sprintf('%s %s.size %s', name, key, mat2str(size(v))), ...
            sprintf('%s %s.sum %.17g', name, key, sum(v(:))), ...
            sprintf('%s %s.sumsq %.17g', name, key, sum(v(:).^2)), ...
            sprintf('%s %s.first %.17g', name, key, v(1)), ...
            sprintf('%s %s.last %.17g', name, key, v(end)), ...
            sprintf('%s %s.checksum %.17g', name, key, bytehash(v))};
      end
   end
   lines = [parts{:}];
end

function h = bytehash(v)
   % Polynomial hash of the IEEE-754 bytes modulo a prime, kept exact in
   % double: h = sum(b(i) * r^(n-i)) mod p. Any changed bit or element order
   % changes h unless the difference polynomial vanishes modulo p, a
   % structured coincidence with probability near 1/p.
   p = 33554393;                   % largest prime below 2^25, so p^2 < 2^53
   r = 257;                        % base above the byte range
   B = 4096;                       % block length: 4096 * 255 * p < 2^53
   w = ones(B, 1);                 % w(k) = r^(B-k) mod p, built by repeated
   for k = B-1:-1:1                % modular multiplication
      w(k) = mod(w(k+1) * r, p);
   end
   rB = mod(w(1) * r, p);          % r^B mod p shifts the running hash
   b = double(typecast(v(:), 'uint8'));
   nb = ceil(numel(b) / B);
   b = reshape([b; zeros(nb*B - numel(b), 1)], B, nb);
   h = 0;
   for k = 1:nb
      h = mod(h * rB + mod(sum(b(:, k) .* w), p), p);
   end
end
