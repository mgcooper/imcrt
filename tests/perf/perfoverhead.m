function o = perfoverhead(ndraw)
   % Measure MATLAB function-call overhead against the photon loop's work.
   %
   %  o = perfoverhead(ndraw) times scalar loops of ndraw iterations and
   %  returns seconds per iteration for each block the kernel could write
   %  out or call, so the call decisions in src/mcrt.m can be reproduced:
   %    call       a loop that calls the one-line function passthrough,
   %               the cost of a bare call;
   %    inline     the direction-cosine update written out in the loop;
   %    chgdir     the same update as a chgdir call;
   %    clamp      the radial index with its two clamps written out;
   %    binindex   the same index as a binindex call;
   %    hginline   the Henyey-Greenstein cosine with coefficients that were
   %               precomputed outside the loop;
   %    hgcos      the same cosine as an hgcos call, which recomputes the
   %               coefficients from g on every call;
   %    rrinline   Russian roulette written out;
   %    roulette   the same roulette as a roulette call;
   %    tally      one indexed accumulate into a 2-d array, as one tally
   %               event costs;
   %    tallyss    the same plus a sum-of-squares accumulate, the extra
   %               cost of a per-bin variance tally.
   %  o.ratio is o.chgdir/o.inline. The kernel calls chgdir, hgcos, and
   %  roulette inside its loop and binindex at exits; the two per-step
   %  index clamps stay written out because two calls per step measured
   %  10 to 13% of a step.
   if nargin < 1
      ndraw = 1e5;
   end
   rng(3, 'twister');
   us = 1 - 2*rand(1, ndraw);
   ps = 2*pi*rand(1, ndraw);
   r = 0.0025*rand(1, ndraw);
   iz = randi(21, 1, ndraw);
   ir = randi(21, 1, ndraw);
   o.call = timeit(@() callloop(ndraw))/ndraw;
   o.inline = timeit(@() inlineloop(us, ps))/ndraw;
   o.chgdir = timeit(@() chgdirloop(us, ps))/ndraw;
   o.clamp = timeit(@() clamploop(r))/ndraw;
   o.binindex = timeit(@() binindexloop(r))/ndraw;
   o.hginline = timeit(@() hginlineloop(ndraw))/ndraw;
   o.hgcos = timeit(@() hgcosloop(ndraw))/ndraw;
   o.rrinline = timeit(@() rrinlineloop(ndraw))/ndraw;
   o.roulette = timeit(@() rouletteloop(ndraw))/ndraw;
   o.tally = timeit(@() tallyloop(iz, ir))/ndraw;
   o.tallyss = timeit(@() tallyssloop(iz, ir))/ndraw;
   o.ratio = o.chgdir/o.inline;
end

function x = passthrough(x)
   % The cheapest possible function: its call cost is the overhead.
end

function x = callloop(ndraw)
   % ndraw calls of passthrough.
   x = 0;
   for n = 1:ndraw
      x = passthrough(x);
   end
end

function [ux, uy, uz] = inlineloop(us, ps)
   % The direction update written out in the loop, as the kernel had it
   % before the update became a chgdir call. Same arithmetic as chgdir.
   ux = 0;
   uy = 0;
   uz = 1;
   for n = 1:numel(us)
      sinth = sqrt(1-uz*uz);
      sinths = sqrt(1-us(n)*us(n));
      cps = cos(ps(n));
      sps = sin(ps(n));
      if sinth < 1e-12
         ux = sinths*cps;
         uy = sinths*sps;
         uz = sign(uz)*us(n);
      else
         uxn = sinths/sinth*(ux*uz*cps-uy*sps)+ux*us(n);
         uyn = sinths/sinth*(uy*uz*cps+ux*sps)+uy*us(n);
         uz = -sinths*sinth*cps+uz*us(n);
         ux = uxn;
         uy = uyn;
      end
   end
end

function [ux, uy, uz] = chgdirloop(us, ps)
   % The same update through the chgdir call the kernel makes.
   ux = 0;
   uy = 0;
   uz = 1;
   for n = 1:numel(us)
      [ux, uy, uz] = chgdir(ux, uy, uz, us(n), ps(n));
   end
end

function s = clamploop(r)
   % The radial index with its two clamps written out, as the kernel's
   % per-step indices are.
   dr = 0.001;
   nr = 20;
   s = 0;
   for n = 1:numel(r)
      ir = ceil(r(n)/dr);
      if ir<1; ir = 1; end
      if ir>nr; ir = nr+1; end
      s = s + ir;
   end
end

function s = binindexloop(r)
   % The same index through binindex.
   dr = 0.001;
   nr = 20;
   s = 0;
   for n = 1:numel(r)
      s = s + binindex(r(n), dr, nr+1);
   end
end

function s = hginlineloop(ndraw)
   % The Henyey-Greenstein cosine with the five coefficients precomputed
   % once, as the kernel had it before the draw became an hgcos call.
   g = 0.75;
   hg1 = 1/(2*g);
   hg2 = (1+g^2);
   hg3 = (1-g^2);
   hg4 = 1+g;
   hg5 = -2*g;
   s = 0;
   for n = 1:ndraw
      s = s + hg1*(hg2-(hg3/(hg4+hg5*rand))^2);
   end
end

function s = hgcosloop(ndraw)
   % The same draw through hgcos, which recomputes the coefficients.
   g = 0.75;
   s = 0;
   for n = 1:ndraw
      s = s + hgcos(g);
   end
end

function s = rrinlineloop(ndraw)
   % Roulette written out, played on every iteration.
   wrr = 10;
   s = 0;
   for n = 1:ndraw
      wt = 5e-5;
      if rand < 1/wrr
         wt = wt*wrr;
      else
         wt = 0;
      end
      s = s + wt;
   end
end

function s = rouletteloop(ndraw)
   % The same roulette through the roulette call.
   wrr = 10;
   s = 0;
   for n = 1:ndraw
      s = s + roulette(5e-5, wrr);
   end
end

function A = tallyloop(iz, ir)
   % One indexed accumulate per iteration, the cost of one tally event.
   A = zeros(21, 21);
   awt = 0.1;
   for n = 1:numel(iz)
      A(iz(n), ir(n)) = A(iz(n), ir(n)) + awt;
   end
end

function A = tallyssloop(iz, ir)
   % The same accumulate plus the sum of squares a variance tally adds.
   A = zeros(21, 21);
   S = zeros(21, 21);
   awt = 0.1;
   for n = 1:numel(iz)
      A(iz(n), ir(n)) = A(iz(n), ir(n)) + awt;
      S(iz(n), ir(n)) = S(iz(n), ir(n)) + awt*awt;
   end
   A = A + S;
end
