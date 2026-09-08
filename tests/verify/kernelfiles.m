function files = kernelfiles()
   % Repository-relative paths of every source file that shapes an mcrt
   % result: the kernel, the functions its photon loop calls, the grid
   % builder, and the derivative the builder calls. The driver's checkpoint
   % fingerprint, the source-clean check, and the impact report's kernel
   % extraction all read this one list, so a new helper is added here once.
   files = {'src/mcrt.m', 'src/chgdir.m', 'src/hgcos.m', 'src/roulette.m', ...
      'src/binindex.m', 'src/buildgrid.m', 'src/derivative/derivative.m'};
end
