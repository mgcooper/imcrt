function files = kernelfiles()
   % Repository-relative paths of every source file that shapes an mcrt
   % result: the kernel, the functions its photon loop calls, the compute
   % functions that normalize its tallies, and the grid builder. The
   % driver's checkpoint fingerprint, the source-clean check, and the
   % impact report's kernel extraction all read this one list, so a new
   % helper is added here once. The third-party derivative is not here:
   % buildgrid no longer calls it, and only historical kernels need it
   % (the impact report adds it for those).
   files = {'src/mcrt.m', 'src/chgdir.m', 'src/hgcos.m', 'src/roulette.m', ...
      'src/binindex.m', 'src/computeReflectance.m', ...
      'src/computeTransmittance.m', 'src/computeAbsorption.m', ...
      'src/mcstderr.m', 'src/buildgrid.m'};
end
