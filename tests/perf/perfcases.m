function cases = perfcases()
   % Timing cases for the perf harness: two of the shared cases of
   % mcrtcases at an N that makes one run take about a tenth of a second,
   % long enough for timeit to average out call overhead and short enough
   % for the fast suite. vdh_reflectance is the few-step slab; mini_fluence
   % has albedo 0.99, so its packets take many more steps and its time is
   % the per-step cost of the loop.
   cases = mcrtcases();
   cases = cases(strcmp({cases.name}, 'vdh_reflectance') ...
      | strcmp({cases.name}, 'mini_fluence'));
   [cases.N] = deal(2e5, 2e4);
end
