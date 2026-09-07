function f = reportname(outdir, stamp)
   % Report path outdir/vdh_report_<stamp>.txt that never overwrites an
   % earlier report: a resumed run can finish within a second, so an
   % existing name gets a numeric suffix. stamp is the caller's timestamp
   % text, so tests can force a collision.
   f = fullfile(outdir, ['vdh_report_' stamp '.txt']);
   k = 1;
   while exist(f, 'file')
      k = k + 1;
      f = fullfile(outdir, sprintf('vdh_report_%s_%d.txt', stamp, k));
   end
end
