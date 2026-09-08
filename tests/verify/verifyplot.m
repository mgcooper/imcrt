function fig = verifyplot(RT, casename)
   %VERIFYPLOT Plot one verification run against its reference.
   %
   %  fig = verifyplot(RT, casename) draws the angular reflectance and
   %  transmittance of a run against the van de Hulst table for case
   %  'reflect', or the fluence depth profile for case 'fluence', and
   %  returns the figure. mcrt_verify shows it; vdhverify saves it to a file
   %  next to its report.

   % Plots: angular tables against the reference, or fluence against depth.
   switch casename

      case 'reflect'

         ref = vdhtable35();

         fig = figure('Units', 'in', 'Position', [3 3 12 5]);
         subplot(1, 2, 1);
         hold on
         box on
         scatter(RT.grid.ai/pi, RT.Rdf_a, 80, 'filled', 's');
         scatter(acos(ref.mu)/pi, ref.R_sr, 80, 'filled');
         xlabel('exiting angle, \alpha [\pi rad]', 'Interpreter', 'tex');
         ylabel('R_d(\alpha) [sr^{-1}]', 'Interpreter', 'tex');
         legend('iMCRT', 'van de Hulst');

         subplot(1, 2, 2);
         hold on
         box on
         scatter(RT.grid.ai/pi, RT.Tdf_a, 80, 'filled', 's');
         scatter(acos(ref.mu)/pi, ref.T_sr, 80, 'filled');
         xlabel('exiting angle, \alpha [\pi rad]', 'Interpreter', 'tex');
         ylabel('T_d(\alpha) [sr^{-1}]', 'Interpreter', 'tex');
         legend('iMCRT', 'van de Hulst');

      case 'fluence'

         % see Fig. 4 in Wang et al. 1995 for comparison
         fig = figure;
         scatter(RT.grid.zi, RT.phi_z);
         set(gca, 'YScale', 'log', 'YLim', [0.5 10], 'XLim', [0 1]);
         xlabel('z (cm)');
         ylabel('Fluence (-)');
   end
end
