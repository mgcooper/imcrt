function Setup()
% Setup.m setup the model paths etc.

% Resolve paths from this file so Setup works from any working directory.
thispath = fileparts(mfilename('fullpath'));

% Add explicit dirs, not genpath over the root: the git-excluded sandbox
% archive holds a script named mcrt.m that would shadow src/mcrt.m (G).
addpath(thispath, fullfile(thispath, 'src'), ...
   fullfile(thispath, 'src', 'derivative'));

% remove paths containing example code
% Only dirs a prior session put on the path: rmpath warns on absent dirs,
% and Octave's warning has no id, so warning('off', id) cannot silence it.
exdirs = strsplit(genpath(fullfile(thispath, 'examples')), pathsep);
exdirs = exdirs(ismember(exdirs, strsplit(path, pathsep)));
if ~isempty(exdirs)
   rmpath(exdirs{:});
end

%try
%   rmpath(genpath(fullfile(thispath,'examples')));
%catch
%end

% display install message
fprintf('\n * ice Monte Carlo Radiative Transfer activated *\n\n')

end
