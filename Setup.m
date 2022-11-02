function Setup()
% Setup.m setup the model paths etc.

% turn off complaints about paths not already being on the path
warning off

% add paths containing source code
addpath(genpath([pwd() filesep 'src']));

% remove git paths
rmpath(genpath([pwd() filesep '.git*']));

% remove paths containing example code
rmpath(genpath([pwd() filesep 'examples']));

%try
%   rmpath([pwd() filesep 'examples']);
%catch
%end

warning on 

% display install message
fprintf('\n * ice-Monte Carlo Radiative Transfer activated *\n\n')
