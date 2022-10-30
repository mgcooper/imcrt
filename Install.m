function Install()
% Setup.m setup the model paths etc.

% add paths containing source code
addpath(genpath([pwd() filesep 'src']));

% remove paths containing example code 
rmpath(genpath([pwd() filesep 'examples']));

% display install message
fprintf('\n * ice-Monte Carlo Radiative Transfer activated *\n\n')