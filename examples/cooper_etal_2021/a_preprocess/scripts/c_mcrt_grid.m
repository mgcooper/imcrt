%==========================================================================
clean
%==========================================================================

% this is an older implementation that saved the output but this is no
% longer used, but is kept here as reference

save_data   = true;
july20      = false;
july21      = true;

p.data      = 'path/to/b_input/';
p.save      = 'path/to/b_input/';
%==========================================================================
%% set paths
%==========================================================================
if july20 == true
    p.data  = [p.data '20july/'];
    p.save  = [p.save '20july/'];
elseif july21 == true
    p.data  = [p.data '21july/'];
    p.save  = [p.save '21july/'];
end
p           = setpath(p);
%==========================================================================
%% build the scoring grid
%==========================================================================

% R           = 200;        % radius of detection [cm]
% A           = 1;          % cosine angular detection radius
% Z           = 100;        % thickness of medium
% nr          = 50;         % number of radial bins
% na          = 20;         % number of angular bins

% radial
nr          = R/dr;
% angular
na          = A/da;     % angular bin width
% vertical
nz          = Z/dz;     % number of vertical layers



