%==========================================================================
clean
%==========================================================================

save_data   = true;
july20      = false;
july21      = true;

p.data      = 'GREENLAND/field/2018/a_submitted/monte_carlo/';
p.save      = 'GREENLAND/field/2018/a_submitted/monte_carlo/';
%==========================================================================
%% set paths
%==========================================================================
if july20 == true
    p.data  = [p.data '20july/b_input/'];
    p.save  = [p.save '20july/b_input/'];
elseif july21 == true
    p.data  = [p.data '21july/b_input/'];
    p.save  = [p.save '21july/b_input/'];
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



