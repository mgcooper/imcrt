function [ri,ai,zi,dr,da,dz] = buildgrid(R,A,Z,dr,da,dz)
%BUILDGRID build grids for scoring observable quantities

% note that U/du are in radians here i.e. theta, but for scoring within the
% mcrt program they need to be in u = cos(theta)

% example:
% R           = 200;        % radius of detection [cm]
% A           = pi/2;       % angular detection radius
% Z           = 100;        % thickness of medium
% nr          = 50;         % number of radial bins
% na          = 25;         % number of angular bins

% number of bins in each dimension
%     nr      = roundn(R/dr,0);     % radial
%     na      = roundn(A/da,0);     % angular
%     nz      = roundn(Z/dz,0);     % number

% grid center coordinates
    ri      = dr/2:dr:R-dr/2;
    zi      = dz/2:dz:Z-dz/2;
    ai      = da/2:da:A-da/2;

% r and z need an extra overflow coordinat
	zi(end+1)   = zi(end)+dz;
   ri(end+1)   = ri(end)+dr;

% optimized to minimize error (https://core.ac.uk/download/pdf/84314628.pdf)
    ri      = ri + (dr*dr/12)./ri;                    % Eq. 8
    ai      = ai + cot(ai).*(1-da/2*cot(da/2));       % Eq. 14
% the optimal z-coordinate is the center of each element i.e. zi

% make ai and zi columns
    ai      = ai(:);
    zi      = zi(:);

% get new da/dr/dz ('derivative' can be downloaded on the file exchange)
    da      = derivative(ai);
    dr      = derivative(ri);
    dz      = derivative(zi);

% see the difference
%     figure;
%     tiledlayout(1,2); nexttile;
%     plot(dr/2:dr:R-dr/2,'o'); hold on; ylabel('radial grid')
%     plot(ri,'o'); legend('grid centers','optimized');
%     nexttile;
%     plot((dr/2:dr:R-dr/2)-ri,'o'); legend('grid centers - optimized');
%
%     figure;
%     tiledlayout(1,2); nexttile;
%     plot(du/2:du:U-du/2,'o'); hold on; ylabel('Angular grid')
%     plot(ui,'o'); legend('grid centers','optimized');
%     nexttile;
%     plot((du/2:du:U-du/2)-ui,'o'); legend('grid centers - optimized');

