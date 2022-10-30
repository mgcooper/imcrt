%--------------------------------------------------------------------------
function [Adf_rz,Adf,phi_rz,phi_z] = scaleA(Adf_rz,Adr_z,kabs,geom,opts)
%--------------------------------------------------------------------------
% scaleA scale absorption (A), also compute fluence
% Matt Cooper, guycooper@ucla.edu, Dec 2020

% Inputs:
%   R       = radial element thickness                              [cm]
%   A       = angular element thickness                             [rad]
%   Z       = vertical element thickness                            [cm]
%   dr      = radial grid spacing                                   [cm]
%   da      = angular grid spacing                                  [rad]
%   dz      = vertical grid spacing                                 [cm]
%   Adf_rz  = diffuse photon absorptance in cylindrical shells      [-]
%   Adr_z   = unscattered photon absorptance along the z-axis       [-]
%   kabs    = absorption coefficient                                [cm-1]
%   N       = number of photons

% Outputs:
%   Adf_rz  = diffuse absorptance (abs. probability / unit volume)  [W/m3]
%   Adf     = diffuse absorptance (total abs. probability)          [W]
%   phi_rz  = internal fluence                                      [W/m2]
%   phi_z   = fluence along the z-axis                              [W]

% number of photons
%--------------------------------------------------------------------------
    N       = opts.N;

% build a grid to calculate observable quantities (eq. 4.1/4.2 Wang)
%--------------------------------------------------------------------------
    [ri,~,~,dr,~,dz] = mcrt_build_grid(geom);
    
% differential surface area and volume of annular rings 
%--------------------------------------------------------------------------
    dA      = 2*pi*dr.*ri;          % 2*pi*r for the entire ring
    dV      = dA.*dz;               % 2*pi*r*dz for the entire ring    

% sum the 2-d arrays into 1-d and 0-d arrays
%--------------------------------------------------------------------------
    Adf_z   = sum(Adf_rz,2);            % Eq. 4.20
    Adf     = sum(Adf_z);               % Eq. 4.22
    
% for absorption, we divide the raw photon weight by the number of photons
% to get power (W) and then by the volume of each elemental control
% volume to get W/m3, and by the surface area through which the power
% passes to get W/m2

% convert raw photon density into power (W) and then into W/m3 and W/m
%--------------------------------------------------------------------------
    Adf_rz  = Adf_rz./dV./N;            % Eq. 4.23      [W/m3]
    Adr_z   = Adr_z./dz./N;             % Eq. 4.24      [W/m3] (area = 1)
    Adf_z   = Adf_z./dz./N;             % Eq. 4.25      [W/m]
    Adf     = Adf./N;                   % Eq. 4.27      [W]

% compute internal fluence (phi_rz) and fluence along the z-axis (phi_z)
%--------------------------------------------------------------------------
    phi_rz  = (Adf_rz+Adr_z)./kabs;       % Eq. 4.28      [W/m2]
    phi_z   = (Adf_z+Adr_z)./kabs;        % Eq. 4.29      [W]
  % phi_rz  = Adf_rz./ka;               % Eq. 4.28
  % phi_z   = Adf_z./ka;                % Eq. 4.29

end

