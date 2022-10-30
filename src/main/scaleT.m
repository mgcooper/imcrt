%--------------------------------------------------------------------------
function [Tdf_ra,Tdf_r,Tdf_a,Tdf,Tdr,Tt] = scaleT(Tdf_ra,Tdr,geom,opts)
%--------------------------------------------------------------------------
% scaleT scale transmittance (T)
% Matt Cooper, guycooper@ucla.edu, Dec 2020

% Inputs:
%   R       = radial element thickness                          [cm]
%   A       = angular element thickness                         [rad]
%   Z       = vertical element thickness                        [cm]
%   dr      = radial grid spacing                               [cm]
%   da      = angular grid spacing                              [rad]
%   dz      = vertical grid spacing                             [cm]
%   Tdf_ra  = diffuse photon density in spherical coord's       [-]
%   Tdr     = unscattered photon density along the z-axis       [-]

% Outputs
%   Tdf_ra  = transmitted diffuse radiance                      [W/m2/sr]
%   Tdf_r   = transmitted diffuse irradiance                    [W/m2]
%   Tdf_a   = transmitted diffuse radiant intensity             [W/sr]
%   Tdf     = transmitted diffuse radiant flux (power)          [W]
%   Tdr     = transmitted direct radiant flux (power)           [W]
%   Tt      = transmitted direct+diffuse radiant flux (power)   [W]

% number of photons
%--------------------------------------------------------------------------
    N       = opts.N;

% build a grid to calculate observable quantities (eq. 4.1/4.2 Wang)
%--------------------------------------------------------------------------
    [ri,ai,~,dr,da,~] = mcrt_build_grid(geom);
    
% differential surface area/steradians of annular rings, projection factor
%--------------------------------------------------------------------------
    dsr     = 2*pi.*sin(ai).*da;
    dA      = 2*pi*dr.*ri;                  % 2*pi*r for the entire ring
    cosa    = cos(ai);                      % projection factor

% sum the 2-d arrays into 1-d and 0-d arrays
%--------------------------------------------------------------------------
    Tdf_r   = sum(Tdf_ra,1);                        % Eq. 4.5
    Tdf_a   = sum(Tdf_ra,2);                        % Eq. 4.6
    Tdf     = sum(Tdf_r);                           % Eq. 4.8

% convert raw photon density into W, W/m2, W/sr, and W/m2/sr
%--------------------------------------------------------------------------
    Tdf_ra  = Tdf_ra./(dA.*dsr.*cosa.*N);           % Eq. 4.10
    Tdf_r   = Tdf_r./(dA.*N);                       % Eq. 4.14
    Tdf_a   = Tdf_a./(dsr.*N);                      % Eq. 4.16
    Tdr     = Tdr/N;
    Tdf     = Tdf/N;                                % Eq. 4.18
    Tt      = Tdf+Tdr;
end

