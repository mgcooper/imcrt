%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [Rdf_ra,Rdf_r,Rdf_a,Rdf,Rdr,Rt] = scaleR(Rdf_ra,Rdr,geom,opts)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% scaleT scale transmittance (T)
% Matt Cooper, guycooper@ucla.edu, Dec 2020

% Inputs:
%   R       = radial element thickness                          [cm]
%   A       = angular element thickness                         [rad]
%   Z       = vertical element thickness                        [cm]
%   dr      = radial grid spacing                               [cm]
%   da      = angular grid spacing                              [rad]
%   dz      = vertical grid spacing                             [cm]
%   Rdf_ra  = diffuse photon density in spherical coord's       [-]
%   Rdr     = unscattered photon density along the z-axis       [-]

% Outputs
%   Rdf_ra  = reflected diffuse radiance                      [W/m2/sr]
%   Rdf_r   = reflected diffuse irradiance                    [W/m2]
%   Rdf_a   = reflected diffuse radiant intensity             [W/sr]
%   Rdf     = reflected diffuse radiant flux (power)          [W]
%   Rdr     = reflected direct radiant flux (power)           [W]
%   Rt      = reflected direct+diffuse radiant flux (power)   [W]

% number of photons
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    N       = opts.N;

% build a grid to calculate observable quantities (eq. 4.1/4.2 Wang)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    [ri,ai,~,dr,da,~] = mcrt_build_grid(geom);
    
% differential surface area/steradians of annular rings, projection factor
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    dsr     = 2*pi.*sin(ai).*da;
    dA      = 2*pi*dr.*ri;                  % 2*pi*r for the entire ring
    cosa    = cos(ai);                      % projection factor

% sum the 2-d arrays into 1-d and 0-d arrays
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    Rdf_r       = sum(Rdf_ra,1);                        % Eq. 4.3
    Rdf_a       = sum(Rdf_ra,2);                        % Eq. 4.4
    Rdf         = sum(Rdf_r);                           % Eq. 4.7

% convert raw photon density into W, W/m2, W/sr, and W/m2/sr
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Rdf_ra      = Rdf_ra./(dA.*dsr.*cosa.*N);           % Eq. 4.9
    Rdf_r       = Rdf_r./(dA.*N);                       % Eq. 4.13
    Rdf_a       = Rdf_a./(dsr.*N);                      % Eq. 4.15
    Rdr         = Rdr/N;                                % + Rsp/N
    Rdf         = Rdf/N;                                % Eq. 4.17
    Rt          = Rdf+Rdr;
end

