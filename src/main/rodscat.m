%--------------------------------------------------------------------------
function [x,y,z,mu_x,mu_y,mu_z] = rodscat(x,y,z,mu_x,mu_y,mu_z,l)
%--------------------------------------------------------------------------
% RODSCAT scatters a "photon" at position x/y/z with ray trajectory mu_x/y/z
% into a new trajectory with pathlength l. Scattering is isotropic. 
% Matt Cooper, guycooper@ucla.edu, Dec 2020

% new scattering angles   
    phi_s=2*pi*(rand);              % azimuthal symmetry
    mu_s=1-2*rand;                  % isotropic scattering
% new direction cosines       
   [mu_x,mu_y,mu_z]=chgdir(mu_x,mu_y,mu_z,mu_s,phi_s);
% new component path lengths
    dx=l*mu_x;
    dy=l*mu_y;
    dz=l*mu_z;
% new position    
    x=x+dx;
    y=y+dy;
    z=z+dz;

