%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [xx,yy,zz,mu_xx,mu_yy,mu_zz] = rodscat(x,y,z,mu_x,mu_y,mu_z,l)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% RODSCAT scatters a "photon" at position x/y/z with ray trajectory
% mu_x/mu_y/mu_z into a new trajectory with pathlength l. Scattering is
% isotropic. Note: assignment in place is not supported, meaning input
% varnames must be different than output i.e. in = mu_x, out = mu_xx
% Matt Cooper, guycooper@ucla.edu, Dec 2020

% new scattering angles   
    phi_s=2*pi*(rand);              % azimuthal symmetry
    mu_s=1-2*rand;                  % isotropic scattering
% new direction cosines       
   [mu_xx,mu_yy,mu_zz]=chgdir(mu_x,mu_y,mu_z,mu_s,phi_s);
% new component path lengths
    dx=l*mu_xx;
    dy=l*mu_yy;
    dz=l*mu_zz;
% new position    
    xx=x+dx;
    yy=y+dy;
    zz=z+dz;
end

