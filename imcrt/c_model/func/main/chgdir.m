%==========================================================================
function [uxx,uyy,uzz] = chgdir(ux,uy,uz,us,phis)
%==========================================================================
% Determines the new direction cosines (ux,uy,uz) from old direction
% cosines (ux, uy, uz), polar scattering cosine (us), and
% scattering azimuthal angle (phis). Robert A Leathers & Trijntje Downes,
% 17 March 2000. Rewritten by Matt Cooper, 6 Dec 2020. guycooper@ucla.edu.

% first two are for just for reference, don't activate them
%costheta   = uz;
%costheta_s = us
sintheta    = sqrt(1-uz*uz);    % sin(theta) - previous photon polar angle
sintheta_s  = sqrt(1-us*us);    % sin(theta_s) - new photon polar angle
cosphi_s    = cos(phis);        % cos(phi_s) - new photon azim. angle
sinphi_s    = sin(phis);        % sin(phi_s) - new photon azim. angle

% if initial direction not straight up or straight down
if sintheta < 1e-12
    s       = sign(uz);
    uxx     = sintheta_s*cosphi_s;
    uyy     = sintheta_s*sinphi_s;
    uzz     = s*us;
% initial direction straight up or down (sin(theta)=0)    
else
    uxx     = sintheta_s/sintheta*(ux*uz*cosphi_s-uy*sinphi_s) + ux*us;
    uyy     = sintheta_s/sintheta*(uy*uz*cosphi_s+ux*sinphi_s) + uy*us;
    uzz     = -sintheta_s*sintheta*cosphi_s + uz*us;
end

% mgc consistency checks, these should = 0
% uscheck1 = us - dot([uxx,uyy,uzz],[ux,uy,uz]); 
% uscheck2 = 1 - norm([ux,uy,uz]);

% mgc NOTE: input varnames must be different than output i.e. assignment in
% place is not supported. I am not sure why, but I think it relates to the
% trigonometric functions. Assignment in place works fine with the original
% Leathers/Downes version. 
% mgc updated notation to ux/uy/uz/us in place of alpha/beta/gamma/gammas
% mgc replaced matrix multiplication with explicit (faster) expressions for
% the updated ux, uy, uz and predefined cosphi_s/sinphi_s
% rewritten code is ~0.1 ms faster than Leathers/Downes version

% mu_x = x-direction cosine
% mu_y = y-direction cosine
% mu_z = z-direction cosine
% mu_s = cosine of polar scattering angle (=cos(Psi)) (or cos(theta))
% phi_s = azimuthal scattering angle (note, the angle (phi), not cos(phi))

end