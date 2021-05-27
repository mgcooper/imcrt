
function [mux,muy,muz,phis]  = mcrt_init_photons(opts,iops)

    N       = opts.N;
    rcr     = iops.rcr;
% ideal case - isotropic source
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if opts.ideal_case == true
    phis    = 2*pi.*rand(N,1);              % phi_s, isotropic (0->2pi)
    muz     = sqrt(rand(N,1));               % cos(theta_s), isotropic 
    mux     = sqrt(1.0-muz.^2).*cos(phis);  % mu_x    
    muy     = sqrt(1.0-muz.^2).*sin(phis);  % mu_y
    return
else
    
% ideal cosine response
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    phis    = 2*pi.*rand(N,1);      % phi_s, isotropic (0->2pi)
    val     = (0.5:1:N-0.5)/N;      % step uniformly over interval (0,1)
    muz     = sqrt(1.0-val');       % mu_z, isotropic (cosine) (0->pi/2)
    muy     = sqrt(1.0-muz.^2).*sin(phis);
    mux     = sqrt(1.0-muz.^2).*cos(phis);
% non-ideal cosine - use empirical response function
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if opts.test_rcr == true
        if length(rcr)~=N
            newfreq = length(rcr)/N;
            newang = newfreq:newfreq:length(rcr);
            rcr = interp1(1:length(rcr),rcr,newang); rcr = rcr(:);
        end
        muz = rcr;
        muy = sqrt(1.0-muz.^2).*sin(phis);
        mux = sqrt(1.0-muz.^2).*cos(phis);
    end
    
        
%     figure; 
%     scatter3(xd,yd,zd,'filled'); hold on;
%     for i = 1:length(muz)
%         xi = xd+mux(i)*rand;
%         yi = yd+muy(i)*rand;
%         zi = zd+muz(i)*rand;
%         plot3([xd,xi],[yd,yi],[zd,zi]);
% %         pause; 
%     end
    
end
