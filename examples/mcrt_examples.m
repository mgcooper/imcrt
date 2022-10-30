% might replace path settings with p = mfilename and then point to relative
% dirs 

% see example at end for intersection of photon with detector rod


opts.test_refl       = true;  % verify total reflectance
opts.test_ang        = true;  % verify angular reflectance
opts.test_fluence    = false; % verify fluence (slower)
opts.save_figs       = false;

opts.path.data  = '/full/path/to/b_input/';
opts.path.save  = '/full/path/to/d_output/verify/';
opts.save_data  = false;

% load the predefined geometry
load([opts.path.data 'mcrt_geometry_input']);
opts.geom = geom;
	
%--------------------------------------------------------------------------
%% experimental setup
%--------------------------------------------------------------------------

% commented values are set conditionally, but would be activated for a new
% custom setup for your own problem, in which case you would remove the 
% conditional blocks for test_refl, test_ang, and test_fluence

% PHOTON SETTINGS
%~~~~~~~~~~~~~~~~~~~~~~~~~~ 
% N       = 5e5;            % number of photons (set below based on the test)
wmin    = 1e-4;             % min photon weight before discarding it
wrr     = 10;               % 1/wrr photons are reinjected (russian roulette)

% Sect. 5.1 Total diffuse reflectance and total transmittance
% Sect. 5.2 Angularly resolved diffuse reflectance and transmittance
if test_refl == true || test_ang == true
    opts.N   = 5e5;
    opts.Z   = 0.02;
    opts.g   = 0.75;
    opts.ka  = 10;
    opts.ks  = 90;
    opts.dz  = 0.001;
end

% Sect. 5.3 Depth resolved internal fluence
if test_fluence == true
    opts.N   = 1e6;
    opts.Z   = 4;
    opts.g   = 0.9;
    opts.ka  = .1;
    opts.ks  = 100;
    opts.dz  = 0.005;
end

% GRID SETTINGS
%~~~~~~~~~~~~~~~~~~~~~~~~~~ 
% Z       = 1;              % thickness of medium                   [cm]
% dz      = 0.005;          % vertical bin width                    [cm]
opts.R       = 2;                % cylindrical detection radius          [cm]
opts.A       = pi/2;             % angular detection radius              [rad]
opts.dr      = 0.001;            % radial bin width                      [cm]
opts.da      = A/30;             % angular bin width                     [rad]
opts.du      = da/(pi/2);        % angular bin width in cos(theta) coordinates
opts.nr      = roundn(R/dr,0);   % radial
opts.na      = roundn(A/da,0);   % angular
opts.nz      = roundn(Z/dz,0);   % vertical

% OPTICAL PROPERTIES
%~~~~~~~~~~~~~~~~~~~~~~~~~~ 
% g       = 0.75;             % asymmetry parameter               [-]
% ka      = 10;               % absorption coefficient            [cm-1]
% ks      = 90;               % scattering coefficient            [cm-1]
opts.w       = ks/(ka+ks);       % single-scattering albedo          [-]
opts.a       = 1-w;              % co-albedo                         [-]
opts.c       = 1/(ka+ks);        % extinction path length            [cm]
	
	

% CALL THE FUNCTION
%~~~~~~~~~~~~~~~~~~~~~~~~~~ 	
mcrt_out = mcrt(mcrt_opts)	
	

% this would go in b_mcrt_geom, after the detector rod is built, but only
% because the 'shp' object is not saved in 'rod', otherwise it could go
% here / at the end of mcrt, the user would pause the model when tfs = true
% and this could be used to plot the prior point, current point, and verify
% the solution agains linexlines2D
%==========================================================================
%% Example of plotting the photon position during the simulation
%==========================================================================

% % if, during the simulation, you want to pause and plot the photon:
% 
% % build a cylinder at a specific depth:
% %Zd = 11; % detector rod z-coordinate, Zd
% zc  = zc-Zd; % offset zc to be centered at Zd
% zd  = zd-Zd; % offset zd (photon launch point)
% zc0 = zc0-Zd;
% shp = alphaShape([xc(:),yc(:),zc(:)]); % rebuild the rod
% 
% % recover the prior photon location:
% x1  = x-ux*l;
% y1  = y-uy*l;
% z1  = z-uz*l;
% 
% % this is the intersection point, returned by rodintersect in mcrt.m
% x0  = x;
% y0  = y;
% z0  = z;
% 
% % plot it, assumes x,y,z are from mcrt.m (i.e., the simulation was paused)
% figure; plot(shp); hold on;
% plot3(xd,yd,zd,'*');                % detector z-coord (rcr launch point)
% plot3(xc0,yc0,zc0,'*');             % center of cylinder end-cap
% scatter3(x0,y0,z0,20,'g','filled')  % current photon location
% scatter3(x1,y1,z1,20,'r','filled')  % prior photon location
% plot3([x1 x0],[y1 y0],[z1 z0],'m'); % ray from prior location to intersection point
% view(90,0); axis tight              % change to 2-d view 
% xlabel('x'); ylabel('y'); zlabel('z');
% 
% % test against linexlines2D (work in the y-z plane)
% ycirc       = -rc:0.0001:rc;
% zcirc       = sqrt(rc^2 - (ycirc).^2);
% ycirc       = yc0+[ycirc -1.*ycirc];
% zcirc       = zc0+[zcirc -1.*zcirc];
% yzshp       = alphaShape(ycirc',zcirc');
% yzpoly      = polyshape(ycirc,zcirc);
% 
% % move the intersection point along the line into the cylinder (necessary
% % because the x,y,z that gave an intersection were overwritten in the main
% % function by rodintersect, so we just make up a new location
% b   = (z0-z1)/(y0-y1); a = z0-b*y0;
% x2  = x0;
% y2  = y0-0.5;
% z2  = a+y2*b;
% 
% % plot it
% figure; plot(yzshp); hold on; 
% scatter(y0,z0,'filled'); 
% scatter(y1,z1,'filled'); 
% plot([y1 y0],[z1 z0]);
% scatter(y2,z2); plot([y y2],[z z2]);    % new photon location
% 
% % test linexlines2D vs rodintersect
% % linexlines2D finds the intersections of the line segment with end
% % points [x1,y1] and [x2,y2] with the boundary of the polyshape object poly.
% 
% for n = 1:1000
%     tic;
%     yz  = linexlines2D(yzpoly,[y1,z1],[y2,z2]);
%     t1(n) = toc;
% end
% 
% for n = 1:1000
%     tic
%     [x,y,z,l,tfs] = rodintersect(x0,y0,z0,ux,uy,uz,l,rod);
%     t2(n) = toc;
% end
% 
% mean(t1)
% mean(t2)
% mean(t1)/mean(t2)
% std(t1)/std(t2)
% 
% % compare with rodintersect solution
% sprintf('%s%.5f%s%.5f','rodintersect: ',y,'  ',z)
% sprintf('%s%.5f%s%.5f','linexlines2D: ',yz(1),'  ',yz(2))
% 
% % plot it, assumes x,y,z are from mcrt.m (i.e., the simulation was paused)
% figure; plot(shp); hold on;
% scatter3(x1,y1,z1,'filled')
% scatter3(x2,y2,z2,'filled')
% scatter3(x,y,z,'filled')
% plot3([x1 x],[y1 y],[z1 z])
% xlabel('x'); ylabel('y'); zlabel('z');
% view(90,0); axis tight
% 
% % reflect the photon 
% us = -sqrt(rand); 
% phis = 2*pi*rand;
% [uxx,uyy,uzz] = chgdir(ux,uy,uz,us,phis)
% 
% M = AxelRot(us,u,x0)
% 
% [XYZnew, R, t] = AxelRot([x1;y1;z1],acosd(us),[uxx uyy uzz],[0 0 0])
% [XYZnew, R, t] = AxelRot([x;y;z],acosd(us),[uxx uyy uzz],[0 0 0])
% 
% scatter3(XYZnew(1),XYZnew(2),XYZnew(3),'filled');