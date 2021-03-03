clean
%==========================================================================

save_data   = true;
july20      = false;
july21      = true;

% need to remove the redundant dates folders 

p.data      = 'GREENLAND/field/2018/a_submitted/monte_carlo/';
p.save      = 'GREENLAND/field/2018/a_submitted/monte_carlo/';
%==========================================================================
%% set paths
%==========================================================================
if july20 == true
    p.data  = [p.data '20july/a_preprocess/data/'];
    p.save  = [p.save '20july/b_input/'];
    p       = setpath(p);
    load([p.data 'IEQ_with_screw']);
    Zd      = 100.*cumsum(IEQ.dIEQ_k20_constant);
elseif july21 == true
    p.data  = [p.data '21july/a_preprocess/data/'];
    p.save  = [p.save '21july/b_input/'];
    load([p.data 'IEQ_with_screw']);
    p       = setpath(p);
    Zd      = 100.*cumsum(IEQ.dIEQ_k21_constant);
end

% static values to reduce overhead
two_pi      = 2*pi;
cos_zero    = 1-1e-12;      % cosine of about 1e-6 rad.
cos_90d     = 1e-6;         % cosine of about 1.57 - 1e-6 rad.
%==========================================================================
%% build the geometry
%==========================================================================

% define the origin, from where photons are launched
x0          = 0.0;              % optional x-y offset, useful for plotting
y0          = 0.0;              
z0          = 0.0;              % z-position (depth) of detector (cm)

% update the origin to be the position of the detector
xd          = x0+0.0;           % x-position of detector (cm)
yd          = y0+0.0;           % y-position of detector (cm)
zd          = z0+0.0;           % z-position of detector (cm)

% detector rod geometry
rc          = 2.065;            % radius of detector cylinder (cm) 
hc          = 200;              % length/height of detector cylinder (cm)

% position of detector rod (center of cylinder base)
xc0         = x0-5.0;           % x-position of detector rod end-center
yc0         = y0;               % y-position of detector rod end-center
zc0         = z0-rc;            % z-position of detector rod end-center

% generate a 3-d cylinder to represent the rod
nc          = 1000;             % # of points to define the circumference
[zc,yc,xc]  = cylinder(rc, nc);
xc          = xc0+xc.*hc;
yc          = yc0+yc;
zc          = zc0+zc;

% generate a shape object to access matlab computational geometry functions
shp         = alphaShape([xc(:),yc(:),zc(:)]);

% detector rod geometry, used as input to ray tracing functions
rod         = [xc0,yc0,zc0,rc,hc];

% plot it
figure;
plot(shp); hold on;
plot3(xd,yd,zd,'*');
plot3(xc0,yc0,zc0,'*');

%==========================================================================
%% put the values into a structure and save it
%==========================================================================

geom.Zd     = Zd;
geom.x0     = x0;
geom.y0     = y0;
geom.z0     = z0;
geom.xd     = xd;
geom.yd     = yd;
geom.zd     = zd;
geom.rod    = rod;

if save_data == true
    save([p.save 'mcrt_geometry_input'],'geom','two_pi');
end

% Z = total depth, dZ = sampling frequency for logging intensity/power
% dz          = 1;                % (cm)
% nz          = Z/dz;             % # of layers (must be a round number)

