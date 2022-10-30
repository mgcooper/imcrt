function RT = mcrt_test(ka,ks,g,Z,dz,N)

% optical coefficients
w       = ks/(ka+ks);       % single-scattering albedo          [-]
a       = 1-w;              % co-albedo                         [-]
c       = 1/(ka+ks);        % extinction path length            [cm]

% henyey-greenstein terms (pre-computed so it runs fast)
hg1     = 1/(2*g);
hg2     = (1+g^2);
hg3     = (1-g^2);
hg4     = 1+g;
hg5     = -2*g;
two_pi  = 2*pi;

% default settings
wmin    = 1e-4;             % min photon weight before discarding it
wrr     = 10;               % 1/wrr photons are reinjected (russian roulette)

% grid settings
R       = 2;                % cylindrical detection radius          [cm]
A       = pi/2;             % angular detection radius              [rad]
dr      = 0.001;            % radial bin width                      [cm]
da      = A/30;             % angular bin width                     [rad]
du      = da/(pi/2);        % angular bin width in cos(theta) coordinates
nr      = round(R/dr,0);    % radial
na      = round(A/da,0);    % angular
nz      = round(Z/dz,0);    % vertical

% initialize output grids with +1 for overflow
Adf_rz  = zeros(nz+1,nr+1);     % absorption, diffuse
Adr_z   = zeros(nz+1,1);        % absorption, direct
Tdf_ra  = zeros(na,nr+1);       % transmittance, diffuse
Rdf_ra  = zeros(na,nr+1);       % reflectance, diffuse 
Tdr     = 0;                    % transmittance, direct (unscattered)
Rdr     = 0;                    % reflectance, direct (unscattered)

% monte carlo
for n=1:N
    wt  = 1;                % new photon, weight = 1
    x   = 0;                % generated at position (0,0,0):
    y   = 0;               
    z   = 0;
    ux  = 0;                % with trajectory 0,0,1 (mu_x, mu_y, mu_z)
    uy  = 0;
    uz  = 1;                % = 1 for vertical, = 1-2*rand for isotropic
    ns  = 0;                % number of scattering events
    
% keep going until min intensity or z>0 is satisfied
    while wt > wmin
        l = -c*log(rand);                           % path length
        x = x+ux*l;                                 % new x-position
        y = y+uy*l;                                 % new y-position
        z = z+uz*l;                                 % new z-position

% grid indices
        ir = ceil(sqrt(x*x+y*y)/dr);                % radial index
        iz = ceil(z/dz);                            % vertical index
        if ir>nr; ir = nr+1; end                    % radial overflow 
        if iz>nz; iz = nz+1; end                    % vertical overflow

% score transmittance / reflectance
        if z>Z                                  	 % transmittance
            iu = ceil(acos(uz)/da);                 % angular index
            if ns==0 || uz < du/2
                Tdr = Tdr+wt;                       % direct
            else 
                Tdf_ra(iu,ir) = Tdf_ra(iu,ir)+wt;   % diffuse
            end
            break                                   % photon escapes
        end
        if z<0                                      % reflection
            iu = ceil(acos(-uz)/da);                % angular index
            if ir>nr; ir = nr+1; end                % overflow
            if ns==0 || -uz < du/2
                Rdr = Rdr+wt;                       % direct
            else
                Rdf_ra(iu,ir) = Rdf_ra(iu,ir)+wt;   % diffuse
            end   
            break                                   % photon escapes
        end

% score absorption/fluence
        awt = a*wt;
        if ns==0
            Adr_z(iz)=Adr_z(iz)+awt;                % direct
        else
            Adf_rz(iz,ir)=Adf_rz(iz,ir)+awt;        % diffuse
        end

% absorption and scattering by ice (henyey-greenstein with azim. symmetry)
        wt = wt*w; ns = ns+1;
        if g == 0
            us = 1-2*rand;
        else
            us = hg1*(hg2-(hg3/(hg4+hg5*rand))^2);
        end
        ps = two_pi*rand;                       % azimuth angle, phi_s
        if sqrt(1-uz*uz) < 1e-12
           % if initial direction not straight up or straight down
            ux = sqrt(1-us*us)*cos(ps);
            uy = sqrt(1-us*us)*sin(ps);
            uz = sign(uz)*us;
        else
            % initial direction straight up or down (sin(theta)=0)
            ux = sqrt(1-us*us)/sqrt(1-uz*uz)*(ux*uz*cos(ps)-uy*sin(ps))+ux*us;
            uy = sqrt(1-us*us)/sqrt(1-uz*uz)*(uy*uz*cos(ps)+ux*sin(ps))+uy*us;
            uz = -sqrt(1-us*us)*sqrt(1-uz*uz)*cos(ps)+uz*us;
        end
% russian roulette, keep one of every wrr photons and multiply by wrr
        if wt<wmin&&rand<(1/wrr)
            wt=wt*wrr;
        end
    end
end

% build a grid to calculate observable quantities (eq. 4.1/4.2 Wang)
[ri,ai,zi,dr,da,dz] = mcrt_build_grid(R,A,Z,dr,da,dz);

dsr         = 2*pi.*sin(ai).*da;                    % steradians
dA          = 2*pi.*ri.*dr;                         % radians
dV          = dA.*dz;                               % m3
cosa        = cos(ai);

% sum the 2-d arrays into 1-d and 0-d arrays (R=reflection, T=transmission)
Rdf_r       = sum(Rdf_ra,1);                        % Eq. 4.3
Rdf_a       = sum(Rdf_ra,2);                        % Eq. 4.4
Tdf_r       = sum(Tdf_ra,1);                        % Eq. 4.5
Tdf_a       = sum(Tdf_ra,2);                        % Eq. 4.6
Rdf         = sum(Rdf_r);                           % Eq. 4.7
Tdf         = sum(Tdf_r);                           % Eq. 4.8

% convert the photon counts to SI units
Rdf_ra      = Rdf_ra./(dA.*dsr.*cosa.*N);           % Eq. 4.9
Tdf_ra      = Tdf_ra./(dA.*dsr.*cosa.*N);           % Eq. 4.10
Rdf_r       = Rdf_r./(dA.*N);                       % Eq. 4.13
Tdf_r       = Tdf_r./(dA.*N);                       % Eq. 4.14
Rdf_a       = Rdf_a./(dsr.*N);                      % Eq. 4.15
Tdf_a       = Tdf_a./(dsr.*N);                      % Eq. 4.16
Rdr         = Rdr/N;                                % + Rsp/N
Tdr         = Tdr/N;
Rdf         = Rdf/N;                                % Eq. 4.17
Tdf         = Tdf/N;                                % Eq. 4.18
Tt          = Tdf+Tdr;

% Absorption
Adf_z       = sum(Adf_rz,2);                        % Eq. 4.20
Adf         = sum(Adf_z);                           % Eq. 4.22
Adf_rz      = Adf_rz./dV./N;                        % Eq. 4.23
Adr_z       = Adr_z./dz./N;                         % Eq. 4.24
Adf_z       = Adf_z./dz./N;                         % Eq. 4.25
Adf         = Adf./N;                               % Eq. 4.27
% phi_rz    = Adf_rz./ka;                           % Eq. 4.28
% phi_z     = Adf_z./ka;                            % Eq. 4.29
    
phi_rz      = (Adf_rz+Adr_z)./ka;                   % Eq. 4.28
phi_z       = (Adf_z+Adr_z)./ka;                    % Eq. 4.29

% arrange the output
RT.Rdf_ra   = Rdf_ra;
RT.Rdf_r    = Rdf_r;
RT.Rdf_a    = Rdf_a;
RT.Rdf      = Rdf;
RT.Rdr      = Rdr;
RT.Tdf_ra   = Tdf_ra;
RT.Tdf_r    = Tdf_r;
RT.Tdf_a    = Tdf_a;
RT.Tdf      = Tdf;
RT.Tdr      = Tdr;
RT.Tt       = Tt;

RT.Adf_z    = Adf_rz;
RT.Adf      = Adf_z;
RT.Adf_rz   = Adf_rz;
RT.Adr_z    = Adr_z;
RT.Adf_z    = Adf_z;
RT.Adf      = Adf;
RT.phi_rz   = phi_rz;
RT.phi_z    = phi_z;


% return the grid
RT.grid.ri  = ri;
RT.grid.ai  = ai;
RT.grid.zi  = zi;
RT.grid.dr  = dr;
RT.grid.da  = da;
RT.grid.dz  = dz;
