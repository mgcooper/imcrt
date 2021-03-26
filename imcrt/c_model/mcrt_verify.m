%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clean %(alias for clc;clear all;close all)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% a stripped down version to verify accuracy by comparison with van de
% Hulst values and values reported in Wang et al. 1995:

% Wang, L., Jacques, S. L. and Zheng, L.: MCML?Monte Carlo modeling of
% light transport in multi layered tissues, Computer Methods and Programs
% in Biomedicine, 47(2), 131?146, https://doi.org/10.1016/0169
% 2607(95)01640 F, 1995.

test_refl       = false; % verify total reflectance
test_ang        = false; % verify angular reflectance
test_fluence    = true;  % verify fluence
save_figs       = false;

opts.path.data  = 'path/to/b_input/';
opts.path.save  = 'path/to/d_output/wang/';
opts.save_data  = false;

% load the predefined geometry
load([opts.path.data 'mcrt_geometry_input']);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% experimental setup
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% PHOTON SETTINGS
%~~~~~~~~~~~~~~~~~~~~~~~~~~ 
% N       = 5e5;              % number of photons (set below)
wmin    = 1e-4;             % min photon weight before discarding it
wrr     = 10;               % 1/wrr photons are reinjected (russian roulette)

% Sect. 5.1 Total diffuse reflectance and total transmittance
% Sect. 5.2 Angularly resolved diffuse reflectance and transmittance
if test_refl == true || test_ang == true
    N   = 5e5;
    Z   = 0.02;
    g   = 0.75;
    ka  = 10;
    ks  = 90;
    dz  = 0.001;
end

% Sect. 5.3 Depth resolved internal fluence
if test_fluence == true
    N   = 1e6;
    Z   = 4;
    g   = 0.9;
    ka  = .1;
    ks  = 100;
    dz  = 0.005;
end

% GRID SETTINGS
%~~~~~~~~~~~~~~~~~~~~~~~~~~ 
% Z       = 1;              % thickness of medium (12,36,58,77)     [cm]
% dz      = 0.005;          % vertical bin width                    [cm]
R       = 2;                % cylindrical radius of detection       [cm]
A       = pi/2;             % angular detection radius              [rad]
dr      = 0.001;            % radial bin width                      [cm]
da      = A/30;             % angular bin width                     [rad]
du      = da/(pi/2);        % angular bin width in cos(theta) coordinates
nr      = roundn(R/dr,0);   % radial
na      = roundn(A/da,0);   % angular
nz      = roundn(Z/dz,0);   % vertical

% OPTICAL PROPERTIES
%~~~~~~~~~~~~~~~~~~~~~~~~~~ 
% g       = 0.75;             % asymmetry parameter               [-]
% ka      = 10;               % absorption coefficient            [cm-1]
% ks      = 90;               % scattering coefficient            [cm-1]
w       = ks/(ka+ks);       % single-scattering albedo          [-]
a       = 1-w;              % co-albedo                         [-]
c       = 1/(ka+ks);        % extinction path length            [cm]
hg1     = 1/(2*g);          % h-g scattering terms
hg2     = (1+g^2);
hg3     = (1-g^2);
hg4     = 1+g;              % 1-g for BMC
hg5     = -2*g;             % 2*g for BMC
two_pi  = 2*pi;             

% this shouldn't be needed
% opts    = mcrt_set_case(opts,1);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% initialize output grids with +1 for overflow
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Adf_rz  = zeros(nz+1,nr+1);     % absorption, diffuse
Adr_z   = zeros(nz+1,1);        % absorption, direct
Tdf_ra  = zeros(na,nr+1);       % transmittance, diffuse
Rdf_ra  = zeros(na,nr+1);       % reflectance, diffuse 
Tdr     = 0;                    % transmittance, direct (unscattered)
Rdr     = 0;                    % reflectance, direct (unscattered)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% monte carlo
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for n=1:N
    wt  = 1;                % new photon, weight = 1
    ns  = 0;                % number of scattering events
    nd  = 0;                % number of rod events    
    x   = 0;                % generated at position (0,0,0):
    y   = 0;               
    z   = 0;
    ux  = 0;                % with trajectory 0,0,1 (mu_x, mu_y, mu_z)
    uy  = 0;
    uz  = 1;                % = 1 for vertical, = 1-2*rand for isotropic

% keep going until min intensity or z>0 is satisfied
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    while wt > wmin
        l = -c*log(rand);                           % path length
        x = x+ux*l;                                 % new x-position
        y = y+uy*l;                                 % new y-position
        z = z+uz*l;                                 % new z-position
% grid indices
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ir = ceil(sqrt(x*x+y*y)/dr);                % radial index
        iz = ceil(z/dz);                            % vertical index
        if ir>nr; ir = nr+1; end                    % radial overflow 
        if iz>nz; iz = nz+1; end                    % vertical overflow

% score transmittance / reflectance
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
        if z>Z                                  	% transmittance
            iu = ceil(acos(uz)/da);                 % angular index
            if ns==0 || uz < du/2
                Tdr = Tdr+wt;                       % unscattered
            else 
                Tdf_ra(iu,ir) = Tdf_ra(iu,ir)+wt;   % diffuse
            end
            break                                   % photon escapes
        end
        if z<0                                      % reflection
            iu = ceil(acos(-uz)/da);                % angular index
            if ir>nr; ir = nr+1; end                % overflow
            if ns==0 || -uz < du/2
                Rdr = Rdr+wt;                       % unscattered
            else
                Rdf_ra(iu,ir) = Rdf_ra(iu,ir)+wt;   % diffuse
            end   
            break                                   % photon escapes
        end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% score internal radiant flux for absorption/fluence
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        awt = a*wt;
        if ns==0
            Adr_z(iz)=Adr_z(iz)+awt;                % unscattered
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
        [ux,uy,uz] = chgdir(ux,uy,uz,us,ps);
% russian roulette, keep one of every wrr photons and multiply by wrr
        if wt<wmin&&rand<(1/wrr)
            wt=wt*wrr;
        end
    end
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% scale R, T, and A
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% build a grid to calculate observable quantities (eq. 4.1/4.2 Wang)
[ri,ai,zi,dr,da,dz] = mcrt_build_grid(R,A,Z,dr,da,dz);
dsr         = 2*pi.*sin(ai).*da;
dA          = 2*pi.*ri.*dr;
dV          = dA.*dz;
cosa        = cos(ai);

% sum the 2-d arrays into 1-d and 0-d arrays
Rdf_r       = sum(Rdf_ra,1);                        % Eq. 4.3
Rdf_a       = sum(Rdf_ra,2);                        % Eq. 4.4
Tdf_r       = sum(Tdf_ra,1);                        % Eq. 4.5
Tdf_a       = sum(Tdf_ra,2);                        % Eq. 4.6
Rdf         = sum(Rdf_r);                           % Eq. 4.7
Tdf         = sum(Tdf_r);                           % Eq. 4.8

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
Adf_z       = sum(Adf_rz,2);            % Eq. 4.20
Adf         = sum(Adf_z);               % Eq. 4.22
Adf_rz      = Adf_rz./dV./N;            % Eq. 4.23
Adr_z       = Adr_z./dz./N;             % Eq. 4.24
Adf_z       = Adf_z./dz./N;             % Eq. 4.25
Adf         = Adf./N;                   % Eq. 4.27
% phi_rz      = Adf_rz./ka;               % Eq. 4.28
% phi_z       = Adf_z./ka;                % Eq. 4.29
    
phi_rz      = (Adf_rz+Adr_z)./ka;       % Eq. 4.28
phi_z       = (Adf_z+Adr_z)./ka;        % Eq. 4.29

% figure; plot(Adf_z); hold on; plot(Adr_z)
% 
% figure; plot(zi,sum(phi_rz,2)); set(gca,'YScale','log')
% figure; plot(zi,phi_z); set(gca,'YScale','log')
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% validate by comparison with van de hulst, Vol. 2, pg. 435, Table 35 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Rdvdh       = 0.09739;                  % diffuse reflectance
Tdvdh       = 0.66096;                  % diffuse transmittance
Tdrs        = exp(-(ka+ks)*Z);          % direct transmittance
Rd_e        = Rdf-Rdvdh;                % error
Td_e        = Tt-Tdvdh;
Tdr_e       = Tdr-Tdrs;

% v.d.h. values
uvdh        = [0 0.1 0.3 0.5 0.7 0.9 1.0];
Rvdh        = [.10641 .13304 .13856 .11956 .09411 .07132 .06180];
Tvdh        = [.11811 .16486 .22473 .27785 .36889 .72155 2.40270];
Tvdh        = fliplr(Tvdh.*uvdh./pi);
Rvdh        = fliplr(Rvdh.*uvdh./pi);
uvdh        = acos(fliplr(uvdh))./pi;

if test_refl == true || test_ang == true
    
    % for comparison with vdh, add back the direct R/T
    Rdf_a(1)    = Rdf_a(1)+Rdr; % total transmittance includes the direct
    Tdf_a(1)    = Tdf_a(1)+Tdr; % total transmittance includes the direct

    
    fprintf(' Rd\n iMCRT: %.5f\n van de Hulst: %.5f\n error: %.5f\n',Rdf,Rdvdh,Rd_e)
    fprintf(' \n Td\n iMCRT: %.5f\n van de Hulst: %.5f\n error: %.5f\n',Tt,Tdvdh,Td_e)
    fprintf(' \n Tdr\n iMCRT: %.5f\n theory: %.5f\n error: %.5f\n',Tdr,Tdrs,Tdr_e)

    figure('Units','in','Position',[3 3 12 5]); 
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact'); nexttile
    scatter(ai./pi,Rdf_a,80,'filled','s'); hold on; box on
    scatter(uvdh,Rvdh,80,'filled'); 
    xlabel('exiting angle, \alpha [\pi rad]');
    ylabel('R_d(\alpha) [sr^{-1}]');
    legend('This study','van de Hulst');
    set(gca,'TickDir','in','XMinorTick','on','YMinorTick','on')
    
    nexttile
    scatter(ai./pi,Tdf_a,80,'filled','s'); hold on; box on
    scatter(uvdh,Tvdh,80,'filled'); 
    xlabel('exiting angle, \alpha [\pi rad]');
    ylabel('T_d(\alpha) [sr^{-1}]');
    legend('This study','van de Hulst');
    set(gca,'TickDir','in','XMinorTick','on','YMinorTick','on')
    
    if save_figs == true
        export_fig([opts.path.save 'mcrt_ang_validate.png'],'-r400');
    end
    
elseif test_fluence == true
    
    figure;
    scatter(zi,phi_z);
    set(gca,'YScale','log','YLim',[0.5 10],'XLim',[0 1]);
    xlabel('z (cm)');
    ylabel('Fluence (-)');
end

% see Vol. 1, pg. 262, Table 12 for isotropic
