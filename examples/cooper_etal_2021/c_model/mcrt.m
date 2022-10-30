%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clean               % v11
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
opts.path.root      = '/full/path/to/imcrt/';
opts.path.data      = [opts.path.root 'b_input/20july/'];
opts.path.save      = [opts.path.root 'd_output/20july/'];
opts.save_data      = true;
opts.use_ssa        = false;    % use specific surface area or mie theory?
opts.use_dE         = true;     
opts.transmittance  = true;
opts.reflectance    = false;
opts.fluence        = false;

%% load preprocessed data
load([opts.path.data 'mcrt_optical_input']);
load([opts.path.data 'mcrt_geometry_input']);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% main settings
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
opts.N      = 1e4;          % number of photons
opts.wmin   = 1e-5;         % min photon weight before termination
opts.wrr    = 10;           % 1/wrr photons are re-launched (russian roulette)
geom.dr     = 5;            % radial bin width                      [cm]
geom.da     = pi/2/30;      % angular bin width                     [rad]
geom.dz     = 1;            % vertical bin width                    [cm]
geom.R      = 400;          % cylindrical detection radius          [cm]
geom.A      = pi/2;         % angular detection radius              [rad]
nZ          = length(geom.Zd);

for zz = 1:nZ               % Loop through each detector position

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% detector geometry and grid system settings
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Zd      = geom.Zd(zz);              % detector position
    geom.Z  = ceil(Zd);                 % detector position (for grid systems)
    geom.zd = -Zd;                      % detector position (for photon release)
    rod     = geom.rod;                 % x,y,z,radius,length of detector rod 
    rod(3)  = rod(3)-Zd;                % reset the detector z-coordinate
    wrod    = iops.wrod;                % rod albedo
    nr      = roundn(geom.R/geom.dr,0); % radial
    na      = roundn(geom.A/geom.da,0); % angular
    nz      = roundn(geom.Z/geom.dz,0); % vertical

for ii = 1:8                % Loop through each experiment

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% case settings, optical settings, initialize photons, set file names
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    opts                        = mcrt_set_case(opts,ii);
    opts                        = mcrt_set_path(opts,geom);
    [wavl,asym,omeg,clen,kabs]  = mcrt_load_optics(opts,iops);
    [mux,muy,muz,~]             = mcrt_init_photons(opts,iops);
    disp([num2str(geom.Z) 'cm set ' int2str(ii)]);

for j = 1:length(wavl)      % monte carlo
    disp([num2str(wavl(j)) 'nm']);
    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
% optical properties at wavelength j
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    g   = asym(j);          % asymmetry parameter
    w   = omeg(j);          % particle albedo
    a   = 1-w;              % co-albedo
    c   = 100*clen(j);      % transport length (convert to cm)
    wr  = wrod(j);          % rod albedo
    hg1 = 1/(2*g);          % h-g scattering terms
    hg2 = (1+g^2);
    hg3 = (1-g^2);
    hg4 = 1-g;
    hg5 = 2*g;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% initialize output grids with +1 for overflow
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    dr  = geom.dr;                      % for faster indexing without 
    da  = geom.da;                      % accessing the geom struct
    dz  = geom.dz;                      % every time, see below
    if opts.transmittance == true
        Tdf_ra  = zeros(na,nr+1);       % transmitted diffuse radiance
        Tdr     = 0;                    % transmitted direct radiance 
    end
    if opts.reflectance == true
        Rdf_ra  = zeros(na,nr+1);       % reflected diffuse radiance
        Rdr     = 0;                    % reflected direct radiance 
    end
    if opts.fluence == true
        ka      = kabs(j)*100;          % absorption coefficient [cm-1]
        Adf_rz  = zeros(nz+1,nr+1);     % absorbed diffuse radiance
        Adr_z   = zeros(nz+1,1);        % absorbed unscattered radiance
    end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
% initialize photon packets one at a time
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for n=1:opts.N
        wt  = 1;                % new packet, weight = 1, launched at:
        x   = geom.xd;          % xd
        y   = geom.yd;          % yd
        z   = geom.zd;          % zd
        ux  = mux(n);           % x-direction
        uy  = muy(n);           % y-direction
        uz  = muz(n);           % z-direction
        ns  = 0;                % number of scattering events
        nd  = 0;                % number of detector interference events
        tfs = false;            % was prior event a scatter off the rod?
        tfd = false;            % assume the photon is not inside the rod
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% keep going until wt<wmin (abs), z>Z (refl), or z<0 (trans) is satisfied
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        while wt > opts.wmin
            l = -c*log(1-rand);         % path length
            x = x+ux*l;                 % new x-position
            y = y+uy*l;                 % new y-position
            z = z+uz*l;                 % new z-position            
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~            
% transmittance
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if z>0 && opts.transmittance == true
                ia = ceil(acos(uz)/da);                 % angular index
                ir = ceil(sqrt(x*x+y*y)/dr);            % radial index
                if ir>nr; ir = nr+1; end                % radial overflow 
                if ns==0 || ia == 0
                    Tdr = Tdr+wt;                       % unscattered
                else
                    Tdf_ra(ia,ir) = Tdf_ra(ia,ir)+wt;   % diffuse
                end
                break                                   % photon escapes
            end 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~            
% detector rod interference
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if opts.test_rod==true
                [x,y,z,l,tfs] = rodintersect(x,y,z,ux,uy,uz,l,rod);
            end
            if tfs == true
                tfr = true;             % assume it scatters into the rod
                while tfr == true
                    [x,y,z,ux,uy,uz] = rodscat(x,y,z,ux,uy,uz,l);
                    tfr = inrod(x,y,z,rod);
                end
                wt = wt*wr; nd = nd+1;                  % rod absorption
            end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% absorption/fluence
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if opts.fluence == true
                ir = ceil(sqrt(x*x+y*y)/dr);            % radial index        
                iz = ceil(z/dz);                        % vertical index
                if ir>nr; ir = nr+1; end                % radial overflow 
                if iz>nz; iz = nz+1; end                % vertical overflow
                awt = a*wt;                             % absorption
                if ns==0 && opts.fluence == true
                    Adr_z(iz) = Adr_z(iz)+awt;          % unscattered
                elseif opts.fluence == true
                    Adf_rz(iz,ir) = Adf_rz(iz,ir)+awt;  % diffuse
                end
            end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% h-g scattering by ice with azimuthal symmetry
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            wt = wt*w; ns = ns+1;                       % reflection
            us = hg1*(hg2-(hg3/(hg4+hg5*rand))^2);      % polar angle
            ps = two_pi*rand;                           % azimuth angle
            [ux,uy,uz] = chgdir(ux,uy,uz,us,ps);        % new direction
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~            
% russian roulette
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if wt<opts.wmin&&rand<(1/opts.wrr)
                wt=wt*opts.wrr;
            end
        end
    end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% compute observable quantities (absorption, fluence, transmittance)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if opts.transmittance == true
        [Tdf_ra,Tdf_r,Tdf_a,Tdf,Tdr,Tt] = scaleT(Tdf_ra,Tdr,geom,opts);
        mcrt_out.Tdf_ra = Tdf_ra;
        mcrt_out.Tdf_r  = Tdf_r;
        mcrt_out.Tdf_a  = Tdf_a;
        mcrt_out.Tdf    = Tdf;
        mcrt_out.Tdr    = Tdr;
        mcrt_out.Tt     = Tt;
    end
    if opts.fluence == true
        [Adf_rz,Adf,phi_rz,phi_z] = scaleA(Adf_rz,Adr_z,ka,geom,opts);
        mcrt_out.Adf_rz = Adf_rz;
        mcrt_out.Adf    = Adf;
        mcrt_out.phi_rz = phi_rz;
        mcrt_out.phi_z  = phi_z;
    end
    if opts.reflectance == true
        [Rdf_ra,Rdf_r,Rdf_a,Rdf,Rdr,Rt] = scaleR(Rdf_ra,Rdr,geom,opts);
        mcrt_out.Rdf_ra = Rdf_ra;
        mcrt_out.Rdf_r  = Rdf_r;
        mcrt_out.Rdf_a  = Rdf_a;
        mcrt_out.Rdf    = Rdf;
        mcrt_out.Rdr    = Rdr;
        mcrt_out.Rt     = Rt;
    end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% generate the scoring grid and save the data
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    [ri,ai,zi,~,~,~]    = buildgrid(geom);
    mcrt_out.nd         = nd;
    mcrt_out.ns         = ns;
    mcrt_out.ri         = ri;
    mcrt_out.ai         = ai;
    mcrt_out.zi         = zi;
    clear ri ai zi
    if opts.save_data == true
        save([opts.fsave.data '_' num2str(wavl(j)) 'nm'],'mcrt_out',    ...
            'opts','geom','iops');
    end
    clearvars -except opts geom iops two_pi wavl asym omeg clen kabs nz nr ...
                na mux muy muz rod wrod nZ
end
end
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% I cut these out
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%                   Tdr_ss(na) = Tdr_ss(na)+(wt*wt);    % variance                    
%                   Tdf_ss(ia,ir) = Tdf_ss(ia,ir)+(wt*wt);  % variance                
%               Adr_ss(iz) = Adr_ss(iz)+(awt*awt);      % variance                
%               Adf_ss(iz,ir) = Adf_ss(iz,ir)+(awt*awt); % variance

%       Adf_ss  = zeros(nz+1,nr+1);     % squared weights - for variance
%       Adr_ss  = zeros(nz+1,1);        % squared weights - for variance        
%       Tdr     = zeros(na,1);          % transmitted direct radiance 
%       Tdf_ss  = zeros(na,nr+1);       % squared weights - for variance        
%       Tdr_ss  = zeros(na,1);          % squared weights - for variance        
%       Rdr     = zeros(na,1);          % reflected direct radiance 
%       Rdf_ss  = zeros(na,nr+1);       % squared weights - for variance        
%       Rdr_ss  = zeros(na,1);          % squared weights - for variance

% % initialize output grids with +1 for overflow
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     if opts.internal_fluence == true
%         Adf_rz  = zeros(nz+1,nr+1);     % absorbed diffuse radiance
%         Adf_ss  = zeros(nz+1,nr+1);     % squared weights - for variance
%         Adr_z   = zeros(nz+1,1);        % absorbed unscattered radiance
%         Adr_ss  = zeros(nz+1,1);        % squared weights - for variance
%     end
%     if opts.internal_radiance == true
%         Adf_ra  = zeros(na,nr+1);       % absorbed diffuse radiance
%         Adf_ss  = zeros(na,nr+1);       % squared weights - for variance
%         Adr_a   = zeros(na,1);          % absorbed unscattered radiance
%         Adr_ss  = zeros(na,1);          % squared weights - for variance
%     end
%     if opts.transmittance == true
%         Tdf_ra  = zeros(na,nr+1);       % transmitted diffuse radiance
%         Tdf_ss  = zeros(na,nr+1);       % squared weights - for variance
%         Tdr     = zeros(na,1);          % transmitted direct radiance 
%         Tdr_ss  = zeros(na,1);          % squared weights - for variance
%     end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% activate this to score angular 'absorption'
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%           ia = ceil(acos(uz)/da);                         % angular index
%           if ns==0
%               Adr_a(ia) = Adr_a(ia)+awt;                  % unscattered
%               Adr_ss(ia) = Adr_ss(ia)+(awt*awt);          % variance
%           else                                        
%               Adf_ra(ia,ir) = Adf_ra(ia,ir)+awt;          % diffuse
%               Adf_ss(ia,ir) = Adf_ss(ia,ir)+(awt*awt);    % variance
%           end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% optical radii to compute the single-scattering values using Mie theory. 

% reff       = [  0.040, 0.050, 0.065, 0.080, 0.100, ...
%                 0.120, 0.140, 0.170, 0.200, 0.240, ...
%                 0.290, 0.350, 0.420, 0.500, 0.570, ...
%                 0.660, 0.760, 0.870, 1.000, 1.100, ...
%                 1.250, 1.400, 1.600, 1.800, 2.000, ...
%                 2.250, 2.500, 2.750, 3.000, 3.500, ...
%                 4.000, 4.500, 5.000, 5.500, 6.000 ];
