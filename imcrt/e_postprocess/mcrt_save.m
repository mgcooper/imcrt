clean

save_data   = true;
fdate       = '20july';
iops        = 'mie_dE';
between     = false;
wavl        = [350;400;450;500;550;600;650;700;750];

p.root      = 'GREENLAND/field/2018/a_submitted/monte_carlo/';

if strcmp(fdate,'20july') == true
    p.save  = [p.root '20july/f_data/' iops '/'];
    p.root  = [p.root '20july/d_output/' iops '/'];
    Z       = [12;33;53;70];
elseif strcmp(fdate,'21july') == true
    Z       = [50;63;77;116];
    p.save  = [p.root '21july/f_data/' iops '/'];
    p.root  = [p.root '21july/d_output/' iops '/'];
end
p        = setpath(p);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% read / save the data
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nwavl       = length(wavl);
nz          = length(Z);
ncases      = 8;

for ii = 1:nz
    p.data  = [p.root int2str(Z(ii)) 'cm/'];

    for jj = 1:ncases
        switch jj
            case 1                                      % ideal kice
                use_kabs    = 'kice';
                test_rod    = 'ideal';
            case 2                                      % norod kice
                use_kabs    = 'kice';
                test_rod    = 'norod';
            case 3                                      % rod kice
                use_kabs    = 'kice';
                test_rod    = 'rod';
            case 4                                      % rcr kice
                use_kabs    = 'kice';
                test_rod    = 'rcr';
            case 5                                      % ideal kice
                use_kabs    = 'kabs';
                test_rod    = 'ideal';
            case 6                                      % norod kabs
                use_kabs    = 'kabs';
                test_rod    = 'norod';
            case 7                                      % rod kabs
                use_kabs    = 'kabs';
                test_rod    = 'rod';
            case 8                                      % rcr kabs
                use_kabs    = 'kabs';
                test_rod    = 'rcr';
        end

      % load each wavelength        
        for kk = 1:nwavl
            swavl = int2str(wavl(kk));
            load([p.data 'data_' test_rod '_' use_kabs '_' swavl 'nm']);
          % total transmittance
            T(ii,kk,jj)  = mcrt_out.Tt;
    
        end     % kk = number of wavelengths
    end         % jj = number of cases
end             % ii = number of z depths
                % T  = [nz x nwavl x ncases]

% save the data
Zd  = geom.Zd;
if save_data == true
    if ~exist(p.save,'dir'); mkdir(p.save); end
    save([p.save 'mcrt_output'],'T','Zd');
end