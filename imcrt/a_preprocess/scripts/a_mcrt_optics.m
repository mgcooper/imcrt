clean
%==========================================================================

save_data   = true;
july20      = true;
july21      = false;
wavl        = [350;400;450;500;550;600;650;700;750];

%==========================================================================
%% set paths
%==========================================================================
p.data      = 'GREENLAND/field/2018/a_submitted/monte_carlo/';
p.save      = 'GREENLAND/field/2018/a_submitted/monte_carlo/';

if july20 == true
    p.data  = [p.data '20july/a_preprocess/data/'];
    p.save  = [p.save '20july/b_input/'];
elseif july21 == true
    p.data  = [p.data '21july/a_preprocess/data/'];
    p.save  = [p.save '21july/b_input/'];
end
p = setpath(p);

%==========================================================================
%% optical properties - mie
%==========================================================================
load([p.data 'mie_iops_dE.mat']);
load([p.data 'rod_refl.mat']);

iwavl           = find(ismember(iops.kice.wavl,wavl));  % wavelength (nm)
mie.kice.wavl   = wavl;
mie.kice.asym   = iops.kice.asym(iwavl);        % asymmetry (-)
mie.kice.omeg   = iops.kice.omeg(iwavl);        % particle albedo (-)
mie.kice.sige   = iops.kice.sige(iwavl);        % particle ext coef. (1/m)
mie.kice.siga   = iops.kice.siga(iwavl);        % particle abs coef. (1/m)
mie.kice.sigs   = iops.kice.sigs(iwavl);        % particle sca coef. (1/m)
mie.kice.clen   = 1./iops.kice.sige(iwavl);     % particle ext path (m)
mie.kice.kabs   = iops.kice.kabs(iwavl);        % abs coef pure ice (1/m)
mie.kice.reff   = iops.kice.reff;
mie.kice.rhos   = iops.kice.rhos;
check1          = mie.kice.siga+mie.kice.sigs-mie.kice.sige % == 0

% repeat for kabs. this assumes sigs is unaffected by kabs
mie.kabs.wavl   = wavl;
mie.kabs.asym   = iops.kabs.asym(iwavl);
mie.kabs.omeg   = iops.kabs.omeg(iwavl);
mie.kabs.sige   = iops.kabs.sige(iwavl);
mie.kabs.siga   = iops.kabs.siga(iwavl);
mie.kabs.sigs   = iops.kabs.sigs(iwavl);
mie.kabs.clen   = 1./iops.kabs.sige(iwavl);
mie.kabs.kabs   = iops.kabs.kabs(iwavl);
mie.kabs.reff   = iops.kabs.reff;
mie.kabs.rhos   = iops.kabs.rhos;
check1          = mie.kabs.siga+mie.kabs.sigs-mie.kabs.sige % == 0

%==========================================================================
%% optical properties - ssa
%==========================================================================
load([p.data 'ssa_iops_dE.mat']);

ssa.kice.wavl   = wavl;
ssa.kice.asym   = iops.kice.asym(iwavl);        % asymmetry (-)
ssa.kice.omeg   = iops.kice.omeg(iwavl);        % particle albedo (-)
ssa.kice.sige   = iops.kice.sige(iwavl);        % particle ext coef. (1/m)
ssa.kice.siga   = iops.kice.siga(iwavl);        % particle abs coef. (1/m)
ssa.kice.sigs   = iops.kice.sigs(iwavl);        % particle sca coef. (1/m)
ssa.kice.clen   = 1./iops.kice.sige(iwavl);     % particle ext path (m)
ssa.kice.kabs   = iops.kice.kabs(iwavl);        % abs coef pure ice (1/m)
ssa.kice.rhos   = iops.kice.rhos;
ssa.kice.ssa    = iops.kice.ssa;
check3          = ssa.kice.siga+ssa.kice.sigs-ssa.kice.sige % == 0

% repeat for kabs. this assumes sigs is unaffected by kabs
ssa.kabs.wavl   = wavl;
ssa.kabs.asym   = iops.kabs.asym(iwavl);
ssa.kabs.omeg   = iops.kabs.omeg(iwavl);
ssa.kabs.sige   = iops.kabs.sige(iwavl);
ssa.kabs.siga   = iops.kabs.siga(iwavl);
ssa.kabs.sigs   = iops.kabs.sigs(iwavl);
ssa.kabs.clen   = 1./iops.kabs.sige(iwavl);
ssa.kabs.kabs   = iops.kabs.kabs(iwavl);
ssa.kabs.rhos   = iops.kabs.rhos;
ssa.kabs.ssa    = iops.kabs.ssa;
check4          = ssa.kabs.siga+ssa.kabs.sigs-ssa.kabs.sige % == 0

ssa_dE          = ssa;
mie_dE          = mie; 
clear ssa mie

%==========================================================================
%% repeat with Eddington inversion
%==========================================================================
load([p.data 'mie_iops.mat']);

iwavl           = find(ismember(iops.kice.wavl,wavl));  % wavelength (nm)
mie.kice.wavl   = wavl;
mie.kice.asym   = iops.kice.asym(iwavl);        % asymmetry (-)
mie.kice.omeg   = iops.kice.omeg(iwavl);        % particle albedo (-)
mie.kice.sige   = iops.kice.sige(iwavl);        % particle ext coef. (1/m)
mie.kice.siga   = iops.kice.siga(iwavl);        % particle abs coef. (1/m)
mie.kice.sigs   = iops.kice.sigs(iwavl);        % particle sca coef. (1/m)
mie.kice.clen   = 1./iops.kice.sige(iwavl);     % particle ext path (m)
mie.kice.kabs   = iops.kice.kabs(iwavl);        % abs coef pure ice (1/m)
mie.kice.reff   = iops.kice.reff;
mie.kice.rhos   = iops.kice.rhos;
check1          = mie.kice.siga+mie.kice.sigs-mie.kice.sige % == 0

% repeat for kabs. this assumes sigs is unaffected by kabs
mie.kabs.wavl   = wavl;
mie.kabs.asym   = iops.kabs.asym(iwavl);
mie.kabs.omeg   = iops.kabs.omeg(iwavl);
mie.kabs.sige   = iops.kabs.sige(iwavl);
mie.kabs.siga   = iops.kabs.siga(iwavl);
mie.kabs.sigs   = iops.kabs.sigs(iwavl);
mie.kabs.clen   = 1./iops.kabs.sige(iwavl);
mie.kabs.kabs   = iops.kabs.kabs(iwavl);
mie.kabs.reff   = iops.kabs.reff;
mie.kabs.rhos   = iops.kabs.rhos;
check1          = mie.kabs.siga+mie.kabs.sigs-mie.kabs.sige % == 0

%==========================================================================
%% optical properties - ssa
%==========================================================================
load([p.data 'ssa_iops.mat']);

ssa.kice.wavl   = wavl;
ssa.kice.asym   = iops.kice.asym(iwavl);        % asymmetry (-)
ssa.kice.omeg   = iops.kice.omeg(iwavl);        % particle albedo (-)
ssa.kice.sige   = iops.kice.sige(iwavl);        % particle ext coef. (1/m)
ssa.kice.siga   = iops.kice.siga(iwavl);        % particle abs coef. (1/m)
ssa.kice.sigs   = iops.kice.sigs(iwavl);        % particle sca coef. (1/m)
ssa.kice.clen   = 1./iops.kice.sige(iwavl);     % particle ext path (m)
ssa.kice.kabs   = iops.kice.kabs(iwavl);        % abs coef pure ice (1/m)
ssa.kice.rhos   = iops.kice.rhos;
ssa.kice.ssa    = iops.kice.ssa;
check3          = ssa.kice.siga+ssa.kice.sigs-ssa.kice.sige % == 0

% repeat for kabs. this assumes sigs is unaffected by kabs
ssa.kabs.wavl   = wavl;
ssa.kabs.asym   = iops.kabs.asym(iwavl);
ssa.kabs.omeg   = iops.kabs.omeg(iwavl);
ssa.kabs.sige   = iops.kabs.sige(iwavl);
ssa.kabs.siga   = iops.kabs.siga(iwavl);
ssa.kabs.sigs   = iops.kabs.sigs(iwavl);
ssa.kabs.clen   = 1./iops.kabs.sige(iwavl);
ssa.kabs.kabs   = iops.kabs.kabs(iwavl);
ssa.kabs.rhos   = iops.kabs.rhos;
ssa.kabs.ssa    = iops.kabs.ssa;
check4          = ssa.kabs.siga+ssa.kabs.sigs-ssa.kabs.sige % == 0

%==========================================================================
%% reassign the values to an 'iops' structure and save it
%==========================================================================
clear iops
iops.mie        = mie;
iops.ssa        = ssa;
iops.mie_dE     = mie_dE;
iops.ssa_dE     = ssa_dE;

% I am putting the rcr and albedo of detector rod here for now
load([p.data 'rcrfunc.mat']);
iops.rcr        = rcr';
iops.wrod       = (rod_refl.wrod(ismember(rod_refl.wavl,wavl)))';

if save_data == true
    save([p.save 'mcrt_optical_input'],'iops')
end