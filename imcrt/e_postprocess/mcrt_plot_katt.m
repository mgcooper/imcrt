clean

save_figs   = true;
save_data   = true;
fdate       = '20july';
iops        = 'mie_dE';
wavl        = [350;400;450;500;550;600;650;700;750];

% kabs w/o rod: 0-12 is nearly identical to 12-77, just slightly higher
% at 500 they are identical, at 600 and beyond 0-12 is higher
% kabs w/ rod: 0-12 is higher than 12-77, and very close to measured in the
% visible, and lower than measured in the NIR
% kice w/o rod: mirrors above
% based on this, we can just show 0-12 kabs and kice w/rod, since w/o rod
% is identical to the others, what we communicate is the impact of kabs on
% 0-12, which lets us make a statement about the IOPS of the 0-12

% one final thing I should do is a 0-77 cm keff and a 0-124 cm keff

%==========================================================================
%% set paths and load data
%==========================================================================
p.root  = 'GREENLAND/field/2018/a_submitted/monte_carlo/';
p.data  = [p.root fdate '/f_data/' iops '/'];
p.save  = [p.root fdate '/g_figs/' iops '/'];
p.kext  = ['GREENLAND/field/2018/data/processed/' fdate '/d_coefficients/'];
p.afec  = ['GREENLAND/field/2018/data/processed/' fdate '/f_iops/'];
p.tran  = ['GREENLAND/field/2018/data/processed/' fdate '/a_irradiance/']; 
p = setpath(p);

load([p.kext 'Kext.mat']);
load([p.kext 'k12.mat']);
load([p.data 'mcrt_output']);

if strcmp(iops,'mie_dE')
    load([p.afec 'mie_iops_dE']);
elseif strcmp(iops,'mie')
    load([p.afec 'mie_iops']);
elseif strcmp(iops,'ssa_dE')
    load([p.afec 'ssa_iops_dE']);
elseif strcmp(iops,'ssa')
    load([p.afec 'ssa_iops']);
end

[nz,nwavl,ncases]   = size(T);
%==========================================================================
%% % pull out kext from my observations and from asymptotic theory
%==========================================================================
lambda  = Kext.interp.wavl;
iwavl   = find(ismember(lambda,wavl));
kobs    = Kext.interp.kext;
kafec   = iops.kice.afec;
keff    = interp1(k12.wavl,k12.k,lambda);

sigX    = repmat(1.2/100*835/917,4,1);
sigY    = Kext.interp.Stau(iwavl,:)';
rxy     = 0;
alpha   = 0.05;

% compute kmcrt using least squares regression
for i = 1:ncases
    Tmcrt       = squeeze(T(:,:,i)); % depth x wavl, T=transmittance
	kmc12(i,:)  = -1.3/Zd(1)*log(Tmcrt(1,:)); % cases x wavl, keff for 0-12
    for j = 1:nwavl
        qi          = Tmcrt(:,j);
        y           = -log(qi);
        stats       = yorkfit(Zd./100,y,sigX,sigY(:,j),rxy,alpha);
        kmcrt(i,j)  = stats.b;  % cases x wavl
    end
end

if save_data == true
    if ~exist(p.save,'dir'); mkdir(p.save); end
    save([p.save 'mcrt_katt'],'kmcrt');
end
%==========================================================================
%% plot settings
%==========================================================================
c       = getdefaultcolors;
i1      = find(lambda == 350);
i2      = find(lambda == 750);

%==========================================================================
%% make the katt figure with theory for kice
%==========================================================================

f1      = figure('Units','inches','Position',[3 3 8 8]);
h1      = plot(lambda(i1:i2),kobs(i1:i2)); hold on; 
h2      = plot(lambda(i1:i2),kafec(i1:i2));
h3      = plot(lambda(i1:i2),keff(i1:i2));
for i = 1:4
    k1 = kmcrt(i,:);            % kice
    k2 = kmcrt(i+4,:);          % kabs
    s1(i) = scatter(wavl,k1,100,'s','LineWidth',1.5,'MarkerEdgeColor',c(i,:));
    s2(i) = scatter(wavl,k2,100,'o','LineWidth',1.5,'MarkerEdgeColor',c(i,:));
    if i == 4
        k3 = kmc12(i,:).*100;       % kice
        k4 = kmc12(i+4,:).*100;     % kabs
        s3 = scatter(wavl,k3,100,'s','LineWidth',1.5,'MarkerEdgeColor',c(i,:));
        s4 = scatter(wavl,k4,100,'o','LineWidth',1.5,'MarkerEdgeColor',c(i,:));
    end
end
ltext   = {'Field estimate (0-12 cm)','Field estimate (12-77 cm)',      ...
            'Asymptotic theory (12-77 cm)',                             ...
            'Monte Carlo, ideal diffusion, no rod',                     ...
            'Monte Carlo, ideal cosine, no rod',                        ...
            'Monte Carlo, ideal cosine, with rod',                      ...
            'Monte Carlo, non-ideal cosine, with rod'};
legend([h3 h1 h2 s1],ltext,'location','best');
xlabel('Wavelength [nm]');
ylabel('k_{att} [m^{-1}]');

set(gca,'XMinorGrid','on','YMinorGrid','on','YLim',[0.1 15],'XLim',[325 775])
set(gca,'XMinorTick','on','YMinorTick','on','MinorGridLineStyle','-');
set(gca,'YScale','log','TickDir','in','MinorGridAlpha',0.05)
%==========================================================================
%% save it
%==========================================================================
if save_figs == true
    export_fig([p.save 'mcrt_katt.png'],'-r400');
end

% %==========================================================================
% %% make the efold figure
% %==========================================================================
% 
% f2      = figure('Units','inches','Position',[3 3 8 8]);
% h3      = plot(lambda(i1:i2),1./kobs(i1:i2)); hold on; 
% h4      = plot(lambda(i1:i2),1./kafec(i1:i2));
% 
% for i = 1:4
%     e1 = 1./kmcrt(i,:);
%     e2 = 1./kmcrt(i+4,:);
%     s3(i) = scatter(wavl,e1,100,'s','LineWidth',1.5,'MarkerEdgeColor',c(i,:));
%     s4(i) = scatter(wavl,e2,100,'o','LineWidth',1.5,'MarkerEdgeColor',c(i,:));
% end
% 
% ltext   = {'Field estimate','Asymptotic theory',                        ...
%             'Monte Carlo, ideal diffusion, no rod',                     ...
%             'Monte Carlo, ideal cosine, no rod',                        ...
%             'Monte Carlo, ideal cosine, with rod',                      ...
%             'Monte Carlo, non-ideal cosine, with rod'};
% legend([h3 h4 s3],ltext,'location','best');
% xlabel('wavelength [nm]');
% ylabel('1/k_{att} [m^-1]');
% 
% set(gca,'XMinorGrid','on','YMinorGrid','on','YLim',[0 6],'XLim',[350 750])
% set(gca,'XMinorTick','on','YMinorTick','on','MinorGridLineStyle','-');
% set(gca,'TickDir','in','MinorGridAlpha',0.05)
% 
% %==========================================================================
% %% save it
% %==========================================================================
% if save_figs == true
%     export_fig([p.save 'mcrt_efold.png'],'-r400');
% end
% 
% % july 21, mie
% % kice: ideal nails it, and intuitive (higher) attenuation for other cases
% % kabs: intuitive at 500 nm, other wavl mixed, but variation is so small it
% % is probably within the range of variability
% 
% % july 20, mie
% % kice: ideal nails it, non-intuitive (lower) attenuation for other cases
% % kabs: ideal is higher, non-intuitive (lower) attenuation for other cases,
% % but non-ideal rcr + rod nails the observations
% 
% 
% % july 21, ssa
% % kice: ideal is too low, intuitive (higher) attenuation for other cases
% % kabs: ideal is too low, intuitive for other cases, and non-ideal + rod is
% % nearly perfect at 400 nm but too low at all others
% 
% % july 20, ssa
% % kice: is too low, non-intuitive (lower) attenuation for other cases
% % kabs: is too low, non-intuitive (lower) attenuation for other cases
% % overall, the july 20 ssa results are way off, indicating their might be a
% % problem somewhere
% %==========================================================================
% %% stuff I cut
% %==========================================================================
% 
% % load([ptran 'irradiance_combined']);        % transmittance
% % trans   = ice.interp.Qi_norm;
% % % trim the data to the case (ii) and wavelength range (i1:i2) of interest
% % ii      = 1;                                    % norod kice
% % ii      = 4;                                    % norod kabs
% % i1      = find(lambda == 350);
% % i2      = find(lambda == 750);
% % 
% % Tobs    = (trans(iwavl,:))';  % depth x wavl
% % Tmcrt   = squeeze(T(:,:,ii)); % depth x wavl
% % 
% % % compute kmcrt using finite-difference at each depth
% % for i = 1:nz
% %     Ti  = squeeze(T(i,:,ii));   % i = depth, : = wavl, ii = case
% %     kmcrtFD(i,:) = -1/Zd(i).*log(Ti);
% % end
% % 
% % % kmcrtFD     = kmcrtFD.*835/917;
% % % 1./kmcrtFD
% % 
% % % plot transmittance
% % if between == false
% % figure;
% % plot(lambda,trans); hold on;
% % for i = 1:nz
% %     scatter(wavl,Tmcrt(i,:),50,c(i,:),'filled');
% % end
% % 
% % % plot log transmittance as in fitting k-att values
% % figure;
% % for i = 1:nwavl
% %     plot(Zd,Tobs(:,i),'-','Color',c(i,:)); hold on;
% %     plot(Zd,Tmcrt(:,i),':','Color',c(i,:));
% % end
% % set(gca,'YScale','log')
% % xlabel('depth'); ylabel('T');
% % legend('400 nm','500 nm','600 nm','700 nm')
% % end
