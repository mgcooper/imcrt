clean

% set options
test_refl      = true;     % verify total reflectance or not
test_ang       = true;     % verify angular reflectance or not
test_fluence   = false;    % verify fluence (slower) or not

% Sect. 5.1 Total diffuse reflectance and total transmittance
% Sect. 5.2 Angularly resolved diffuse reflectance and transmittance
if test_refl == true || test_ang == true
    N   = 1e6;       % number of photon packets
    Z   = 0.02;      % total thickness of medium
    g   = 0.75;      % asymmetry parameter
    ka  = 10;        % absorption coefficient (m-1)
    ks  = 90;        % scattering coefficient (m-1)
    dz  = 0.001;
end

% Sect. 5.3 internal fluence
if test_fluence == true
    N   = 1e6;
    Z   = 4;
    g   = 0.9;
    ka  = .1;
    ks  = 100;
    dz  = 0.005;
end

%--------------------------------------------------------------------------
% Run the model
%--------------------------------------------------------------------------
RT = mcrt_test(ka,ks,g,Z,dz,N);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% pull out the values to plot them
Rdr   = RT.Rdr;
Rdf   = RT.Rdf;
Tt    = RT.Tt;
Tdr   = RT.Tdr;
Rdf_a = RT.Rdf_a;
Tdf_a = RT.Tdf_a;

% compare with van de hulst, Vol. 2, pg. 435, Table 35 (vdh=van de hulst)
Rd_vdh  = 0.09739;                  % diffuse reflectance
Td_vdh  = 0.66096;                  % diffuse transmittance
Tdrs    = exp(-(ka+ks)*Z);          % direct transmittance
Rd_err  = Rdf-Rd_vdh;               % error
Td_err  = Tt-Td_vdh;
Tdr_err = Tdr-Tdrs;

% v.d.h. values
uvdh  = [0 0.1 0.3 0.5 0.7 0.9 1.0];
Rvdh  = [.10641 .13304 .13856 .11956 .09411 .07132 .06180];
Tvdh  = [.11811 .16486 .22473 .27785 .36889 .72155 2.40270];
Tvdh  = fliplr(Tvdh.*uvdh./pi);
Rvdh  = fliplr(Rvdh.*uvdh./pi);
uvdh  = acos(fliplr(uvdh))./pi;

%% plot the results
if test_refl == true || test_ang == true
    
    % for comparison with vdh, add back the direct R/T
    Rdf_a(1)    = Rdf_a(1)+Rdr; % total transmittance includes the direct
    Tdf_a(1)    = Tdf_a(1)+Tdr; % total transmittance includes the direct
    
    fprintf(' Rd\n iMCRT: %.5f\n van de Hulst: %.5f\n error: %.5f\n',Rdf,Rd_vdh,Rd_err)
    fprintf(' \n Td\n iMCRT: %.5f\n van de Hulst: %.5f\n error: %.5f\n',Tt,Td_vdh,Td_err)
    fprintf(' \n Tdr\n iMCRT: %.5f\n theory: %.5f\n error: %.5f\n',Tdr,Tdrs,Tdr_err)

    figure('Units','in','Position',[3 3 12 5]);
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact'); nexttile
    scatter(RT.grid.ai./pi,Rdf_a,80,'filled','s'); hold on; box on
    scatter(uvdh,Rvdh,80,'filled'); 
    xlabel('exiting angle, \alpha [\pi rad]','Interpreter','tex');
    ylabel('R_d(\alpha) [sr^{-1}]','Interpreter','tex');
    legend('This study','van de Hulst');
    set(gca,'TickDir','in','XMinorTick','on','YMinorTick','on')
    
    nexttile
    scatter(RT.grid.ai./pi,Tdf_a,80,'filled','s'); hold on; box on
    scatter(uvdh,Tvdh,80,'filled'); 
    xlabel('exiting angle, \alpha [\pi rad]','Interpreter','tex');
    ylabel('T_d(\alpha) [sr^{-1}]','Interpreter','tex');
    legend('This study','van de Hulst');
    set(gca,'TickDir','in','XMinorTick','on','YMinorTick','on')
    
elseif test_fluence == true
    
    % see Fig. 4 in Wang et al. 1995 for comparison
    figure;
    scatter(RT.grid.zi,RT.phi_z);
    set(gca,'YScale','log','YLim',[0.5 10],'XLim',[0 1]);
    xlabel('z (cm)');
    ylabel('Fluence (-)');
end

% see Vol. 1, pg. 262, Table 12 for isotropic

% figure; plot(Adf_z); hold on; plot(Adr_z) 
% figure; plot(zi,sum(phi_rz,2)); set(gca,'YScale','log')
% figure; plot(zi,phi_z); set(gca,'YScale','log')