function opts = mcrt_set_path(opts,geom)

if opts.use_ssa == true
    if opts.use_dE == true
        opts.path.case_save  = [opts.path.save 'ssa_dE/' int2str(geom.Z) 'cm/'];
    else
        opts.path.case_save  = [opts.path.save 'ssa/' int2str(geom.Z) 'cm/'];
    end
else
    if opts.use_dE == true
        opts.path.case_save  = [opts.path.save 'mie_dE/' int2str(geom.Z) 'cm/'];
    else
        opts.path.case_save  = [opts.path.save 'mie/' int2str(geom.Z) 'cm/'];
    end
end
        
if ~exist(opts.path.case_save,'dir')
    mkdir(opts.path.case_save); addpath(opts.path.case_save);
end
    
if opts.save_data == true
    if opts.ideal_case == true
        opts.fsave.opts = [opts.path.case_save 'opts_ideal_' opts.use_kabs];
        opts.fsave.data = [opts.path.case_save 'data_ideal_' opts.use_kabs];
        return        
    elseif opts.test_rod == true
        if opts.test_rcr == true
            opts.fsave.opts = [opts.path.case_save 'opts_rcr_' opts.use_kabs];
            opts.fsave.data = [opts.path.case_save 'data_rcr_' opts.use_kabs];
        else
            opts.fsave.opts = [opts.path.case_save 'opts_rod_' opts.use_kabs];
            opts.fsave.data = [opts.path.case_save 'data_rod_' opts.use_kabs];
        end
    else
        opts.fsave.opts = [opts.path.case_save 'opts_norod_' opts.use_kabs];
        opts.fsave.data = [opts.path.case_save 'data_norod_' opts.use_kabs];
    end
end