function [wavl,asym,omeg,clen,kabs] = mcrt_load_optics(opts,iops)
% need to add dE options 
if opts.use_dE == true
    use_iops.mie = iops.mie_dE;
    use_iops.ssa = iops.ssa_dE;
else
    use_iops.mie = iops.mie;
    use_iops.ssa = iops.ssa;
end
    
switch opts.use_kabs
    case 'kice'                                 % norod kice
        if opts.use_ssa == true
            wavl    = use_iops.ssa.kice.wavl;
            asym    = use_iops.ssa.kice.asym;
            omeg    = use_iops.ssa.kice.omeg;
            clen    = use_iops.ssa.kice.clen;
            kabs    = use_iops.ssa.kice.kabs;
        else
            wavl    = use_iops.mie.kice.wavl;
            asym    = use_iops.mie.kice.asym;
            omeg    = use_iops.mie.kice.omeg;
            clen    = use_iops.mie.kice.clen;
            kabs    = use_iops.mie.kice.kabs;
        end
    case 'kabs'
        if opts.use_ssa == true
            wavl    = use_iops.ssa.kabs.wavl;
            asym    = use_iops.ssa.kabs.asym;
            omeg    = use_iops.ssa.kabs.omeg;
            clen    = use_iops.ssa.kabs.clen;
            kabs    = use_iops.ssa.kabs.kabs;
        else
            wavl    = use_iops.mie.kabs.wavl;
            asym    = use_iops.mie.kabs.asym;
            omeg    = use_iops.mie.kabs.omeg;
            clen    = use_iops.mie.kabs.clen;
            kabs    = use_iops.mie.kabs.kabs;
        end
end