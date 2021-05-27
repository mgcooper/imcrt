function opts = mcrt_set_case(opts,case_flag)

% 'kice' is pure ice absorptivity from Warren et al. 2006
% 'kabs' is the in-situ ice absorptivity from Cooper et al. 2021
% the options below setup eight experiments that test the two
% k options, with and without detector rod interference, and 
% with and without the empirical angular response of the rcr.
% see Cooper et al. 2021 for a full description.

switch case_flag
    case 1                              % no rod, no rcr, kice (ideal case)
        opts.use_kabs   = 'kice';
        opts.test_rod   = false;
        opts.test_rcr   = false;
        opts.ideal_case = true;
    case 2                              % no rod, ideal rcr, kice
        opts.use_kabs   = 'kice';
        opts.test_rod   = false;
        opts.test_rcr   = false;
        opts.ideal_case = false;
    case 3
        opts.use_kabs   = 'kice';       % rod, ideal rcr, kice
        opts.test_rod   = true;
        opts.test_rcr   = false;
        opts.ideal_case = false;
    case 4                              % rod, non-ideal rcr, kice
        opts.use_kabs   = 'kice';
        opts.test_rod   = true;
        opts.test_rcr   = true;
        opts.ideal_case = false;
	case 5                              % no rod, no rcr, kabs (ideal case)
        opts.use_kabs   = 'kabs';
        opts.test_rod   = false;
        opts.test_rcr   = false;
        opts.ideal_case = true;
    case 6                              % no rod, ideal rcr, kabs
        opts.use_kabs   = 'kabs';
        opts.test_rod   = false;
        opts.test_rcr   = false;
        opts.ideal_case = false;
    case 7
        opts.use_kabs   = 'kabs';       % rod, ideal rcr, kabs
        opts.test_rod   = true;
        opts.test_rcr   = false;
        opts.ideal_case = false;        
    case 8                              % rod, non-ideal rcr, kabs
        opts.use_kabs   = 'kabs';
        opts.test_rod   = true;
        opts.test_rcr   = true;
        opts.ideal_case = false;
end
