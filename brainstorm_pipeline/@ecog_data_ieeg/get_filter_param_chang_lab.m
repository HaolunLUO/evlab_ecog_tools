function [cfs,sds]=get_filter_param_chang_lab(obj,f_low,f_high)
    % TODO - description

    % Standard signal bands for neuroscience
    %   bands = ['theta', 'alpha', 'beta', 'gamma', 'high gamma']
    %   min_freqs = [4., 8., 15., 30., 70.]
    %   max_freqs = [7., 14., 29., 59., 150.]
    %   HG_freq = 200.

    cfs_round_factor = 1;
    cfs_round_val = 10^(cfs_round_factor);

    sds_round_factor = 2;
    sds_round_val = 10^(sds_round_factor);
    
    fq_min = 4.0749286538265;
    fq_max = 200.;
    scale = 7.;

    cfs = 2.^((log2(fq_min) * scale:1: log2(fq_max) * scale) / scale);
    sds = 10.^( log10(.39) + .5 * (log10(cfs)));
    sds = (sds) * sqrt(2.);

    sds=round(sds*sds_round_val)/sds_round_val;
    cfs=round(cfs*cfs_round_val)/cfs_round_val;

    if nargin<2
        sds=sds;
        cfs=cfs;
    else 
        index=(cfs<f_high) & (cfs>f_low);
        cfs=cfs(index);
        sds=sds(index);
    end 

end