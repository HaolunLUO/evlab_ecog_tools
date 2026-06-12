function [envelopes_rjct,outliers_idx,outlier_prcnt,ignored_channels]=envelope_outliers(obj,envelopes,ops,varargin)
    % Detects outliers in signal envelopes. The difference between the median
    % and the Qth percentile (default: 90%) of the envelope distribution for
    % each electrode is measured. Envelopes are considered outliers if they
    % fall a certain number of units above this value.
    % 
    % 2016-08-15 - Created, Sam NH
    % 
    % 2019-01-21 - Cosmetic changes, Sam NH
    %
    % 2021-11-29 - update to fit with Evlab ecog pipeline 

    n_channels = size(envelopes,2);

    % make threshold parameters electrode specific if not already
    if isvector(ops.threshold)
        ops.threshold = repmat(ops.threshold(:), 1, n_channels);
    else
        assert(size(ops.threshold) == n_channels);
    end

    % % trim first and last bit of envelope
    % samples_to_remove = obj.sample_freq * 1;
    % envelopes = envelopes(samples_to_remove:size(envelopes,1)-samples_to_remove,:);

    % normalize by 90% of the distribution
    s = quantile(envelopes, ops.percentile);
    Z = envelopes ./ repmat(s, size(envelopes,1),1);
    clear s;

    outliers = Z > ops.threshold;
    envelopes_rjct = nan*envelopes;
    outlier_prcnt = nan*zeros(size(envelopes,2),1);
    outliers_idx = {};
    
    ignored_channels = [];

    % loop through channels
    for q=1:n_channels
        outlie = outliers(:,q);
        outlier_prcnt(q) = sum(outlie)./numel(outlie)*1e2;
        
        envl = envelopes(:,q);
        envl_rejct = envl;
        
        rise = find(diff(outlie)==1).';
        fall = find(diff(outlie)==-1).';

        try 
            assert(numel(rise)==numel(fall));
            process = true;
        catch err 
            ignored_channels = [ignored_channels, q];
            process = false;
        end 

        outlie_idx ={ };
        if numel(rise)>0 & process
            for r=rise % go through each outlier
                
                f = fall(find(fall>r,1,'first'));
                
                % find window to look at for interpolation
                interp_win = r:f;
                interp_val = envl(interp_win).';
                lookup_window = max([r-ops.buffer,1]):min([(f+ops.buffer),numel(outlie)]);
                main_sig = envl(lookup_window);
                [C,ia,ib] = intersect(lookup_window,interp_win);

                aux = main_sig;
                fixed_sig = main_sig;
                
                % interpolate
                aux(ia) = [];
                intrp_sig = interp1(setdiff(lookup_window,interp_win,'stable'),aux,interp_win,ops.interpMethod);
                fixed_sig(ia) = intrp_sig;

                % add interpolated values to signal
                envl_rejct(lookup_window) = fixed_sig;
                outlie_idx = [outlie_idx,[interp_win;interp_val]];
            end 
        end

        envelopes_rjct(:,q) = envl_rejct;
        outliers_idx = [outliers_idx;{outlie_idx}];
        
        fprintf(1,'.');

    end 

end