function remove_outliers(obj)
    % Removes outliers from envelope.
    %
    % CANNOT be done as the first preprocessing step. MUST be done 
    % after high gamma extraction (need envelopes).

    ops.trimmed      = obj.for_preproc.outlier.trimmed;
    ops.threshold    = obj.for_preproc.outlier.threshold;
    ops.percentile   = obj.for_preproc.outlier.percentile;
    ops.buffer       = obj.for_preproc.outlier.buffer;
    ops.interpMethod = obj.for_preproc.outlier.interpMethod;

    samples_to_remove = ops.trimmed * obj.sample_freq;

    signal = obj.elec_data';
    signal = obj.zero_out_signal(signal);

    outlierRemoval_results.idxs    = [];
    outlierRemoval_results.prcnts  = [];
    outlierRemoval_results.ignored = [];

    if ~isempty(obj.bip_elec_data)
        signal_bipolar = obj.bip_elec_data';
        signal_bipolar = obj.zero_out_signal(signal_bipolar);

        outlierRemoval_results.idxs_bipolar    = [];
        outlierRemoval_results.prcnts_bipolar  = [];
        outlierRemoval_results.ignored_bipolar = [];
    end
    

    for k=1:length(obj.stitch_index) % number of separate data files with signal
        fprintf(1, '\n> Removing outliers from file %d of %d ... \n',k,length(obj.stitch_index));
        fprintf(1,'[');

        if k == length(obj.stitch_index) % signal for file stops at end of matrix
            stop = size(signal,1);
        else % signal for file stops before stitch index of next file
            stop = obj.stitch_index(k+1)-1;
        end

        % --- UNIPOLAR ---
        signal_ = signal(obj.stitch_index(k):stop,:);
        [signal_,idxs,prcnts,ignored] = obj.envelope_outliers(signal_,ops);
        signal(obj.stitch_index(k):stop,:) = signal_;
        fprintf(1,'] done\n');

        outlierRemoval_results.idxs = [outlierRemoval_results.idxs, idxs];
        outlierRemoval_results.prcnts = [outlierRemoval_results.prcnts, prcnts];
        outlierRemoval_results.ignored = [outlierRemoval_results.ignored; ismember(obj.elec_ch,ignored)'];
        
        ignored_not_noisy = intersect(obj.elec_ch_clean,ignored);
        fprintf(1,'Unequal rise and fall of outliers, ignored unipolar channels: ');
        fprintf(1,'%d ',ignored_not_noisy(:)); fprintf('\n');

        % --- BIPOLAR ---
        if ~isempty(obj.bip_elec_data)
            fprintf(1,'[');
            signal_bipolar_ = signal_bipolar(obj.stitch_index(k):stop,:);
            [signal_bipolar_,idxs_bipolar,prcnts_bipolar,ignored_bipolar] = obj.envelope_outliers(signal_bipolar_,ops);
            signal_bipolar(obj.stitch_index(k):stop,:) = signal_bipolar_;
            fprintf(1,'] done\n');

            outlierRemoval_results.idxs_bipolar = [outlierRemoval_results.idxs_bipolar, idxs_bipolar];
            outlierRemoval_results.prcnts_bipolar = [outlierRemoval_results.prcnts_bipolar, prcnts_bipolar];
            outlierRemoval_results.ignored_bipolar = [outlierRemoval_results.ignored_bipolar, ismember(obj.bip_ch,ignored_bipolar)'];

            ignored_bipolar_not_noisy = intersect(obj.bip_ch,ignored_bipolar);
            fprintf(1,'Unequal rise and fall of outliers, ignored bipolar channels: ');
            fprintf(1,'%d ',ignored_bipolar_not_noisy(:)); fprintf('\n');
        end

    end

    obj.elec_data = signal';

    if ~isempty(obj.bip_elec_data)
        obj.bip_elec_data = signal_bipolar';
    end

    obj.for_preproc.outlierRemoval_results = outlierRemoval_results;

    if obj.for_preproc.isPlotVisible
        obj.plot_channels(signal,...
                        obj.elec_ch_label,...
                        obj.elec_ch_clean,...
                        obj.elec_ch_valid,...
                        'stitch_index',obj.stitch_index,...
                        't_len',60,...
                        'sample_freq',obj.sample_freq...
        );
    end

end