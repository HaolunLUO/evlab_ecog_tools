function zscore_signal(obj)
    % Zscores signal using default MATLAB zscore function.

    signal = obj.elec_data';
    if ~isempty(obj.bip_elec_data)
        signal_bipolar = obj.bip_elec_data';
    end

    for k=1:length(obj.stitch_index) % number of separate data files with signal
        fprintf(1, '\n> Computing z-score of signal from file %d of %d ... \n',k,length(obj.stitch_index));


        if k == length(obj.stitch_index) % signal for file stops at end of matrix
            stop = size(signal,1);
        else % signal for file stops before stitch index of next file
            stop = obj.stitch_index(k+1)-1;
        end

        % --- UNIPOLAR ---
        signal_ = signal(obj.stitch_index(k):stop,:);
        signal_ = zscore(signal_);
        signal(obj.stitch_index(k):stop,:) = signal_;
    
        % --- BIPOLAR ---
        if ~isempty(obj.bip_elec_data)
            signal_bipolar_ = signal_bipolar(obj.stitch_index(k):stop,:);
            signal_bipolar_ = zscore(signal_bipolar_);
            signal_bipolar(obj.stitch_index(k):stop,:) = signal_bipolar_;
        end

    end
    
    obj.elec_data = signal';

    if ~isempty(obj.bip_elec_data)
        obj.bip_elec_data = signal_bipolar';
    end


end