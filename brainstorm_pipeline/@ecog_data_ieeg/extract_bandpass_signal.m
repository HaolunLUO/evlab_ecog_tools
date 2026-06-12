function [filtered_signal, filtered_signal_bipolar] = extract_bandpass_signal(obj, low_cutoff, high_cutoff)
    % EXTRACT_BANDPASS_SIGNAL Bandpass filter (using eegfilt) for stitched ecog_data
    % Inputs:
    %   obj         - ecog_data object
    %   low_cutoff  - low cutoff frequency (Hz)
    %   high_cutoff - high cutoff frequency (Hz)
    % Outputs:
    %   filtered_signal         - bandpass filtered unipolar signal (time x channels)
    %   filtered_signal_bipolar - bandpass filtered bipolar signal (time x channels), if present

    signal = double(obj.elec_data'); % time x channels
    if ~isempty(obj.bip_elec_data)
        signal_bipolar = obj.bip_elec_data'; % time x channels
    else
        signal_bipolar = [];
    end

    sample_freq = obj.sample_freq;
    stitch_index = obj.stitch_index(:); % ensure column vector

    filtered_signal = nan(size(signal));
    if ~isempty(signal_bipolar)
        filtered_signal_bipolar = nan(size(signal_bipolar));
    else
        filtered_signal_bipolar = [];
    end

    for k = 1:length(stitch_index)
        fprintf(1, '\n> Bandpass filtering segment %d of %d ... \n', k, length(stitch_index));

        if k == length(stitch_index)
            stop = size(signal,1);
        else
            stop = stitch_index(k+1)-1;
        end
        seg_idx = stitch_index(k):stop;

        % --- UNIPOLAR ---
        segment = signal(seg_idx, :); % time x channels
        filtered_segment = eegfilt(segment', sample_freq, low_cutoff, high_cutoff, 0, 200)'; % keep time x channels
        filtered_signal(seg_idx, :) = filtered_segment;

        % --- BIPOLAR ---
        if ~isempty(signal_bipolar)
            segment_bip = signal_bipolar(seg_idx, :);
            filtered_segment_bip = eegfilt(segment_bip', sample_freq, low_cutoff, high_cutoff,0, 200)';
            filtered_signal_bipolar(seg_idx, :) = filtered_segment_bip;
        end
    end

    % Optionally, assign back to object if desired:
    obj.elec_data = filtered_signal';
    if ~isempty(filtered_signal_bipolar)
        obj.bip_elec_data = filtered_signal_bipolar';
    end
end
