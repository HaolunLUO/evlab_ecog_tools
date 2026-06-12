function notch_filter(obj)
    % Notch filters signal and marks channels with significant noise

    signal = obj.elec_data(obj.elec_ch_clean,:)';

    notch = obj.for_preproc.notch;

    % measure line noise prior to notch filtering
    signal_noise_before = obj.measure_line_noise(obj.elec_data');

    for k=1:length(obj.stitch_index) % number of separate data files with signal
        fprintf(1, '\n> Notch filtering signal from file %d of %d ... \n',k,length(obj.stitch_index));
        fprintf(1,'[');

        if k == length(obj.stitch_index) % signal for file stops at end of matrix
            stop = size(signal,1);
        else % signal for file stops before stitch index of next file
            stop = obj.stitch_index(k+1)-1;
        end

        signal_ = signal(obj.stitch_index(k):stop,:);

        % notch filter
        for idx_channel=1:size(signal_,2)
            signal_preliminary = double(signal_(:,idx_channel));
            % remove all harmonics of line-noise
            for idx = 1:length(obj.for_preproc.filter_params.notch.fcenter) %#ok<PFBNS>
                signal_preliminary = filtfilt(notch{idx}.b,notch{idx}.a,signal_preliminary); %#ok<PFBNS>
            end 
            signal_(:,idx_channel) = signal_preliminary;
            fprintf(1,'.');
        end

        signal(obj.stitch_index(k):stop,:) = signal_;

        fprintf(1,'] done\n');
    end

    obj.elec_data(obj.elec_ch_clean,:) = signal';

    % measure line noise after notch filtering
    signal_noise_after = obj.measure_line_noise(obj.elec_data');

    % % mark channels with sig line noise after notch filtering
    % obj.elec_ch_with_noise = obj.elec_ch(signal_noise_after(:,2) > (mean(signal_noise_after(:,2))+5*std(signal_noise_after(:,2))));
    % obj.elec_ch_with_noise = intersect(obj.elec_ch_clean,obj.elec_ch_with_noise); % don't mark already noisy electrodes
    % obj.define_clean_channels();

    obj.for_preproc.notchFilter_results.signal_noise_before_notch = signal_noise_before;
    obj.for_preproc.notchFilter_results.mean_signal_noise_before_notch = mean(signal_noise_before(obj.elec_ch_clean,2));
    obj.for_preproc.notchFilter_results.signal_noise_after_notch = signal_noise_after;
    obj.for_preproc.notchFilter_results.mean_signal_noise_after_notch = mean(signal_noise_after(obj.elec_ch_clean,2));

    fprintf(1,'\nReduced 60 Hz noise from %.2f to %.2f uV\n',mean(signal_noise_before(obj.elec_ch_clean,2)),mean(signal_noise_after(obj.elec_ch_clean,2)));
    % fprintf(1,'Electrodes with significant line noise: ');
    % fprintf(1,'%d ', obj.elec_ch_with_noise(:)); fprintf('\n');

    % % plot line noise to select additional channels to remove
    % if obj.for_preproc.isPlotVisible
    %     size(signal_noise_before)
    %     size(signal_noise_after)
    %     f = obj.plot_line_noise(signal_noise_before,signal_noise_after);
    %     prompt1 = '\nUSER INPUT REQUIRED: \nAdditional channels to remove due to significant line noise? (format: [1,2]) - ';
    % 
    %     more_line_noise = input(prompt1)';
    %     obj.elec_ch_with_noise = union(obj.elec_ch_with_noise,more_line_noise);
    %     obj.define_clean_channels();
    %     close(f);
    % end
end