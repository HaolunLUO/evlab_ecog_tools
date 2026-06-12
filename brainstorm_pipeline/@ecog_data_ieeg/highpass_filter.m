function highpass_filter(obj)
    % Highpass filters signal
    %
    %

    signal = obj.elec_data(obj.elec_ch_clean,:)';

    highpass = obj.for_preproc.highpass;

    for k=1:length(obj.stitch_index) % number of separate data files with signal
        fprintf(1, '\n> Highpass filtering signal from file %d of %d ... \n',k,length(obj.stitch_index));
        fprintf(1,'[');

        if k == length(obj.stitch_index) % signal for file stops at end of matrix
            stop = size(signal,1);
        else % signal for file stops before stitch index of next file
            stop = obj.stitch_index(k+1)-1;
        end

        signal_ = signal(obj.stitch_index(k):stop,:);

        % highpass filter signal
        parfor idx_channel=1:size(signal_,2)
            warning('off', 'signal:filtfilt:ParseSOS');
            signal_(:,idx_channel) = filtfilt(highpass.sos,highpass.g,double(signal_(:,idx_channel)));
            fprintf(1,'.');
        end

        signal(obj.stitch_index(k):stop,:) = signal_;

        fprintf(1,'] done\n');
    end

    obj.elec_data(obj.elec_ch_clean,:) = signal';

end