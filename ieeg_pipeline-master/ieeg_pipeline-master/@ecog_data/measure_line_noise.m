function signal_noise=measure_line_noise(obj,signal)
    % 
    % 
    % Returns : measurement of signal noise

    fprintf(1, '\n> Measuring 60Hz noise power ...\n');
    fprintf(1,'[');

    peak = obj.for_preproc.peak;

    signal_noise = zeros(size(signal,2),length(peak));

    % calculate average root-mean-square of the line-noise
    for idx_channel=1:size(signal,2)
        % signal_noise(idx_channel,1) = mean(sqrt(filtfilt(peak.b,peak.a,double(signal(:,idx_channel))).^2),1); %#ok<PFBNS>
        for idx_filter=1:length(peak)
            signal_noise(idx_channel,idx_filter) = mean(abs(filter(peak{idx_filter}.b,peak{idx_filter}.a,signal(:,idx_channel))));
        end
        fprintf(1,'.');
    end

    fprintf(1,'] done\n');

end
