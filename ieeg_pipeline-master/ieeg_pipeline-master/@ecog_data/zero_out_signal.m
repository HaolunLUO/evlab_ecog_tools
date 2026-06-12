%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function signal=zero_out_signal(obj,signal)
    % Zeros out first couple seconds of each file. 
    % 
    % Returns : signal

    samples_to_remove = obj.sample_freq*obj.for_preproc.zero;

    for k=1:length(obj.stitch_index) % number of separate data files with signal

        if k == length(obj.stitch_index) % signal for file stops at end of matrix
            stop = size(signal,1);
        else % signal for file stops before stitch index of next file
            stop = obj.stitch_index(k+1)-1;
        end

        signal_ = signal(obj.stitch_index(k):stop,:);
        mean_signal = mean(signal_,1);
        signal_(1:samples_to_remove,:) = zeros(samples_to_remove,size(signal_,2));
        signal_((end-samples_to_remove+1):end,:) = zeros(samples_to_remove,size(signal_,2));
        signal(obj.stitch_index(k):stop,:) = signal_;

    end

end