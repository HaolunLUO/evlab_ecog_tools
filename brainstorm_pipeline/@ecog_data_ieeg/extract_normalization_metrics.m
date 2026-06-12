%% Modified extract_normalization_metrics
function extract_normalization_metrics(obj, epoch_args)
    arguments
        obj ecog_data_ieeg
        epoch_args.baseTimeRange = [-0.5 0];
        epoch_args.timepad = 0.5;
    end

    [baseData, baseData_bip] = obj.extract_trial_epochs(epoch_tw=epoch_args.baseTimeRange);
    
    % Handle low trial count for unipolar data
    if size(baseData, 2) < 11
        [all_trials, total_len] = get_all_trials(obj);
        baseData = sample_windows(obj.elec_data, all_trials, total_len, obj.sample_freq, 1000, 'outside');
    end

    fprintf(1, '\n>> Extracting normalization metrics for unipolar high gamma envelope\n');
    normFactor = zeros(size(baseData, 1), 2);
    fprintf(1, '[');
    for iChan = 1:size(baseData, 1)
        normFactor(iChan, :) = [mean2(squeeze(baseData(iChan, :, :))), std2(squeeze(baseData(iChan, :, :)))];
        fprintf(1, '.');
    end
    fprintf(1, '] done\n');
    normMetrics.normFactor = normFactor;

    % Handle bipolar data
    if ~isempty(baseData_bip)
        if size(baseData_bip, 2) < 10
            [all_trials, total_len] = get_all_trials(obj);
            baseData_bip = sample_windows(obj.bip_elec_data, all_trials, total_len, obj.sample_freq, 1000, 'outside');
        end

        fprintf(1, '\n>> Extracting normalization metrics for bipolar high gamma envelope\n');
        normFactor = zeros(size(baseData_bip, 1), 2);
        [~, goodtrials] = remove_bad_trials(baseData_bip);
        goodTrialsCommon = extractCommonTrials(goodtrials);
        fprintf(1, '[');
        for iChan = 1:size(baseData_bip, 1)
            normFactor(iChan, :) = [mean2(squeeze(baseData_bip(iChan, goodTrialsCommon, :))), std2(squeeze(baseData_bip(iChan, goodTrialsCommon, :)))];
            fprintf(1, '.');
        end
        fprintf(1, '] done\n');
        normMetrics.normFactor_bip = normFactor;
    end

    obj.stats.normMetrics = normMetrics;
end

%% New Helper Functions
function [all_trials, total_len] = get_all_trials(obj)
    % Extract all trial boundaries from trial_timing
    all_trials = [];
    for k = 1:numel(obj.trial_timing)
        trial_table = obj.trial_timing{k};
        for t = 1:height(trial_table)
            new_trial.start = trial_table.start(t);
            new_trial.end = trial_table.end(t);
            all_trials = [all_trials; new_trial];
        end
    end
    total_len = size(obj.elec_data, 2);
end

function data_out = sample_windows(data, trials, total_len, sample_freq, n_windows, mode)
    window_samples = round(0.5 * sample_freq);
    n_channels = size(data, 1);
    
    % Get valid start indices
    valid_starts = get_valid_starts(trials, total_len, window_samples, mode);
    
    if numel(valid_starts) < n_windows
        error('Insufficient data for %d windows (available: %d)', n_windows, numel(valid_starts));
    end
    
    start_indices = randsample(valid_starts, n_windows);
    data_out = zeros(n_channels, n_windows, window_samples);
    for i = 1:n_windows
        data_out(:, i, :) = data(:, start_indices(i):start_indices(i)+window_samples-1);
    end
end

function valid_starts = get_valid_starts(trials, total_len, window_samples, mode)
    valid_starts = [];
    
    if strcmp(mode, 'outside')
        % Sort trials chronologically
        [~, order] = sort([trials.start]);
        trials = trials(order);
        
        % Before first trial
        if trials(1).start > 1
            valid_starts = [valid_starts, 1:(trials(1).start - window_samples)];
        end
        
        % Between trials
        for i = 1:length(trials)-1
            gap_start = trials(i).end + 1;
            gap_end = trials(i+1).start - 1;
            if gap_end >= gap_start + window_samples - 1
                valid_starts = [valid_starts, gap_start:(gap_end - window_samples + 1)];
            end
        end
        
        % After last trial
        if trials(end).end < total_len - window_samples + 1
            valid_starts = [valid_starts, (trials(end).end + 1):(total_len - window_samples + 1)];
        end
    else % 'inside' mode
        for t = 1:length(trials)
            trial = trials(t);
            if trial.end - trial.start + 1 >= window_samples
                valid_starts = [valid_starts, trial.start:(trial.end - window_samples + 1)];
            end
        end
    end
end