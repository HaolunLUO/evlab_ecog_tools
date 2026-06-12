%% Modified extract_significant_channel
function extract_significant_channel(obj, args)
    arguments
        obj ecog_data
        args.baseTime = [-0.5 0];
        args.epochTime = [0 0.5];
        args.numPerm = 10000;
        args.p_val = 0.05;
    end

    [baseData, baseData_bip] = obj.extract_trial_epochs(epoch_tw=args.baseTime);
    [epochData, epochData_bip] = obj.extract_trial_epochs(epoch_tw=args.epochTime);

    % Handle low trial count for unipolar data
    if size(baseData, 2) < 11 || size(epochData, 2) < 11
        [all_trials, total_len] = get_all_trials(obj);
        baseData = sample_windows(obj.elec_data, all_trials, total_len, obj.sample_freq, 1000, 'outside');
        epochData = sample_windows(obj.elec_data, all_trials, total_len, obj.sample_freq, 1000, 'inside');
    end

    basePower = mean(baseData.^2, 3);
    epochPower = mean(epochData.^2, 3);

    fprintf(1, '\n>> Extracting significant channels for unipolar high gamma envelope\n');
    fprintf(1, '[');
    [~, goodtrials] = remove_bad_trials(epochData);

    pSig = struct();
    for iChan = 1:size(basePower, 1)
        pSig.pChan(iChan) = permtest(...
            epochPower(iChan, goodtrials(iChan, :)), ...
            basePower(iChan, goodtrials(iChan, :)), ...
            args.numPerm);
        fprintf(1, '.');
    end

    pSig.h_fdr_05 = fdr_bh(pSig.pChan, 0.05);
    pSig.h_fdr_01 = fdr_bh(pSig.pChan, 0.01);
    fprintf(1, '] done\n');

    % Handle bipolar data
    if ~isempty(baseData_bip)
        if size(baseData_bip, 2) < 11 || size(epochData_bip, 2) < 11
            [all_trials, total_len] = get_all_trials(obj);
            baseData_bip = sample_windows(obj.bip_elec_data, all_trials, total_len, obj.sample_freq, 1000, 'outside');
            epochData_bip = sample_windows(obj.bip_elec_data, all_trials, total_len, obj.sample_freq, 1000, 'inside');
        end

        basePower_bip = mean(baseData_bip.^2, 3);
        epochPower_bip = mean(epochData_bip.^2, 3);

        fprintf(1, '\n>> Extracting significant channels for bipolar high gamma envelope\n');
        fprintf(1, '[');
        [~, goodtrials_bip] = remove_bad_trials(epochData_bip);

        for iChan = 1:size(basePower_bip, 1)
            pSig.pChan_bip(iChan) = permtest(...
                epochPower_bip(iChan, goodtrials_bip(iChan, :)), ...
                basePower_bip(iChan, goodtrials_bip(iChan, :)), ...
                args.numPerm);
            fprintf(1, '.');
        end

        pSig.h_bip_fdr_05 = fdr_bh(pSig.pChan_bip, 0.05);
        pSig.h_bip_fdr_01 = fdr_bh(pSig.pChan_bip, 0.01);
        fprintf(1, '] done\n');
    end

    obj.stats.sig_hg_channel = pSig;
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



function p = permtest(sample1, sample2, numperm)
% permtest - Perform one-sided permutation test to compare the means of two samples.
%
% Syntax: p = permtest(sample1, sample2, numperm)
%
% Inputs:
%   sample1     - First sample data (1 x n1) array
%   sample2     - Second sample data (1 x n2) array
%   numperm     - Number of permutations to perform
%
% Outputs:
%   p           - p-value indicating the significance of the difference between the means
%
% Example:
%   sample1 = [1, 2, 3, 4, 5]; % Example first sample
%   sample2 = [6, 7, 8, 9, 10]; % Example second sample
%   p = permtest(sample1, sample2, 1000); % Perform one-sided permutation test with 1000 permutations
% Author - Kumar Duraivel
samples = [sample1 sample2]; % Combine the two samples
samplediff = mean(sample1) - mean(sample2); % Calculate the difference between the means of the samples
sampdiffshuff = zeros(1, numperm); % Initialize an array to store shuffled sample differences

for n = 1:numperm
    sampshuff = samples(randperm(length(samples))); % Shuffle the combined samples
    sampdiffshuff(n) = mean(sampshuff(1:length(sampshuff)/2)) - mean(sampshuff(length(sampshuff)/2+1:end)); % Calculate the difference between means for the shuffled samples
end

p = length(find(sampdiffshuff > samplediff)) / numperm; % Calculate the p-value as the proportion of shuffled sample differences greater than the observed difference

end
