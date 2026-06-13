
function test_s_vs_n_wordboundaries(obj, varargin)
    % WORD BOUNDARIES TIME-SERIES PERMUTATION TEST to identify S > N electrodes
    % Based on the word boundaries analysis in generateReportLangloc
    % Performs timePermCluster test for sentence vs nonword comparison across word positions

    p = inputParser();
    addParameter(p, 'S_condition_flag', 'sentence');
    addParameter(p, 'N_condition_flag', 'nonword');
    addParameter(p, 'n_rep', 1000); % numTail parameter for timePermCluster
    addParameter(p, 'threshold', 0.05);
    addParameter(p, 'sessions', []);
    addParameter(p, 'do_plot', false);
    addParameter(p, 'epoch_range', [-0.5, 0.5]); % Time range around each word
    addParameter(p, 'num_words', 12); % Number of word positions
    addParameter(p,'use_odd_for_inference',true); % NEW
    parse(p, varargin{:});
    ops = p.Results;

    obj.s_vs_n_wordboundaries_ops = ops;

    fprintf(1, '\n> Running word boundaries time-series permutation test for %s vs %s ...\n', ...
            ops.S_condition_flag, ops.N_condition_flag);

    % Validate that required conditions exist
    availableConditions = unique(obj.condition);
    if ~ismember(ops.S_condition_flag, availableConditions) || ~ismember(ops.N_condition_flag, availableConditions)
        error('Required conditions "%s" and "%s" not found. Available conditions: %s', ...
              ops.S_condition_flag, ops.N_condition_flag, strjoin(availableConditions, ', '));
    end

    % Filter sessions if specified
    if ~isempty(ops.sessions)
        keep_trials = any(obj.session == ops.sessions, 2);
    else
        keep_trials = true(size(obj.condition));
    end

    % Initialize statistical analysis storage
    pSig = cell(size(obj.elec_ch_label, 1), ops.num_words);
    pSig_bip = cell(size(obj.bip_ch_label, 1), ops.num_words);

    % Initialize concatenated epochs storage
    concatenatedEpochsSentence = [];
    concatenatedEpochsNonword = [];
    concatenatedEpochsSentence_bip = [];
    concatenatedEpochsNonword_bip = [];

    fprintf(1, '> Processing word boundaries for %s experiment...\n', obj.experiment);

    % Loop through each word position to extract and analyze data
    for wordPos = 1:ops.num_words
        fprintf(1, '  Processing word position %d/%d...\n', wordPos, ops.num_words);

        % Extract epochs for current word position
        [epochData, epochData_bip] = obj.extract_trial_epochs('epoch_tw', ops.epoch_range, ...
                                                              'probe_key', wordPos+1);

        % Apply session filter and extract sentence/nonword trials
        sentence_trials_idx = strcmp(obj.condition, ops.S_condition_flag) & keep_trials;
        nonword_trials_idx = strcmp(obj.condition, ops.N_condition_flag) & keep_trials;

        % --- ODD/EVEN SPLIT ---
        if ops.use_odd_for_inference
            % Convert logical index -> integer indices, keep only odd ones
            s_idx = find(sentence_trials_idx);
            n_idx = find(nonword_trials_idx);
            sentence_trials_idx = false(size(sentence_trials_idx));
            nonword_trials_idx  = false(size(nonword_trials_idx));
            sentence_trials_idx(s_idx(1:2:end)) = true;
            nonword_trials_idx(n_idx(1:2:end))  = true;
        end

        sentenceTrials = epochData(:, sentence_trials_idx, :);
        nonwordTrials = epochData(:, nonword_trials_idx, :);

        % Concatenate along time axis (dimension 3)
        if isempty(concatenatedEpochsSentence)
            concatenatedEpochsSentence = sentenceTrials;
            concatenatedEpochsNonword = nonwordTrials;
        else
            concatenatedEpochsSentence = cat(3, concatenatedEpochsSentence, sentenceTrials);
            concatenatedEpochsNonword = cat(3, concatenatedEpochsNonword, nonwordTrials);
        end

        % Perform timePermCluster test for significance (unipolar)
        fprintf(1, '    Running permutation tests for unipolar channels...\n');
        parfor iChan = 1:size(epochData, 1)
            aTrialData = squeeze(epochData(iChan, sentence_trials_idx, :));
            bTrialData = squeeze(epochData(iChan, nonword_trials_idx, :));

            if size(aTrialData, 1) > 1 && size(bTrialData, 1) > 1
                try
                    pSig{iChan, wordPos} = timePermCluster(aTrialData, bTrialData, 'nPerm', ops.n_rep, 'statstype','corr','pThresh',ops.threshold);
                catch ME
                    warning('TimePermCluster failed for unipolar channel %d, word %d: %s', ...
                           iChan, wordPos, ME.message);
                    % Handle insufficient trials or other errors
                    pSig{iChan, wordPos} = struct('h_sig_05', zeros(1, size(epochData, 3)), ...
                                                 'p_val', ones(1, size(epochData, 3)));
                end
            else
                % Handle insufficient trials
                pSig{iChan, wordPos} = struct('h_sig_05', zeros(1, size(epochData, 3)), ...
                                             'p_val', ones(1, size(epochData, 3)));
            end
        end

        % Perform timePermCluster test for bipolar data if available
        if ~isempty(epochData_bip)
            fprintf(1, '    Running permutation tests for bipolar channels...\n');

            % Handle bipolar concatenation
            sentenceTrials_bip = epochData_bip(:, sentence_trials_idx, :);
            nonwordTrials_bip = epochData_bip(:, nonword_trials_idx, :);

            if isempty(concatenatedEpochsSentence_bip)
                concatenatedEpochsSentence_bip = sentenceTrials_bip;
                concatenatedEpochsNonword_bip = nonwordTrials_bip;
            else
                concatenatedEpochsSentence_bip = cat(3, concatenatedEpochsSentence_bip, sentenceTrials_bip);
                concatenatedEpochsNonword_bip = cat(3, concatenatedEpochsNonword_bip, nonwordTrials_bip);
            end

            parfor iChan = 1:size(epochData_bip, 1)
                aTrialData = squeeze(epochData_bip(iChan, sentence_trials_idx, :));
                bTrialData = squeeze(epochData_bip(iChan, nonword_trials_idx, :));

                if size(aTrialData, 1) > 1 && size(bTrialData, 1) > 1
                    try
                        pSig_bip{iChan, wordPos} = timePermCluster(aTrialData, bTrialData, 'nPerm', ops.n_rep, 'statstype','corr','pThresh',ops.threshold);
                    catch ME
                        warning('TimePermCluster failed for bipolar channel %d, word %d: %s', ...
                               iChan, wordPos, ME.message);
                        pSig_bip{iChan, wordPos} = struct('h_sig_05', zeros(1, size(epochData_bip, 3)), ...
                                                         'p_val', ones(1, size(epochData_bip, 3)));
                    end
                else
                    % Handle insufficient trials
                    pSig_bip{iChan, wordPos} = struct('h_sig_05', zeros(1, size(epochData_bip, 3)), ...
                                                     'p_val', ones(1, size(epochData_bip, 3)));
                end
            end
        end
    end

    % Calculate word boundary parameters
    timePointsPerWord = size(concatenatedEpochsSentence, 3) / ops.num_words;
    totalTimePoints = timePointsPerWord * ops.num_words;
    wordBoundaries = 0:timePointsPerWord:totalTimePoints;

    % Store statistical results in object
    obj.s_vs_n_wordboundaries_pSig = pSig;
    obj.s_vs_n_wordboundaries_wordBoundaries = wordBoundaries;
    obj.s_vs_n_wordboundaries_timePointsPerWord = timePointsPerWord;
    obj.s_vs_n_wordboundaries_totalTimePoints = totalTimePoints;

    if ~isempty(concatenatedEpochsSentence_bip)
        obj.s_vs_n_wordboundaries_pSig_bip = pSig_bip;
    end

    % Create significance summary for unipolar channels
    sigChannels = sum(cellfun(@(x) any(cell2mat(arrayfun(@(y) any(y.h_sig_05), x, 'UniformOutput', false))), pSig'));
    sigUnipolarChannels = find(sigChannels > 0);

    fprintf(1, '> Found %d/%d unipolar channels with significant word boundary clusters\n', ...
            length(sigUnipolarChannels), size(concatenatedEpochsSentence, 1));

    % Create significance summary for bipolar channels if available
    if ~isempty(concatenatedEpochsSentence_bip)
        sigChannels_bip = sum(cellfun(@(x) any(cell2mat(arrayfun(@(y) any(y.h_sig_05), x, 'UniformOutput', false))), pSig_bip'));
        sigBipolarChannels = find(sigChannels_bip > 0);

        fprintf(1, '> Found %d/%d bipolar channels with significant word boundary clusters\n', ...
                length(sigBipolarChannels), size(concatenatedEpochsSentence_bip, 1));

        obj.s_vs_n_wordboundaries_sigBipolarChannels = sigBipolarChannels;
    end

    obj.s_vs_n_wordboundaries_sigUnipolarChannels = sigUnipolarChannels;

    % Store concatenated data for plotting
    obj.s_vs_n_wordboundaries_concatenatedEpochsSentence = concatenatedEpochsSentence;
    obj.s_vs_n_wordboundaries_concatenatedEpochsNonword = concatenatedEpochsNonword;

    if ~isempty(concatenatedEpochsSentence_bip)
        obj.s_vs_n_wordboundaries_concatenatedEpochsSentence_bip = concatenatedEpochsSentence_bip;
        obj.s_vs_n_wordboundaries_concatenatedEpochsNonword_bip = concatenatedEpochsNonword_bip;
    end

    % Generate plots if requested
    if ops.do_plot
        obj.plot_s_vs_n_wordboundaries(ops);
    end

    fprintf(1, '> Word boundaries time-series permutation test completed\n');
end
