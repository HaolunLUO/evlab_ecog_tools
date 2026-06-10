function brainstorm_to_mit_crunched_new(allDataFiles, params)
    % Converts Brainstorm data to MIT's "crunched" format
    % Output can be directly used with analyze_MITLangloc.m
    
    %% =====================================================================
    % LOAD BRAINSTORM DATA
    % ======================================================================
    fprintf('\n=== LOADING BRAINSTORM DATA ===\n');
    
    % Initialize
    allSignal = [];
    allTrialTiming = {};
    allConditions = {};
    allSessions = [];
    stitch_index = [1];
    targetFs = 500
    
    % Get channel info from first file
    DataMat1 = in_bst_data(allDataFiles{1});
    [sStudy, ~] = bst_get('DataFile', allDataFiles{1});
    ChannelMat = in_bst_channel(sStudy.Channel.FileName);
    
    chanLabels = {ChannelMat.Channel.Name}';
    chanTypes = {ChannelMat.Channel.Type}';
    nChannels = length(chanLabels);
    
    % Identify valid SEEG channels
    validChanMask = contains(chanTypes, 'SEEG') & ~contains(chanTypes, 'NO_LOC');
    validChanIdx = find(validChanMask);
    
    fprintf('Found %d channels, %d are valid SEEG\n', nChannels, sum(validChanMask));
    
    % Load all runs
    for iRun = 1:length(allDataFiles)
        fprintf('Loading run %d/%d...\n', iRun, length(allDataFiles));
        
        DataMat = in_bst_data(allDataFiles{iRun});
        runSignal = DataMat.F;
        
        if iRun == 1
            sRate = round(1 / (DataMat.Time(2) - DataMat.Time(1)));
            fprintf('Sampling rate: %.2f Hz (rounded to integer)\n', sRate);
        end
        
        % Extract trials
        [timing, conds] = extract_trials_from_bst(DataMat, params, size(allSignal, 2), sRate);
        
        allSignal = [allSignal, runSignal];
        allTrialTiming = [allTrialTiming; timing];
        allConditions = [allConditions; conds];
        allSessions = [allSessions; repmat(iRun, length(conds), 1)];
        
        if iRun < length(allDataFiles)
            stitch_index = [stitch_index; size(allSignal, 2) + 1];
        end
    end
    
    fprintf('Total: [%d × %d], %d trials\n', size(allSignal), length(allTrialTiming));
    
    
    %% =====================================================================
    % BUILD MIT-COMPATIBLE OBJECT
    % ======================================================================
    fprintf('\n=== CREATING MIT ecog_data_seeg OBJECT ===\n');
    
    % Compute downsampled stitch_index (per-run length mapping)
    runStarts = stitch_index(:);
    runStops  = [stitch_index(2:end)-1; size(allSignal,2)];
    runLensRaw = runStops - runStarts + 1;
    runLensDec = ceil(runLensRaw * targetFs / sRate);
    stitch_index_dec = [1; 1 + cumsum(runLensDec(1:end-1))];
    
    % Compute downsampled trial timing (per run, to keep offsets consistent)
    trialTiming_dec = allTrialTiming;
    for i = 1:numel(allTrialTiming)
        run = allSessions(i);              % which run this trial came from
        T = allTrialTiming{i};
    
        rawOffset = stitch_index(run) - 1;       % raw samples before this run
        decOffset = stitch_index_dec(run) - 1;   % dec samples before this run
    
        % map sample indices using zero-based time mapping
        T.start = decOffset + round((T.start - rawOffset - 1) * targetFs / sRate) + 1;
        T.end   = decOffset + round((T.end   - rawOffset - 1) * targetFs / sRate) + 1;
    
        trialTiming_dec{i} = T;
    end
    
    % Create for_preproc structure
    for_preproc = struct();
    for_preproc.elec_data_raw      = allSignal;
    for_preproc.stitch_index_raw   = stitch_index;
    for_preproc.stitch_index_dec   = stitch_index_dec;
    
    for_preproc.sample_freq_raw    = sRate;
    for_preproc.decimation_freq    = targetFs;          % <<< CHANGED
    for_preproc.decimation_factor  = sRate / targetFs;  % <<< CHANGED
    
    for_preproc.elecs_per_amp      = 64;
    for_preproc.trial_timing_raw   = allTrialTiming;
    for_preproc.trial_timing_dec   = trialTiming_dec;   % <<< CHANGED
    
    % Create ecog_data_seeg object
    obj = ecog_data_seeg(...
        for_preproc, ...
        params.SubjectName, ...
        params.ProtocolName, ...
        fullfile(params.outputPath, [params.SubjectName '_MITLangloc_crunched.mat']), ...
        params.outputPath, ...
        allDataFiles, ...
        '', ...  % Raw path (not applicable)
        chanLabels, ...
        (1:nChannels)', ...
        find(~validChanMask), ...  % Prelim deselect
        chanTypes);
    
    % Set processed data
    obj.elec_data = allSignal;
    obj.stitch_index = stitch_index;
    obj.sample_freq = sRate;
    
    % Set trial info
    obj.trial_timing = allTrialTiming;
    obj.condition = allConditions;
    obj.session = allSessions;
    
    % Mark clean channels
    obj.elec_ch_clean = validChanIdx;
    obj.elec_ch_valid = validChanMask;
    
    fprintf('Created ecog_data_seeg object\n');
    
    
    %% =====================================================================
    % SAVE AS "CRUNCHED" FILE
    % ======================================================================
    savePath = fullfile(params.outputPath, [params.SubjectName '_MITLangloc_crunched.mat']);
    
    fprintf('\n=== SAVING CRUNCHED FILE ===\n');
    fprintf('Saving to: %s\n', savePath);
    
    if ~exist(params.outputPath, 'dir')
        mkdir(params.outputPath);
    end
    
    save(savePath, 'obj', '-v7.3');
    
    fprintf('\n=== DONE ===\n');
    fprintf('You can now run: analyze_MITLangloc(''doOneSub'', ''%s'', ''fromScratch'', true)\n', ...
        params.SubjectName);
end


%% =========================================================================
% HELPER: Extract trials from Brainstorm events
% ==========================================================================
function [trialTiming, conditions] = extract_trials_from_bst(DataMat, params, offsetSamples, sRate)
    
    events = DataMat.Events;
    blockStartTime = DataMat.Time(1);
    blockEndTime = DataMat.Time(end);
    
    fprintf('  Block time range: %.2f - %.2f seconds\n', blockStartTime, blockEndTime);
    
    % Get event pattern from task config
    if isfield(params, 'taskConfig') && isfield(params.taskConfig, 'eventPattern')
        eventPattern = params.taskConfig.eventPattern;
    else
        eventPattern = '%s_tp%d';  % Default pattern
    end
    
    uniqueConditions = {};
    trialsByCondition = {};
    
    for iEvt = 1:length(events)
        label = events(iEvt).label;
        
        % Parse event label based on pattern
        [condName, wordPos] = parseEventLabel(label, eventPattern, params.nWordPositions);
        
        if isempty(condName) || isnan(wordPos)
            continue;
        end
        
        % Get event times
        eventTimes = events(iEvt).times;
        validEvents = (eventTimes >= blockStartTime) & (eventTimes <= blockEndTime);
        eventTimes = eventTimes(validEvents);
        
        if isempty(eventTimes)
            continue;
        end
        
        eventSamplesRelative = round((eventTimes - blockStartTime) * sRate) + 1;
        eventSamplesAbsolute = eventSamplesRelative + offsetSamples;
        
        condIdx = find(strcmp(uniqueConditions, condName));
        if isempty(condIdx)
            uniqueConditions{end+1} = condName;
            trialsByCondition{end+1} = {};
            condIdx = length(uniqueConditions);
        end
        
        for iTrial = 1:length(eventTimes)
            if length(trialsByCondition{condIdx}) < iTrial
                trialsByCondition{condIdx}{iTrial}.wordOnsets = nan(params.nWordPositions, 1);
            end
            trialsByCondition{condIdx}{iTrial}.wordOnsets(wordPos) = eventSamplesAbsolute(iTrial);
        end
    end
    
    % Convert to MIT trial_timing format
    trialTiming = {};
    conditions = {};
    wordDurationSamples = round(params.wordDuration * sRate);
    
    for iCond = 1:length(uniqueConditions)
        condName = uniqueConditions{iCond};
        trials = trialsByCondition{iCond};
        
        for iTrial = 1:length(trials)
            if any(isnan(trials{iTrial}.wordOnsets))
                continue;
            end
            
            keys = arrayfun(@(i) sprintf('word_%d', i), 1:params.nWordPositions, 'UniformOutput', false)';
            strings = arrayfun(@(i) sprintf('word%d', i), 1:params.nWordPositions, 'UniformOutput', false)';
            starts = trials{iTrial}.wordOnsets;
            ends = starts + wordDurationSamples - 1;
            
            timingTable = table(keys, strings, starts, ends, ...
                'VariableNames', {'key', 'string', 'start', 'end'});
            
            trialTiming{end+1, 1} = timingTable;
            conditions{end+1, 1} = condName;
        end
    end
    
    fprintf('  Extracted %d trials from this block\n', length(trialTiming));
end

% Helper function to parse event labels
% Helper function to parse event labels
function [condName, wordPos] = parseEventLabel(label, pattern, maxWords)
    condName = '';
    wordPos = NaN;
    
    % === NEW: Handle simple '%s' pattern (no word position) ===
    if strcmp(pattern, '%s')
        condName = label;
        wordPos = 1;  % Default to position 1 for single-event conditions
        return;
    end
    % === END NEW ===
    
    % Try different parsing strategies based on pattern
    if contains(pattern, '_tp')
        % Pattern like 'sentence_tp1'
        parts = strsplit(label, '_tp');
        if length(parts) == 2
            condName = parts{1};
            wordPos = str2double(parts{2});
        end
    elseif contains(pattern, '_')
        % Pattern like 'condition_1'
        parts = strsplit(label, '_');
        if length(parts) >= 2
            condName = strjoin(parts(1:end-1), '_');
            wordPos = str2double(parts{end});
        end
    else
        % Try regex for general pattern
        tokens = regexp(label, '([a-zA-Z_]+)(\d+)', 'tokens');
        if ~isempty(tokens)
            condName = tokens{1}{1};
            wordPos = str2double(tokens{1}{2});
        end
    end
    
    % Validate word position
    if ~isnan(wordPos) && (wordPos < 1 || wordPos > maxWords)
        wordPos = NaN;
    end
end