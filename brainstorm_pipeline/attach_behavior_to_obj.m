function obj = attach_behavior_to_obj(obj, behaviorFile, varargin)
% ATTACH_BEHAVIOR_TO_OBJ  Load behavioral log and set obj.events_table.
%
%   obj = attach_behavior_to_obj(obj, behaviorFile)
%   obj = attach_behavior_to_obj(obj, behaviorFile, 'alignBy', 'order')
%
%   Required columns (directly or derivable):
%     accuracy  - 0/1, 1 = correct probe response
%     RT        - reaction time in SECONDS (report displays ms)
%
%   If accuracy is missing, it is computed as (response == probe) when both
%   columns exist. RT can be supplied as RT, rt, or reaction_time.
%
%   Alignment (alignBy):
%     'order'         - row i in behavior matches trial i (default)
%     'session_trial' - merge on session + trial columns

p = inputParser();
addParameter(p, 'alignBy', 'order');
addParameter(p, 'sessionCol', 'session');
addParameter(p, 'trialCol', 'trial');
addParameter(p, 'conditionCol', 'condition');
addParameter(p, 'excludeFixation', true);
addParameter(p, 'requireCompleted', false);
addParameter(p, 'completedCol', 'trial_completed');
parse(p, varargin{:});
opts = p.Results;

if isempty(behaviorFile) || ~(ischar(behaviorFile) || isstring(behaviorFile))
    return;
end
behaviorFile = char(behaviorFile);
if ~isfile(behaviorFile)
    error('attach_behavior_to_obj: behavior file not found:\n  %s', behaviorFile);
end

beh = readtable(behaviorFile);
beh = normalize_behavior_table(beh, opts);

nTrials = numel(obj.condition);
switch opts.alignBy
    case 'order'
        events_table = align_behavior_by_order(beh, nTrials, opts);
    case 'session_trial'
        events_table = align_behavior_by_session_trial(beh, obj, opts);
    otherwise
        error('Unknown alignBy: %s', opts.alignBy);
end

obj.events_table = events_table;
fprintf('Attached behavior: %d trials | accuracy %.1f%% | mean RT %.0f ms (correct trials)\n', ...
    height(events_table), 100 * mean(events_table.accuracy), ...
    1000 * mean(events_table.RT(events_table.accuracy == 1), 'omitnan'));
end


function beh = normalize_behavior_table(beh, opts)

if opts.excludeFixation && ismember(opts.conditionCol, beh.Properties.VariableNames)
    condCol = beh.(opts.conditionCol);
    if iscell(condCol)
        isFix = cellfun(@(x) strcmpi(strtrim(x), 'F') || strcmpi(strtrim(x), 'fixation'), condCol);
    else
        isFix = strcmpi(strtrim(string(condCol)), 'F') | strcmpi(strtrim(string(condCol)), 'fixation');
    end
    beh = beh(~isFix, :);
end

if opts.requireCompleted && ismember(opts.completedCol, beh.Properties.VariableNames)
    beh = beh(beh.(opts.completedCol) == 1, :);
end

if ~ismember('accuracy', beh.Properties.VariableNames)
    if ismember('response', beh.Properties.VariableNames) && ismember('probe', beh.Properties.VariableNames)
        beh.accuracy = double(beh.response == beh.probe);
    else
        error(['Behavior table must contain ''accuracy'' or both ''response'' and ''probe''. ' ...
            'Columns: %s'], strjoin(beh.Properties.VariableNames, ', '));
    end
end

if ~ismember('RT', beh.Properties.VariableNames)
    rtCol = find_rt_column(beh.Properties.VariableNames);
    if isempty(rtCol)
        error('Behavior table must contain ''RT'' (seconds) or a recognized RT column.');
    end
    beh.RT = beh.(rtCol);
end

beh.accuracy = double(beh.accuracy ~= 0);
beh.RT = double(beh.RT);
end


function rtCol = find_rt_column(varNames)
candidates = {'RT', 'rt', 'reaction_time', 'ReactionTime', 'react_time'};
rtCol = '';
for i = 1:numel(candidates)
    if ismember(candidates{i}, varNames)
        rtCol = candidates{i};
        return;
    end
end
end


function events_table = align_behavior_by_order(beh, nTrials, opts)

if height(beh) < nTrials
    error(['Behavior file has %d rows after filtering but crunched obj has %d trials. ' ...
        'Check alignment or use alignBy=''session_trial''.'], height(beh), nTrials);
end
if height(beh) > nTrials
    warning('attach_behavior_to_obj: behavior has %d rows; using first %d to match neural trials.', ...
        height(beh), nTrials);
    beh = beh(1:nTrials, :);
end

events_table = table(beh.accuracy, beh.RT, 'VariableNames', {'accuracy', 'RT'});

% Carry optional metadata columns through when present
optionalCols = {'trial', 'list', 'session', 'probe', 'response', 'condition', ...
    'trial_onset', 'trial_completed'};
for i = 1:numel(optionalCols)
    col = optionalCols{i};
    if ismember(col, beh.Properties.VariableNames)
        events_table.(col) = beh.(col);
    end
end
end


function events_table = align_behavior_by_session_trial(beh, obj, opts)

if ~ismember(opts.sessionCol, beh.Properties.VariableNames)
    error('alignBy session_trial requires column ''%s'' in behavior file.', opts.sessionCol);
end
if ~ismember(opts.trialCol, beh.Properties.VariableNames)
    error('alignBy session_trial requires column ''%s'' in behavior file.', opts.trialCol);
end
if isempty(obj.session)
    error('alignBy session_trial requires obj.session to be populated.');
end

nTrials = numel(obj.condition);
accuracy = nan(nTrials, 1);
RT = nan(nTrials, 1);

behSession = beh.(opts.sessionCol);
behTrial = beh.(opts.trialCol);

for i = 1:nTrials
    trialIdx = sum(obj.session(1:i) == obj.session(i));
    match = find(behSession == obj.session(i) & behTrial == trialIdx, 1);
    if isempty(match)
        error('No behavior row for neural trial %d (session=%d, trial=%d).', ...
            i, obj.session(i), trialIdx);
    end
    accuracy(i) = beh.accuracy(match);
    RT(i) = beh.RT(match);
end

events_table = table(accuracy, RT, 'VariableNames', {'accuracy', 'RT'});
end
